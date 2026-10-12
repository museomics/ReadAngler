#!/bin/bash
#SBATCH --job-name=read_angler
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=read_angler_%A_%a.out
#SBATCH --error=read_angler_%A_%a.err

# ---------------------------------------------------------------------------
# read_angler.sh
#
# Extract FASTQ records whose headers are present in a FASTA alignment.
#
# Matching is STRICT by default: only read IDs that appear in the FASTA are
# pulled, so R1 and R2 outputs are NOT guaranteed to be in sync.
#
# With -k (pair rescue), matching is done on the FRAGMENT key (the header with
# the trailing _1 / _2 read-number stripped). Any fragment matched by either
# mate is pulled from both FASTQs, so R1 and R2 are equal in size and in sync.
#
# Read-pair assignment is taken from the read-number field of the AVITI /
# bases2fastq header (e.g. ...:0526_1:N:0:IDX+IDX -> R1). If no such field is
# found, all IDs are treated as single-end.
#
# Requires: seqkit (>=2.0), GNU coreutils, awk
#
# Usage (single sample):
#   sbatch read_angler.sh -f aln.fasta -1 sample_R1.fastq.gz \
#       [-2 sample_R2.fastq.gz] -o out_dir
#     (-2 is optional; if omitted the mate is auto-detected from the -1
#      filename, and the run falls back to single-end if none is found)
#
# Usage (batch / array over many samples):
#   sbatch --array=0-$((N-1)) read_angler.sh \
#       -F fasta_dir -Q fastq_dir -o out_dir
#     Samples are matched by the FASTA basename (minus extension and any
#     trailing suffix given with -x) against FASTQ files in -Q.
#     Without --array, batch mode loops over every sample serially.
#
# Options:
#   -f FILE   input FASTA alignment (single-sample mode)
#   -1 FILE   input R1 FASTQ(.gz) (single-sample mode)
#   -2 FILE   input R2 FASTQ(.gz) (optional)
#   -F DIR    directory of FASTA alignments (batch mode)
#   -Q DIR    directory of FASTQ files (batch mode)
#   -x STR    suffix to strip from FASTA basename to get sample name
#             (e.g. -x '_aligned') [default: none]
#   -o DIR    output directory [default: ./extracted_reads]
#   -t INT    threads [default: SLURM_CPUS_PER_TASK or 4]
#   -e STR    conda environment to activate [default: read_qc]
#   -k        pair rescue: keep complete pairs even when only one mate is
#             present in the alignment [default: strict, off]
#   -i        keep intermediate ID list files
#   -v        verbose: log a per-fragment status table and the full read
#             headers written to each output file
#   -h        show this help
# ---------------------------------------------------------------------------

set -euo pipefail

FASTA=""; R1=""; R2=""
FASTA_DIR=""; FASTQ_DIR=""; SAMPLE_SUFFIX=""
OUTDIR="./extracted_reads"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
CONDA_ENV="read_qc"
KEEP_IDS=0
KEEP_PAIRS=0
VERBOSE=0

usage() { sed -n '/^# ----/,/^# ----$/p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

while getopts "f:1:2:F:Q:x:o:t:e:kivh" opt; do
    case "$opt" in
        f) FASTA="$OPTARG" ;;
        1) R1="$OPTARG" ;;
        2) R2="$OPTARG" ;;
        F) FASTA_DIR="$OPTARG" ;;
        Q) FASTQ_DIR="$OPTARG" ;;
        x) SAMPLE_SUFFIX="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        e) CONDA_ENV="$OPTARG" ;;
        k) KEEP_PAIRS=1 ;;
        i) KEEP_IDS=1 ;;
        v) VERBOSE=1 ;;
        h) usage 0 ;;
        *) usage 1 ;;
    esac
done

log()  { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
die()  { printf '[ERROR] %s\n' "$*" >&2; exit 1; }
warn() { printf '[WARN]  %s\n' "$*" >&2; }

# --- environment -----------------------------------------------------------
if [[ -n "$CONDA_ENV" ]]; then
    # shellcheck disable=SC1091
    if [[ -n "${CONDA_EXE:-}" ]]; then
        source "$(dirname "$(dirname "$CONDA_EXE")")/etc/profile.d/conda.sh"
    elif command -v conda >/dev/null 2>&1; then
        source "$(conda info --base)/etc/profile.d/conda.sh"
    else
        die "conda not found but -e '$CONDA_ENV' was requested"
    fi
    conda activate "$CONDA_ENV" || die "could not activate conda env '$CONDA_ENV'"
fi

command -v seqkit >/dev/null 2>&1 || die "seqkit not found in PATH (use -e to activate a conda env)"

# --- mode resolution -------------------------------------------------------
MODE=""
if [[ -n "$FASTA" || -n "$R1" ]]; then
    [[ -n "$FASTA" && -n "$R1" ]] || die "single-sample mode needs both -f and -1"
    [[ -n "$FASTA_DIR" || -n "$FASTQ_DIR" ]] && die "-f/-1 cannot be combined with -F/-Q"
    MODE="single"
elif [[ -n "$FASTA_DIR" || -n "$FASTQ_DIR" ]]; then
    [[ -n "$FASTA_DIR" && -n "$FASTQ_DIR" ]] || die "batch mode needs both -F and -Q"
    MODE="batch"
else
    usage 1
fi

mkdir -p "$OUTDIR" logs

# --- helpers ---------------------------------------------------------------

# awk function shared by the ID-splitting and pair-rescue steps. Returns the
# 'fragment key' for a read ID (the header up to and including the read-number
# field, with the trailing _1 / _2 removed).
AWK_FRAG='function frag(id,   a, n, i, key) {
    n = split(id, a, ":")
    if (n >= 4) { key = a[1]; for (i = 2; i <= n - 3; i++) key = key ":" a[i] }
    else        { key = id }
    sub(/_[12]$/, "", key)
    return key
}'

# Find the R2 partner of an R1 filename. Returns empty string if none exists.
find_mate() {
    local r1="$1" dir base cand
    dir="$(dirname "$r1")"; base="$(basename "$r1")"
    for pat in '_R1_:_R2_' '_R1.:_R2.' '.R1.:.R2.' '_R1_001:_R2_001' '_1.:_2.' '_r1_:_r2_' '_r1.:_r2.'; do
        local from="${pat%%:*}" to="${pat##*:}"
        if [[ "$base" == *"$from"* ]]; then
            cand="$dir/${base/$from/$to}"
            [[ -f "$cand" ]] && { printf '%s' "$cand"; return 0; }
        fi
    done
    printf ''
}

# Pull IDs from a FASTA and split them into R1 / R2 / unpaired lists.
# The read-number field is the 4th-from-last colon-delimited field of the
# bases2fastq header (…:0526_1:N:0:IDX+IDX). Duplicate IDs are collapsed.
split_ids() {
    local fasta="$1" prefix="$2"
    local all="${prefix}.all.ids"

    if [[ "$fasta" == *.gz ]]; then gzip -cd -- "$fasta"; else cat -- "$fasta"; fi \
        | awk '/^>/ { sub(/^>/, ""); sub(/[ \t].*$/, ""); if (length($0)) print }' \
        | LC_ALL=C sort -u > "$all"

    [[ -s "$all" ]] || die "no sequence headers found in $fasta"

    awk -F: -v r1="${prefix}.R1.ids" -v r2="${prefix}.R2.ids" -v un="${prefix}.SE.ids" '
        {
            n = (NF >= 4) ? $(NF-3) : $0
            if      (n ~ /_1$/) print > r1
            else if (n ~ /_2$/) print > r2
            else                print > un
        }' "$all"

    for f in "${prefix}.R1.ids" "${prefix}.R2.ids" "${prefix}.SE.ids"; do
        [[ -f "$f" ]] || : > "$f"
    done

    # No read-number field anywhere -> treat the whole set as single-end.
    if [[ ! -s "${prefix}.R1.ids" && ! -s "${prefix}.R2.ids" ]]; then
        cp "$all" "${prefix}.SE.ids"
        : > "${prefix}.R1.ids"; : > "${prefix}.R2.ids"
    fi
}

# Pair rescue. Collapse the alignment IDs to 'fragment keys', then scan the real
# headers of each FASTQ and keep every header whose fragment key was matched.
# This gives exact full-header ID lists for seqkit grep without having to guess
# what an absent mate's header looked like.
rescue_ids() {
    local prefix="$1" r1="$2" r2="$3"
    local keys="${prefix}.frag.keys"

    awk "$AWK_FRAG"' { print frag($0) }' "${prefix}.all.ids" | LC_ALL=C sort -u > "$keys"
    log "pair rescue: $(wc -l < "$keys") unique fragment key(s) in alignment"

    local f
    for f in R1 R2; do
        local fq; [[ "$f" == R1 ]] && fq="$r1" || fq="$r2"
        [[ -n "$fq" && -f "$fq" ]] || continue
        seqkit seq -n -j "$THREADS" "$fq" \
            | awk "$AWK_FRAG"'
                NR == FNR { want[$0]; next }
                {
                    id = $0; sub(/[ \t].*$/, "", id)
                    if (frag(id) in want) print $0
                }' "$keys" - > "${prefix}.${f}.ids"
    done
}

# Verbose: per-fragment status table built from the alignment IDs.
status_table() {
    local prefix="$1"
    log "  fragment status table (alignment: $(wc -l < "${prefix}.all.ids") header(s)):"
    awk -v rescue="$KEEP_PAIRS" "$AWK_FRAG"'
        {
            id = $0
            n = split(id, a, ":")
            num = (n >= 4) ? a[n-3] : id
            k = frag(id)
            if (!(k in seen)) { seen[k] = 1; order[++c] = k }
            if      (num ~ /_1$/) r1[k] = 1
            else if (num ~ /_2$/) r2[k] = 1
            else                  un[k] = 1
        }
        END {
            for (i = 1; i <= c; i++) {
                k = order[i]
                if      (k in un)                 s = "unassigned (no read-number field)"
                else if ((k in r1) && (k in r2)) {s = "R1+R2";   both++}
                else if (k in r1)                {s = rescue ? "R1-only -> R2 rescued" : "R1-only (dropped from R2)"; only1++}
                else                             {s = rescue ? "R2-only -> R1 rescued" : "R2-only (dropped from R1)"; only2++}
                printf "    %-58s %s\n", k, s
            }
            printf "    ---- %d fragment(s): %d complete, %d R1-only, %d R2-only\n", c, both, only1, only2
        }' "${prefix}.all.ids"
}

# Verbose: dump the read headers actually written to an output file.
log_headers() {
    local out="$1" label="$2"
    log "  $label: headers written to $(basename "$out"):"
    seqkit seq -n "$out" | sed 's/^/    /'
}

# seqkit grep + record-count sanity check.
extract() {
    local fastq="$1" idlist="$2" out="$3" label="$4"
    local wanted found
    wanted=$(wc -l < "$idlist")

    if [[ "$wanted" -eq 0 ]]; then
        warn "$label: no IDs to extract, skipping"
        return 0
    fi

    log "$label: extracting $wanted read(s) from $(basename "$fastq")"
    seqkit grep -n -f "$idlist" -j "$THREADS" "$fastq" -o "$out"

    found=$(seqkit stats -T "$out" | awk 'NR==2 {print $4}')
    log "$label: wrote $found record(s) -> $out"
    if [[ "$found" -ne "$wanted" ]]; then
        warn "$label: $wanted ID(s) requested but $found record(s) found — check that the FASTA headers came from this FASTQ (and were not renamed by the aligner)"
    fi

    [[ "$VERBOSE" -eq 1 ]] && log_headers "$out" "$label"
    return 0
}

process_sample() {
    local sample="$1" fasta="$2" r1="$3" r2="$4"
    local prefix="$OUTDIR/${sample}.ids"

    log "=== $sample ==="
    [[ -f "$fasta" ]] || die "$sample: FASTA not found: $fasta"
    [[ -f "$r1" ]]    || die "$sample: R1 FASTQ not found: $r1"

    split_ids "$fasta" "$prefix"
    log "$sample: $(wc -l < "${prefix}.R1.ids") R1 / $(wc -l < "${prefix}.R2.ids") R2 / $(wc -l < "${prefix}.SE.ids") unassigned (no read-number field) ID(s) in alignment"
    [[ "$VERBOSE" -eq 1 ]] && status_table "$prefix"

    if [[ "$KEEP_PAIRS" -eq 1 ]]; then
        if [[ -s "${prefix}.SE.ids" ]]; then
            warn "$sample: -k has no effect on unpaired/single-end IDs"
        fi
        if [[ -z "$r2" || ! -f "$r2" ]]; then
            warn "$sample: -k requested but no R2 FASTQ found — falling back to strict"
        else
            rescue_ids "$prefix" "$r1" "$r2"
            local n1 n2
            n1=$(wc -l < "${prefix}.R1.ids"); n2=$(wc -l < "${prefix}.R2.ids")
            log "$sample: after rescue, $n1 R1 / $n2 R2 ID(s)"
            [[ "$n1" -eq "$n2" ]] || warn "$sample: rescued R1 ($n1) and R2 ($n2) counts differ — some mates are missing from the FASTQ files themselves"
        fi
    fi

    if [[ -s "${prefix}.SE.ids" ]]; then
        extract "$r1" "${prefix}.SE.ids" "$OUTDIR/${sample}.extracted.fastq.gz" "$sample SE"
    fi
    if [[ -s "${prefix}.R1.ids" ]]; then
        extract "$r1" "${prefix}.R1.ids" "$OUTDIR/${sample}_R1.extracted.fastq.gz" "$sample R1"
    fi
    if [[ -s "${prefix}.R2.ids" ]]; then
        if [[ -n "$r2" && -f "$r2" ]]; then
            extract "$r2" "${prefix}.R2.ids" "$OUTDIR/${sample}_R2.extracted.fastq.gz" "$sample R2"
        else
            warn "$sample: alignment contains R2 IDs but no R2 FASTQ was found — skipping"
        fi
    fi

    if [[ "$KEEP_IDS" -eq 0 ]]; then
        rm -f "${prefix}".all.ids "${prefix}".R1.ids "${prefix}".R2.ids "${prefix}".SE.ids "${prefix}".frag.keys
    fi
}

# --- run -------------------------------------------------------------------
log "seqkit: $(seqkit version 2>&1 | head -1); threads: $THREADS; mode: $MODE; pair rescue: $([[ "$KEEP_PAIRS" -eq 1 ]] && echo on || echo off)"

if [[ "$MODE" == "single" ]]; then
    if [[ -z "$R2" ]]; then
        R2="$(find_mate "$R1")"
        [[ -n "$R2" ]] && log "auto-detected mate: $R2" || log "no mate found — treating input as single-end"
    fi
    sample="$(basename "$FASTA")"; sample="${sample%.gz}"
    sample="${sample%.*}"; sample="${sample%$SAMPLE_SUFFIX}"
    process_sample "$sample" "$FASTA" "$R1" "$R2"
else
    mapfile -t FASTAS < <(find "$FASTA_DIR" -maxdepth 1 -type f \
        \( -name '*.fa' -o -name '*.fas' -o -name '*.fasta' -o -name '*.fna' \
           -o -name '*.fa.gz' -o -name '*.fas.gz' -o -name '*.fasta.gz' \) | LC_ALL=C sort)
    [[ "${#FASTAS[@]}" -gt 0 ]] || die "no FASTA files found in $FASTA_DIR"
    log "found ${#FASTAS[@]} FASTA file(s) in $FASTA_DIR"

    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        idx="$SLURM_ARRAY_TASK_ID"
        [[ "$idx" -lt "${#FASTAS[@]}" ]] || die "array index $idx exceeds sample count ${#FASTAS[@]}"
        TARGETS=("${FASTAS[$idx]}")
    else
        TARGETS=("${FASTAS[@]}")
        warn "not a SLURM array job — processing all ${#FASTAS[@]} sample(s) serially"
    fi

    for fa in "${TARGETS[@]}"; do
        sample="$(basename "$fa")"; sample="${sample%.gz}"
        sample="${sample%.*}"; sample="${sample%$SAMPLE_SUFFIX}"

        mapfile -t hits < <(find "$FASTQ_DIR" -maxdepth 1 -type f \
            \( -name "${sample}*_R1*.f*q*" -o -name "${sample}*_1.f*q*" \) | LC_ALL=C sort)
        if [[ "${#hits[@]}" -eq 0 ]]; then
            mapfile -t hits < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name "${sample}*.f*q*" | LC_ALL=C sort)
        fi
        [[ "${#hits[@]}" -gt 0 ]] || { warn "$sample: no FASTQ found in $FASTQ_DIR — skipping"; continue; }
        [[ "${#hits[@]}" -gt 1 ]] && warn "$sample: multiple R1 candidates, using ${hits[0]}"

        r1="${hits[0]}"
        r2="$(find_mate "$r1")"
        process_sample "$sample" "$fa" "$r1" "$r2"
    done
fi

log "done. outputs in $OUTDIR"
