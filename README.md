# ReadAngler
Extract FASTQ records whose headers are present in an (MGE) FASTA alignment output during a BeeGees pipeline run.

## Installation

read_angler requires [seqkit](https://bioinf.shenwei.me/seqkit/) >= 2.0, GNU coreutils, and awk.

To install Seqkit:
```bash
conda install bioconda::seqkit
```

GNU coreutils and awk are pre-installed on most Linux systems and macOS.

### Getting started

Clone or download `ReadAngler.sh` and run it directly.
```bash
git clone https://github.com/museomics/ReadAngler.git
cd ReadAngler
bash ReadAngler.sh -h
```

## Overview

`read_angler.sh` filters FASTQ files to retain only reads whose IDs appear in a FASTA alignment. Two matching modes are available:

- **Strict mode** (default): only reads whose IDs are present in the FASTA are extracted. R1 and R2 outputs are not guaranteed to be in sync.
- **Pair rescue mode** (`-k`): matching is performed on the FRAGMENT key (the read header with the trailing `_1`/`_2` read-number stripped). Any fragment matched by either mate is pulled from both FASTQs, so R1 and R2 are equal in size and in sync.

Read-pair assignment uses the read-number field of the AVITI/bases2fastq header (e.g. `...:0526_1:N:0:IDX+IDX` -> R1). If no such field is found, all IDs are treated as single-end.

> [!NOTE]
> Matching is STRICT by default

## Usage

### Single sample

```bash
sbatch read_angler.sh -f aln.fasta -1 sample_R1.fastq.gz \
    [-2 sample_R2.fastq.gz] -o out_dir
```

The `-2` argument is optional; if omitted, the mate FASTQ is auto-detected from the `-1` filename, and the run falls back to single-end if none is found.

### Batch mode (SLURM array)

```bash
sbatch --array=0-$((N-1)) read_angler.sh \
    -F fasta_dir -Q fastq_dir -o out_dir
```

Samples are matched by FASTA basename (minus extension and any trailing suffix given with `-x`) against FASTQ files in `-Q`. Without `--array`, batch mode loops over every sample serially.

## Options

| Option | Description |
|--------|-------------|
| `-f FILE` | Input FASTA alignment (single-sample mode) |
| `-1 FILE` | Input R1 FASTQ(.gz) (single-sample mode) |
| `-2 FILE` | Input R2 FASTQ(.gz) (optional) |
| `-F DIR` | Directory of FASTA alignments (batch mode) |
| `-Q DIR` | Directory of FASTQ files (batch mode) |
| `-x STR` | Suffix to strip from FASTA basename to get sample name (e.g. `-x '_aligned'`) |
| `-o DIR` | Output directory (default: `./extracted_reads`) |
| `-t INT` | Threads (default: `$SLURM_CPUS_PER_TASK` or 4) |
| `-e STR` | Conda environment to activate (default: `read_qc`) |
| `-k` | Pair rescue: keep complete pairs even when only one mate is in the alignment |
| `-i` | Keep intermediate ID list files |
| `-v` | Verbose: log a per-fragment status table and the full read headers written to each output file |
| `-h` | Show help |