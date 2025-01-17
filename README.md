# CDSmsa
Coding sequence MSA (Multiple Sequence Alignment) with optional HMM profile building

## Contents Index
* [MSA example](#msa-example)
* [Usage](#usage)
  * [Single FASTA file](#single-fasta-file-use)
  * [Multiple FASTA files](#multiple-fasta-file-use)

## MSA example
![workflow](https://github.com/ToHanwei/CDSmsa/blob/master/imgs/CDSmsa.jpg)

## Usage
### Single FASTA file use
```python
python single_msa.py --help

usage: single_msa.py [-h] [--msafile] [--cdsfile] [--transcode] [--output] [--build_hmm]

CDS MSA tool (single)

optional arguments:
  -h, --help         show this help message and exit
  --msafile , -m     input file path, it is a MSA output file (FASTA format)
  --cdsfile , -c     nucleotide sequence file (FASTA format)
  --transcode , -t   translate code file path, default='code.csv'
  --output , -o      output file path, CDS alignment file
  --build_hmm        build HMM profile from MSA and save to output directory
```

### Multiple FASTA file use
```python
python batch_msa.py --help

usage: batch_msa.py [-h] [--indir] [--cdsdir] [--transcode] [--outdir] [--build_hmm]

CDS MSA tool (batch)

optional arguments:
  -h, --help         show this help message and exit
  --indir , -i       input folder path, there are MSA output files (FASTA format)
  --cdsdir , -c      save nucleotide sequence files (FASTA format), have the
                     same name with input files
  --transcode , -t   translate code file path
  --outdir , -o      output dir path, save CDS alignment files
  --build_hmm        build HMM profile from MSA and save to output directory
```

### Output Files
The program generates the following output files:
1. MSA files (`.fa` or `.fasta`)
2. HMM profile files (`.hmm`) - when `--build_hmm` option is used
   - HMM profiles are named using the input filename as prefix
   - Saved in the same output directory as MSA files

### Example Commands
```bash
# Single file with HMM profile
python single_msa.py -c input.fa -o ./output --build_hmm

# Batch processing with HMM profiles
python batch_msa.py -i ./input_dir -o ./output_dir --build_hmm
```
