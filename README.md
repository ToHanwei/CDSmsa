# CDSmsa
coding sequence MSA（multiple sequence alignment）

## Contents Index
* [MSA example](#msa-example)
* [Usage](#usage)

## MSA example
![workfolw](https://github.com/ToHanwei/CDSmsa/blob/master/imgs/CDSmsa.jpg)

## Usage
### Single FASTA file use
```python
python single_msa.py --help

usage: single_msa.py [-h] [--msafile] [--cdsfile] [--transcode] [--output]

CDS MSA tool (sigle)

optional arguments:
  -h, --help         show this help message and exit
  --msafile , -m     inputfile path, it is a MSA output file (FASTA format)
  --cdsfile , -c     uncleotide sequence file (FASTA formfat)
  --transcode , -t   translate code file path, default='code.csv'
  --output , -o      output file path, CDS alignment file
```
### Multiple FASTA file use
```python
python batch_msa.py --help

usage: batch_msa.py [-h] [--indir] [--cdsdir] [--transcode] [--outdir]

CDS MSA tool (batch)

optional arguments:
  -h, --help         show this help message and exit
  --indir , -i       input folder path, there are MSA output files (FASTA format)
  --cdsdir , -c      save uncleotide sequence files (FASTA formfat), have the
                     same name with input files
  --transcode , -t   translate code file path
  --outdir , -o      output dir path, save CDS alignment files
```
