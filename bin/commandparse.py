#!coding: utf-8

import argparse


class Argparse():
    def __init__(self):
        self.single_args = ''
        self.batch_args = ''

    def single_parse(self):
        """
        single file parse
        run a MSA only one time
        """
        parser = argparse.ArgumentParser(description="CDS MSA tool (sigle)")
        parser.add_argument('--cdsfile', '-c', 
                            type=str, metavar='', 
                            help="uncleotide sequence file (FASTA formfat)")
        parser.add_argument('--outdir', '-o', 
                            type=str, 
                            metavar='', 
                            help="output dir path")
        parser.add_argument('--keepfile', '-k',
                            action='count',
                            help="Keep file",
                            default=False,)
        self.single_args = parser.parse_args()

    def batch_parse(self):
        """
        batch runing
        input folder and batch run ever file
        """
        parser = argparse.ArgumentParser(description="CDS MSA tool (batch)")
        parser.add_argument('--cdsdir', '-c', 
                            type=str, metavar='', 
                            help=" save uncleotide sequence files (FASTA formfat), have the same name with input files")
        parser.add_argument('--outdir', '-o', 
                            type=str, metavar='', 
                            help="output dir path, save CDS alignment files")
        parser.add_argument('--keepfile', '-k',
                            action='count',
                            help="Keep file",
                            default=False,)
        self.batch_args = parser.parse_args()
