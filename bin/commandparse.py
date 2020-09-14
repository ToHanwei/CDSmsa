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
        parser = argparse.ArgumentParser(description="CDS MSA tool")
        parser.add_argument('--msafile', '-m', type=str, metavar='', help="inputfile path, it is a MSA output file (FASTA format)")
        parser.add_argument('--cdsfile', '-c', type=str, metavar='', help="uncleotide sequence file (FASTA formfat)")
        parser.add_argument('--transcode', '-t', type=str, metavar='', help="translate code file path", default="./code.csv")
        parser.add_argument('--output', '-o', type=str, metavar='', help="output file path, CDS alignment file", default="output.fas")
        self.single_args = parser.parse_args()

    def batch_parse(self):
        """
        batch runing
        """
        parser = argparse.ArgumentParser(description="CDS MSA tool")
        parser.add_argument('--msafile', '-m', type=str, metavar='', help="inputfile path, it is a MSA output file (FASTA format)")
        parser.add_argument('--cdsfile', '-c', type=str, metavar='', help="uncleotide sequence file (FASTA formfat)")
        parser.add_argument('--transcode', '-t', type=str, metavar='', help="translate code file path", default="code.csv")
        self.batch_args = parser.parse_args()

        

