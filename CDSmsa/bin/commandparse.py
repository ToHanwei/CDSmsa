#!coding: utf-8

import argparse


class Argparse():
    def __init__(self):
        self.single_args = ''
        self.batch_args = ''
        self.single_parser = argparse.ArgumentParser(
            description='Single sequence MSA'
        )

    def single_parse(self):
        """
        single file parse
        run a MSA only one time
        """
        self.single_parser.add_argument('-c', '--cdsfile',
                                      required=True,
                                      help='Input CDS file')
        self.single_parser.add_argument('-o', '--outdir',
                                      required=True,
                                      help='Output directory')
        self.single_parser.add_argument('-k', '--keepfile',
                                      action='store_true',
                                      help='Keep temporary files')
        self.single_parser.add_argument('--build_hmm', 
                                      action='store_true',
                                      default=False,
                                      help='Build HMM profile from MSA and save to output directory')
        
        self.single_args = self.single_parser.parse_args()

    def batch_parse(self):
        """Parse arguments for batch processing"""
        self.batch_parser = argparse.ArgumentParser(
            description='Batch sequence MSA'
        )
        self.batch_parser.add_argument('-i', '--indir',
                                      required=True,
                                      help='Input directory containing CDS files')
        self.batch_parser.add_argument('-o', '--outdir',
                                      required=True,
                                      help='Output directory')
        self.batch_parser.add_argument('-k', '--keepfile',
                                      action='store_true',
                                      help='Keep temporary files')
        self.batch_parser.add_argument('--build_hmm',
                                      action='store_true',
                                      default=False,
                                      help='Build HMM profile from MSA and save to output directory')
        
        self.batch_args = self.batch_parser.parse_args()
