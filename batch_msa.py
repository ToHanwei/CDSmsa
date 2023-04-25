#!coding: utf-8

import os
import sys

from bin.commandparse import Argparse
from bin.CDSalignment import splitseq
from bin.CDSalignment import mapseq
from bin.CDSalignment import write2file
from bin.CDSalignment import run_ptot_msa


def commands():
    """
    command lin eparser
    """
    if len(sys.argv) < 2:
        # no command parma
        os.system("python batch_msa.py -h")
        sys.exit()
    ArgParse = Argparse()
    ArgParse.batch_parse()
    args = ArgParse.batch_args
    cdsdir = args.cdsdir
    outdir = args.outdir
    keepfile = args.keepfile
    return cdsdir, outdir, keepfile


def run_single(cdsfile, outdir, keepfile):
    """
    run single file
    """
    msas_dict, uncl_dict, pfile, mfile = run_ptot_msa(cdsfile, outdir)
    uncl_dict = splitseq(uncl_dict)
    resu_dict = mapseq(msas_dict, uncl_dict)
    outfile = os.path.join(outdir, os.path.basename(cdsfile))
    write2file(resu_dict, outfile)
    if not keepfile:
        os.remove(pfile)
        os.remove(mfile)


def main():
    cdsdir, outdir, keepfile = commands()
    cdsfiles = os.listdir(cdsdir)
    for _file in cdsfiles:
        cdsfile = os.path.join(cdsdir, _file)
        run_single(cdsfile, outdir, keepfile)


if __name__ == "__main__":
    main()

