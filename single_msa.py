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
        os.system("python single_msa.py -h")
        sys.exit()
    ArgParse = Argparse()
    ArgParse.single_parse()
    args = ArgParse.single_args
    # msafile = args.msafile
    cdsfile = args.cdsfile
    keepfile = args.keepfile
    outdir = args.outdir
    return cdsfile, outdir, keepfile


def main():
    cdsfile, outdir, keepfile = commands()
    msas_dict, uncl_dict, pfile, mfile = run_ptot_msa(cdsfile, outdir)
    uncl_dict = splitseq(uncl_dict)
    resu_dict = mapseq(msas_dict, uncl_dict)
    outfile = os.path.join(outdir, os.path.basename(cdsfile))
    write2file(resu_dict, outfile)
    if not keepfile:
        os.remove(pfile)
        os.remove(mfile)


if __name__ == "__main__":
    main()


