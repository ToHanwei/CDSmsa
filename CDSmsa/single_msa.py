#!coding: utf-8

import os
import sys

from bin.commandparse import Argparse
from bin.CDSalignment import splitseq
from bin.CDSalignment import mapseq
from bin.CDSalignment import write2file
from bin.CDSalignment import run_ptot_msa
from bin.CDSalignment import build_hmm_profile


def commands() -> tuple[str, str, bool, bool]:
    """
    Parse command line arguments
    
    Returns:
        tuple: (cdsfile, outdir, keepfile, build_hmm)
    """
    if len(sys.argv) < 2:
        os.system("python single_msa.py -h")
        sys.exit()
    ArgParse = Argparse()
    ArgParse.single_parse()
    args = ArgParse.single_args
    cdsfile = args.cdsfile
    keepfile = args.keepfile
    outdir = args.outdir
    build_hmm = args.build_hmm
    return cdsfile, outdir, keepfile, build_hmm


def main():
    # Parse command line arguments
    cdsfile, outdir, keepfile, build_hmm = commands()
    
    # Run existing MSA pipeline
    msas_dict, uncl_dict, pfile, mfile = run_ptot_msa(cdsfile, outdir)
    uncl_dict = splitseq(uncl_dict)
    resu_dict = mapseq(msas_dict, uncl_dict)
    outfile = os.path.join(outdir, os.path.basename(cdsfile))
    write2file(resu_dict, outfile)
    
    # Build HMM profile if requested
    if build_hmm:
        hmm_out = os.path.join(outdir, f"{os.path.splitext(os.path.basename(cdsfile))[0]}.hmm")
        build_hmm_profile(outfile, hmm_out)
    
    # Clean up temporary files
    if not keepfile:
        os.remove(pfile)
        os.remove(mfile)


if __name__ == "__main__":
    main()


