#!coding: utf-8

import os
import sys

from bin.commandparse import Argparse
from bin.CDSalignment import makedict
from bin.CDSalignment import buildcode
from bin.CDSalignment import splitseq
from bin.CDSalignment import mapseq
from bin.CDSalignment import write2file


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
    indir = args.indir
    cdsdir = args.cdsdir
    codefile = args.transcode
    outdir = args.outdir
    return indir, cdsdir, codefile, outdir

def run_single(msafile, cdsfile, codefile, outfile):
    """
    run single file
    """
    prot_dict = makedict(msafile)
    uncl_dict = makedict(cdsfile)
    code_dict = buildcode(codefile)
    uncl_dict = splitseq(uncl_dict)
    resu_dict = mapseq(prot_dict, uncl_dict, code_dict)
    write2file(resu_dict, outfile)


def main():
    indir, cdsdir, codefile, outdir = commands()
    infiles = os.listdir(indir)
    cdsfiles = os.listdir(cdsdir)
    assert set(infiles) == set(cdsfiles)
    for _file in infiles:
        msafile = os.path.join(indir, _file)
        cdsfile = os.path.join(cdsdir, _file)
        outfile = os.path.join(outdir, _file)
        run_single(msafile, cdsfile, codefile, outfile)


if __name__ == "__main__":
    main()

