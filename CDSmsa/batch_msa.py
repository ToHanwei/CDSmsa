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
    Parse command line arguments for batch processing
    
    Returns:
        tuple: (indir, outdir, keepfile, build_hmm)
    """
    if len(sys.argv) < 2:
        os.system("python batch_msa.py -h")
        sys.exit()
    ArgParse = Argparse()
    ArgParse.batch_parse()
    args = ArgParse.batch_args
    indir = args.indir
    outdir = args.outdir
    keepfile = args.keepfile
    build_hmm = args.build_hmm
    return indir, outdir, keepfile, build_hmm

def main():
    # Parse command line arguments
    indir, outdir, keepfile, build_hmm = commands()
    
    # Get list of valid input files
    input_files = [f for f in os.listdir(indir) if f.endswith(('.fa', '.fasta'))]
    total_files = len(input_files)
    
    print(f"Found {total_files} files to process")
    for idx, cdsfile in enumerate(input_files, 1):
        print(f"\nProcessing file {idx}/{total_files}: {cdsfile}")
        
        cds_path = os.path.join(indir, cdsfile)
        
        # Step 1: Run MSA pipeline
        print("Step 1: Running multiple sequence alignment...")
        msas_dict, uncl_dict, pfile, mfile = run_ptot_msa(cds_path, outdir)
        
        # Step 2: Split sequences
        print("Step 2: Splitting sequences...")
        uncl_dict = splitseq(uncl_dict)
        
        # Step 3: Map sequences
        print("Step 3: Mapping sequences...")
        resu_dict = mapseq(msas_dict, uncl_dict)
        
        # Step 4: Write results
        print("Step 4: Writing results to file...")
        outfile = os.path.join(outdir, cdsfile)
        write2file(resu_dict, outfile)
        
        # Step 5: Build HMM profile if requested
        if build_hmm:
            print("Step 5: Building HMM profile...")
            hmm_out = os.path.join(outdir, f"{os.path.splitext(cdsfile)[0]}.hmm")
            build_hmm_profile(outfile, hmm_out)
        
        # Clean up temporary files
        if not keepfile:
            print("Cleaning up temporary files...")
            os.remove(pfile)
            os.remove(mfile)
        
        print(f"Completed processing {cdsfile}")

    print("\nAll files processed successfully!")

if __name__ == "__main__":
    main()

