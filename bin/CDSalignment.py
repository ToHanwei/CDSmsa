import os
import re
import sys

from subprocess import Popen


def aad2codes():
    """get uncleotide code table"""
    codes = {
                "A": ["GCT", "GCC", "GCA", "GCG"], 
                "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"], 
                "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"], 
                "K": ["AAA", "AAG"], 
                "N": ["AAT", "AAC"], 
                "M": ["ATG"], 
                "D": ["GAT", "GAC"], 
                "F": ["TTT", "TTC"], 
                "C": ["TGT", "TGC"], 
                "P": ["CCT", "CCC", "CCA", "CCG"], 
                "Q": ["CAA", "CAG"], 
                "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"], 
                "E": ["GAA", "GAG"], 
                "T": ["ACT", "ACC", "ACA", "ACG"], 
                "G": ["GGT", "GGC", "GGA", "GGG"], 
                "W": ["TGG"], 
                "H": ["CAT", "CAC"], 
                "Y": ["TAT", "TAC"], 
                "I": ["ATT", "ATC", "ATA"], 
                "V": ["GTT", "GTC", "GTA", "GTG"], 
            }
    return codes


def get_codon_table():
    """get codes, translate code to amino acid"""
    CODON_TABLE = {
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'CGT': 'R',
        'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S',
        'AGC': 'S', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'TTA': 'L',
        'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'GTT': 'V',
        'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'ACT': 'T', 'ACC': 'T',
        'ACA': 'T', 'ACG': 'T', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
        'CCG': 'P', 'AAT': 'N', 'AAC': 'N', 'GAT': 'D', 'GAC': 'D',
        'TGT': 'C', 'TGC': 'C', 'CAA': 'Q', 'CAG': 'Q', 'GAA': 'E',
        'GAG': 'E', 'CAT': 'H', 'CAC': 'H', 'AAA': 'K', 'AAG': 'K',
        'TTT': 'F', 'TTC': 'F', 'TAT': 'Y', 'TAC': 'Y', 'ATG': 'M',
        'TGG': 'W', 'TAG': '*', 'TGA': '*', 'TAA': '*'
    }
    return CODON_TABLE


def translate(cdsdict:dict) -> dict:
    """traslate dna to protein"""
    codon_table = get_codon_table()
    prot_dict = {}
    for name, seq in cdsdict.items():
        codes = re.findall("...", seq)
        prots = [codon_table.get(c, "X") for c in codes]
        prots = "".join(prots)
        prot_dict[name] = prots
    return prot_dict    


def sequence_align(seqf, alignf):
    """
    Function:
        sequence_align
        call MAFFT(linsi) multiple sequence alignment
    Parameter:
        seqf, sequence file name (FASTA format)
        alignf, Multiple sequence alignment results
    Return: None
    """

    command = ('mafft'
               # run linsi
               #+ " --localpair "
               #+ "--maxiterate 1000 "
               + " --thread -1 "
               + "--quiet "
               + seqf
               + " > "
               + alignf
    )
    p = Popen(command, shell=True)
    p.wait()


def run_ptot_msa(cdsfile:str, outdir:str) -> dict:
    """run MSA for protein sequences"""
    uncl_dict = makedict(cdsfile)
    prot_dict = translate(uncl_dict)
    # reformat FASTA
    prots = [f">{name}\n{seq}\n" for name, seq in prot_dict.items()]
    prefix = os.path.basename(cdsfile)
    prefix = os.path.splitext(prefix)[0]
    protfile = os.path.join(outdir, f"{prefix}_ptot.fasta")
    msasfile = os.path.join(outdir, f"{prefix}_msas.fasta")
    # write protein sequences to file
    with open(protfile, "w") as protf:
        protf.writelines(prots)
    # run msa
    sequence_align(protfile, msasfile)
    msas_dict = makedict(msasfile)

    return msas_dict, uncl_dict, protfile, msasfile

def makedict(infile):
    """build dict from sequences"""
    seq_dict = {}
    with open(infile) as seqf:
        seqs = seqf.read().split(">")[1:]
    for seq in seqs:
        lines = seq.split("\n")
        name = lines[0]
        seq = "".join(lines[1:])
        seq_dict[name] = seq

    return seq_dict


def splitseq(indict, num=3):
    """split sequence by step"""
    names = tuple(indict.keys())
    outdict = {}
    for name in names:
        seq = indict[name]
        Length = len(seq)
        assert Length % 3 == 0
        #indict[name] = [seq[i:i+3] for i in range(0, Length, 3)]
        codes = re.findall("...", seq)
        if codes[-1] in ["TAG", "TGA", "TAA"]:
            codes.remove(codes[-1])
        outdict[name] = codes
    return outdict


def mapseq(pdict, udict):
    """map MSA protein sequence to uncleotide"""
    cdict = aad2codes()
    names = tuple(pdict.keys())
    mapdict = {}
    for name in names:
        i, mseq, pseq, useq = 0, "", pdict[name], udict[name]
        for aad in pseq:
            if aad == "-":
                uncl = "-" * 3
            else:
                uncl = useq[i]
                assert uncl in cdict[aad]
                i += 1
            mseq += uncl
        mapdict[name] = mseq
    return mapdict


def write2file(indict, outfile):
    """write map sequence to file"""
    with open(outfile, "w") as out:
        for name in indict.keys():
            seq = ">" + name + "\n" + indict[name]+ "\n"
            out.write(seq)


#def main():
#    # protein sequence MSA file
#    prot_seqfile = sys.argv[1]
#    # uncleotide sequence file
#    uncl_seqfile = sys.argv[2]
#    # code file
#    codefile = sys.argv[3]
#    
#    filepath,tempfilename = os.path.split(uncl_seqfile)
#    outfile,extension = os.path.splitext(tempfilename)
#    outfile += ".fas"
#
#    prot_dict = makedict(prot_seqfile)
#    uncl_dict = makedict(uncl_seqfile)
#    code_dict = buildcode(codefile)
#    uncl_dict = splitseq(uncl_dict)
#    resu_dict = mapseq(prot_dict, uncl_dict, code_dict)
#    write2file(resu_dict, outfile)
#
#if __name__ == "__main__":
#    main()
