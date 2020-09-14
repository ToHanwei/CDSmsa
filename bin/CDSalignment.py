import os
import re
import sys


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

def buildcode(infile):
	"""build uncleotide code table"""
	with open(infile) as codef:
		codes = codef.readlines()
	codes = [line.strip().split(",") for line in codes]
	codes = {lines[0]:lines[1:] for lines in codes}
	return codes

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

def mapseq(pdict, udict, cdict):
	"""map MSA protein sequence to uncleotide"""
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
#	# protein sequence MSA file
#	prot_seqfile = sys.argv[1]
#	# uncleotide sequence file
#	uncl_seqfile = sys.argv[2]
#	# code file
#	codefile = sys.argv[3]
#	
#	filepath,tempfilename = os.path.split(uncl_seqfile)
#	outfile,extension = os.path.splitext(tempfilename)
#	outfile += ".fas"
#
#	prot_dict = makedict(prot_seqfile)
#	uncl_dict = makedict(uncl_seqfile)
#	code_dict = buildcode(codefile)
#	uncl_dict = splitseq(uncl_dict)
#	resu_dict = mapseq(prot_dict, uncl_dict, code_dict)
#	write2file(resu_dict, outfile)
#
#if __name__ == "__main__":
#	main()
