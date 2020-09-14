#!coding:utf-8

import os
import re
import sys
import tempfile


# Codon mapping table
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

def ReadSampleFasta(seqfile):
	"""
	Function:
		ReadSampleFasta
		read small FASTA format file.
		Note,large files consume a lot of memory.
	Parameter:
		seqfile, input file name, FASTA format file
	Return:
		return type -> List
		seq_list, [(name1, seq1), (name2, seq2), ...]
	"""

	seq_list = []
	with open(seqfile) as seqf:
		text = seqf.read().replace('\r', '')
		seqs = text.split('>')[1:]
		for seq in seqs:
			lines = seq.split('\n')
			name = lines[0]
			aads = ''.join(lines[1:])
			aads = aads.replace('*', '')
			seq_list.append((name, aads))
	return seq_list


def find_all(substring, string):
	"""
	Function:
		find_all
		find all indexes of a substring in a string.
		Overlapping is considered.
	Parameter:
		substring,
		string,
	Return:
		return type -> List
		substring starting position
	"""

	indexes = [m.start() for m in re.finditer(substring, string)]
	return indexes


def find_stop_codons(seq, SeqLengthLimit):
	"""
	Function:
		find_stop_codons
		Find stop codons in seq, no matter what reading frame is
	Parameter:
		seq, DNA sequence without header
	Return:
		return type -> List
		All 'stop codon' of seq(input DNA sequence)
	"""

	stop_tag = find_all('(?=(TAG))', seq)
	stop_tga = find_all('(?=(TGA))', seq)
	stop_taa = find_all('(?=(TAA))', seq)
	stops = stop_tag + stop_tga + stop_taa
	# Filter stop long enough
	stops = sorted(i for i in stops if i > SeqLengthLimit)
	return stops


def dna_translation(dna_seq):
	"""
	Function:
		dna_translation
		Translate DNA sequence to protein sequence
	Parameter:
		dna_seq, DNA sequence without header
	Return:
		return type -> String
		Translated protein sequence
	"""

	protein = ''
	i = 0
	len_dna = len(dna_seq)
	if len_dna % 3 != 0:
		logging.error("The length of CDS({0}) is not divisible by 3."
					  .format(len_dna))
		raise LengthError(len_dna)
	plen = len(dna_seq) / 3
	while i < plen:
		n = dna_seq[3 * i: 3 * (i + 1)]
		if 'N' in n:
			r = '*'
		else:
			r = CODON_TABLE[n]
		i += 1
		protein += r
	return protein


def find_cds(seqs, SeqLengthLimit):
	"""
	Function:
		find_cds
		try find ATG and STOP codons for each seq,
		put good cds info into fun = {},
		put others into outliers = {} for later processing
	Parameter:
		hmmout, funtion 'proc_nhmmer_out' output
		hmmout_seq, function 'extract_cds' outpuy
		SeqLengthLimit, artificially set OR's sequence length threshold
	Return:
	return type(fun) -> Dictory
	key   -> hit gene name
	value -> function OR CDS
	return type(outliers) -> Dictory
	key   -> hit gene name
	value -> outlier OR CDS
	"""

	hit = 1
	finals = []
	for seq in seqs:
		gname, raw_cds = seq
		raw_cds += 'TAG'
		len_cds = len(raw_cds)
		# The CDS length is too short
		if len_cds < SeqLengthLimit: continue
		atgs = find_all('(?=(ATG))', raw_cds)
		starts = sorted(i for i in atgs if i < len_cds - SeqLengthLimit)
		stops = find_stop_codons(raw_cds, SeqLengthLimit)
		iso = 1
		for iatg in starts:
			for istop in stops:
				klen = istop - iatg
				if (klen > SeqLengthLimit) and (klen % 3 == 0):
					fine_cds = raw_cds[iatg:istop]
					a_list = re.findall('...', fine_cds)
					# interrupting stop codon
					if 'TAG' in a_list or 'TGA' in a_list or 'TAA' in a_list:
						continue
					if 'N' not in fine_cds:
						hits = "hit" + str(hit)
						seqn = hits + '_' + gname + '_' \
							   + str(iatg) + '_' \
							   + str(istop + 3) + '_' \
							   + 'iso' + str(iso)
						final_cds = raw_cds[iatg:(istop + 3)]
						iso += 1
						finals.append((seqn, final_cds))
		hit += 1
	return finals


def unredundancy(hits):
	temp = tempfile.NamedTemporaryFile(mode='w+b', prefix='redun')
	fname = temp.name
	oname = fname + '_redun'
	lines = ['>'+h[0]+'\n'+h[1]+'\n' for h in hits]
	tempf = open(fname, 'w')
	tempf.writelines(lines)
	tempf.seek(0)
	os.system('cd-hit -i ' + fname + ' -c 1.0 -o ' + oname)
	redun_hits = ReadSampleFasta(oname)
	os.system('rm ' + oname +"*")
	tempf.close()
	return redun_hits


def write2file(seq_list, mapdict, profile, dnafile):
	"""write map sequence to file"""
	with open(dnafile, "w") as dnaf, open(profile, 'w') as prof:
		for name, seq in seq_list:
			if seq[-3:] in ['TAG', 'TGA', 'TAA']:
				seq = seq[:-3]
			dnaline = ">" + name + "\n" + seq + "\n"
			dnaf.write(dnaline)
			pro = dna_translation(seq)
			proline = ">"  + name + "\n" + pro + "\n"
			prof.write(proline)


infile = sys.argv[1]
prefix = infile.split('.')[0]
profile = prefix + "_pro.fasta"
dnafile = prefix + "_dna.fasta"
seqs = ReadSampleFasta(infile)
hits = find_cds(seqs, 870)
redun_hits = unredundancy(hits)
write2file(redun_hits, CODON_TABLE, profile, dnafile)

