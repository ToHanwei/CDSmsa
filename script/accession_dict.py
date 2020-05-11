import sys
import pickle

dbfile = sys.argv[1]
outfile = sys.argv[2]

lineage_dict = {}
with open(dbfile) as dbf:
	for line in dbf:
		lines = line.rstrip().split('\t')
		key = lines[0]
		lineage = '\t'.join(lines[2:])
		lineage_dict[key] = lineage

with open(outfile, 'wb') as outf:
	pickle.dump(lineage_dict, outf)
