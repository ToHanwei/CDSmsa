import sys
import pickle

dbfile = sys.argv[1]
infile = sys.argv[2]
outfile = sys.argv[3]

with open(dbfile, 'rb') as pkl:
	mapdict = pickle.load(pkl)

with open(infile) as namef:
	names = [name.rstrip() for name in namef]

with open(outfile, 'w') as outf:
	for name in names:
		line = name + '\t' + mapdict[name]+'\n'
		outf.write(line)
