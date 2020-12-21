import sys

dbfile = sys.argv[1]
namefile = sys.argv[2]
outfile = sys.argv[3]


dbdict = {}
with open(dbfile) as dbf:
	seqs = dbf.read().split('>')[1:]
for seq in seqs:
	lines = seq.split('\n')
	key = '>'+lines[0]+'\n'
	seq = ''.join(lines[1:])+'\n'
	dbdict[key] = seq

with open(namefile) as namef:
	#names = [name.strip() for name in namef]
	names = namef.readlines()

for i, name in enumerate(names):
    names[i] = name if name[0] == '>' else '>'+name

#namemap = {}
#for name in names:
#	for key in dbdict:
#		if name in key:
#			namemap[name] = key

with open(outfile, 'w') as outf:
	for key in names:
		try:
			#seq = dbdict[namemap[key]]
			#line = namemap[key] + seq
			seq = dbdict[key]
			line = key + seq
			outf.write(line)
		except KeyError:
			print(key)
