import random

file = "data/test_9Degree10AS"
prob = 2
outputFile = file + "_biodirection"

with open(file, 'r') as fin:
	inLines = fin.readlines()

with open(outputFile, 'w') as fout:
	fout.write(inLines[0])
	edges = set()
	for line in inLines[1:-1]:
		src, dst, weight, bw, srcAS, dstAS = line.rstrip().split("\t")
		if (src, dst) not in edges:
			fout.write(line)
			edges.add((src, dst))
		if (dst, src) not in edges and random.uniform(0,1) < prob:
		    fout.write(dst + '\t' + src + '\t' + weight + '\t' + bw + '\t' + dstAS + '\t' + srcAS + '\n')
		    edges.add((dst, src))
