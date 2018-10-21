file = "data/test_6Degree10AS"
outputFile = file + "_undirected"

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
		if (dst, src) not in edges:
		    fout.write(dst + '\t' + src + '\t' + weight + '\t' + bw + '\t' + dstAS + '\t' + srcAS + '\n')
		    edges.add((dst, src))
