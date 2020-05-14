#!/usr/bin/env python3

file = open(str(sys.argv[1])+'HGDP_Map.txt', 'r')
file1 = open(str(sys.argv[1])+'HGDP.map', 'w')

for lines in file.readlines():
    line = lines.strip().split("\t")
    id = line[0]
    chr = line[1]
    basepaire = line[2]

    file1.write(str(chr) + '\t' + str(id) + '\t0\t' + str(basepaire) + '\n')

file.close()
file1.close()
