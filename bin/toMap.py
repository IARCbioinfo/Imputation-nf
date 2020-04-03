#!/usr/bin/env python3

file = open('/data/gep/MR_Signatures/work/Boris/protocol_min/data/HGDP_Map.txt', 'r')
file1 = open('/data/gep/MR_Signatures/work/Boris/protocol_min/data/HGDP.map', 'w')

for lines in file.readlines():
    line = lines.strip().split("\t")
    id = line[0]
    chr = line[1]
    basepaire = line[2]

    file1.write(str(chr) + '\t' + str(id) + '\t0\t' + str(basepaire) + '\n')

file.close()
file1.close()
