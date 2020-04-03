#!/usr/bin/env python3

import numpy as np

file = open('/data/gep/MR_Signatures/work/Boris/protocol_min/data/HGDP_FinalReport_Forward.txt', 'r')
list_list = [file.readline().split()]
for lines in file.readlines():
    line = lines.split()
    l = line[1:len(line)]
    list_list.append(l)

matrix = np.array(list_list)
np.savetxt('/data/gep/MR_Signatures/work/Boris/protocol_min/data/matrice.txt',matrix.T, fmt='%s')
file.close()

mat = open('/data/gep/MR_Signatures/work/Boris/protocol_min/data/matrice.txt', 'r')
file2 = open('/data/gep/MR_Signatures/work/Boris/protocol_min/data/HGDP.ped', 'w')

for lines in mat.readlines():
    line = lines.strip().split()
    id = line[0]

    infos = line[1:len(line)]
    infos = "".join(infos)
    infos = " ".join(list(infos))

    ped = str(id) + ' ' + str(id) + ' 0 0 0 0 ' + str(infos.replace('-','0')) + "\n"

    file2.write(ped)

mat.close()
file2.close()
