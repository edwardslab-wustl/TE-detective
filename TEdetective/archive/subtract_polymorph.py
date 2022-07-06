#!/home/mksingh/anaconda3/bin/python

import sys
import os
from io_functions import eprint
#inp_isz = 362
#
inp_ref_file = sys.argv[1]
inp_pred_file = sys.argv[2]
try:
    inp_isz = int( sys.argv[3] )
except IndexError:
    inp_isz = 50

eprint("ref_file: " + inp_ref_file)
eprint("pred_file: " + inp_pred_file)
eprint("isz: " + str(inp_isz))

#
with open(inp_ref_file, 'r') as ref_file:
    ref_file_lines = ref_file.readlines()
ref_file.close()
#
with open(inp_pred_file, 'r') as pred_file:
    pred_file_lines = pred_file.readlines()
pred_file.close()
#
dict_ref = {}
ref_list = []
line_count = 0
header = ''
for line in ref_file_lines:
    line_count += 1
    if line_count > 1:
        items = line.strip().split()
        #chrom = items[0]
        #ip = int(items[1])
        chrom = items[1]
        ip = int(items[2])
        ref_list.append((chrom, ip))
        try:
            dict_ref[chrom].append(ip)
        except KeyError:
            dict_ref[chrom] = []
            dict_ref[chrom].append(ip)
    else:
        header = line.strip()
     
#
tp = []
fp = []
ref_list_pred = []
line_count = 0
for line in pred_file_lines:
    line_count +=1
    if line_count > 1:
        items = line.strip().split()
        #chrom = items[0]
        #ip = int(items[1])
        chrom = items[1]
        ip = int(items[2])
        tp_flg = 0
        try:
            for i in dict_ref[chrom]:
                if ip in range(i-inp_isz, i+inp_isz+1):
                    tp.append((chrom, ip, abs(i-ip) ))
                    ref_list_pred.append((chrom, i))
                    tp_flg = 1
                    break
        except KeyError:
            pass
        #
        if tp_flg == 0:
            fp.append((chrom, ip))

#
# Write output
for item in fp:
    sys.stdout.write(str(item[0])+'\t'+str(item[1])+'\n')
#
#
