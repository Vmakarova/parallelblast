#!/usr/bin/env python
import os, sys
from re import split

use_blast_res = False
if len(sys.argv) > 1:
  use_blast_res = True

f_right = open("Octopus_bimaculoides.og.tsv", 'r')
f_my = open("orthologGroups", "r")
f_results_all = open('checking_results', 'w')
f_wr_gr = open('wrong_group', 'w')
f_wr_match = open('wrong_match', 'w')
f_wr_gr_match = open('wrong_group_and_match', 'w')
f_wr_no_gr = open('wrong_no_group', 'w')
f_wr_choosen = open("wrong_choosen", 'w')

#original OrthoMCL output
protdb = {}
for line in f_right:
  if not len(line):
    continue;
  prot = line.split()
  protdb[ prot[0] ] = [prot[1], prot[2], prot[6], prot[3], prot[4], prot[5] ]


#blast result output data
blast_result = {}
if use_blast_res:
  for filename in os.listdir("tmp"):
    f = open ("tmp/" + filename, 'r')
    for line in f:
      if len(line) > 0 :
        blast_result[line.split()[0]] = filename
        break
    f.close()


str_res = ''
count = 0
count_no =0
count_wr_gr =0
count_wr_match = 0
count_wr_all = 0
count_blast_wr = 0
this_evalue_and_pident_exist = 0
all_coincides = 0 


for line in f_my:
  str_res = ''
  if not len(line):
    continue;
  prot = line.split()
  if not blast_result.has_key(prot[0]): 
    blast_result[prot[0] ] = ''
  if (protdb.has_key(prot[0])):
    right_prot = protdb[prot[0]]
    if (prot[1] != right_prot[0]):
  		count += 1
  		str_res = "wrong_group " + prot[0] + " " + blast_result[prot[0] ] + " " \
                               + prot[1]+" / " + right_prot[0]+" "+ prot[2]+' / ' \
  		                         + right_prot[1] + " " + prot[6] + " / " + right_prot[2] + '\n'
  		if (prot[6] != right_prot[2]):
  			count_wr_all += 1
  			f_wr_gr_match.write(str_res)
  		else:
  			count_wr_gr += 1
  			f_wr_gr.write(str_res)
    elif (prot[6] != right_prot[2]): # and (int(right_prot[2]) - int(prot[6])) != 1 and (int(prot[6]) - int(right_prot[2])) != 1):
        count += 1
        str_res = "wrong_match " + prot[0]  + " " + blast_result[prot[0] ] \
                           + " " + prot[6] + "/" + right_prot[2] + "\n"
        count_wr_match += 1
        #f_wr_match.write(str_res)
        if use_blast_res:
          filename = "tmp/"+blast_result[prot[0] ]
          f = open(filename, 'r')
          for line in f:
            #h[1] sseqid h[2] qstart h[3] qend h[4] sstart h[5] send
            #h[6] qlen h[7] slen h[8] evalue h[9] pident h[10] qcovs
            blast_result_hit = line.split()
            evalue = split(r'[e.]',blast_result_hit[8])
            if (right_prot[3] == evalue[0]  and right_prot[4] == evalue[1]) and \
                right_prot[5] == str(int(round(float(blast_result_hit[9])))):
                this_evalue_and_pident_exist += 1
                if (blast_result_hit[10] == right_prot[2]):
                  all_coincides +=1
                  f_wr_choosen.write(str_res)
                  f_wr_choosen.write(line)
  else:
    count_no += 1
    str_res = "NO_PROTEIN " + prot[0]  + " " + blast_result[prot[0] ] + "\n"
    f_wr_no_gr.write(str_res)
  if (len(str_res) != 0):
	f_results_all.write(str_res)


print "All mistake", count, "Group", count_wr_gr, "Match", count_wr_match, "Gr+Match", count_wr_all, "No_prot", count_no 
print "Wrong choosen", all_coincides, "Exist evalue and pident", this_evalue_and_pident_exist



