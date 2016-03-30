#!/usr/bin/python
# -*- coding: utf-8
import argparse
from textwrap import dedent
from uuid import uuid4 as uuid
from collections import deque
from ruffus import *
import subprocess
from threading import Thread
from optparse import OptionParser
import multiprocessing 
from math import frexp
import os, sys
exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
sys.path.insert(0,os.path.abspath(os.path.join(exe_path,"..", "..")))


import logging
logger = logging.getLogger("run_parallel_blast")

import time
test_time = time.time()

MATCH_PART = 50

parser = OptionParser(version="%prog 1.0", usage = "\n\n    %prog --input_file QUERY_FASTA --database_file FASTA_DATABASE --omcl_file OrthoMCL_Groups [more_options]")
"""	                  formatter_class=argparse.RawDescriptionHelpFormatter, 
                                 description=dedent('''\
Parallel blastp + OrthoMCL 
Author: Makarova Valentina, makarovavs07@gmail.com, 2015
Use: http://www.ruffus.org.uk/examples/bioinformatics/
In: set of proteins in fasta-file
Out: TSV-table, orthologGroups'''))"""

parser.add_option("-i", "--input_file", dest="input_file",
                  metavar="FILE",
                  type="string",
                  help="Name and path of query sequence file in FASTA format. ")
parser.add_option("-d", "--database_file", dest="database_file",
                  metavar="FILE",
                  type="string",
                  help="Name and path of FASTA database to search. ")
parser.add_option("-o", "--out_file", dest="out_file",
                  metavar="FILE",
                  type="string",
                  default="orthologGroups",
                  help="Name and path to out tsv-table of ortholog groups")
parser.add_option("-t", "--temp_directory", dest="temp_directory",
                  metavar="PATH",
                  type="string",
                  default="tmp",
                  help="Name and path of temporary directory where calculations "
                            "should take place. ")
parser.add_option('-g', "--groups_omcl_file", dest="omcl_file", 
                   metavar='PATH', 
                   default="groups_OrthoMCL-5.txt", 
                   type=str, nargs=1,
                   help='Name and path to OrthoMCL groups file .fasta')
parser.add_option('-b', "--blastp_exe", dest="blastp", 
                   metavar='PATH', 
                   default="blastp", 
                   type=str, 
                   help='Name and path to blastp.exe file')

#
#   general options: verbosity / logging
#
parser.add_option("-v", "--verbose", dest = "verbose",
                  action="count", default=0,
                  help="Print more detailed messages for each additional verbose level."
                       " E.g. run_parallel_blast --verbose --verbose --verbose ... (or -vvv)")

#
#   pipeline
#
parser.add_option("-j", "--jobs", dest="jobs",
                  default=1,
                  metavar="jobs",
                  type="int",
                  help="Specifies the number of jobs (operations) to run in parallel.")
parser.add_option("--flowchart", dest="flowchart",
                  metavar="FILE",
                  type="string",
                  help="Print flowchart of the pipeline to FILE. Flowchart format "
                       "depends on extension. Alternatives include ('.dot', '.jpg', "
                       "'*.svg', '*.png' etc). Formats other than '.dot' require "
                       "the dot program to be installed (http://www.graphviz.org/).")
parser.add_option("-n", "--just_print", dest="just_print",
                    action="store_true", default=False,
                    help="Only print a trace (description) of the pipeline. "
                         " The level of detail is set by --verbose.")
parser.add_option("-p", '--pack_size', dest="pack_size", metavar='INT', default=1, type=int, nargs=1,
                   help='Count of protein for one run blast')


(options, remaining_args) = parser.parse_args()


if not options.flowchart:
    if not options.database_file:
        parser.error("\n\n\tMissing parameter --database_file FILE\n\n")
    if not options.input_file:
        parser.error("\n\n\tMissing parameter --input_file FILE\n\n")
    if not options.omcl_file:
        parser.error("\n\n\tMissing parameter --omcl_file FILE\n\n")

if options.verbose:
    logger.setLevel(logging.DEBUG)
    stderrhandler = logging.StreamHandler(sys.stderr)
    stderrhandler.setFormatter(logging.Formatter("    %(message)s"))
    stderrhandler.setLevel(logging.DEBUG)
    logger.addHandler(stderrhandler)

original_fasta = options.input_file
database_file  = options.database_file
temp_directory = options.temp_directory
out_file    = options.out_file
pack_size = options.pack_size
output_file = open(options.out_file, 'w', 0)
blastp = options.blastp
file_queue = []
if __name__ == '__main__':
  m = multiprocessing.Manager()
  file_queue = m.Queue()



#groups_OrthoMCL_file= open(options.omcl_file, 'r')
groups_OrthoMCL = {};
for line in open(options.omcl_file, 'r'):
	group_proteins = line.split()
	group_name = group_proteins.pop(0)[:-1]
	for protein in group_proteins:
	  temp = protein.split('|')
	  if temp[0] not in groups_OrthoMCL:
	  	groups_OrthoMCL[temp[0]] = {}
	  groups_OrthoMCL[ temp[0] ][temp[1]] = group_name


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Functions


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


def run_cmd(cmd_str):
  
    #Throw exception if run command fails
    iteration = 0
    while (iteration < 10):
      try:
        process = subprocess.Popen(cmd_str, stdout = subprocess.PIPE,
                                stderr = subprocess.PIPE, shell = True)
        #stdout_str, stderr_str = process.communicate()
        stdout_str, stderr_str = process.communicate()
        if process.returncode != 0:
          print "Failed to run '%s'\nNon-zero exit status %s" %(cmd_str, process.returncode)
        break
      except: 
        continue


def orthomcl_daemon():
   print "daemon start"
   while not os.path.exists(os.path.join(temp_directory, "end")):
      time.sleep(2)
   f = open(os.path.join(temp_directory, "end"), 'r')
   pr_q = int(f.read())
   f.close()
   os.remove(os.path.join(temp_directory, "end"))
   completed_file_quantity = 0
   while (completed_file_quantity < pr_q or os.listdir(temp_directory) > 0):
     time.sleep(2)
     while not file_queue.empty():
        file_ind = file_queue.get()

        #learn blast file
        completed_file_quantity +=1
        print completed_file_quantity, "/", pr_q
        protein_hits = {}
        hit_priority = deque()
        hit_q = ''
        hit_q_len = 0
        f =  open(file_ind, 'r')
        for line in f:
          if not len(line):
            continue;
          hit = line.split()
          #h[0] qseqid 
          #h[1] sseqid
          #h[2] qstart
          #h[3] qend
          #h[4] sstart
          #h[5] send
          #h[6] qlen
          #h[7] slen
          #h[8] evalue
          #h[9] pident
          #h[10] qcovs
          if not hit_q_len:
            hit_q_len = float(hit[6])
            hit_q = hit[0]
          hit_priority.append(hit[1])
          if hit[1] not in protein_hits:
            protein_hits[hit[1]] = {"evalue": hit[8], "pident": hit[9]}
            if (hit_q_len < float(hit[7])):
              protein_hits[hit[1]]["is_shorter"] = False
              protein_hits[hit[1]]["match"] = hit[10]
            else:
              protein_hits[hit[1]]["is_shorter"] = True
              protein_hits[hit[1]]["len"] = hit[7]
              protein_hits[hit[1]]["sstart_send"] = [(hit[4], hit[5])] 
          elif protein_hits[hit[1]]["is_shorter"]:
              protein_hits[hit[1]]["sstart_send"].append((hit[4], hit[5]))


        #choose best hit
        while (len(hit_priority) != 0):
            hit_s = hit_priority.popleft()
            match = 0
            if not (protein_hits[hit_s]["is_shorter"]):
                if (protein_hits[hit_s]["match"] < MATCH_PART):
                    continue;
                match =  protein_hits[hit_s]["match"]
            else:
                protein_hits[hit_s]["sstart_send"].sort(key=(lambda x: x[0]))
                match_union = protein_hits[hit_s]["sstart_send"][1:]
                match_len = 0;
                st = int(protein_hits[hit_s]["sstart_send"][0][0])
                en = int(protein_hits[hit_s]["sstart_send"][0][1])

                for st_end_pair in match_union:
                    if (int(st_end_pair[0]) <= en):
                        en = int(st_end_pair[1])
                    else:
                        match_len += (en - st + 1)
                        st = int(st_end_pair[0])
                        en = int(st_end_pair[1])
                match_len += (en - st + 1)
                match = int(round(match_len * 100 / float(protein_hits[hit_s]["len"])))

                print "hit_q", hit_q

            #write best hit into file
            if match > MATCH_PART:
                tmp = hit_s.split("|") 
              
                gr = "NO_GROUP"
                if (tmp[0] in groups_OrthoMCL) and (tmp[1] in groups_OrthoMCL[tmp[0]]):
                    gr = groups_OrthoMCL[tmp[0]][tmp[1]]

                #evalue parsing
                ev = (protein_hits[hit_s]["evalue"]).split("e")
                if (ev[0] == "0.0") or (ev[0] == "0") or (int(ev[1]) < -180):
                   ev = [0, "-181"]

                output_file.write(str(hit_q + " "+ gr + " " + hit_s + " " + str(ev[0]) + " " + str(int(ev[1]))+ \
                                     " " + str(int(round(float(protein_hits[hit_s]["pident"]))))+ " " + str(match) +"\n"))
                break

        #close all one protein file and remove them       
        f.close()
        #os.remove(file_ind)
        tmp = file_ind.split(".")
        tmp[-1] = "segment"
        os.remove(".".join(tmp))
        tmp[-1] =  "blastSuccess"
        if os.path.exists(".".join(tmp)):
          for i in range(1,10):
            try:
              os.remove(".".join(tmp))
              break
            except:
               continue
   os.rmdir(temp_directory)
   print (time.time() - test_time)


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Pipeline tasks


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

@follows(mkdir(temp_directory))

@split(original_fasta, os.path.join(temp_directory, "*.segment"))
def splitFasta (seqFile, segments):
    #Split sequence file into
    #  as many fragments as appropriate
    #   depending on the size of original_fasta
    #
    #   Clean up any segment files from previous runs before creating new one
    #
    for i in segments:
        os.unlink(i)
    current_file_index = 0
    #current_size = 0
    for line in open(original_fasta):
        #
        # start a new file for each accession line
        #
        if line[0] == '>':
            #current_size +=1
            current_file_index += 1
            #if (current_size >= pack_size)
            file_name = "%d.segment" % current_file_index
            file_path = os.path.join(temp_directory, file_name)
            current_file = open(file_path, "w")
            #current_size = 0
        if current_file_index:
           current_file.write(line) 
    end_file = open(os.path.join(temp_directory, "end_tmp"), 'w')
    end_file.write(str(current_file_index))
    end_file.close()
    os.rename(os.path.join(temp_directory, "end_tmp"), os.path.join(temp_directory, "end"))

  
@transform(splitFasta, suffix(".segment"), [".blastResult", ".blastSuccess"], file_queue)
def runBlast(seqFile,  output_files, file_queue):
    #
    blastResultFile, flag_file = output_files
    cmd_str = blastp + " -db %s -query %s -out %s -evalue 1e-5 -outfmt \"6 qseqid sseqid qstart qend sstart send qlen slen evalue pident qcovs\""
    run_cmd(cmd_str % (database_file, seqFile, blastResultFile))
    file_queue.put(blastResultFile)
    open(flag_file, "w")
    time.sleep(5)
    f =  open(blastResultFile, 'r')
    file_len = 0
    for line in f:
      if not len(line):
         continue;
      file_len += 1
    print blastResultFile, " ", file_len


if __name__ == '__main__':
  if options.just_print:
      pipeline_printout(sys.stdout, [runBlast], verbose=options.verbose)
  elif options.flowchart:
      # use file extension for output format
      output_format = os.path.splitext(options.flowchart)[1][1:]
      pipeline_printout_graph (open(options.flowchart, "w"),
                               output_format,
                               [combineBlastResults],
                               no_key_legend = True)
  else:
      """result_daemon = Thread(target=orthomcl_daemon)
      result_daemon.setDaemon(True)
      result_daemon.start() """
      pipeline_run([runBlast],  multiprocess = options.jobs,
                          logger = logger, verbose=options.verbose)

      
      """if (result_daemon.isAlive()):
          result_daemon.join() 
      """ 
      print "without daemon"