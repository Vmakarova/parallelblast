#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
WRONG_MAP = {
"not_exist_such_taxon_in_my_res":0,"wrong_match_soft":0,"wrong_match":0,"wrong_pident":0,"wrong_group":0, \
"is_shorter":0,"number_wrong_taxon":0, "no_taxon_in_right_res":0, "empty": 0}
dist_map = {a: 0 for a in range(1000)}
f_not_exist_such_taxon_in_my_res= open ("not_exist_such_taxon_in_my_res", 'w')
f_wrong_match_soft = open ("wrong_match_soft", 'w')
f_wrong_pident = open ("wrong_pident", 'w')
f_wrong_group = open ("wrong_group", 'w') 
f_no_taxon_in_right_res = open ("no_taxon_in_right_res", 'w')

parser = OptionParser(version="%prog 1.0", usage = "\n\n    %prog --input_file QUERY_FASTA --database_file \
FASTA_DATABASE --groups_omcl_file OrthoMCL_Groups [more_options]")



"""formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=dedent('''\
Parallel blastp + OrthoMCL 
Author: Makarova Valentina, makarovavs704@gmail.com, 2015 
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
                  help="Name of output")
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

Right_Results = {}
r_res_file = open("Octopus_bimaculoides.og.tsv", 'r')
for line in r_res_file:
  qid,gr,sid,evalue_mant,evalue_exp,pident,pmatch = line.split()
  Right_Results[qid] = [sid,evalue_mant,evalue_exp,pident,pmatch, gr] 


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Functions


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


#  flip orientation if nec.
def getStartEnd(h): 
    hspStart = h[0]
    hspEnd = h[1]
    if (hspStart > hspEnd):
      hspEnd = h[0]
      hspStart = h[1]
    return hspStart,hspEnd


def computeNonOverlappingMatchLength(subject):
    hsps = sorted(subject["hspspans"], key=lambda x: x[0])
    first = hsps[0]
    hsps = hsps[1::]
    if len(first) == 0:
      return 0 
    start, end = getStartEnd(first)
    match_len = 0
    for h in hsps:
      hspStart, hspEnd = getStartEnd(h)
      if hspEnd <= end: ##does not extend
        continue
      if hspStart <= end:  ##overlaps
        end = hspEnd #extend end ... already dealt with if new end is less
      else:   ##there is a gap in between ..
        match_len += end - start + 1
        start = hspStart
        end = hspEnd
    match_len += end - start + 1 # deal with the last one 
    return match_len


def formatEvalue(evalue):
    if evalue[0] == 'e':
      evalue = '1' + evalue
    evalue = "%.3e" % float(evalue)
    evalue_mant,evalue_exp = evalue.split("e")
    evalue_mant = "%.2f" % float(evalue_mant)
    if float(evalue_mant) == int(float(evalue_mant)):
      evalue_mant = str(int(float(evalue_mant)))
    if (evalue_exp[0] == '+'):
      evalue_exp = evalue_exp[1::]
    if evalue_exp == "00":
      evalue_exp = '0'
    return evalue_mant, evalue_exp


def printPreviousSubject(subject, file_ind):
    if  subject["find_right"]:
      return
    nonOverlapMatchLen = computeNonOverlappingMatchLength(subject)
    percentIdent =int(round((float(subject["totalIdentities"])/ float(subject["totalLength"]) * 10 + 0.5)/10))
    shorterLength = subject["queryLength"] if subject["queryShorter"] else subject["subjectLength"]
    percentMatch = int(round((float(nonOverlapMatchLen) / float(shorterLength)* 1000 + 0.5) / 10))
    if percentMatch > 50:
      ans= str(subject["queryId"]+"\t"+subject["omcl_group"].rjust(12)+" "+subject["subjectId"].ljust(26) \
           + subject["evalueMant"].rjust(5)+" "+subject["evalueExp"].rjust(5)+" "\
           + str(percentIdent).rjust(5)+" "+ str(percentMatch).rjust(5) + "\n") 
      output_file.write(ans)
      #print "PERCENT MAtch ", percentMatch
      if subject["queryId"] in Right_Results:
        right_subject = Right_Results[subject["queryId"]]
        if subject["subjectId"] == right_subject[0]:
          subject["find_right"] = True
          if (subject["wrong_dist"]) != 0:
            WRONG_MAP["number_wrong_taxon"] +=1 
            if subject["queryShorter"]:
              WRONG_MAP["is_shorter"] +=1
          dist_map[subject["wrong_dist"]] += 1
          right_match = int(right_subject[4])
          if right_match - percentMatch > 1 or right_match < -1:
            WRONG_MAP["wrong_match_soft"] += 1
            f_wrong_match_soft.write(subject["queryId"] + " "+file_ind + " " + str(subject["queryShorter"]) + "\n")
          elif right_match != percentMatch:
            WRONG_MAP["wrong_match"] += 1
          if right_subject[5] != subject["omcl_group"]:
            WRONG_MAP["wrong_group"] += 1
            f_wrong_group(subject["queryId"] + " "+file_ind + " " + str(subject["queryShorter"]) + "\n")
          if right_subject[3] != percentIdent:
            WRONG_MAP["wrong_pident"] +=1
            f_wrong_pident.write(subject["queryId"] + " "+file_ind + " " + str(subject["queryShorter"]) + "\n")
      else:
        subject["find_right"] = True
        WRONG_MAP["no_taxon_in_right_res"] +=1
        WRONG_MAP["number_wrong_taxon"] +=1 
        f_no_taxon_in_right_res.write(subject["queryId"] + " " +file_ind + " " + str(subject["queryShorter"]) + "\n")
        if subject["queryShorter"]:
            WRONG_MAP["is_shorter"] +=1
      subject["wrong_dist"] += 1



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
   #groups_OrthoMCL_file= open(options.omcl_file, 'r')
   groups_OrthoMCL = {}
   for line in open(options.omcl_file, 'r'):
    group_proteins = line.split()
    group_name = group_proteins.pop(0)[:-1]
    for protein in group_proteins:
      temp = protein.split('|')
      if temp[0] not in groups_OrthoMCL:
        groups_OrthoMCL[temp[0]] = {}
      groups_OrthoMCL[ temp[0] ][temp[1]] = group_name
   while not os.path.exists(os.path.join(temp_directory, "end")):
      time.sleep(3)
      print "wait blast"
   f = open(os.path.join(temp_directory, "end"), 'r')
   pr_q = int(f.read())
   f.close()
   os.remove(os.path.join(temp_directory, "end"))
   completed_file_quantity = 0
   while (completed_file_quantity < pr_q): 
     while not file_queue.empty():
        #try:
        file_ind = file_queue.get()
        #learn blast file
        completed_file_quantity +=1
        f =  open(file_ind, 'r')
        if os.path.getsize(file_ind) < 10:
          print completed_file_quantity, "/", pr_q
          WRONG_MAP["empty"] +=1
          continue
        prevSubjectId = 'blah'
        prevQueryId = 'blah'
        subject  = {} # hash to hold subject info
        queryShorter = False
        subject["find_right"] = False
        subject["wrong_dist"] = 0
        for line in f:
          if subject["find_right"]:
            break
          if len(line) == 0:
            continue
          #qseqid sseqid pident length mismatch  qstart qend sstart send evalue bitscore qlen slen
          queryId, subjectId, percentIdentity, length, mismatches, \
          queryStart, queryEnd, subjectStart, subjectEnd, evalue, bist, qlen, slen = line.split()
          if queryId != prevQueryId or subjectId  != prevSubjectId:
            # print previous subject
            if len(subject) > 2:
              printPreviousSubject(subject, file_ind)
            # initialize new one from first HSP
            prevSubjectId = subjectId
            prevQueryId = queryId
            subject["hspspans"] = []
            subject["totalIdentities"] = 0
            subject["totalLength"] = 0
            subject["queryId"] = queryId
            subject["subjectId"] = subjectId
            tmp = subject["subjectId"].split("|")           
            subject["omcl_group"] = "NO_GROUP"
            if (tmp[0] in groups_OrthoMCL) and (tmp[1] in groups_OrthoMCL[tmp[0]]):
              subject["omcl_group"]  = groups_OrthoMCL[tmp[0]][tmp[1]]
            subject["queryLength"] = int(qlen)
            subject["subjectLength"] = int(slen)
            subject["queryShorter"] = subject["queryLength"] < subject["subjectLength"]
            subject["evalueMant"], subject["evalueExp"] = formatEvalue(evalue) # from first hsp
          # get additional info from subsequent HSPs
          hspspan = [int(subjectStart), int(subjectEnd)]
          if subject["queryShorter"]:
            hspspan = [int(queryStart), int(queryEnd)]
          subject["hspspans"].append(hspspan)
          subject["totalIdentities"] += float(percentIdentity) * float(length)
          subject["totalLength"] += int(length)
        #if first:
        if len(subject) > 2:
          printPreviousSubject(subject, file_ind)
        if not subject["find_right"]:
         if subject["queryId"] in Right_Results:
          WRONG_MAP["not_exist_such_taxon_in_my_res"] += 1
          WRONG_MAP["number_wrong_taxon"] +=1 
          f_not_exist_such_taxon_in_my_res.write(subject["queryId"] + " "+ file_ind + " " + str(subject["queryShorter"]) + "\n")
          if subject["queryShorter"]:
            WRONG_MAP["is_shorter"] +=1
         else:
          dist_map[0] += 1
        print completed_file_quantity, "/", pr_q
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
      #except Exception as e:
      #  print e
      #  continue
   #os.rmdir(temp_directory)
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
    current_file_name = 0
    #current_size = 0
    for line in open(original_fasta):
        #
        # start a new file for each accession line
        #
        if line[0] == '>':
            current_file_name +=1
            file_name = "%d.segment" % current_file_name
            file_path = os.path.join(temp_directory, file_name)
            current_file = open(file_path, "w")
            #current_size = 0
        if current_file_name != "":
           current_file.write(line) 
    end_file = open(os.path.join(temp_directory, "end_tmp"), 'w')
    end_file.write(str(current_file_name))
    end_file.close()
    os.rename(os.path.join(temp_directory, "end_tmp"), os.path.join(temp_directory, "end"))


  
@transform(splitFasta, suffix(".segment"), [".blastResult", ".blastSuccess"], file_queue)
def runBlast(seqFile,  output_files, file_queue):
    #
    blastResultFile, flag_file = output_files
    cmd_str = blastp + " -db %s -query %s -out %s -evalue 1e-5 -max_target_seqs 10000 -outfmt \"6 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore qlen slen\""
    run_cmd(cmd_str % (database_file, seqFile, blastResultFile))
    file_queue.put(blastResultFile)
    open(flag_file, "w")
    #time.sleep(5)
    f =  open(blastResultFile, 'r')
    file_len = 0
    for line in f:
      if not len(line):
         continue
      file_len += 1


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
      result_daemon = Thread(target=orthomcl_daemon)
      result_daemon.setDaemon(True)
      result_daemon.start() 
      pipeline_run([runBlast],  multiprocess = options.jobs,
                          logger = logger, verbose=options.verbose)

      
      if (result_daemon.isAlive()):
          result_daemon.join() 
      print "wrong_match_soft ", WRONG_MAP["wrong_match_soft"], "wrong_match ", WRONG_MAP["wrong_match"],\
            "wrong_group ", WRONG_MAP["wrong_group"], "not_exist_such_taxon ", WRONG_MAP["not_exist_such_taxon_in_my_res"], \
            "no_taxon_in_right_res", WRONG_MAP["no_taxon_in_right_res"],"wrong_taxon_number", \
            WRONG_MAP["number_wrong_taxon"], "is_shorter", WRONG_MAP["is_shorter"], "empty", WRONG_MAP["empty"], "wrong_pident", WRONG_MAP["wrong_pident"]

      string_map = ""
      for i in range(1000):
        if  dist_map[i] != 0:
          string_map +=  str(i) +"\t" + str(dist_map[i]) + "\n"
      print string_map
