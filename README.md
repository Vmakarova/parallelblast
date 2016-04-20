# parallelblast
This is a python program for parallel run of ncbi-blast and assign your proteins to OrthoMCL Groups 

Parallel blastp + OrthoMCL 
Author: Makarova Valentina, makarovavs07@gmail.com, 2015
Use: http://www.ruffus.org.uk/examples/bioinformatics/
     blastp Protein-Protein BLAST 2.2.30+  http://blast.ncbi.nlm.nih.gov/
     algorithm from OrthoMCL's service Assign your proteins to OrthoMCL Groups 
     http://www.orthomcl.org/orthomcl/proteomeUpload.do
In: set of proteins in fasta-file
Out: TSV-table, orthologGroups


##Build and run

For run this use Python2.7. Use --help for getting more info

Usage: prog --input_file QUERY_FASTA --database_file FASTA_DATABASE --omcl_file OrthoMCL_Groups [more_options]")
Example of running: C:\Python27\python.exe parallelblast.py -i GNAT.fasta  -d myblastdb -o groups_OrthoMCL-5.txt -j4 -b blast-2.2.31+\bin\blastp.exe
One else example: parallelblast -i octotest.fas -d omcl5 -j 4

##work data
parallelblastwithcheck.py - переписанная версия c относительно правильным match и подсчетом расстояния между ответом моей программки и оригинальной tsv таблицей(имя таблицы нужно вручную менять в коде, если используете сравнение с другим организмом)
wrong_dist_results - результат работы программы для 926 белков с расстоянием больше 1. 
