#!/bin/bash -login

#### start of configuration

###### Adjust values for these parameters ####
#       SEQFILE, WORKDIR, REF_DIR, HMMSEARCH
#       ARG, E_VALUE, THREADS
##############################################

## THIS SECTION MUST BE MODIFIED FOR YOUR FILE SYSTEM. MUST BE ABSOLUTE PATH

#path to refsoil sequence files
SEQDIR=/mnt/research/ShadeLab/WorkingSpace/Dunivin/as_ncbi/refsoil

#path to working directory
WORKDIR=/mnt/research/ShadeLab/Dunivin/hmmsearch/beta_lactam

#path to hmms
REF_DIR=/mnt/research/ShadeLab/Dunivin/gene_targeted_assembly/RDPTools/Xander_assembler

#path to hmmsearch
HMMSEARCH=/opt/software/HMMER/3.1b2--GCC-4.8.3/bin/hmmsearch

## THIS SECTION NEED TO BE MODIFIED BASED ON EXPERIMENTAL QUESTIONS
#Which antibiotic resistance genes would you like to search?
ARG="ClassA ClassB ClassC AAC6-Ia ANT3 ANT6 ANT9 strA strB CAT CEP ermB ermC qnr adeB mexC mexE tolC sul2 tetA tetD tetM tetQ tetW tetX dfra1 dfra12 vanA vanC vanH vanT vanW vanX vanY vanZ intI repA"

#what evalue cutoff (for hmmsearch) would you like to use?
evalue=0.0000000001

## THIS SECTION NEEDS TO BE MODIFIED BASED ON RESOURCES
THREADS=8


