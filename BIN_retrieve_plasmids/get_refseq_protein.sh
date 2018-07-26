#!/bin/bash

#downlaod all of refseq & plasmids
cd /mnt/scratch/dunivint/refseq_2018/
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/*plasmid.faa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/*plasmid.faa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/archaea/*plasmid.faa.gz

#unzip files
gunzip * 

#move grouped files to clean directory
cd /mnt/research/ShadeLab/Dunivin/hmmsearch/refseq
for i in plasmid bacteria archaea; do cat /mnt/scratch/dunivint/refseq_2018/*${i}* > protein_files/${i}.refseq.protein.fa; done

#linearize fasta files
for i in plasmid bacteria archaea; do awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ${i}.refseq.protein.fa > ${i}.refseq.protein.long.fa; done

#clean up fasta headers
for i in plasmid bacteria archaea; do awk '/^>/{print $1;next}{print}' ${i}.refseq.protein.long.fa > ${i}.protein.fa; done
