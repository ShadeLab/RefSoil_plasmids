#!/bin/bash

#downlaod all of refseq & plasmids
cd /mnt/scratch/dunivint/refseq_2018/
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/*.faa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/archaea/*.faa.gz

#unzip files
gunzip * 

#move grouped files to clean directory
cd /mnt/research/ShadeLab/Dunivin/hmmsearch/refseq
for i in bacteria archaea; do cat /mnt/scratch/dunivint/refseq_2018/*${i}* > protein_files/${i}.refseq.protein.fa; done

#linearize fasta files
for i in bacteria archaea; do awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ${i}.refseq.protein.fa > ${i}.refseq.protein.long.fa; done

#clean up fasta headers
for i in bacteria archaea; do awk '/^>/{print $1;next}{print}' ${i}.refseq.protein.long.fa > ${i}.protein.fa; done


#download ALL plasmid genbank files for refseq
cd /mnt/scratch/dunivint/refseq_2018/
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/*genomic.gbff.gz

#unzip files
gunzip *gbff.gz

#grab number of bp
grep '^LOCUS' *gbff > plasmid_size_refseq_tab.txt
grep -A1 -e '^DEFINITION' -e '^ACCESSION' -e '^ORGANISM' *gbff > plasmid_identifier2_refseq_tab.txt

