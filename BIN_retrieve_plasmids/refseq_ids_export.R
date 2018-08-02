#load dependencies or install 
library(tidyverse)
library(reshape2)
library(taxize)

#print working directory for future references
setwd("/Users/dunivint/Documents/GitHubRepos/RefSoil_plasmids/ANALYSIS_antibiotic_resistance/data/refseq/genomes/")
wd <- print(getwd())

#read in refseq plasmid quality data
#ultimately to remove incomplete plasmids
#get data for full refseq
refseq.quality <- read_delim("organism_identifier_refseq_tab.txt",
                             delim = "$", col_names = FALSE)


##EXTRACT COMPLETE BACTERIA/ARCHAEA PLASMIDS ONLY
#tidy refseq quality information
refseq.quality.tidy <- refseq.quality %>%
  mutate(X1 = gsub(".*.genomic.gbff:", "", X1),
         X1 = gsub(".*.genomic.gbff-", "", X1),
         X1 = trimws(X1)) %>%
  subset(X1 !="--") %>%
  subset(X1 !="ACCESSION") %>%
  mutate(Addition = ifelse(X1 !="LOCUS", X1, NA),
         Addition = ifelse(Addition == "DEFINITION", NA, Addition),
         Addition = lead(Addition, n = 1)) %>%
  subset(X1 %in% c("LOCUS", "DEFINITION")) %>%
  mutate(Label = paste(X2, ifelse(is.na(Addition), "", Addition), sep = " "),
         Label = str_squish(Label)) %>%
  select(X1, Label)
  

#subset and join definition and accession information
refseq.quality.def <- refseq.quality.tidy %>%
  subset(X1 == "DEFINITION")
refseq.quality.accno <- refseq.quality.tidy %>%
  subset(X1 == "LOCUS") %>%
  separate(Label, into = c("Label", "size"), sep = " ") %>%
  mutate(size = as.numeric(size))
refsoi.quality.complete <- cbind(refseq.quality.def, refseq.quality.accno)

#fix col names
colnames(refsoi.quality.complete) <- c("d", "Definition", "a", "Accession", "Size")

#read in refsoil data to remove
refsoil <- read_delim("../../../output/refsoil_metadata_long.csv", delim = ",")

#subset refseq based on refsoil
refseq <- refsoi.quality.complete %>%
  subset(!Accession %in% refsoil$NCBI.ID)

#extract plasmids
plasmids <- refseq %>%
  select(Definition, Accession, Size) %>%
  filter(grepl("lasmid",Definition))

chromosomes <- refseq %>%
  select(Definition, Accession) %>%
  filter(!grepl("lasmid",Definition))

#export accession numbers
write.table(plasmids$Accession, file = "../../../output/refseq_plasmid.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(chromosomes$Accession, file = "../../../output/refseq_chromosomes.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

#export plasmid sizes
write.table(select(plasmids, Accession, Size), file = "../../../data/refseq/plasmid_size_refseq_tab.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
