#load dependencies or install 
library(tidyverse)
library(reshape2)


#print working directory for future references
setwd("/Users/dunivint/Documents/GitHubRepos/RefSoil_plasmids/ANALYSIS_antibiotic_resistance/data/refseq/")
wd <- print(getwd())

#read in refseq plasmid quality data
#ultimately to remove incomplete plasmids
#get data for full refseq
refseq.quality <- read_delim("organism_identifier_refseq_tab.txt",
                             delim = "$", col_names = FALSE)


##EXTRACT COMPLETE BACTERIA/ARCHAEA PLASMIDS ONLY
#tidy refseq quality information
refseq.quality.tidy <- refseq.quality %>%
  separate(X1, into = c("Assembly", "Type"), sep = "genomic.gbff") %>%
  subset(Assembly !="--") %>%
  mutate(Type = gsub("-", "", Type),
         Type = gsub(":", "", Type)) %>%
  subset(Type !="ACCESSION") %>%
  mutate(Definition = ifelse(Type !="LOCUS", X2, NA),
         Definition = ifelse(Definition == "DEFINITION", NA, Definition),
         Definition = lead(Definition, n =1)) %>%
  subset(Type !="DEFINITION") %>%
  mutate(Addition = ifelse(Type == "LOCUS", " ", Type),
         Addition = lead(Addition, n = 1),
         Definition = paste(Definition, Addition, sep = " ")) %>%
  drop_na(X2) %>%
  select(-Type) %>%
  mutate(X2 = str_squish(X2)) %>%
  separate(col = X2, into = c("Accession", "size"), sep = " ") %>%
  mutate(size = as.numeric(size))


#read in refsoil data to remove
refsoil <- read_delim("../../output/refsoil_metadata_long.csv", delim = ",")

#subset refseq based on refsoil
refseq <- refseq.quality.tidy %>%
  subset(!Accession %in% refsoil$NCBI.ID) 

#extract plasmids
plasmids <- refseq %>%
  filter(grepl("lasmid",Definition))

chromosomes <- refseq %>%
  filter(!grepl("lasmid",Definition)) %>%
  filter(!grepl("extrachromosom",Definition))

#export accession numbers
write.table(select(plasmids, Accession), file = "../../output/refseq_plasmid.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(select(chromosomes, Accession), file = "../../output/refseq_chromosomes.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

#export plasmid sizes
write.table(select(plasmids, Assembly, Accession, size), file = "../../data/refseq/plasmid_size_refseq_tab.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

write.table(select(chromosomes, Assembly, Accession, size), file = "../../data/refseq/chromosome_size_refseq_tab.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
