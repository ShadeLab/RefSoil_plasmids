#load dependencies or install 
library(tidyverse)
library(reshape2)
library(taxize)

#print working directory for future references
setwd("/Users/dunivint/Documents/GitHubRepos/RefSoil_plasmids/ANALYSIS_antibiotic_resistance/data/refseq/")
wd <- print(getwd())

#read in refseq plasmid quality data
#ultimately to remove incomplete plasmids
#get data for full refseq
refseq.quality <- read_delim("plasmid_identifier_refseq_tab.txt",
                             delim = "$", col_names = FALSE)


##EXTRACT COMPLETE BACTERIA/ARCHAEA PLASMIDS ONLY
#tidy refseq quality information
refseq.quality.tidy <- refseq.quality %>%
  mutate(X1 = gsub("plasmid.*.genomic.gbff:", "", X1),
         X1 = gsub("plasmid.*.genomic.gbff-", "", X1)) %>%
  subset(X1 !="--") %>%
  subset(X1 !="VERSION")

#subset and join definition and accession information
refseq.quality.def <- refseq.quality.tidy %>%
  subset(X1 == "DEFINITION")
refseq.quality.accno <- refseq.quality.tidy %>%
  subset(X1 == "ACCESSION")
refsoi.quality.complete <- cbind(refseq.quality.def, refseq.quality.accno)

#fix col names
colnames(refsoi.quality.complete) <- c("d", "Definition", "a", "Accession")

#tidy refseq information
refsoi.q.c.tidy <- refsoi.quality.complete %>%
  select(-c(a, d)) %>%
  dplyr::filter(grepl("complete sequence", Definition)) %>%
  separate(Definition, into = c("Taxonomy1", "Taxonomy2"), sep = " ", remove = FALSE) %>%
  unite(col = "query", c("Taxonomy1", "Taxonomy2"), 
        sep = " ", remove = TRUE) %>%
  mutate(query = gsub("\\[", "", query),
         query = gsub("\\]", "", query),
         query = gsub("\\'", "", query))

#pull taxonomy... only want bacteria/archaea plasmids
#this step takes a long time
taxize_results <- tax_name(query = unique(refsoi.q.c.tidy$query), db = "ncbi", get = c("kingdom", "phylum"))
#need to manually examine NA entries

#remove all non-prokaryotic lineages
#then subset all refseq data based on this
ok <- c("Microscilla sp.", "Clostridium chauvoei,", "Cylindrospermum sp.", 'Alpha proteobacterium', 'Methylibium sp.', "Erwinina amylovora", "Gloeocapsa sp.", "Geobacillus genomosp.", "Marivivens sp.", "Cyanothece sp.", "Persicobacter sp.", "Halorientalis sp.", "Nostoc azollae", "Frankia symbiont", "Chelativorans sp.", "Rivularia sp.", "Sulfuriferula sp.", "Tenericutes bacterium")

#get list of actual plasmid accession numbers
refsoi.q.c.tidy.tax <- refsoi.q.c.tidy %>%
  left_join(taxize_results, by = "query") %>%
  dplyr::filter(!grepl("Fungi|Metazoa|Viridiplantae", kingdom)) %>%
  dplyr::filter(!grepl("Arthropoda|Streptophyta|Ascomycota|Basidiomycota", phylum)) %>%
  mutate(phylum2 = ifelse(is.na(phylum), ifelse(query %in% ok, "ok", "bad"), phylum)) %>%
  dplyr::filter(!grepl("bad", phylum2)) %>%
  separate(col = Accession, into = "Accession", sep = " ")

#write table
write.table(refsoi.q.c.tidy.tax$Accession, file = "../../output/refseq.plasmid.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
