######################################
#READ IN AND PREPARE DATA FOR ANALYSIS#
#######################################

#load dependencies or install 
library(tidyverse)
library(reshape2)
library(forcats)
library(scales)

#print working directory for future references
wd <- print(getwd())

#############################################
#PREPARE ARG HMM SEARCH RESULTS FOR ANALYSIS#
#############################################
#temporarily change working directory to data to bulk load files
setwd("data")
  
#read in abundance data
names=list.files(pattern="*.tbl.txt")
data <- do.call(rbind, lapply(names, function(X) {
  data.frame(id = basename(X), read.table(X))}))

#fix working directory
setwd(wd)

#remove unnecessary columns
data <- data %>%
  mutate(id = gsub(".0.0000000001.tbl.txt", "", id)) %>%
  separate(col = id, into = c("Gene", "Sample"), sep = "[.]", extra = "merge") %>%
  select(-c(V2, V5, V23))

#add column names
#known from http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf
colnames(data) <- c("Gene", "Sample", "t.name", "t.length", "q.name", "q.length", "e.value", "score1", "bias1", "#", "of", "c.evalue", "i.evalue", "score2", "bias2", "from.hmm", "to.hmm", "from.ali", "to.ali", "from.env", "to.env", "acc", "NCBI.ID")


#Calculate the length of the alignment
data.l <- data %>%
  mutate(length = to.ali - from.ali) %>%
  mutate(perc.ali = length / q.length)

#plot data quality distribution
(quality <- ggplot(data.l, aes(x = perc.ali, y = score1)) +
    geom_point(alpha = 0.5, shape = 1) +
    facet_wrap(~Gene, scales = "free_y", ncol = 5) +
    ylab("Score") +
    xlab("Percent alignment") +
    theme_bw(base_size = 8))

ggsave(quality, filename = paste(wd, "/figures/search.quality.png", sep = ""), width = 7, height = 7, units = "in")

#calculate score cutoff (1/3 of max score)
data.summary <- data.l %>%
  group_by(Gene) %>%
  summarise(max.score = max(score1)) %>%
  mutate(calc.min = max.score*0.4)

#remove hits with score less than 30% of maximum score
data.quality <- data.l %>%
  left_join(data.summary, by = "Gene") %>%
  mutate(NCBI.ID = as.character(NCBI.ID),
         t.name = as.character(t.name))

data.quality <- data.quality[which(data.quality$score1 > data.quality$calc.min),]

#remove rows of low quality (based on figure above)
#of hmm length (std)
data.quality.f <- data.quality[which(data.quality$perc.ali > 0.90 & data.quality$perc.ali < 1.1 & data.quality$acc > 0.70),]

#add more stringent filtering based on figure
data.quality.f <- data.quality.f[-which(data.quality.f$Gene == "vanA" & data.quality.f$score1 < 400),]

data.quality.f <- data.quality.f[-which(data.quality.f$Gene == "CEP" & data.quality.f$score1 < 500),]

data.quality.f <- data.quality.f[-which(data.quality.f$Gene == "adeB" & data.quality.f$score1 < 1000),]

data.quality.f <- data.quality.f[-which(data.quality.f$Gene == "vanC" & data.quality.f$score1 < 240),]

data.quality.f <- data.quality.f[-which(data.quality.f$Gene == "vanH" & data.quality.f$score1 < 400),]

data.quality.f <- data.quality.f[-which(data.quality.f$Gene == "vanX" & data.quality.f$score1 < 300),]

data.quality.f <- data.quality.f[-which(data.quality.f$Gene == "vanW" & data.quality.f$score1 < 250),]


#examine if any HMM hits apply to two genes 
#on the same organism
duplicates <- data.quality.f[duplicated(data.quality.f[,c(3,23)]),]
#88 duplicates detected


#fix any of the duplicates
#we will accept the one with a higher score
data.quality.f <- data.quality.f[order(data.quality.f$t.name, abs(data.quality.f$score1), decreasing = TRUE), ] 
data.quality.f <- data.quality.f[!duplicated(data.quality.f$t.name),]

(quality.f <- ggplot(data.quality.f, aes(x = perc.ali, y = score1)) +
    geom_point(alpha = 0.5) +
    facet_wrap(~Gene, scales = "free_y") +
    ylab("Score") +
    xlab("Percent alignment") +
    theme_bw(base_size = 10))

#save table to output
write.table(data.quality.f, paste(wd, "/output/ARG_summary.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)

#prep table for supplemental material
data.quality.f.pub <- data.quality.f %>%
  select(-c(t.length:q.length, bias1:acc, max.score,calc.min)) %>%
  mutate(Sample=gsub("bacteria.protein.fa","Bacterial genome",Sample),
         Sample=gsub("archaea.protein.fa","Archaeal genome", Sample),
         Sample=gsub("plasmid.protein.fa", "Plasmid", Sample),
         Gene = as.factor(Gene),
         Gene = fct_relevel(Gene)) %>%
  separate(t.name, into = c("p1", "p2", "p3"), sep = "_") %>%
  unite(col = `Protein accession`, p1, p2, sep = "_") %>%
  rename(Score = score1,
         `Alignment length` = length,
         `Percent alignment` = perc.ali) %>%
  select(-p3)

#save table to output
write.table(data.quality.f.pub, paste(wd, "/output/ARG_summary_clean.csv", sep = ""), sep = ",", quote = FALSE, row.names = FALSE)

#######################################
#PREPARE REFSOIL METADATA FOR ANALYSIS#
#######################################
#read in gene classifications
classification <- read_delim(file = paste(wd, "/data/gene_classifications.txt", sep = ""), delim = "\t")

#read in RefSoil data
#read in taxanomic information
ncbi <- read_delim(file = paste(wd, "/data/ismej2016168x6_plasmid.csv", sep = ""), delim = ",")

#separate out extra NCBI ID's and
#remove version number on NCBI ID (.##)
ncbi.tidy <- ncbi %>%
  mutate(Phylum = ifelse(Phylum == "Proteobacteria", Class, Phylum)) %>%
  group_by(Phylum) %>%
  mutate(length.phy = length(Phylum)) %>%
  separate(col = NCBI.ID, into = c("NCBI.ID", "NCBI.ID2", "NCBI.ID3"), sep = ",") %>%
  mutate(Contains.plasmid = !is.na(Plasmid.ID)) %>%
  separate(col = Plasmid.ID, into = c("plasmid1", "plasmid2","plasmid3","plasmid4","plasmid5","plasmid6","plasmid7","plasmid8","plasmid9","plasmid10","plasmid11","plasmid12","plasmid13","plasmid14"), sep = ",") %>%
  mutate(NCBI.ID = gsub("\\..*","",NCBI.ID),
         NCBI.ID2 = gsub("\\..*","",NCBI.ID2),
         NCBI.ID3 = gsub("\\..*","",NCBI.ID3)) %>%
  select(`RefSoil ID`, `Taxon ID`, Kingdom, Phylum, Organism, Contains.plasmid, NCBI.ID, NCBI.ID2, NCBI.ID3, plasmid1:plasmid14, length.phy) %>%
  melt(id = c("RefSoil ID", "Taxon ID","Kingdom", "Phylum", "Organism", "Contains.plasmid", "length.phy"), variable.name = "Source", value.name = "NCBI.ID") 


#remove all rows where value = NA
ncbi.tidy <- ncbi.tidy[!is.na(ncbi.tidy$NCBI.ID),]

#########################################
#MERGE DATASETS AND ANALYZE ARG LOCATION#
#########################################
#calculate n plasmids by taxonomy
tax.summary <- ncbi.tidy %>%
  mutate(Source = gsub("plasmid.", "plasmid", Source),
         Source = gsub("plasmid.", "plasmid", Source),
         Source = gsub("NCBI.ID.", "chromosome", Source),
         Source = gsub("NCBI.ID", "chromosome", Source)) %>%
  group_by(Phylum) %>%
  mutate(totalElements = length(`RefSoil ID`)) %>%
  ungroup() %>%
  group_by(Source, Phylum, totalElements) %>%
  summarise(nElements = length(Phylum)) %>%
  dcast(Phylum+totalElements ~ Source) %>%
  mutate(plasmid = ifelse(is.na(plasmid), 0, plasmid))
  
#join AsRG data with RefSoi IDs
data.tax <- ncbi.tidy %>%
  mutate(Source = gsub("plasmid.", "plasmid", Source),
         Source = gsub("plasmid.", "plasmid", Source),
         Source = gsub("NCBI.ID.", "chromosome", Source),
         Source = gsub("NCBI.ID", "chromosome", Source)) %>%
  group_by(`RefSoil ID`) %>%
  mutate(Elements = length(`RefSoil ID`)) %>%
  left_join(data.quality.f, by = "NCBI.ID") %>%
  left_join(classification, by = "Gene") 

#change NA gene to "None"
data.tax$Gene[is.na(data.tax$Gene)] <- "None"

#prepare data for plotting
data.tax.cast <- data.tax %>%
  select(`RefSoil ID`, Kingdom, Phylum, length.phy, Source, NCBI.ID, Contains.plasmid, Elements, Gene, Description) %>%
  dcast(`RefSoil ID`+Phylum+length.phy+Description+Gene~Source) %>%
  group_by(Description) %>%
  mutate(group = ifelse(plasmid == 0 & chromosome > 0, "Genome only", "Plasmid only"),
         group = ifelse(plasmid > 0 & chromosome > 0, "Plasmid and genome", group),
         group = as.factor(group), 
         Location = fct_relevel(group, "Plasmid and genome", "Plasmid only", "Genome only"),
         N = length(unique(Gene))) %>%
  ungroup() %>%
  mutate(Number.genes = chromosome+plasmid,
    Description = as.factor(Description), 
         Description = paste(Description, " (", N, ")", sep = ""),
         Description = fct_relevel(Description, "Trimethoprim (2)",  "Macrolide (3)", "Sulfonamide (1)","Plasmid (2)","Vancomycin (7)","Chloramphenicol (2)","Aminoglycoside (6)", "Tetracycline (6)", "Beta-Lactam (3)", "MultiDrugEfflux (4)")) %>%
  mutate(Phylum = as.factor(Phylum)) %>%
  group_by(Phylum) %>%
  subset(Gene !="None") %>%
  mutate(phy.count = sum(chromosome)+sum(plasmid)) %>%
  ungroup() %>%
  left_join(tax.summary, by = "Phylum") %>%
  mutate(Phylum = paste(Phylum, " (", length.phy, ")", sep = ""),
         Phylum = fct_reorder(Phylum, phy.count)) 
         
  
#plot genes on different locations
(gene_location <- ggplot(data.tax.cast, aes(x = Description, y = Number.genes, fill = Location)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Paired") +
  ylab("Number of genes") +
  coord_flip() +
  theme_bw(base_size = 10))

#save plot
ggsave(gene_location, filename = "figures/RefSoil_desc_location.eps", width = 5.5, height = 2.5, units = "in")

#make tidy table for plasmid-only sequences
plasmid_only <- data.tax.cast %>%
  subset(Location == "Plasmid only") %>%
  group_by(Gene) %>%
  summarise(Number = sum(Number.genes)) %>%
  mutate(Gene = as.factor(Gene)) %>%
  arrange(-Number)

#plot tidy plasmid-only table
library(ggpubr)
(plasmid_only_table <- ggtexttable(plasmid_only, 
                                   rows = NULL, 
                                   theme = ttheme(base_size = 8,
                                   tbody.style = tbody_style(hjust = 0, x=0.1, fill = "white", size = 8),
                                   colnames.style = colnames_style(color = "black", face = "bold", size = 9, linewidth = 1, linecolor = "white", fill = "white", hjust = 0, x = 0.1))))

#save tidy plasmid-only table
ggsave(plasmid_only_table, filename = "output/plasmid_only_table.png", width = 1.1, height = 5.2)


#plot genes on different locations
(phy_location <- ggplot(data.tax.cast, aes(x = Phylum, y = Number.genes, fill = Location)) +
    geom_bar(stat = "identity") +
    scale_fill_brewer(palette = "Paired") +
    ylab("Number of genes") +
    coord_flip() +
    theme_bw(base_size = 10))

#save plot
ggsave(phy_location, filename = "figures/RefSoil_phy_location.eps", width = 5.4, height = 3, units = "in")

#look at phylum-level plasmid &
#ARG-plasmid distributions
element.summary <- tax.summary %>%
  group_by(Phylum) %>%
  melt(id.vars = "Phylum", measure.vars = c("chromosome", "plasmid"), value.name = "N.element", variable.name = "Location.element")  %>%
  subset(Phylum %in% c( "Gammaproteobacteria", "Betaproteobacteria", "Alphaproteobacteria", "Firmicutes", "Deltaproteobacteria", "Actinobacteria", "Cyanobacteria", "Bacteroidetes", "Epsilonproteobacteria", "Acidobacteria", "Deinococcus-Thermus", "Verrucomicrobia", "Spirochaetes", "Planctomycetes", "Chlorobi", "Nitrospirae", "Fusobacteria", "Deferribacteres", "Chlamydiae", "Gemmatimonadetes", "Armatimonadetes")) %>%
  mutate(Phylum = fct_relevel(Phylum, "Armatimonadetes", "Gemmatimonadetes","Chlamydiae", "Deferribacteres","Fusobacteria","Nitrospirae","Chlorobi","Planctomycetes", "Spirochaetes","Verrucomicrobia","Deinococcus-Thermus","Acidobacteria","Epsilonproteobacteria", "Bacteroidetes","Cyanobacteria","Actinobacteria", "Deltaproteobacteria","Firmicutes","Alphaproteobacteria", "Betaproteobacteria",
"Gammaproteobacteria"), 
      Location.element = as.factor(Location.element),
Location.element = fct_relevel(Location.element, "plasmid", "chromosome")) 

#plot genes on different locations
(phy_location_element <- ggplot(element.summary, aes(x = Phylum, y = N.element, fill = Location.element)) +
    geom_bar(stat = "identity") +
    scale_fill_brewer(palette = "Paired") +
    ylab("Number of genetic elements") +
    coord_flip() +
    theme_bw(base_size = 10))

#save plot
ggsave(phy_location_element, filename = "figures/RefSoil_phy_location_element.eps", width = 5.4, height = 3, units = "in")

#######################################
#ANALYZE REFSOIL GENOME/PLASMID MAKEUP#
#######################################
#calculate refsoil organisms with 
#and without plasmids
refsoil.elements <- ncbi.tidy %>%
  group_by(`RefSoil ID`) %>%
  mutate(Source = gsub("plasmid.", "plasmid", Source),
         Source = gsub("plasmid.", "plasmid", Source),
         Source = gsub("NCBI.ID.", "chromosome", Source),
         Source = gsub("NCBI.ID", "chromosome", Source)) %>%
  dcast(`RefSoil ID`+Phylum~Source) %>%
  mutate(Location = ifelse(plasmid == 0 & chromosome > 0, "Genome", "Genome and plasmid"),
         Location = as.factor(Location), 
         Location = fct_relevel(Location, "Genome and plasmid", "Genome"),
         RefSoil = "RefSoil") %>%
  group_by(Location) %>%
  summarise(Proportion=length(Location)) 

#plot refsoil proportions
(refsoil.plasmid <- ggplot(refsoil.elements, aes(x = "", y = Proportion, fill = Proportion)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_gradient(low = "#1F78B4", high = "#B2DF8A") +
  ylab("Number of genes") +
  coord_polar("y",start = 0) +
    theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold")
    ) +
  theme(axis.text.x=element_blank()) +
  geom_text(aes(y = Proportion/2 + c(0, cumsum(Proportion)[-length(Proportion)]), label = percent(Proportion/922)), size=5))
ggsave(refsoil.plasmid, filename = "figures/refsoil.pie.eps")

#calculate plasmid number
#calculate refsoil organisms with 
#and without plasmids
refsoil.plasmids <- ncbi.tidy %>%
  group_by(`RefSoil ID`) %>%
  mutate(Source = gsub("plasmid.", "plasmid", Source),
         Source = gsub("plasmid.", "plasmid", Source),
         Source = gsub("NCBI.ID.", "chromosome", Source),
         Source = gsub("NCBI.ID", "chromosome", Source)) %>%
  dcast(`RefSoil ID`+Phylum~Source)

(plasmid.hist <- ggplot(subset(refsoil.plasmids, plasmid > 0), aes(x = plasmid)) +
  geom_histogram(stat = "count", color = "black") +
    scale_x_continuous(breaks = c(1,3,5,7,9,11,13,15)) +
    theme_classic(base_size = 10)+
    ylab("RefSoil organisms") +
    xlab("Plasmid number"))
ggsave(plasmid.hist, filename = "figures/plasmid.hist.eps", units = "in", height = 1.8, width = 2.5)

#calculate mean plasmid number
mean(refsoil.plasmids$plasmid)
median(refsoil.plasmids$plasmid)

(number.dist <- refsoil.plasmids %>%
  melt() %>%
ggplot(aes(x = variable, y = value)) +
  geom_boxplot() +
  ylab("Number of genetic elements") +
  xlab("Type of genetic element") +
  theme_bw())
ggsave(number.dist, filename = "figures/plasmid_number_dist.eps", units = "in", width = 3.5, height = 3)

#######################################
#ANALYZE REFSOIL GENOME/PLASMID MAKEUP#
#######################################
#temporarily change working directory to data to bulk load files
setwd("data")

#read in abundance data
names=list.files(pattern="*_size_refsoil_tab.txt")
size <- do.call(rbind, lapply(names, function(X) {
  data.frame(id = basename(X), read.table(X))}))

#fix working directory
setwd(wd)

#calculate bp in plasmids and genomes
size.groups <- size %>%
  mutate(id = gsub("archaea", "genome", id),
         id = gsub("bacteria", "genome", id)) %>%
  group_by(id) %>%
  summarise(total.Kbp = sum(as.numeric(V3)/1000),
            N.group = length(as.numeric(V3))) %>%
  mutate(percent = total.Kbp/sum(total.Kbp))

#plot refsoil bp proportions
(refsoil.bp.plasmid <- ggplot(size.groups, aes(x = "", y = total.Kbp, fill = total.Kbp)) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_gradient(low = "#1F78B4", high = "#B2DF8A") +
    ylab("Number of genes") +
    coord_polar("y",start = 0) +
    theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold")
    ) +
    theme(axis.text.x=element_blank()) +
    geom_text(aes(y = total.Kbp/2 + c(0, cumsum(total.Kbp)[-length(total.Kbp)]), label = percent(total.Kbp/sum(total.Kbp))), size=5))
ggsave(refsoil.bp.plasmid, filename = "figures/refsoil.bp.pie.eps")

#extract plasmids only 
size.p <- size %>%
  subset(id == "plasmid_size_refsoil_tab.txt") %>%
  mutate(NCBI.ID = as.character(V2))

#plot plasmid sizes
(size.plasmid <- ggplot(size.p, aes(x = V3/1000)) +
  geom_histogram(bins = 30) +
  scale_x_log10(breaks = c(1,5,10,35,100,1000)) +
  ylab("Number of plasmids") +
  xlab("Plasmid size (Kbps)") +
  theme_classic(base_size = 10))
ggsave(size.plasmid, filename = "figures/plasmid.size.hist.eps", units = "in", height = 1.8, width = 2.5)

###############################
#PLASMID DIFFER BTWN ORGANISMS#
###############################
duplicates <- ncbi[which(duplicated(ncbi$Organism)),]

#extract duplicated organisms
dup.df <- ncbi.tidy %>%
  subset(Organism %in% duplicates$Organism) %>%
  subset(Contains.plasmid !=FALSE) %>%
  mutate(NCBI.ID = gsub("\\..*", "", NCBI.ID)) %>%
  left_join(size.p, by = "NCBI.ID") %>%
  mutate(Source = ifelse(Contains.plasmid == FALSE, "plasmid1", Source)) %>%
  subset(Source !="NCBI.ID") %>%
  subset(Source !="NCBI.ID2") %>%
  mutate(Organism = as.factor(Organism))

#order by organism
dup.df$`RefSoil ID` <- factor(dup.df$`RefSoil ID`, 
                              dup.df$`RefSoil ID`[order(dup.df$Organism)])

(duplicate.plot <- ggplot(dup.df, aes(x = `RefSoil ID`, y = Source,
                                      size = V3/1000,
                                      color = Organism)) +
    geom_point() +
    scale_color_brewer(palette = "Dark2") +
    coord_flip() +
    theme_bw(base_size = 10) +
    ylab("Plasmid") +
    ylim(2,NA) +
    labs(size = "Plasmid size (kbp)") +
    theme(axis.text.x = element_blank()))

#############################
#PLASMID SIZE VS GENOME SIZE#
#############################
#prep taxonomy for size annotation
ncbi.tidy.size <- ncbi.tidy %>%
  mutate(NCBI.ID = gsub("\\..*", "", NCBI.ID),
         NCBI.ID = trimws(NCBI.ID))

#annotate element size by taxonomy
size.annotated <- size %>%
  rename(NCBI.ID = V2) %>%
  left_join(ncbi.tidy.size, by = "NCBI.ID") %>%
  mutate(id = gsub("plasmid_size_refsoil_tab.txt", "plasmid", id),
         id = gsub(".*_size_refsoil_tab.txt", "genome", id)) %>%
  group_by(Kingdom, Phylum, Organism, `RefSoil ID`, id) %>%
  summarise(mean_size = mean(V3),
            tot_size = sum(V3),
            n_size = length(V3)) %>%
  melt(measure.vars = c("mean_size", "tot_size", "n_size"), variable.name = "Measurement") %>%
  dcast(Kingdom+Phylum+Organism+`RefSoil ID`~id+Measurement, value.var = "value") %>%
  mutate(plasmid_n_size = ifelse(is.na(plasmid_n_size), 0, plasmid_n_size),
         plasmid_tot_size = ifelse(is.na(plasmid_tot_size), 0, plasmid_tot_size),
         plasmid_mean_size = ifelse(is.na(plasmid_mean_size), 0, plasmid_mean_size))

(plas_v_genome <- ggplot(size.annotated, aes(x = genome_tot_size/1000, y = plasmid_tot_size/1000)) +
    geom_point(alpha = 0.4, size = 2) +
    ylab("Total plasmid size (kbp)") +
    xlab("Total genome size (kbp)") +  
    scale_y_log10() +
    labs(fill = "Number of
plasmids") +
    theme_bw(base_size = 12) +
  theme(legend.position = c(0.85, 0.6), 
        legend.background = element_rect(color = "black", 
                                         fill = "white", 
                                         size = 0.05, 
                                         linetype = "solid"), 
        axis.title.x = element_blank(), 
        legend.title=element_text(size=9), 
        legend.text=element_text(size=8),
        axis.text.x = element_blank()))

ggsave(plas_v_genome, filename = "figures/plasmid_v_genome.png", units = "in", width = 3.5, height = 3, dpi = 300)


(plas.n_v_genome <- ggplot(size.annotated, aes(x = genome_tot_size/1000, y = plasmid_n_size, size = plasmid_mean_size/1000, alpha = plasmid_mean_size/1000)) +
    geom_jitter(height = 0.2, width = 0) +
    scale_fill_gradientn(colours = heat.colors(400)) +
    ylab("Number of plasmids") +
    xlab("Total genome size (kbp)") +
    labs(size = "Mean plasmid
size (kbp)", 
         alpha = "Mean plasmid
size (kbp)") +
    theme_bw(base_size = 12) +
    theme(legend.position = c(0.85, 0.6), 
          legend.background = element_rect(color = "black", 
                                           fill = "white",
                                           size = 0.05, 
                                           linetype = "solid"), 
          legend.title=element_text(size=9), 
          legend.text=element_text(size=8)))


ggsave(plas.n_v_genome, filename = "figures/plasmid.number_v_genome.png", units = "in", width = 3.5, height = 3, dpi = 300)

(density.genome <- size.annotated %>%
  mutate(plasmid.logical = ifelse(plasmid_n_size == 0, "None", ifelse(plasmid_n_size == 1, "One", "Multiple")),
         plasmid.logical = as.factor(plasmid.logical),
         plasmid.logical = fct_relevel(plasmid.logical, "None", "One", "Multiple")) %>%
 ggplot(aes(x = genome_tot_size/1000, fill = plasmid.logical)) +
  geom_density(alpha = 0.6) +
  labs(fill = "Number of
plasmids") +
  scale_fill_manual(values = c("#B2DF8A", "#1F78B4", "#2007ff")) +
  theme_void() +
  theme(legend.position = c(0.85, 0.6), legend.background = element_rect(color = "black", fill = "white", size = 0.05, linetype = "solid"), axis.title.x = element_blank(), legend.title=element_text(size=9), 
        legend.text=element_text(size=8)))

(density.plasmid <- size.annotated %>%
    mutate(plasmid.logical = ifelse(plasmid_n_size == 0, "None", ifelse(plasmid_n_size == 1, "One", "Multiple")),
           plasmid.logical = as.factor(plasmid.logical),
           plasmid.logical = fct_relevel(plasmid.logical, "None", "One", "Multiple")) %>%
    ggplot(aes(x = plasmid_tot_size/1000, fill = plasmid.logical)) +
    geom_density(alpha = 0.6) +
    labs(fill = "Number of
plasmids") +
    scale_x_log10() +
    scale_fill_manual(values = c("#1F78B4", "#2007ff")) +
    theme_void() +
    coord_flip() +
    theme(legend.position = "none",
          axis.title.x = element_blank()))

library(ggpubr)
(Fig_S3 <- ggarrange(density.genome,NULL, plas_v_genome, density.plasmid, plas.n_v_genome, 
          ncol = 2, nrow = 3,  align = "hv",
          widths = c(2, 0.9), heights = c(0.9, 1.5, 1.5)))
ggsave(Fig_S3, filename = "figures/figure_s3.png", height = 9, width = 6, units = "in", dpi = 300)

#set up df for correlation stats
size.stats <- size.annotated %>%
  select(genome_tot_size, plasmid_tot_size, genome_n_size, plasmid_n_size)

#test correlations
library(psych)
corr.stats <- corr.test(size.stats)
corr.stats.r <- corr.stats$r
corr.stats.p <- corr.stats$p
