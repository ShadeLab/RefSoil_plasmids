######################################
#READ IN AND PREPARE DATA FOR ANALYSIS#
#######################################

#load dependencies or install 
library(tidyverse)
library(reshape2)
library(forcats)
library(scales)
library(gridExtra)
library(ggpubr)
library(psych)

#print working directory for future references
wd <- print(getwd())

#temporarily change working directory to data to bulk load files
setwd("data")

#read in abundance data for refsoil
names=list.files(pattern="*.tbl.txt")
data <- do.call(rbind, lapply(names, function(X) {
  data.frame(id = basename(X), database = "RefSoil", read.table(X))}))

#get data for full refseq
setwd("refseq/")
refseq.names=list.files(pattern="*.tbl.txt")
refseq <- do.call(rbind, lapply(refseq.names, function(X) {
  data.frame(id = basename(X), database = "RefSeq", read.table(X))}))

#fix working directory
setwd(wd)

#remove residual refsoil accnos from refseq results
refseq <- refseq %>%
  subset(V24 !="NC_002679.1") %>%
  subset(V24 !="NC_008697.1") %>%
  subset(V24 !="NC_014000.1") %>%
  subset(V24 !="NZ_HG938354.1") %>%
  subset(V24 !="NZ_HG938356.1") 

#join datasets together
data <- rbind(data, refseq)
data <- data %>%
  mutate(id = gsub(".tbl.txt", "", id)) %>%
  separate(col = id, into = c("Gene", "Sample"), sep = "[.]", extra = "merge") %>%
  select(-c(V2, V5, V23))

#add column names
#known from http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf
colnames(data) <- c("group", "extra", "database","Query", "t.length", "Gene", "q.length", "e.value", "score1", "bias1", "#", "of", "c.evalue", "i.evalue", "score2", "bias2", "from.hmm", "to.hmm", "from.ali", "to.ali", "from.env", "to.env", "acc", "NCBI.ID")

#Calculate the length of the alignment
data.l <- data %>%
  mutate(length = to.ali - from.ali) %>%
  mutate(perc.ali = length / q.length)

#plot data quality distribution
(quality <- data.l %>%
    subset(database == "RefSoil") %>%
    ggplot(aes(x = perc.ali, y = score1)) +
    geom_point(alpha = 0.5, shape = 1) +
    facet_wrap(~Gene, scales = "free_y", ncol = 7) +
    ylab("Score") +
    xlab("Percent alignment") +
    theme_bw(base_size = 8))

ggsave(quality, filename = paste(wd, "/figures/figure_s1.png", sep = ""), width = 18, height = 18, units = "in")

#calculate score cutoff
data.summary <- data.l %>%
  group_by(Gene) %>%
  summarise(max.score = max(score1)) %>%
  mutate(calc.min = max.score*0.3)

#remove hits with score less than 30% of maximum score
data.quality <- data.l %>%
  left_join(data.summary, by = "Gene") %>%
  mutate(NCBI.ID = as.character(NCBI.ID),
         Gene = as.character(Gene))

data.quality <- data.quality[which(data.quality$score2 > data.quality$calc.min),]

#remove rows of low quality (based on figure above)
#of hmm length (std)
data.quality.f <- data.quality[which(data.quality$perc.ali > 0.90 & data.quality$perc.ali < 1.1 & data.quality$acc > 0.80),]

#examine if any HMM hits apply to two genes 
#on the same organism
duplicates <- data.quality.f[duplicated(data.quality.f[,c(4,24)]),]
#152025 duplicates detected

#fix any of the duplicates
#we will accept the one with a higher score
data.quality.f.d <- data.quality.f[order(data.quality.f$Gene, abs(data.quality.f$score2), decreasing = TRUE),] 

data.quality.final <- data.quality.f.d[!duplicated(data.quality.f.d[,c(4,24)]),]

#check that duplicates were removed
duplicates.f.check <- data.quality.final[duplicated(data.quality.final[,c(4,24)]),]

(quality.f <- data.quality.final %>%
    ggplot(aes(x = perc.ali, y = score1)) +
    geom_point(alpha = 0.5) +
    facet_wrap(~Gene, scales = "free_y") +
    ylab("Score") +
    xlab("Percent alignment") +
    theme_bw(base_size = 10))
ggsave(quality.f, filename = "lala.png", width = 18, height = 18)

#save table to output
write.table(subset(data.quality.final, database == "RefSeq"), paste(wd, "/output/ARG_summary.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)

#prep table for supplemental material
data.quality.f.pub <- data.quality.f %>%
  select(-c(extra, t.length,q.length, bias1:acc, max.score,calc.min)) %>%
  mutate(Sample = group,
         Sample=gsub("bacteria","Chromosome",Sample),
         Sample=gsub("archaea","Chromosome", Sample),
         Sample=gsub("chromosome","Chromosome", Sample),
         Sample=gsub("plasmid", "Plasmid", Sample),
         Gene = as.factor(Gene),
         Gene = fct_relevel(Gene)) %>%
  arrange(Gene) %>%
  # separate(t.name, into = c("p1", "p2", "p3"), sep = "_") %>%
  #unite(col = `Protein accession`, p1, p2, sep = "_") %>%
  rename(Score = score1,
         `Alignment length` = length,
         `Percent alignment` = perc.ali) 

#save table to output
write.table(subset(data.quality.f.pub, database == "RefSoil"), paste(wd, "/output/ARG_summary_clean.csv", sep = ""), sep = ",", quote = FALSE, row.names = FALSE)


#######################################
#PREPARE REFSOIL METADATA FOR ANALYSIS#
#######################################
#read in gene classifications
classification <- read_delim(file = "data/180102_resfams_metadata_updated_v1.2.2.csv", delim = ",")

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

#save ncbi table
write.table(ncbi.tidy, file = "output/refsoil_metadata_long.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

##################################
#ANALYZE ARG LOCATION for REFSOIL#
##################################
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
  left_join(subset(data.quality.f, database == "RefSoil"), by = "NCBI.ID") %>%
  mutate(Gene = gsub("efflux_EmrB", "emrB", Gene)) %>%
  left_join(classification, by = "Gene") 

#change NA gene to "None"
data.tax$Gene[is.na(data.tax$Gene)] <- "None"

#calculate n organisms per phylum
n.org.tax <- ncbi %>%
  mutate(Phylum = ifelse(Phylum == "Proteobacteria", Class, Phylum)) %>%
  group_by(Phylum) %>%
  summarise(n.org = length(Phylum))

#prepare data for plotting
data.tax.cast <- data.tax %>%
  select(`RefSoil ID`, Phylum, Classification, length.phy, Source, NCBI.ID, Contains.plasmid, Elements, Gene, Description) %>%
  dcast(`RefSoil ID`+Phylum+length.phy+Description+Classification+Gene~Source) %>%
  group_by(Classification) %>%
  mutate(group = ifelse(plasmid == 0 & chromosome > 0, "Chromosome only", "Plasmid only"),
         group = ifelse(plasmid > 0 & chromosome > 0, "Plasmid and chromosome", group),
         N = length(unique(Gene))) %>%
  ungroup() %>%
  mutate(group = as.factor(group),
         Location = fct_relevel(group,"Plasmid only", "Plasmid and chromosome")) %>%
  ungroup() %>%
  mutate(Number.genes = chromosome+plasmid,
         Classification = as.factor(Classification), 
         Classification = paste(Classification, " (", N, ")", sep = "")) %>%
  mutate(Phylum = as.factor(Phylum),
         Classification = fct_relevel(Classification,"ABC Transporter (7)","Aminoglycoside (24)","Gene Modulating Resistance (15)", "Other (7)", "MFS Transporter (6)", "Beta-Lactam (32)", "Other Efflux (6)", "RND Antibiotic Efflux (10)", "Gylcopeptide (7)", "Trimethoprim (2)", "Tetracycline (9)", "Chloramphenicol (4)", "rRNA Methyltransferase (7)", "Quinolone (2)", "Macrolide (1)")) %>%
  group_by(Phylum) %>%
  subset(Gene !="None") %>%
  mutate(phy.count = sum(chromosome)+sum(plasmid)) 

#plot genes on different locations
(gene_location <- data.tax.cast %>%
    ggplot(aes(x = Classification, y = Number.genes, fill = Location)) +
    geom_bar(stat = "identity", color = NA) +
    scale_fill_brewer(palette = "Paired") +
    ylab("") +
    labs(fill = "") +
    theme_bw(base_size = 10) +
    coord_flip())

(gene_location.prop <- ggplot(data.tax.cast, aes(x = Classification, y = Number.genes, fill = Location)) +
    geom_bar(stat = "identity", color = NA, position = "fill") +
    scale_fill_brewer(palette = "Paired") +
    labs(fill = "") +
    ylab("") +
    xlab("") +
    theme_bw(base_size = 10) +
    coord_flip() +
    theme(axis.text.y = element_blank()))


#make tidy table for plasmid-only sequences
plasmid_only <- data.tax.cast %>%
  subset(Location == "Plasmid only") %>%
  group_by(Gene) %>%
  summarise(Number = sum(Number.genes)) %>%
  mutate(Gene = as.factor(Gene)) %>%
  arrange(-Number)

#get plasmid only gene info
no.plas <- data.tax.cast %>% subset(group != "Plasmid only")
yes.plas <- data.tax.cast %>% subset(group == "Plasmid only")

yes.plas <- yes.plas[-which(yes.plas$Gene %in% no.plas$Gene),]

yes.plas.tidy <- yes.plas %>%
  group_by(Gene) %>%
  summarise(Number = sum(plasmid))

#plot tidy plasmid-only table
library(ggpubr)
(plasmid_only_table <- ggtexttable(yes.plas.tidy, 
                                   rows = NULL, 
                                   theme = ttheme(base_size = 8,
                                                  tbody.style = tbody_style(hjust = 0, x=0.1, fill = "white", size = 8),
                                                  colnames.style = colnames_style(color = "black", face = "bold", size = 9, linewidth = 1, linecolor = "white", fill = "white", hjust = 0, x = 0.1))))

#save tidy plasmid-only table
ggsave(plasmid_only_table, filename = "figures/figure_4_c.eps", width = 2, height = 4)

#set up data for phylum comparison
data.tax.cast.phy <- data.tax.cast %>%
  ungroup() %>%
  left_join(n.org.tax, by = "Phylum") %>%
  mutate(Phylum = as.factor(Phylum)) %>%
  group_by(Phylum, n.org) %>%
  mutate(gene.total = sum(Number.genes), 
         n.plasmid = sum(plasmid), 
         n.chrom = sum(chromosome),
         norm.gene.total = gene.total/n.org) %>%
  ungroup() %>%
  mutate(Phylum = paste(Phylum, " (", n.org, ")", sep = ""),
         Phylum = fct_reorder(as.factor(Phylum), -gene.total))

#plot genes on different locations
(figure_4c <- ggplot(data.tax.cast.phy, aes(x = Phylum, y = Number.genes, fill = Location)) +
    geom_bar(stat = "identity") +
    scale_fill_brewer(palette = "Paired") +
    ylab("Number of gene hits") +
    labs(fill = "ARG location") +
    coord_flip() +
    theme_bw(base_size = 10))

(figure_4d <- ggplot(data.tax.cast.phy, aes(x = Phylum, y = Number.genes, fill = Location)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_brewer(palette = "Paired") +
    xlab("") +
    ylab("Proportion") +
    labs(fill = "ARG location") +
    coord_flip() +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_blank()))

#save figure 4
#save plot
(figure_4 <- ggarrange(gene_location, gene_location.prop, figure_4c, figure_4d, labels = c("A","B","C", "D"), common.legend = TRUE, ncol = 2, nrow = 2, heights = c(0.6, 1.1), align = "v"))

ggsave(figure_4, filename = "figures/Figure_4.png", width = 8, height = 7.5, units = "in", dpi = 300)


#######################################
#COMPARE REFSOIL ARGS WITH REFSEQ ARGS#
#######################################
#RefSeq release 89
db.numbers <- read_delim("data/database_numbers.txt", delim = "\t", col_names = TRUE)

#tidy data for comparison
data.full <- data.quality.final
data.full.tidy <- data.full %>%
  left_join(classification, by = "Gene") %>%
  mutate(Sample = group,
         Sample = gsub("archaea", "Chromosome", Sample),
         Sample = gsub("bacteria", "Chromosome", Sample),
         Sample = gsub("chromosome", "Chromosome", Sample),
         Sample = gsub("plasmid", "Plasmid", Sample)) %>%
  group_by(Classification) %>%
  mutate(N = length(unique(Gene))) %>%
  group_by(database, Sample, Classification, N, Gene) %>%
  summarise(ARG = length(Gene)) %>%
  left_join(db.numbers, by = c("database", "Sample")) %>%
  ungroup() %>%
  mutate(p.ARG = ARG/Number) 

#make table of results
options(digits = 1)
proportion.table <- data.full.tidy %>%
  select(database,Sample, Classification, Gene, p.ARG) %>%
  rename(Database = database, `Proportion with ARG` = p.ARG) %>%
  dcast(Classification+Gene~Database+Sample)
proportion.table[is.na(proportion.table)] <- 0


data.full.tidy.s4 <- data.full.tidy %>%
  ungroup() %>%
  select(-c(ARG, Number)) %>%
  dcast(database+Sample~Gene, value.var = "p.ARG") %>%
  mutate_if(is.numeric , replace_na, replace = 0) %>%
  melt(variable.name = "Gene") %>%
  left_join(classification, by = "Gene") %>%
  unique() %>%
  mutate(axis = paste(database, Sample, sep = " "),
         value = value+0.00001) 

f5_stats <- compare_means(value ~ axis,  data = data.full.tidy.s4, paired = FALSE, method = "wilcox")

f5_comparisons <- list(c("RefSeq Chromosome", "RefSeq Plasmid"), c("RefSoil Chromosome", "RefSeq Plasmid"), c("RefSoil Chromosome", "RefSoil Plasmid"), c("RefSeq Chromosome", "RefSoil Chromosome"),  c("RefSeq Plasmid", "RefSoil Plasmid"), c("RefSeq Chromosome", "RefSoil Plasmid"))
options(digits = 3)
(figure_5 <- data.full.tidy.s4 %>%
    ggplot(aes(x = axis, y = value)) +
    geom_boxplot(aes(color = database), outlier.alpha = 0) +
    geom_jitter(aes(fill = Classification),
                color = "black", 
                size = 2,
                shape = 21,
                width = 0.2) +
    scale_fill_manual(values = c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#42d4f4", "#f032e6", "#fabebe", "#469990", "#e6beff", "#9A6324", "#fffac8", "#800000", "#aaffc3", "#000075", "#a9a9a9", "#ffffff", "#000000")) +
    #facet_wrap(~Classification, scales = "free_y") +
    ylab("ARGs per element") +
    xlab("Database") +
    labs(color = "Database") +
    theme_bw(base_size = 8) +
    scale_color_manual(values = c("grey20", "#7fb8d7")) +
    scale_y_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100)) +
    theme(strip.background = element_rect(fill = "white")) +
    stat_compare_means(comparisons = f5_comparisons, method = "wilcox", paired = FALSE, label = "p.signif", size = 2.5) +
    stat_compare_means(method = "kruskal.test", paired = FALSE, label.x.npc = "left", inherit.aes = TRUE, label.y = 5, size = 2.5) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10)))

ggsave(figure_5, filename = "figures/figure_5.eps", units = "in", width = 6.8, height = 5)


(figure_s4 <- data.full.tidy.s4 %>%
    ggplot(aes(x = axis, y = value)) +
    geom_boxplot(outlier.alpha = 0) +
    geom_jitter(aes(fill = Classification),
                color = "black", 
                size = 2,
                shape = 21,
                width = 0.2) +
    scale_fill_manual(values = c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#42d4f4", "#f032e6", "#fabebe", "#469990", "#e6beff", "#9A6324", "#fffac8", "#800000", "#aaffc3", "#000075", "#a9a9a9", "#ffffff", "#000000")) +
    facet_wrap(~Classification,  ncol = 3) +
    ylab("ARGs per element") +
    xlab("Database") +
    labs(color = "Database") +
    theme_bw(base_size = 8) +
    scale_y_log10() +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      legend.position = "bottom"
    ) +
    stat_compare_means(method = "kruskal.test", paired = FALSE, label.x.npc = "left", inherit.aes = TRUE,  label.y.npc = "top", size = 2.5) +
    stat_compare_means(comparisons = f5_comparisons, method = "wilcox", paired = FALSE, label = "p.signif", size = 2.5) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10)))

ggsave(figure_s4, filename = "figures/figure_s4.eps", width = 6.8, height = 8, units = "in")

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
    geom_bar(stat = "identity", width = 1, color = "black") +
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
      plot.title=element_text(size=8, face="bold")
    ) +
    theme(axis.text.x=element_blank(),
          legend.position = "none") +
    geom_text(aes(y = Proportion/2 + c(0, cumsum(Proportion)[-length(Proportion)]), label = percent(Proportion/922)), size=4))

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
    geom_histogram(stat = "count", color = "black", fill = "#A6CEE3") +
    scale_x_continuous(breaks = c(1,3,5,7,9,11,13,15)) +
    theme_classic(base_size = 10)+
    ylab("RefSoil organisms") +
    xlab("Plasmid number"))

#plot and save figure 1
(figure_1 <- ggarrange(refsoil.plasmid, plasmid.hist, labels = c("A", "B"), widths = c(0.7, 1)))
ggsave(figure_1, filename = "figures/figure_1.eps", height = 2.5, width = 5.5, units = "in")

#calculate mean plasmid number
options(digits = 5)
mean(refsoil.plasmids$plasmid)
#1.0065

#calculate variance of plasmid number
options(digits = 5)
var(refsoil.plasmids$plasmid)
#3.1987



#calculate RefSeq plasmid number/ genome
refseq_assembly_chromosome <- data.frame(element = "chromosome", read_delim(file = "data/refseq/chromosome_size_refseq_tab.txt", delim = " "))
refseq_assembly_plasmid <- data.frame(element = "plasmid", read_delim(file = "data/refseq/plasmid_size_refseq_tab.txt", delim = " "))

refseq_assembly <- rbind(refseq_assembly_chromosome, refseq_assembly_plasmid)

refseq_assembly_tidy <- refseq_assembly %>%
  group_by(Assembly, element) %>%
  summarise(Length = length(element)) %>%
  dcast(Assembly ~ element) %>%
  mutate(plasmid = ifelse(is.na(plasmid), 0, plasmid))

var(refseq_assembly_tidy$plasmid)
#var = 2.92

mean(refseq_assembly_tidy$plasmid)
#mean = 0.856

wilcox.test(refseq_assembly_tidy$plasmid, refsoil.plasmids$plasmid, paired = FALSE)
#p = 0.0014

#test again with yes/no
wilcox.test(as.numeric(refseq_assembly_tidy$plasmid>0), as.numeric(refsoil.plasmids$plasmid>0), paired = FALSE)
mean(as.numeric(refseq_assembly_tidy$plasmid>0))
mean(as.numeric(refsoil.plasmids$plasmid>0))

#percentage w plasmids
refseq_percent <- as.numeric(refseq_assembly_tidy$plasmid>0)
length(refseq_percent[refseq_percent> 0])/length(refseq_percent)

#######################################
#ANALYZE REFSOIL GENOME/PLASMID MAKEUP#
#######################################
#temporarily change working directory to data to bulk load files
setwd("data")

#read in size data for refsoil
names=list.files(pattern="*_size_refsoil_tab.txt")
size <- do.call(rbind, lapply(names, function(X) {
  data.frame(id = basename(X), read.table(X))}))

#read in size data for refseq
setwd("refseq")
refseq.size <- read_delim("plasmid_size_refseq_tab.txt", delim = " ")

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

#extract plasmids only 
size.p <- size %>%
  subset(id == "plasmid_size_refsoil_tab.txt") %>%
  mutate(NCBI.ID = as.character(V2))

#plot plasmid sizes
(size.plasmid <- ggplot(size.p, aes(x = V3/1000)) +
    geom_histogram(bins = 30, fill = "#A6CEE3") +
    scale_x_log10(breaks = c(1,5,10,35,100,1000)) +
    ylab("Number of plasmids") +
    xlab("") +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_blank()))

#plot size distribution between refsoil and refseq
(size.plasmid.refseq <- ggplot() +
    #geom_density(data = size.tot, aes(x = Size/1000, fill = "All"), alpha = 0.75) +
    geom_density(data = refseq.size, aes(x = size/1000, fill = "RefSeq"), alpha = 0.75) +
    geom_density(data = size.p, aes(x = V3/1000, fill = "RefSoil"), alpha = 0.75) +
    scale_x_log10(breaks = c(1,5,10,35,100,1000)) +
    scale_fill_manual(name="Database", values=c(All = "black", RefSeq ="grey20", RefSoil="#A6CEE3")) +
    ylab("Density") +
    xlab("Plasmid size (kbp)") +
    theme_classic(base_size = 10))

#plot and save figure 2
(figure_2 <- ggarrange(size.plasmid, size.plasmid.refseq, labels = c("A", "B"),  common.legend = TRUE, ncol = 1, nrow = 2, align = "v"))
ggsave(figure_2, filename = "figures/figure_2.png", height = 4.5, width = 3, units = "in", dpi = 300)

#test for difference 
wilcox.test(refseq.size$size, size.p$V3, paired = FALSE)
size.35 <- size.p[size.p$V3 > 35000,]
median(size.p$V3)

#test for bimodality
library(diptest)
dip.test(refseq.size$size)
#refseq plasmid size is at least bimodal

dip.test(size.p$V3)
#refsoil plasmid size is at least bimodal

#check bimodality coefficients
library(modes)
bimodality_coefficient(size.p$V3)
bimodality_coefficient(refseq.size$size)

lala <- t(modes(size.p$V3,  nmore = 8))
alal <- t(modes(refseq.size$size, type = 2, digits = 4, nmore = 8))

la <- t(modes(refseq.size$size, nmore = 13))

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
    xlab("Genome Size (kbp)") +
    theme_bw(base_size = 12) +
    theme(legend.position = c(0.85, 0.6), 
          legend.background = element_rect(color = "black", 
                                           fill = "white", 
                                           size = 0.05, 
                                           linetype = "solid"),
          legend.title=element_text(size=9), 
          legend.text=element_text(size=8)))

#test for correlation between genome and plasmid size
stats:::cor.test.default(x = size.annotated$genome_tot_size, y = size.annotated$plasmid_mean_size, method = "spearman")

#plasmid density
density.data <- size.annotated %>%
  mutate(plasmid.logical = ifelse(plasmid_n_size == 0, "None", ifelse(plasmid_n_size == 1, "One", "Multiple")),
         plasmid.logical = as.factor(plasmid.logical),
         plasmid.logical = fct_relevel(plasmid.logical, "None", "One", "Multiple")) 

(density.genome <- density.data %>%
    ggplot(aes(x = genome_tot_size/1000, fill = plasmid.logical)) +
    geom_density(alpha = 0.6) +
    labs(fill = "Number of plasmids") +
    scale_fill_manual(values = c("#B2DF8A", "#1F78B4", "#2007ff")) +
    theme_void() +
    theme(legend.position = c(0.85, 0.6), legend.background = element_rect(color = "black", fill = "white", size = 0.05, linetype = "solid"), axis.title.x = element_blank(), legend.title=element_text(size=9), 
          legend.text=element_text(size=8)))

(density.plasmid <- density.data %>%
    ggplot(aes(x = plasmid_tot_size/1000, fill = plasmid.logical)) +
    geom_density(alpha = 0.6) +
    labs(fill = "Number of
         plasmids") +
    scale_x_log10() +
    scale_fill_manual(values = c("#1F78B4", "#2007ff")) +
    theme_void() +
    coord_flip() +
    theme(legend.position = "none", axis.title.x = element_blank()))

library(ggpubr)
(Fig_3 <- ggarrange(density.genome, NULL, plas_v_genome, density.plasmid, common.legend = TRUE, ncol = 2, nrow = 2, align = "hv", heights = c(0.5, 1), widths = c(1, 0.5)))

ggsave(Fig_3, filename = "figures/figure_3.png", height = 4, width = 5, units = "in", dpi = 300)

#set up df for correlation stats
size.stats <- size.annotated %>%
  select(genome_tot_size, plasmid_tot_size, genome_n_size, plasmid_n_size)

#test correlations
library(psych)
corr.stats <- corr.test(size.stats)
corr.stats.r <- corr.stats$r
corr.stats.p <- corr.stats$p

#test difference between genome size based on 
#plasmid number groups
library(broom)
t.test <- density.data %>%
  group_by(plasmid.logical) %>%
  do(tidy(t.test(.$genome_tot_size)))

a <- aov(genome_tot_size ~ plasmid.logical, density.data)
summary(a)

(genome.size.anova <- ggplot(density.data, aes(x = plasmid.logical, y = genome_tot_size/1000)) +
    geom_boxplot() +
    stat_compare_means(method = "anova", label.x = 2.5) +
    ylab("Genome size (kbp)") +
    xlab("Number of plasmids") +
    theme_bw())

ggsave(genome.size.anova, filename = "figures/figure_s2.eps", width = 4, height = 3, units = "in")

###################################
#PLASMID SIZE IN THE ENVIRONMENT#
###################################

#read in env info from RefSoil ISME paper
type <- read_delim("data/ismej2016168x12_soiltype.csv", delim = ",", n_max = 381,
                   col_types = list(col_character(),
                                    col_character(),
                                    col_character(),
                                    col_number(),
                                    col_number(),
                                    col_number(),
                                    col_number(),
                                    col_number(),
                                    col_number(),
                                    col_number(),
                                    col_number(),
                                    col_number()))

#tidy soil type information
#tidy type data
type.tidy <- type %>%
  separate(`RefSoil ID`, into = c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12", "R13", "R14", "R15", "R16", "R17", "R18", "R19", "R20", "R21", "R22", "R23", "R24", "R25", "R26", "R27", "R28", "R29", "R30", "R31", "R32", "R33", "R34", "R35", "R36", "R37", "R38", "R39", "R40", "R41", "R42", "R43"), sep = ",") %>%
  melt(id.vars = c("OTU ID in (Rideout et al., 2014)", "Phylum", "Andisols", "Gelisols", "Vertisols", "Mollisols", "Inceptisols","Alfisols","Ultisols","Sand,Rock,Ice","Entisols"), value.name = "RefSoil ID", na.rm = TRUE) %>%
  select(-variable)


qual <- c("Mollisols", "Alfisols", "Vertisols")

size.type <- ncbi.tidy %>%
  separate(NCBI.ID, into = "NCBI.ID", sep = "\\.") %>%
  right_join(size.p, by = "NCBI.ID") %>%
  rename(Size = V3) %>%
  left_join(type.tidy, by = "RefSoil ID") %>%
  select(`RefSoil ID`, Phylum.x, Organism, NCBI.ID, Size, Andisols:Entisols) %>%
  replace(is.na(.), 0) %>%
  group_by(`RefSoil ID`, Phylum.x, Organism, Size) %>%
  melt(id.vars = c("RefSoil ID", "Phylum.x", "Organism","NCBI.ID", "Size"), variable.name = "Order", value.name = "Count") %>%
  subset(Order %in% qual) 

size.type.final <- size.type[rep(seq_len(dim(size.type)[1]), size.type$Count),]

(figure_2c <- size.type.final %>%
    ggplot(aes(x = Size/1000, fill = Order)) +
    geom_histogram(binwidth = 0.1) +
    facet_wrap(~Order, scales = "free_y", ncol = 1) +
    scale_x_log10(breaks = c(1, 5, 10, 35, 100, 1000)) +
    xlab("Size (kbp)") +
    ylab("Number of plasmids") +
    scale_fill_brewer(palette = "Set2") +
    theme_classic(base_size = 10) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = "top"))
ggsave(figure_2c, filename = "figures/figure_2c.eps", width = 3.8, height = 4.5)

size.type.median <- size.type.final %>%
  group_by(Order) %>%
  summarise(Mean = mean(Size),
            Median = median(Size),
            Variance = var(Size),
            sd = sd(Size),
            N = length(Size))

kruskal.test(Size ~ Order, data = size.type.final)
options(digits = 20)

fx_stats <- compare_means(Size ~ Order,  data = size.type.final, paired = FALSE, method = "wilcox")

tidy(pairwise.wilcox.test(size.type.final$Size, size.type.final$Order, p.adjust.method = "fdr", paired = FALSE))

sx_comparisons <- list(c("Vertisols", "Mollisols"), c("Mollisols", "Alfisols"), c("Vertisols", "Alfisols"))
size.type.final %>%
  ggplot(aes(x = Order, y = Size/1000)) +
  geom_boxplot() +
  stat_compare_means(comparisons = sx_comparisons, method = "wilcox", paired = FALSE,  size = 2.5) +
  scale_y_log10() 

#plasmid size variance
stats::var(size.p$V3)
sd(size.p$V3/1000)
mean(size.p$V3/1000)
ggplot(size.p, aes(x = "a", y = V3)) +
  geom_boxplot() +
  scale_y_log10()





###########################
#PLASMID ARGS NORM TO SIZE#
###########################
data.full.tidy.size <- data.full.tidy %>%
  dcast(Classification+Gene~database+Sample, value.var = "ARG") 
data.full.tidy.size[is.na(data.full.tidy.size)] <- 0

data.full.tidy.size2 <- data.full.tidy.size %>%
  melt() %>%
  separate(variable, into = c("database", "Sample"), sep = "_") %>%
  left_join(db.numbers, by = c("database", "Sample")) %>%
  mutate(size.norm = value/Size,
         axis = paste(database, Sample, sep = " "),
         size.norm = ifelse(size.norm > 0, size.norm, size.norm + 0.00000000001))


options(digits = 3)
(figure_s5 <- data.full.tidy.size2 %>%
    ggplot(aes(x = axis, y = size.norm)) +
    geom_boxplot(aes(color = database), outlier.alpha = 0) +
    geom_jitter(aes(fill = Classification),
                color = "black", 
                size = 2,
                shape = 21,
                width = 0.2) +
    scale_fill_manual(values = c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#42d4f4", "#f032e6", "#fabebe", "#469990", "#e6beff", "#9A6324", "#fffac8", "#800000", "#aaffc3", "#000075", "#a9a9a9", "#ffffff", "#000000")) +
    #facet_wrap(~Classification, scales = "free_y") +
    ylab("ARGs per element") +
    xlab("Database") +
    labs(color = "Database") +
    theme_bw(base_size = 10) +
    scale_color_manual(values = c("grey20", "#7fb8d7")) +
    scale_y_log10() +
    theme(strip.background = element_rect(fill = "white")) +
    stat_compare_means(comparisons = f5_comparisons, method = "wilcox", paired = FALSE, label = "p.signif", size = 2.5) +
    stat_compare_means(method = "kruskal.test", paired = FALSE, label.x.npc = "left", inherit.aes = TRUE, label.y = 1, size = 2.5) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10)))

ggsave(figure_s5, filename = "figures/figure_s5.eps", units = "in", width = 6.8, height = 5)
