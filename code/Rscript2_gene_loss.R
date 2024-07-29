library(ggplot2)
library(RColorBrewer)
library(gplots)
library(ggpubr)
library(ggthemes)
library(ggsci)
library(doBy)
library(reshape2)
library(ggExtra)
library(data.table)
library(fmsb)
library(ggridges)
library(ggrepel)
library(scales)
library(psych)
library(viridis)
library(stringr)
library(readxl)
library(plyr)
library(rotl)
library(phytools)
#library(EBImage)
library(ggtree)
library(RRphylo)
library(tidyr)
library(lattice)
library(dplyr)
library(MCMCtreeR, quietly = TRUE, warn.conflicts = FALSE)
library(lme4)

options(scipen=999)
setwd("~/work/w-degeneration-scirep/results/")


text_size_colour = list(theme_bw(base_family="Helvetica", base_size = 12) + 
                          theme(axis.text.x= element_text(colour="black", size=12)) +
                          theme(axis.text.y= element_text(colour="black", size=12)) +
                          theme(axis.title.x = element_text(colour="black",size=14)) + 
                          theme(axis.title.y = element_text(colour="black",size=14)) + 
                          theme(axis.ticks = element_line(colour = "black", size = 0.5)) +
                          theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)) +
                          theme(plot.title=element_text(family="Helvetica", size=14, colour="black", hjust = 0.5)) +
                          theme(panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                                                colour = "lightgrey"), 
                                panel.grid.minor = element_line(size = 0.15, linetype = 'solid',
                                                                colour = "lightgrey")))

theme_title <- function(...) {
  theme_gray(base_family = "Helvetica") + 
    theme(plot.title = element_text(face = "bold"))
}
title_theme <- calc_element("plot.title", theme_title())

theme_description <- function(...) {
  theme_gray(base_family = "Courier") + 
    theme(plot.title = element_text(face = "plain",size = 12, hjust = 0))
}
theme_description <- calc_element("plot.title", theme_description())

# ==================================== # 
# ==================================== # 
# SECTION 1 - LOAD AND MANIPULATE DATA #
# ==================================== # 
# ==================================== # 

# ==================================== # 
# ======= 1.1 Load data tables ======= #
# ==================================== # 

# Load ZF gene info
zfPep <- read.table("../data/meta/zfPep.index")
zfPep <- plyr::rename(zfPep, c("V5"="Gene","V2"="Chr", "V3"="Start", "V4"="End", "V6"="Trans", "V7"="protein_coding", "V1"="Protein"))
zfPep$Protein <- sub("\\.1", "", zfPep$Protein)
zfPep$Gene <- sub("\\.1", "", zfPep$Gene)
zfPep$Trans <- sub("\\.1", "", zfPep$Trans)
zfPep$Protein <- sub("\\.2", "", zfPep$Protein)
zfPep$Gene <- sub("\\.2", "", zfPep$Gene)
zfPep$Trans <- sub("\\.2", "", zfPep$Trans)
zfPep <- subset(zfPep, select = c(Chr, Start, End, Trans, Gene))
length(zfPep$Trans) # 18204
zfPep <- unique(subset(zfPep, select = -c(Start, End) ))

perGeneData.IQR.noShort <- read.delim("perGeneData.IQR.noShort.tsv")

# ==================================== # 
# ==== 1.2 Summarize perGeneData ===== #
# ==================================== # 

# Calculate min and max values of total number of genes in each species
perGeneData.IQR.noShort %>% count(species) %>% summarize(Min_Value = min(n), Max_Value = max(n)) # 10537 10562

# Calculate total number of genes auto vs. sex-linked
totalGenes <- perGeneData.IQR.noShort %>% count(species, sex.linkage) %>% 
  rename(total_genes = n) %>% pivot_wider(names_from = sex.linkage, values_from = total_genes) %>%
  rename(sex_linked_total_genes = "sex-linked") %>% rename(auto_total_genes = "autosomal") 

# Calculate min and max values of autosomal vs. sex-linked genes in each species (incl "unsure")
totalGenes %>% summarize(Min_Value = min(auto_total_genes), Max_Value = max(auto_total_genes)) # 8920      9802
totalGenes %>% summarize(Min_Value = min(sex_linked_total_genes), Max_Value = max(sex_linked_total_genes)) # 755      1642

# Calculate counts of each lost_status category, for autosomal vs. sex-linked genes in each species
summaryGenes <- perGeneData.IQR.noShort %>% count(species, sex.linkage, lost_status)
summaryGenes.table <- summaryGenes %>% pivot_wider(names_from = lost_status, values_from = n)
summaryGenes.table$total_nr_genes <- summaryGenes.table$diploid_female + summaryGenes.table$haploid_female + summaryGenes.table$unsure

# Supplementary Table
outname <- sprintf("tables/summaryGenes.table.tsv")
write.table(summaryGenes.table, file = outname, sep = "\t", quote = FALSE, row.names = F)

# Perc. of haploid genes in autosomal vs. sex-linked regions (proxy for error rate)
summaryGenes.table %>% group_by(sex.linkage) %>% summarize(
  Mean = mean(haploid_female/total_nr_genes*100), 
  Min = min(haploid_female/total_nr_genes*100), 
  Max = max(haploid_female/total_nr_genes*100))
# sex.linkage   Mean    Min   Max
#1 autosomal    0.692  0.437  1.09
#2 sex-linked  59.4   35.6   71.0 

# Perc. of unsure genes in autosomal vs. sex-linked regions <-- GENES WERE REMOVED
summaryGenes.table %>% group_by(sex.linkage) %>% summarize(
  Mean = mean(unsure/total_nr_genes*100), 
  Min = min(unsure/total_nr_genes*100), 
  Max = max(unsure/total_nr_genes*100))
#sex.linkage  Mean    Min   Max
#1 autosomal   0.201 0.113 0.332
#2 sex-linked  1.94  0.974 3.01 


ggplot(perGeneData.IQR.noShort, aes(x = Gene.cov_ratio_scaled, y = Gene.diff.hetSites.prop.mean, colour = lost_status, fill = lost_status)) +
  geom_point(alpha = 0.8) + scale_colour_jco() +
  text_size_colour + 
  ylab("female-male heterozygosity") +
  xlab("normalized female/male genome coverage") +
  facet_wrap(~species + sex.linkage, ncol = 2)

# ======================================================================== # 
# = 1.3 Extracting only genes in sex-linked strata and no "unsure" genes = #
# ======================================================================== # 

df <- perGeneData.IQR.noShort
df %>% count(lost_status) # Count rows with different categories of lost_status
#1 diploid_female 11991
#2 haploid_female  4547


df <- subset(df, c(df$lost_status!="unsure" & df$sc.region!="Autosomal"))
df_new <- subset(df, select=c(Trans, species,sc.region,lost_status,sex.linkage))

# Calculate counts of each lost_status category, for autosomal vs. sex-linked genes in each species (Figure 4c stats)
summaryGenes <- df_new %>% count(species, sex.linkage, lost_status)
summaryGenes.table <- summaryGenes %>% pivot_wider(names_from = lost_status, values_from = n)
summaryGenes.table$total_nr_genes <- summaryGenes.table$diploid_female + summaryGenes.table$haploid_female 

wide_df <- reshape2::dcast(df_new, Trans + sc.region ~ species, value.var = "lost_status")
complete_rows <- complete.cases(wide_df)
wide_df <- wide_df[complete_rows, ]
wide_df %>% count()
wide_df %>% count(sc.region) # Count number of genes remaining on each strata
#sc.region   n
#1    chr3_a  89
#2    chr3_b  29
#3    chr3_c 479
#4      chr4 262
#5   chr4A_a  64
#6   chr4A_b  32
#7      chr5 345
#8      chr8 163
#9      chrZ 468

# ============================================== # 
# = 1.4 Naming phylogenetic branches and nodes = #
# ============================================== #

### This code generates a phylogenetic tree where each node is labeled by the branch number leading up to the node
animals = c("Alauda arvensis", "Alauda razae", "Panurus biarmicus", "Cisticola juncidis", 
            "Sylvietta_brachyura", "Calandrella_cinerea", "Camaroptera_brevicaudata", 
            "Eremophilis_alpestris", "Taeniopygia guttata")
(resolved_names <- tnrs_match_names(animals, context_name = "Animals"))
my_tree <- tol_induced_subtree(ott_ids=ott_id(resolved_names))

my_tree$tip.label <- sub("_ott[0-9]+", "", my_tree$tip.label)
my_tree$tip.label <- my_tree$tip.label

my_tree <- rotateNodes(my_tree, 15)

pdf("phyloTree_fully_named.pdf", width = 10, height = 5)
par(family = "Helvetica")
node_nr <- c("NA","n1","n6","n7","n2","n3","n4","n5")
plotTree(my_tree, offset=3)
nodelabels(node_nr, cex=0.8)
Midbranch_nr <- c("NA","b1","b6","b7","bCbre","bCjun","bSbra","b2","b3","b5","bAarv","bAraz","b4","bEalp","bCcin","bPbia")
edgelabels(Midbranch_nr, cex=0.8)
tiplabels(c("NA", "tCbre", "tCjun", "tSbra", "tAarv", "tAraz", "tEalp","tCcin","tPbia"), cex=0.8)
dev.off()


phy.MCMC <- readMCMCtree("FigTree.tre")

MCMC.tree.plot(phy.MCMC, analysis.type = "MCMCtree", cex.tips = 2, cex.age = 1.5, cex.labels = 2,
               time.correction = 100, plot.type = "phylogram", lwd.bar = 10, add.abs.time = TRUE, abs.age.lwd = TRUE,
               scale.res = c("Eon", "Period"), node.method = "bar", col.age = "#02315E90", no.margin = TRUE, label.offset = 4)


phy <- read.tree("FigTree_noNodeAges.tre")
plot(phy)


tip_species <- c("Dromaius_novaehollandiae", "Apteryx_haastii","Lepidothrix_coronata", "Alauda_razae", "Alauda_arvensis", "Eremophilis_alpestris", "Calandrella_cinerea", "Panurus_biarmicus", "Camaroptera_brevicaudata",  "Cisticola_juncidis", "Sylvietta_brachyura", "Taeniopygia_guttata", "Ficedula_albicollis", "Melopsittacus_undulatus", "Gallus_gallus", "Anas_platyrhynchos", "Anolis_carolinensis")


phy$tip.label <- tip_species

plot(phy)
nodelabels()
phy <- rotateNodes(phy, 31)

pdf("plots/phylogeny_node_labels_ages_noGT.pdf", height = 15, width = 25)
par(family = "Helvetica")
plot.phylo(phy, type = "tidy", show.tip.label = TRUE, label.offset = 0.01)
nodelabels(round(phy.MCMC$nodeAges[,1]*100, digits = 1), bg = "lightgray")
dev.off()

pdf("plots/phylogeny_node_labels_ages_noGT_cladogram.pdf", height = 15, width = 15)
par(family = "Helvetica")
plot.phylo(phy, type = "cladogram", use.edge.length = FALSE,  show.tip.label = TRUE, label.offset = 0.01)
nodelabels(round(phy.MCMC$nodeAges[,1]*100, digits = 1), bg = "lightgray")
dev.off()

# Node ages (Values from the MCMCTree file): 
node0Age <- 29.7
node1Age <- 26.4
node2Age <- 22.2
node3Age <- 10.4
node4Age <- 6.7
node5Age <- 4.5
node6Age <- 23.6
node7Age <- 12.3

node0AgeLabel <- paste0("n0:", node0Age)
node1AgeLabel <- paste0("n1:", node1Age)
node2AgeLabel <- paste0("n2:", node2Age)
node3AgeLabel <- paste0("n3:", node3Age)
node4AgeLabel <- paste0("n4:", node4Age)
node5AgeLabel <- paste0("n5:", node5Age)
node6AgeLabel <- paste0("n6:", node6Age)
node7AgeLabel <- paste0("n7:", node7Age)


# Calculating fusion ages (midpoint of branch on which it appears)
chrZFusionAge <- 100
chr4A_aFusionAge <- node1Age+(node0Age-node1Age)/2
chr4A_bFusionAge_n2 <- node2Age+(node1Age-node2Age)/2
chr3_aFusionAge <- node2Age+(node1Age-node2Age)/2
chr3_bFusionAge <- node3Age+(node2Age-node3Age)/2
chr3_cFusionAge <- node5Age+(node3Age-node5Age)/2
chr5FusionAge <- node5Age+(node3Age-node5Age)/2
chr4FusionAge <- node7Age+(node6Age-node7Age)/2
chr4A_bFusionAge_n7 <- node7Age+(node6Age-node7Age)/2
chr8FusionAge_n7 <- 0+(node6Age-0)/2

chrZFusionAgeLabel <- paste0("chrZ:", chrZFusionAge)
chr4A_aFusionAgeLabel <- paste0("chr4A_a:", chr4A_aFusionAge)
chr4A_bFusionAgeLabel_n2 <- paste0("chr4A_b:", chr4A_bFusionAge_n2)
chr3_aFusionAgeLabel <- paste0("chr3_a:", chr3_aFusionAge)
chr3_bFusionAgeLabel <- paste0("chr3_b:", chr3_bFusionAge)
chr3_cFusionAgeLabel <- paste0("chr3_c:", chr3_cFusionAge)
chr5FusionAgeLabel <- paste0("chr5:", chr5FusionAge)
chr4FusionAgeLabel <- paste0("chr4:", chr4FusionAge)
chr4A_bFusionAgeLabel_n7 <- paste0("chr4A_b:", chr4A_bFusionAge_n7)
chr8FusionAgeLabel <- paste0("chr8:", chr8FusionAge_n7)


chrZ_4AaFusionAgeLabel <- paste0(chrZFusionAgeLabel, " \n ", chr4A_aFusionAgeLabel)
chr3a_4AbFusionAgeLabel <- paste0(chr3_aFusionAgeLabel, " \n ", chr4A_bFusionAgeLabel_n2)
chr4_4AbFusionAgeLabel <- paste0(chr4FusionAgeLabel, " \n ", chr4A_bFusionAgeLabel_n7)
chr3c_5FusionAgeLabel <- paste0(chr3_cFusionAgeLabel, " \n ", chr5FusionAgeLabel)


# Specify the tips to be dropped

tips_to_drop <- c("Apteryx_haastii", "Dromaius_novaehollandiae", "Gallus_gallus", "Anas_platyrhynchos", 
                  "Melopsittacus_undulatus", "Ficedula_albicollis", "Lepidothrix_coronata", "Anolis_carolinensis")

# Drop the specified tips from the tree
tree_dropped <- drop.tip(phy, tips_to_drop)

plot(tree_dropped)
# Original order of tips
original_order <- tree_dropped$tip.label
print(original_order)

# Swap the order of tips
new_order <- c("Panurus_biarmicus", "Calandrella_cinerea", "Eremophilis_alpestris", "Alauda_arvensis", "Alauda_razae",       
               "Sylvietta_brachyura", "Cisticola_juncidis", "Camaroptera_brevicaudata", "Taeniopygia_guttata" )  

tree <- rotateConstr(tree_dropped, rev(new_order))
pdf("plots/studySpecies_phylogeny_node_branch_labels.pdf", height = 7, width = 4)
par(family = "Helvetica")
node_nr <- c("NA","n1","n6","n7","n2","n3","n5","n4")
#node_nr <- c("NA","n1","n2","n3","n5","n4","n6","n7")
plotTree(tree, use.edge.length = TRUE, offset=3)
nodelabels(node_nr, cex=0.8, )
Midbranch_nr <- c("NA","b1","b6","b7","bCbre","bCjun","bSbra","b2","b3","b5","bAarv","bAraz","b4","bEalp","bCcin","bPbia")
#Midbranch_nr <- c("NA","b1","b6","b7","bCbre","bCjun","bSbra","b2","b3","b4","bEalp","bCcin","b5","bAarv","bAraz","bPbia")
edgelabels(Midbranch_nr, cex=0.8)
tiplabels(c("tCbre", "tCjun", "tSbra", "tAarv", "tAraz", "tCcin", "tEalp", "tPbia", "tTgut"), cex=0.8)
dev.off()

pdf("plots/studySpecies_phylogeny_node_strata_ages.pdf", height = 7, width = 10)
par(family = "Helvetica")
plotTree(tree, use.edge.length = TRUE, offset=3)
#Midbranch_nr <- c("NA",chrZ_4AaFusionAgeLabel,NA,chr4_4AbFusionAgeLabel,NA,NA,chr8FusionAgeLabel,chr3a_4AbFusionAgeLabel,chr3_bFusionAgeLabel,NA,NA,NA,chr3c_5FusionAgeLabel,NA,NA,NA)
Midbranch_nr <- c("NA",chrZ_4AaFusionAgeLabel,NA,chr4_4AbFusionAgeLabel,NA,NA,chr8FusionAgeLabel,chr3a_4AbFusionAgeLabel,chr3_bFusionAgeLabel,chr3c_5FusionAgeLabel,NA,NA,NA,NA,NA,NA)
edgelabels(Midbranch_nr, cex=0.8, adj = c(0.8,-1), frame = "none")
node_ages <- c(node0Age,node1AgeLabel, node6AgeLabel, node7AgeLabel, node2AgeLabel, node3AgeLabel, node4AgeLabel, node5AgeLabel)
nodelabels(node_ages, cex=0.8, bg = "lightgray")
edgelabels(pch = 21, cex=0.8)
dev.off()

# ============================================= # 
# ==== 1.5 Evaluating timing of W gene loss === #
# ============================================= #

wide_df <- wide_df %>% mutate(node1 = 0, node2 = 0, node3 = 0, node4 = 0, node5 = 0, node6 = 0, node7 = 0)

matching_rows <- wide_df$Alauda_arvensis == "diploid_female" | wide_df$Alauda_razae == "diploid_female"
wide_df$node5[matching_rows] <- 1

matching_rows <- wide_df$Eremophilis_alpestris =="diploid_female" | wide_df$Calandrella_cinerea =="diploid_female"
wide_df$node4[matching_rows] <- 1

matching_rows <- wide_df$Camaroptera_brevicaudata == "diploid_female" | wide_df$Cisticola_juncidis == "diploid_female"
wide_df$node7[matching_rows] <- 1

matching_rows <- wide_df$node5 == 1 | wide_df$node4 == 1
wide_df$node3[matching_rows] <- 1

matching_rows <- wide_df$node3 == 1 | wide_df$Panurus_biarmicus == "diploid_female"
wide_df$node2[matching_rows] <- 1

matching_rows <- wide_df$node7 == 1 | wide_df$Sylvietta_brachyura == "diploid_female"
wide_df$node6[matching_rows] <- 1

matching_rows <- wide_df$node6 == 1 | wide_df$node2 == 1
wide_df$node1[matching_rows] <- 1

wide_df %>% count(node6) 




outname <- sprintf("tables/allGenesBranchLeafStatus.tsv")
write.table(wide_df, file = outname, sep = "\t", quote = FALSE, row.names = F)



long_df <- melt(wide_df, id.vars = c("Trans", "sc.region", "node1", "node2","node3", "node4","node5", "node6","node7"))
long_df <- long_df %>% rename(species = variable)
long_df$value[which(c(long_df$value=="diploid_female"))] <- 1
long_df$value[which(c(long_df$value=="haploid_female"))] <- 0


# Add info on whether region is sex-linked or not
df_sc <- unique(subset(df_new, select=c(sc.region, species, sex.linkage)))
long_df.sc <- merge(long_df, df_sc, by = c("sc.region", "species"))

outname <- sprintf("tables/gene_long_table.tsv")
write.table(long_df.sc, file = outname, sep = "\t", quote = FALSE, row.names = F)


# Count number of genes being lost on each branch
column_names <- colnames(wide_df)[-c(1, 2)]
column_names <- gsub("node", "", column_names)

all_combinations <- expand.grid(sc.region = unique(wide_df$sc.region),
                                node_nr = column_names)


b1_loss <- wide_df[which(wide_df$node1== 0),] %>% count(sc.region) %>% mutate(node_nr = 1)
b2_loss <- wide_df[which(c(wide_df$node1== 1 & wide_df$node2== 0)),] %>% count(sc.region) %>% mutate(node_nr = 2)
b3_loss <- wide_df[which(c(wide_df$node2== 1 & wide_df$node3== 0)),] %>% count(sc.region) %>% mutate(node_nr = 3)
b4_loss <- wide_df[which(c(wide_df$node3== 1 & wide_df$node4== 0)),] %>% count(sc.region) %>% mutate(node_nr = 4)
b5_loss <- wide_df[which(c(wide_df$node3== 1 & wide_df$node5== 0)),] %>% count(sc.region) %>% mutate(node_nr = 5)
b6_loss <- wide_df[which(c(wide_df$node1== 1 & wide_df$node6== 0)),]  %>% count(sc.region) %>% mutate(node_nr = 6)
b7_loss <- wide_df[which(c(wide_df$node6== 1 & wide_df$node7== 0)),]  %>% count(sc.region) %>% mutate(node_nr = 7)

Aarv_tip_loss <- long_df[which(c(long_df$species=="Alauda_arvensis" & long_df$node5== 1 & long_df$value== 0)),] %>% count(sc.region) %>% mutate(node_nr = "Alauda_arvensis")
Araz_tip_loss <- long_df[which(c(long_df$species=="Alauda_razae" & long_df$node5== 1 & long_df$value== 0)),] %>% count(sc.region) %>% mutate(node_nr = "Alauda_razae")
Pbia_tip_loss <- long_df[which(c(long_df$species=="Panurus_biarmicus" & long_df$node2== 1 & long_df$value== 0)),] %>% count(sc.region) %>% mutate(node_nr = "Panurus_biarmicus")
Ealp_tip_loss <- long_df[which(c(long_df$species=="Eremophilis_alpestris" & long_df$node4== 1 & long_df$value== 0)),] %>% count(sc.region) %>% mutate(node_nr = "Eremophilis_alpestris")
Ccin_tip_loss <- long_df[which(c(long_df$species=="Calandrella_cinerea" & long_df$node4== 1 & long_df$value== 0)),] %>% count(sc.region) %>% mutate(node_nr = "Calandrella_cinerea")
Cjun_tip_loss <- long_df[which(c(long_df$species=="Cisticola_juncidis" & long_df$node7== 1 & long_df$value== 0)),] %>% count(sc.region) %>% mutate(node_nr = "Cisticola_juncidis")
Cbre_tip_loss <- long_df[which(c(long_df$species=="Camaroptera_brevicaudata" & long_df$node7== 1 & long_df$value== 0)),] %>% count(sc.region) %>% mutate(node_nr = "Camaroptera_brevicaudata")
Sbra_tip_loss <- long_df[which(c(long_df$species=="Sylvietta_brachyura" & long_df$node6== 1 & long_df$value== 0)),] %>% count(sc.region) %>% mutate(node_nr = "Sylvietta_brachyura")

lossTable <- rbind(b1_loss, b6_loss, b7_loss, Cjun_tip_loss, Cbre_tip_loss, Sbra_tip_loss, b2_loss, b3_loss, b4_loss, 
                   Ealp_tip_loss, Ccin_tip_loss, b5_loss, Aarv_tip_loss, Araz_tip_loss, Pbia_tip_loss)

# Complete the dataframe with missing combinations and set Value to zero
lossTable <- all_combinations %>%
  left_join(lossTable, by = c("sc.region", "node_nr")) %>%
  mutate(Value = replace_na(n, 0))

lossTable <- select(lossTable, -3)

outname <- sprintf("tables/lossTable.table.tsv")
write.table(lossTable, file = outname, sep = "\t", quote = FALSE, row.names = F)


# ========================================= # 
# === Counting present and lost W genes === #
# ========================================= #


### Count number of original genes per strata
wide_df %>% count(sc.region) 
orig_genes <- wide_df %>% count(sc.region) %>% rename(orig_nr_genes = n)

#sc.region   n
#1    chr3_a  89
#2    chr3_b  29
#3    chr3_c 479
#4      chr4 262
#5   chr4A_a  64
#6   chr4A_b  32
#7      chr5 345
#8      chr8 163
#9      chrZ 468


## These lines calculates how many unique chrZ genes lost their W gametolog within Sylvioidea
long_df.node1.1 <- subset(long_df, c(long_df$node1==1 & long_df$sc.region=="chrZ")) 
long_df.node1.1 <- long_df.node1.1 %>% filter(across(3:11, ~ . == 0)%>% rowSums() > 0)
length(unique(long_df.node1.1$Trans)) # 17

long_df.count <- long_df %>% count(sc.region, species, value) 
long_df.count <- merge(long_df.count, orig_genes, by = "sc.region")
long_df.count <- merge(long_df.count, df_sc, by = c("sc.region", "species"))

long_df.count.1 <-subset(long_df.count, long_df.count$value==1) %>% rename(remaining_nr_genes = n)
long_df.count.1$prop.remaining <- long_df.count.1$remaining_nr_genes / long_df.count.1$orig_nr_genes
long_df.count.1$freq.int <- round(long_df.count.1$prop.remaining, digits = 2)
long_df.count.1$perc.remaining <- long_df.count.1$prop.remaining * 100

#subset(long_df.count, long_df.count$sex.linkage=="sex-linked") 


# Summary table: How many sex-linked W genes remain per species and strata?
gene_count <- subset(long_df.count.1, long_df.count.1$sex.linkage=="sex-linked") %>%
  group_by(species) %>%
  summarize(Sum_orig = sum(orig_nr_genes), Sum_remain = sum(remaining_nr_genes))
gene_count$gene_lost <- gene_count$Sum_orig - gene_count$Sum_remain
gene_count

#species                  Sum_orig Sum_remain gene_lost
#<fct>                       <int>      <int>     <int>
#1 Alauda_arvensis              1506        969       537
#2 Alauda_razae                 1506        977       529
#3 Calandrella_cinerea           682        175       507
#4 Camaroptera_brevicaudata      826        297       529
#5 Cisticola_juncidis            826        277       549
#6 Eremophilis_alpestris         682        179       503
#7 Panurus_biarmicus             653        182       471
#8 Sylvietta_brachyura           695        222       473

# Mean orig. genes values Alauda species
subset(gene_count, c(gene_count$species=="Alauda_arvensis" | gene_count$species=="Alauda_razae")) %>%
  summarize(Mean = mean(Sum_orig)) # 1506

# Mean orig. genes NOT values Alauda species
subset(gene_count, c(gene_count$species!="Alauda_arvensis" & gene_count$species!="Alauda_razae")) %>%
  summarize(Mean = mean(Sum_orig)) # 727

# Mean gene loss across all species and strata
gene_count %>%
  summarize(Mean = mean(gene_lost)) # 512

# Mean, Min and Max prop remaining genes per strata
subset(long_df.count.1, long_df.count.1$sex.linkage=="sex-linked") %>%
  group_by(c(sc.region)) %>%
  summarize(Mean = mean(perc.remaining), Min = min(perc.remaining), Max = max(perc.remaining))

#`c(sc.region)`  Mean   Min    Max
#<chr>          <dbl> <dbl>  <dbl>
#1 chr3_a         68.3  62.9   79.8 
#2 chr3_b         80.2  79.3   82.8 
#3 chr3_c         97.7  97.5   97.9 
#4 chr4           77.7  74.4   80.9 
#5 chr4A_a        78.5  68.8   89.1 
#6 chr4A_b        76.3  65.6  100   
#7 chr5           95.8  95.4   96.2 
#8 chr8           90.8  90.8   90.8 
#9 chrZ            4.17  3.42   4.70


4.17/100

# ChrZ
perc_geneloss <- 100-4.17
chrZ_percMyrGeneLoss <- perc_geneloss / 100

# Chr4a-a
perc_geneloss <- 100-78.5
chrZ_percMyrGeneLoss <- perc_geneloss / 100


# Mean remaining genes for mean values of all added strata (i.e., not chrZ)
subset(long_df.count.1, c(long_df.count.1$sex.linkage=="sex-linked" & long_df.count.1$sc.region!="chrZ")) %>%
  group_by(c(sc.region)) %>% summarize(Mean = mean(perc.remaining)) %>% summarize(Mean = mean(Mean)) 

# Mean, Min and Max prop of "false Wlost genes" (haploid in female despite not being in a sex-linked region)
subset(long_df.count.1, long_df.count.1$sex.linkage!="sex-linked") %>%
  group_by(c(sc.region)) %>%
  summarize(Mean = mean(perc.remaining), Min = min(perc.remaining), Max = max(perc.remaining))

# Figure 2b
Fig2a_data <- as.data.frame(gene_count) #%>% rename("species" = "c(species)")
p1 <- ggplot(Fig2a_data, aes(x = species)) +
  geom_bar(aes(y = Sum_orig, fill = "Sum_orig"), color = "black", stat = "identity") +
  geom_bar(aes(y = gene_lost, fill = "gene_lost"), color = "black", stat = "identity") +
  theme_minimal() +
  coord_flip() +
  scale_x_discrete(limits=rev(c("Panurus_biarmicus", "Calandrella_cinerea", "Eremophilis_alpestris", "Alauda_arvensis", "Alauda_razae","Sylvietta_brachyura", "Cisticola_juncidis","Camaroptera_brevicaudata"))) +
  #scale_x_discrete(limits=rev(c("Panurus_biarmicus","Alauda_arvensis", "Alauda_razae", "Calandrella_cinerea", #"Eremophilis_alpestris", "Sylvietta_brachyura", "Cisticola_juncidis","Camaroptera_brevicaudata"))) +
  scale_fill_manual(values = c("#FEF0D9", "#B30000")) + text_size_colour + theme(axis.text.y = element_blank(),
                                                                                 axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position="bottom")

col.number <- colorRampPalette(brewer.pal(6, 'OrRd'))(10)
# Figure 2a - heatmap perc genes remaining
p2 <- ggplot(long_df.count.1, aes(sc.region, species, fill= prop.remaining*100)) + 
  geom_tile() + geom_tile(color="white", linewidth=0.7) + 
  # geom_label(aes(label=freq.int), size = 5, fill = "white", position = position_nudge(y = 0.15)) +
  geom_label(aes(label=paste0(freq.int*100, "%")), size = 5, fill = "white", position = position_nudge(y = 0.15)) +
  #geom_text(stat = "prop", size = 5, fill = "white", position = position_nudge(y = 0.15)) +
  geom_point(aes(sc.region, species),data = . %>% filter(sex.linkage=="sex-linked"), pch = 21, color = "white", fill = "black", size = 5, position = position_nudge(y = -0.15) ) +
  labs(x = "Sex-linked region", y =" ", fill = "% remaining W genes") +  
  scale_y_discrete(limits=rev(c("Panurus_biarmicus", "Calandrella_cinerea", "Eremophilis_alpestris", "Alauda_arvensis", "Alauda_razae","Sylvietta_brachyura", "Cisticola_juncidis","Camaroptera_brevicaudata"))) +
  scale_x_discrete(limits=c("chrZ", "chr4A_a", "chr4A_b", "chr3_a", "chr3_b", "chr3_c", "chr5","chr8", "chr4")) +
  text_size_colour + scale_fill_gradient(high = "#FEF0D9", low = "#B30000") + theme(legend.position="bottom")

ggplot(long_df.count.1, aes(sc.region, species, fill= prop.remaining*100)) + 
  geom_tile() + geom_tile(color="white", linewidth=0.7) + 
  # geom_label(aes(label=freq.int), size = 5, fill = "white", position = position_nudge(y = 0.15)) +
  geom_label(aes(label=paste0(freq.int*100, "%")), size = 5, fill = "white", position = position_nudge(y = 0.15)) +
  #geom_text(stat = "prop", size = 5, fill = "white", position = position_nudge(y = 0.15)) +
  geom_point(aes(sc.region, species),data = . %>% filter(sex.linkage=="sex-linked"), pch = 21, color = "white", fill = "black", size = 5, position = position_nudge(y = -0.15) ) +
  labs(x = "Sex-linked region", y =" ", fill = "% remaining W genes") +  
  scale_y_discrete(limits=rev(c("Panurus_biarmicus", "Calandrella_cinerea", "Eremophilis_alpestris", "Alauda_arvensis", "Alauda_razae","Sylvietta_brachyura", "Cisticola_juncidis","Camaroptera_brevicaudata"))) +
  scale_x_discrete(limits=c("chrZ", "chr4A_a", "chr4A_b", "chr3_a", "chr3_b", "chr3_c", "chr5","chr8", "chr4")) +
  text_size_colour + scale_fill_gradient(high = "#FEF0D9", low = "#B30000") + theme(legend.position="bottom")

pdf("plots/Fig2a_decimals.pdf", height = 10, width = 15)
ggplot(long_df.count.1, aes(sc.region, species, fill= prop.remaining*100)) + 
  geom_tile() + geom_tile(color="white", linewidth=0.7) + 
  # geom_label(aes(label=freq.int), size = 5, fill = "white", position = position_nudge(y = 0.15)) +
  geom_label(aes(label=round(prop.remaining, digits = 3)), size = 5, fill = "white", position = position_nudge(y = 0.15)) +
  #geom_text(stat = "prop", size = 5, fill = "white", position = position_nudge(y = 0.15)) +
  geom_point(aes(sc.region, species),data = . %>% filter(sex.linkage=="sex-linked"), pch = 21, color = "white", fill = "black", size = 5, position = position_nudge(y = -0.15) ) +
  labs(x = "Sex-linked region", y =" ", fill = "% remaining W genes") +  
  scale_y_discrete(limits=rev(c("Panurus_biarmicus", "Calandrella_cinerea", "Eremophilis_alpestris", "Alauda_arvensis", "Alauda_razae","Sylvietta_brachyura", "Cisticola_juncidis","Camaroptera_brevicaudata"))) +
  scale_x_discrete(limits=c("chrZ", "chr4A_a", "chr4A_b", "chr3_a", "chr3_b", "chr3_c", "chr5","chr8", "chr4")) +
  text_size_colour + scale_fill_gradient(high = "#FEF0D9", low = "#B30000") + theme(legend.position="bottom") 
dev.off()

tree_noZF <-drop.tip(tree, "ZF")
tree_plot <- ggtree(tree_noZF, ladderize = F) 
legend_1 <- get_legend(p1)
legend_2 <- get_legend(p2)
legends <- ggarrange(legend_2, legend_1, nrow=1)
rm_legend <- function(p){p + theme(legend.position = "none")}
plots <- ggarrange(tree_plot, rm_legend(p2), rm_legend(p1), widths =  c(0.15,1, 0.35), nrow = 1, labels = c(" ", "a", "b"))


pdf("plots/Figure2_tree.pdf", height = 10, width = 15)
#ggarrange(tree_plot, p2,p1, widths = c(0.15, 1, 0.35), heights = c(0.5, 1,1), common.legend = TRUE, nrow=1, labels = c(" ", "A", "B"))
ggarrange(plots, legends, nrow=2, heights = c(1,0.1))
dev.off()


pdf("plots/Figure2.pdf", height = 10, width = 15)
cowplot::plot_grid(p2 + theme(legend.position="bottom"), 
                   p1 + theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_blank(), 
                              legend.position="bottom"), 
                   nrow = 1,
                   rel_widths = c(1,0.33),
                   labels = "auto", align = "h", vjust = 0.5, hjust = 0.5) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
dev.off()



# Plot gene loss on each branch
col.prop <- colorRampPalette(c("black", "white"))
col.prop <- colorRampPalette(brewer.pal(5, 'Spectral'))(10)
col.number <- colorRampPalette(rev(brewer.pal(5, 'Spectral')))(10)
show_col(colorRampPalette(brewer.pal(6, 'OrRd'))(10))
col.number <- colorRampPalette(brewer.pal(6, 'OrRd'))(10)
# Define a vector to specify the sorting order
sorting_vector <- c("1", "6", "7", "Camaroptera_brevicaudata","Cisticola_juncidis",
                    "Sylvietta_brachyura", "2", "3", "5", "Alauda_razae", "Alauda_arvensis",
                    "4", "Eremophilis_alpestris","Calandrella_cinerea", 
                    "Panurus_biarmicus")

# Sort the dataframe based on the "ID" column using the sorting vector
df_sorted <- lossTable[order(match(lossTable$node_nr, sorting_vector)), ] 
df_sorted <- df_sorted %>% replace_na(list(n = 0))

layout(matrix(1:3,1,3),widths=c(0.39,0.22,0.39))
plot(tree)
edgelabels(bg = "white", frame="rect", cex = 2)
plot.new()

col.number <- colorRampPalette(brewer.pal(9, 'OrRd'))(10)
col.number <- col.number[-1]
col.number <- col.number[-1]
## Chr Z
lost <- c(0, subset(df_sorted$Value, df_sorted$sc.region=="chrZ"))
thing <- rep(16, x = NA)
thing[2] <- "*"
pdf("plots/chrZ_geneloss.pdf", width = 6, height = 8)
par(family = "Helvetica")
plotBranchbyTrait(tree, lost, mode="edge", y.lim = c(0,10), xlims = c(min(lost),max(lost)), show.tip.label = FALSE, palette = colorRampPalette(col.number, bias = 1), outline = TRUE, edge.width = 10, cex= 1, font = 1, legend=TRUE)
edgelabels(lost, bg = "white", frame = "rect", cex = 2)
edgelabels(thing ,adj = c(0.5,-0.5), frame="none", cex = 4)
title(main="chrZ (n = 468)", line = -3, cex.main = 2)
dev.off()

## Chr 4A_a
pdf("plots/chr4A_a_geneloss.pdf", width = 6, height = 8)
par(family = "Helvetica")
lost <- c(0, subset(df_sorted$Value, df_sorted$sc.region=="chr4A_a"))
thing <- rep(16, x = NA)
thing[2] <- "*"
plotBranchbyTrait(tree, lost, mode="edge", y.lim = c(0,10), xlims = c(min(lost),max(lost)), show.tip.label = FALSE, palette = colorRampPalette(col.number, bias = 1), outline = TRUE, edge.width = 10, cex= 1, font = 1, legend=TRUE)
edgelabels(lost, bg = "white", frame="rect", cex = 2)
edgelabels(thing ,adj = c(0.5,-0.5), frame="none", cex = 4)
title(main="chr4A_a (n = 64)", line = -3, cex.main = 2)
dev.off()

## Chr 4A_b
pdf("plots/chr4A_b_geneloss.pdf", width = 6, height = 8)
thing <- rep(16, x = NA)
thing[4] <- "*"
thing[8] <- "*"
par(family = "Helvetica")
lost <- c(0, subset(df_sorted$Value, df_sorted$sc.region=="chr4A_b"))
plotBranchbyTrait(tree, lost, mode="edge", y.lim = c(0,10), xlims = c(min(lost),max(lost)), show.tip.label = FALSE, palette = colorRampPalette(col.number, bias = 1), outline = TRUE, edge.width = 10, cex= 1, font = 1, legend=TRUE)
edgelabels(lost, bg = "white", frame="rect", cex = 2)
edgelabels(thing ,adj = c(0.5,-0.5), frame="none", cex = 4)
title(main="chr4A_b (n = 32)", line = -3, cex.main = 2)
dev.off()

## Chr chr3_a
pdf("plots/chr3_a_geneloss.pdf", width = 6, height = 8)
par(family = "Helvetica")
thing <- rep(16, x = NA)
thing[8] <- "*"
lost <- c(0, subset(df_sorted$Value, df_sorted$sc.region=="chr3_a"))
plotBranchbyTrait(tree, lost, mode="edge", y.lim = c(0,10), xlims = c(min(lost),max(lost)), show.tip.label = FALSE, palette = colorRampPalette(col.number, bias = 1), outline = TRUE, edge.width = 10, cex= 1, font = 1, legend=TRUE)
edgelabels(lost, bg = "white", frame="rect", cex = 2)
edgelabels(thing ,adj = c(0.5,-0.5), frame="none", cex = 4)
title(main="chr3_a (n = 89)", line = -3, cex.main = 2)
dev.off()

## Chr chr3_b
pdf("plots/chr3_b_geneloss.pdf", width = 6, height = 8)
par(family = "Helvetica")
thing <- rep(16, x = NA)
thing[9] <- "*"
lost <- c(0, subset(df_sorted$Value, df_sorted$sc.region=="chr3_b"))
plotBranchbyTrait(tree, lost, mode="edge", y.lim = c(0,10), xlims = c(min(lost),max(lost)), show.tip.label = FALSE, palette = colorRampPalette(col.number, bias = 1), outline = TRUE, edge.width = 10, cex= 1, font = 1, legend=TRUE)
edgelabels(lost, bg = "white", frame="rect", cex = 2)
edgelabels(thing ,adj = c(0.5,-0.5), frame="none", cex = 4)
title(main="chr3_b (n = 29)", line = -3, cex.main = 2)
dev.off()

## Chr chr3_c
pdf("plots/chr3_c_geneloss.pdf", width = 6, height = 8)
par(family = "Helvetica")
thing <- rep(16, x = NA)
thing[10] <- "*"
lost <- c(0, subset(df_sorted$Value, df_sorted$sc.region=="chr3_c"))
plotBranchbyTrait(tree, lost, mode="edge", y.lim = c(0,10), xlims = c(min(lost),max(lost)), show.tip.label = FALSE, palette = colorRampPalette(col.number, bias = 1), outline = TRUE, edge.width = 10, cex= 1, font = 1, legend=TRUE)
edgelabels(lost, bg = "white", frame="rect", cex = 2)
edgelabels(thing ,adj = c(0.5,-0.5), frame="none", cex = 4)
title(main="chr3_c (n = 479)", line = -3, cex.main = 2)
dev.off()

## Chr chr5
pdf("plots/chr5_geneloss.pdf", width = 6, height = 8)
par(family = "Helvetica")
thing <- rep(16, x = NA)
thing[10] <- "*"
lost <- c(0, subset(df_sorted$Value, df_sorted$sc.region=="chr5"))
plotBranchbyTrait(tree, lost, mode="edge", y.lim = c(0,10), xlims = c(min(lost),max(lost)), show.tip.label = FALSE, palette = colorRampPalette(col.number, bias = 1), outline = TRUE, edge.width = 10, cex= 1, font = 1, legend=TRUE)
edgelabels(lost, bg = "white", frame="rect", cex = 2)
edgelabels(thing ,adj = c(0.5,-0.5), frame="none", cex = 4)
title(main="chr5 (n = 345)", line = -3, cex.main = 2)
dev.off()

## Chr chr4
pdf("plots/chr4_geneloss.pdf", width = 6, height = 8)
par(family = "Helvetica")
thing <- rep(16, x = NA)
thing[4] <- "*"
lost <- c(0, subset(df_sorted$Value, df_sorted$sc.region=="chr4"))
plotBranchbyTrait(tree, lost, mode="edge", y.lim = c(0,10), xlims = c(min(lost),max(lost)), show.tip.label = FALSE, palette = colorRampPalette(col.number, bias = 1), outline = TRUE, edge.width = 10, cex= 1, font = 1, legend=TRUE)
edgelabels(lost, bg = "white", frame="rect", cex = 2)
edgelabels(thing ,adj = c(0.5,-0.5), frame="none", cex = 4)
title(main="chr4 (n = 262)", line = -3, cex.main = 2)
dev.off()

## Chr chr8
pdf("plots/chr8_geneloss.pdf", width = 6, height = 8)
par(family = "Helvetica")
thing <- rep(16, x = NA)
thing[7] <- "*"
lost <- c(0, subset(df_sorted$Value, df_sorted$sc.region=="chr8"))
plotBranchbyTrait(tree, lost, mode="edge", y.lim = c(0,10), xlims = c(min(lost),max(lost)), show.tip.label = FALSE, palette = colorRampPalette(col.number, bias = 1), outline = TRUE, edge.width = 10, cex= 1, font = 1, legend=TRUE)
edgelabels(lost, bg = "white", frame="rect", cex = 2)
edgelabels(thing ,adj = c(0.5,-0.5), frame="none", cex = 4)
title(main="chr8 (n = 163)", line = -3, cex.main = 2)
dev.off()

# Combine in terminal using pdfjam: 
# pdfjam chrZ_geneloss.pdf chr4A_a_geneloss.pdf chr4A_b_geneloss.pdf chr3_* chr5_geneloss.pdf chr4_geneloss.pdf chr8_geneloss.pdf --nup 3x3 --outfile test.pdf


chrZ <- rep(16, x = NA)
chrZ[2] <- chrZFusionAgeLabel

chr4A_a <- rep(16, x = NA)
chr4A_a[2] <- chr4A_aFusionAgeLabel

chr4A_b <- rep(16, x = NA)
chr4A_b[8] <- chr4A_bFusionAgeLabel_n2
chr4A_b[4] <- chr4A_bFusionAgeLabel_n7

chr4 <- rep(16, x = NA)
chr4[4] <- chr4FusionAgeLabel

chr8 <- rep(16, x = NA)
chr8[7] <- chr8FusionAgeLabel

chr5 <- rep(16, x = NA)
chr5[10] <- chr5FusionAgeLabel

chr3_a <- rep(16, x = NA)
chr3_a[8] <- chr3_aFusionAgeLabel

chr3_b <- rep(16, x = NA)
chr3_b[9] <- chr3_bFusionAgeLabel

chr3_c <- rep(16, x = NA)
chr3_c[10] <- chr3_cFusionAgeLabel

agesNodes <- c(node0Age,node1AgeLabel, node6AgeLabel, node7AgeLabel, node2AgeLabel, node3AgeLabel, node5AgeLabel, node4AgeLabel)

node2tip <- rep(16, x = 0)
node2tip[16] <- 1 # n2_Pbia
node2tip[13] <- 2 # n4_Ccin
node2tip[15] <- 2 # n4_Ccin
node2tip[12] <- 3 # n5_Aarv
node2tip[6] <- 4 # n7_Cjun

show_col(pal_jco("default")(10))
colours_plot <- pal_jco("default")(10)[1:5]

pdf("plots/rate_dndsHI_node2tip_tree.pdf", height = 8, width = 8)
plotBranchbyTrait(tree, node2tip, mode="edge", y.lim = c(0,10), show.tip.label = TRUE, label.offset = 0.01, palette = colorRampPalette(colours_plot), outline = TRUE, edge.width = 5, cex= 1, font = 1, legend=FALSE)
edgelabels(chrZ, adj = c(0.5,-3), frame="none", cex = 1)
edgelabels(chr4A_a, adj = c(0.5,-1), frame="none", cex = 1)
edgelabels(bg = "gray", pch = 21, frame="none", cex = 2)
edgelabels(chr4A_b, adj = c(0.5,-3), frame="none", cex = 1)
edgelabels(chr4, adj = c(0.5,-0.5), frame="none", cex = 1)
edgelabels(chr3_a, adj = c(0.5,-1), frame="none", cex = 1)
edgelabels(chr3_b, adj = c(0.5,-1), frame="none", cex = 1)
edgelabels(chr3_c, adj = c(0.5,-3), frame="none", cex = 1)
edgelabels(chr5, adj = c(0.5,-1), frame="none", cex = 1)
edgelabels(chr8, adj = c(0.5,-1), frame="none", cex = 1)
#node_ages <- c("n0: 25.31","n1: 21.5","n2: 17.8","n3: 8.9","n5: 4.3","n4: 5.8","n6: 18.3","n7: 9.2")
node_ages <- c("n0: 25.31","n1: 21.5","n6: 18.3","n7: 9.2","n2: 17.8","n3: 8.9","n5: 4.3","n4: 5.8")
nodelabels(node_ages, bg = "lightgray", cex=0.8)
dev.off()

# ====================================================================== # 

# ========================== HI and dnds  ============================== # 

# ====================================================================== # 


### Read in HI scores
HI_scores <- read.delim("HI_Predictions_Version3_edit.bed", header = F)
HI_scores <- subset(HI_scores, select = c(V4, V6))
HI_scores <- plyr::rename(HI_scores, c("V4"="Gene_name","V6"="HI_score"))

### Merge HI scores for human genes with one2one zebra finch orthologs
zf.human.orthologs <- read.table("zf_human_orthologs.txt", sep="\t", header = T)
zf.human.orthologs.one2one <- subset(zf.human.orthologs, zf.human.orthologs$Human.homology.type =="ortholog_one2one")
zf.HI.scores <- merge(zf.human.orthologs.one2one, HI_scores, by.x="Human.gene.name", by.y="Gene_name")
zf.HI.scores <- subset(zf.HI.scores, select = c(Transcript.stable.ID, HI_score))
zf.HI.scores <- plyr::rename(zf.HI.scores, c("Transcript.stable.ID"="Trans"))
length(zf.HI.scores$Trans) # 9819

### Read in dnds values
bird.dnds <- read.delim("bird_dnds_to_zebra_finch.txt", sep="\t", header = T, na.strings=c("","NA"))
bird.dnds <- unique(subset(bird.dnds, Chicken.homology.type=="ortholog_one2one"))
bird.dnds <- subset(bird.dnds, dS.with.Chicken<=3)

bird.dnds$dnds <- bird.dnds$dN.with.Chicken/bird.dnds$dS.with.Chicken
bird.dnds <- subset(bird.dnds, select = c(Transcript.stable.ID, dnds))
bird.dnds <- plyr::rename(bird.dnds, c("Transcript.stable.ID"="Trans"))
bird.dnds <- subset(bird.dnds, !is.na(bird.dnds$dnds))
length(bird.dnds$Trans) # 9838

### Merge datasets
dnds.HI <- merge(bird.dnds, zf.HI.scores, by = "Trans", all = TRUE)
long_df.dnds.HI <- merge(long_df, dnds.HI, by = "Trans", all = TRUE)
long_df.dnds.HI <- merge(long_df.dnds.HI, df_sc, by = c("sc.region", "species"))

length(unique(long_df.dnds.HI$Trans)) # 1931

options(scipen=999)



# Extracting genes present or lost from n2 to Panurus biarmicus 
sp = "Panurus_biarmicus"
age = node2Age*-1
selectNode = "node2"

dnds.HI.1species <- subset(long_df.dnds.HI, c(long_df.dnds.HI$species== sp & long_df.dnds.HI$sex.linkage=="sex-linked"))
dnds.HI.1species$value <- as.numeric(dnds.HI.1species$value)
matching_rows <- which(dnds.HI.1species[selectNode] == 1)
dnds.HI.1species.selectedNode <- dnds.HI.1species[matching_rows,]

PanBia.data <- dnds.HI.1species.selectedNode
PanBia.data.dnds <- PanBia.data[!is.na(PanBia.data$dnds), ]
PanBia.data.HI <- PanBia.data[!is.na(PanBia.data$HI_score), ]

# How many genes are in the different tests?
length(PanBia.data$Trans) # 207
length(subset(PanBia.data$Trans, PanBia.data$value==1)) # 190
length(PanBia.data.dnds$Trans) # 157
length(PanBia.data.HI$Trans) # 144


summary(glm(as.numeric(value) ~ dnds, data=PanBia.data.dnds, family="binomial"))
summary(glm(as.numeric(value) ~ HI_score, data=PanBia.data.HI, family="binomial"))


dnds.HI.1species.selectedNode.count <- dnds.HI.1species.selectedNode %>% count(value) 
dnds.HI.1species.selectedNode.count <- pivot_wider(dnds.HI.1species.selectedNode.count, names_from = value, values_from = n) %>% rename(now = "0", orig = "1") 
dnds.HI.1species.selectedNode.count$now <- dnds.HI.1species.selectedNode.count$orig - dnds.HI.1species.selectedNode.count$now
dnds.HI.1species.selectedNode.count <- pivot_longer(dnds.HI.1species.selectedNode.count, everything(), names_to = "Category", values_to = "Value")
dnds.HI.1species.selectedNode.count$Age <- NA
dnds.HI.1species.selectedNode.count$Age[1] <- 0
dnds.HI.1species.selectedNode.count$Age[2] <- age
dnds.HI.1species.selectedNode.count$species <- sp
dnds.HI.1species.selectedNode.count$per_cent <- (dnds.HI.1species.selectedNode.count$Value / dnds.HI.1species.selectedNode.count$Value[2]) * 100
PanBia.count.data <- dnds.HI.1species.selectedNode.count


sp = "Alauda_arvensis"
age = node5Age*-1
selectNode = "node5"

dnds.HI.1species <- subset(long_df.dnds.HI, c(long_df.dnds.HI$species== sp & long_df.dnds.HI$sex.linkage=="sex-linked"))
dnds.HI.1species$value <- as.numeric(dnds.HI.1species$value)
matching_rows <- which(dnds.HI.1species[selectNode] == 1)
dnds.HI.1species.selectedNode <- dnds.HI.1species[matching_rows,]

AlaArv.data <- dnds.HI.1species.selectedNode
AlaArv.data.dnds <- AlaArv.data[!is.na(AlaArv.data$dnds), ]
AlaArv.data.HI <- AlaArv.data[!is.na(AlaArv.data$HI_score), ]

# How many genes are in the different tests?
length(AlaArv.data$Trans) # 834
length(subset(AlaArv.data$Trans, AlaArv.data$value==1)) # 810
length(AlaArv.data.dnds$Trans) # 683
length(AlaArv.data.HI$Trans) # 636


summary(glm(as.numeric(value) ~ dnds, data=AlaArv.data.dnds, family="binomial"))
summary(glm(as.numeric(value) ~ HI_score, data=AlaArv.data.HI, family="binomial"))


dnds.HI.1species.selectedNode.count <- dnds.HI.1species.selectedNode %>% count(value) 
dnds.HI.1species.selectedNode.count <- pivot_wider(dnds.HI.1species.selectedNode.count, names_from = value, values_from = n) %>% rename(now = "0", orig = "1") 
dnds.HI.1species.selectedNode.count$now <- dnds.HI.1species.selectedNode.count$orig - dnds.HI.1species.selectedNode.count$now
dnds.HI.1species.selectedNode.count <- pivot_longer(dnds.HI.1species.selectedNode.count, everything(), names_to = "Category", values_to = "Value")
dnds.HI.1species.selectedNode.count$Age <- NA
dnds.HI.1species.selectedNode.count$Age[1] <- 0
dnds.HI.1species.selectedNode.count$Age[2] <- age
dnds.HI.1species.selectedNode.count$species <- sp
dnds.HI.1species.selectedNode.count$per_cent <- (dnds.HI.1species.selectedNode.count$Value / dnds.HI.1species.selectedNode.count$Value[2]) * 100
AlaArv.count.data <- dnds.HI.1species.selectedNode.count


sp = "Calandrella_cinerea"
age = node4Age*-1
selectNode = "node4"

sp = "Calandrella_cinerea"
age = node3Age*-1
selectNode = "node3"

dnds.HI.1species <- subset(long_df.dnds.HI, c(long_df.dnds.HI$species== sp & long_df.dnds.HI$sex.linkage=="sex-linked"))
dnds.HI.1species$value <- as.numeric(dnds.HI.1species$value)
matching_rows <- which(dnds.HI.1species[selectNode] == 1)
dnds.HI.1species.selectedNode <- dnds.HI.1species[matching_rows,]

CalCin.data <- dnds.HI.1species.selectedNode
CalCin.data.dnds <- CalCin.data[!is.na(CalCin.data$dnds), ]
CalCin.data.HI <- CalCin.data[!is.na(CalCin.data$HI_score), ]

# How many genes are in the different tests?
length(CalCin.data$Trans) # 153
length(subset(CalCin.data$Trans, CalCin.data$value==1)) # 134
length(CalCin.data.dnds$Trans) # 118
length(CalCin.data.HI$Trans) # 108

summary(glm(as.numeric(value) ~ dnds, data=CalCin.data.dnds, family="binomial"))
summary(glm(as.numeric(value) ~ HI_score, data=CalCin.data.HI, family="binomial"))


dnds.HI.1species.selectedNode.count <- dnds.HI.1species.selectedNode %>% count(value) 
dnds.HI.1species.selectedNode.count <- pivot_wider(dnds.HI.1species.selectedNode.count, names_from = value, values_from = n) %>% rename(now = "0", orig = "1") 
dnds.HI.1species.selectedNode.count$now <- dnds.HI.1species.selectedNode.count$orig - dnds.HI.1species.selectedNode.count$now
dnds.HI.1species.selectedNode.count <- pivot_longer(dnds.HI.1species.selectedNode.count, everything(), names_to = "Category", values_to = "Value")
dnds.HI.1species.selectedNode.count$Age <- NA
dnds.HI.1species.selectedNode.count$Age[1] <- 0
dnds.HI.1species.selectedNode.count$Age[2] <- age
dnds.HI.1species.selectedNode.count$species <- sp
dnds.HI.1species.selectedNode.count$per_cent <- (dnds.HI.1species.selectedNode.count$Value / dnds.HI.1species.selectedNode.count$Value[2]) * 100
CalCin.count.data <- dnds.HI.1species.selectedNode.count


sp = "Cisticola_juncidis"
age = node7Age*-1
selectNode = "node7"

dnds.HI.1species <- subset(long_df.dnds.HI, c(long_df.dnds.HI$species== sp & long_df.dnds.HI$sex.linkage=="sex-linked"))
dnds.HI.1species$value <- as.numeric(dnds.HI.1species$value)
matching_rows <- which(dnds.HI.1species[selectNode] == 1)
dnds.HI.1species.selectedNode <- dnds.HI.1species[matching_rows,]

CisJun.data <- dnds.HI.1species.selectedNode
CisJun.data.dnds <- CisJun.data[!is.na(CisJun.data$dnds), ]
CisJun.data.HI <- CisJun.data[!is.na(CisJun.data$HI_score), ]

# How many genes are in the different tests?
length(CisJun.data$Trans) # 255
length(subset(CisJun.data$Trans, CisJun.data$value==1)) # 215
length(CisJun.data.dnds$Trans) # 185
length(CisJun.data.HI$Trans) # 179


summary(glm(as.numeric(value) ~ dnds, data=CisJun.data.dnds, family="binomial"))
summary(glm(as.numeric(value) ~ HI_score, data=CisJun.data.HI, family="binomial"))



dnds.HI.1species.selectedNode.count <- dnds.HI.1species.selectedNode %>% count(value) 
dnds.HI.1species.selectedNode.count <- pivot_wider(dnds.HI.1species.selectedNode.count, names_from = value, values_from = n) %>% rename(now = "0", orig = "1") 
dnds.HI.1species.selectedNode.count$now <- dnds.HI.1species.selectedNode.count$orig - dnds.HI.1species.selectedNode.count$now
dnds.HI.1species.selectedNode.count <- pivot_longer(dnds.HI.1species.selectedNode.count, everything(), names_to = "Category", values_to = "Value")
dnds.HI.1species.selectedNode.count$Age <- NA
dnds.HI.1species.selectedNode.count$Age[1] <- 0
dnds.HI.1species.selectedNode.count$Age[2] <- age
dnds.HI.1species.selectedNode.count$species <- sp
dnds.HI.1species.selectedNode.count$per_cent <- (dnds.HI.1species.selectedNode.count$Value / dnds.HI.1species.selectedNode.count$Value[2]) * 100
CisJun.count.data <- dnds.HI.1species.selectedNode.count


### Combine datasets
HI.dnds.rate.data <- rbind(PanBia.data, AlaArv.data, CalCin.data, CisJun.data)
count.data <- rbind(PanBia.count.data, AlaArv.count.data, CalCin.count.data, CisJun.count.data)
HI.dnds.rate.data <- subset(HI.dnds.rate.data, select = c(species, value, dnds, HI_score))
HI.dnds.rate.data$species <- factor(HI.dnds.rate.data$species)
HI.dnds.rate.data.dnds <- HI.dnds.rate.data[!is.na(HI.dnds.rate.data$dnds), ] # Make datasets without dnds NA's
HI.dnds.rate.data.HI <- HI.dnds.rate.data[!is.na(HI.dnds.rate.data$HI_score), ] # Make datasets without HI_score NA's

### How many genes present and lost in each species?
HI.dnds.rate.data %>% count(species, value) 
#species value   n
#1     Alauda_arvensis     0  24
#2     Alauda_arvensis     1 810
#3 Calandrella_cinerea     0  19
#4 Calandrella_cinerea     1 134
#5  Cisticola_juncidis     0  40
#6  Cisticola_juncidis     1 215
#7   Panurus_biarmicus     0  10
#8   Panurus_biarmicus     1 143

### How many genes present and lost in each species? - dnds dataset
HI.dnds.rate.data.dnds %>% count(species, value) 
#species value   n
#1     Alauda_arvensis     0  17
#2     Alauda_arvensis     1 666
#3 Calandrella_cinerea     0  15
#4 Calandrella_cinerea     1 103
#5  Cisticola_juncidis     0  29
#6  Cisticola_juncidis     1 156
#7   Panurus_biarmicus     0   9
#8   Panurus_biarmicus     1 109

### How many genes present and lost in each species? - HI dataset
HI.dnds.rate.data.HI %>% count(species, value) 
#species value   n
#1     Alauda_arvensis     0  16
#2     Alauda_arvensis     1 620
#3 Calandrella_cinerea     0  13
#4 Calandrella_cinerea     1  95
#5  Cisticola_juncidis     0  27
#6  Cisticola_juncidis     1 152
#7   Panurus_biarmicus     0   7
#8   Panurus_biarmicus     1  99


### Linear models per species

#### Pbia
summary(glm(as.numeric(value) ~ dnds, data=subset(HI.dnds.rate.data.dnds, HI.dnds.rate.data.dnds$species=="Panurus_biarmicus"), family="binomial"))
# dnds         -8.0886     4.0556  -1.994      0.0461 *  

summary(glm(as.numeric(value) ~ HI_score, data=subset(HI.dnds.rate.data.HI, HI.dnds.rate.data.HI$species=="Panurus_biarmicus"), family="binomial"))
# HI_score    -0.04588    0.01806  -2.540     0.0111 *  

#### Aarv
summary(glm(as.numeric(value) ~ dnds, data=subset(HI.dnds.rate.data.dnds, HI.dnds.rate.data.dnds$species=="Alauda_arvensis"), family="binomial"))
# dnds         -4.8306     1.6362  -2.952              0.00315 ** 

summary(glm(as.numeric(value) ~ HI_score, data=subset(HI.dnds.rate.data.HI, HI.dnds.rate.data.HI$species=="Alauda_arvensis"), family="binomial"))
# HI_score    -0.03550    0.00973  -3.648             0.000264 ***

#### Cjun
summary(glm(as.numeric(value) ~ dnds, data=subset(HI.dnds.rate.data.dnds, HI.dnds.rate.data.dnds$species=="Cisticola_juncidis"), family="binomial"))
# dnds         -4.4581     1.3833  -3.223           0.00127 ** 

summary(glm(as.numeric(value) ~ HI_score, data=subset(HI.dnds.rate.data.HI, HI.dnds.rate.data.HI$species=="Cisticola_juncidis"), family="binomial"))
# HI_score    -0.04383    0.01075  -4.078 0.0000454415429 ***

#### Ccal
summary(glm(as.numeric(value) ~ dnds, data=subset(HI.dnds.rate.data.dnds, HI.dnds.rate.data.dnds$species=="Calandrella_cinerea"), family="binomial"))
# dnds         -6.3956     3.0098  -2.125       0.0336 *  

summary(glm(as.numeric(value) ~ HI_score, data=subset(HI.dnds.rate.data.HI, HI.dnds.rate.data.HI$species=="Calandrella_cinerea"), family="binomial"))
# HI_score    -0.06073    0.01732  -3.506    0.000454 ***


# Test all species together
summary(glmer(as.numeric(value) ~ dnds + (1 + dnds | species), data=HI.dnds.rate.data.dnds, family="binomial"))
#summary(glm(as.numeric(value) ~ dnds, data=HI.dnds.rate.data, family="binomial"))

summary(glmer(as.numeric(value) ~ HI_score + (1 + HI_score | species), data=HI.dnds.rate.data.HI, family="binomial"))
#summary(glm(as.numeric(value) ~ HI_score, data=HI.dnds.rate.data.HInoNA, family="binomial"))


#============= Manual things (update) SEPT 2023 =============== # 

# % genes lost /Myr per species
summary(lm(per_cent ~ Age, subset(count.data, count.data$species=="Panurus_biarmicus")))
summary(lm(per_cent ~ Age, subset(count.data, count.data$species=="Alauda_arvensis")))
summary(lm(per_cent ~ Age, subset(count.data, count.data$species=="Calandrella_cinerea")))
summary(lm(per_cent ~ Age, subset(count.data, count.data$species=="Cisticola_juncidis")))

# % genes lost /Myr
summary(lm(per_cent ~ Age, count.data))
#Age          -0.6820     0.3426   -1.99       0.0937 .  

