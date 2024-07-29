library(ggplot2)
library(plyr)
library(dplyr)
library(doBy)
library(reshape2)
library(psych)
library(ggsci)
library(data.table)

# ============================================== # 
# ============================================== # 
# SETUP: WORKING DIRECTORY AND PLOTTING SETTINGS #
# ============================================== # 
# ============================================== # 

### Edit path to where the results from the snakemake pipeline are
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
  theme_gray(base_family = "Helvetica") + 
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

#setwd("~/work/lost_w/results/")
# Load data tables
nr_species <- 8 
zfExons <- read.table("gene_coord_filenames_trans_exon_intron.bed") # Genome coordinates for zebra finch (ZF) exons and introns
hetSites.data <- read.table("allSp.nrHetSites.ExonIntron.8.out", header = F) # Per-exon heterozygosity data (n = 8 species)
LOFSites.data <- read.table("allSp.nrLOFSites.ExonIntron.8.out", header = F) # Per-exon loss-of-function mutation data (n = 8 species)
gencov.data <- read.table("allSp.genCov.ExonIntron.8.out", header = F) # Per-exon genome coverage data (n = 8 species)
callable_sites <- read.table("callableSites.tsv", header = F) # Per-exon number of callable sites 

# ======================================== # 
#  1.2 Create exon data table to populate  #
# ======================================== # 

# Edit ZF exon info, and add exon numbers 
zfExons <- plyr::rename(zfExons, c("V5"="ExonIntron.type","V2"="ExonIntron.start", "V3"="ExonIntron.end", "V4"="Trans", "V1"="Chr"))
length(zfExons$Trans) # 313334 exons (n = 160397) and introns (n = 152937) 
zfExons <- zfExons[with(zfExons, order(Trans, ExonIntron.start)), ]
zfExons$number <- 1
zfExons <- zfExons %>%
  group_by(Trans) %>%
  mutate(ExonIntron.nr = cumsum(number))
zfExons <- subset(zfExons, select = -c(number) )

# Calculate start and end of gene
exonList.data.TransStartMin <- aggregate(ExonIntron.start ~ Trans, zfExons, function(x) min(x))
exonList.data.TransEndMax <- aggregate(ExonIntron.end ~ Trans, zfExons, function(x) max(x))
exonStartEnd <- merge(exonList.data.TransStartMin, exonList.data.TransEndMax , by = ("Trans"))
exonStartEnd <- plyr::rename(exonStartEnd, c("ExonIntron.start"="Gene.start", "ExonIntron.end"="Gene.end"))
zfExons <- merge(zfExons, exonStartEnd , by = ("Trans"))

# Make list with all exons/introns, one observation per sex
exonList.base <- zfExons
exonList.base$sex <- "male"
exonList.base2 <- zfExons
exonList.base2$sex <- "female"
exonList.data <- rbind(exonList.base, exonList.base2) 
length(exonList.data$Trans) # Length: 626668
length(exonList.data$Trans)/ 2 # Length: 313334

# Expand data table to one row per exon, sex and species (86930 is the length of the dataset)
exonList.data.8 <- exonList.data[rep(seq_len(nrow(exonList.data)), each = nr_species), ]
exonList.data.8$species <- rep(c("Alauda_arvensis", "Alauda_razae", "Calandrella_cinerea", "Camaroptera_brevicaudata", 
                                 "Cisticola_juncidis", "Eremophila_alpestris", "Panurus_biarmicus", "Sylvietta_brachyura"), length(exonList.data$Trans))

length(exonList.data.8$Trans)/(nr_species*2) # Nr of exons: 313334

# ===================================== # 
# 1.3 Populate with heterozygosity data #
# ===================================== # 

# Rename columns
hetSites.data <- plyr::rename(hetSites.data, c("V5"="ExonIntron.type","V2"="ExonIntron.start", "V3"="ExonIntron.end", "V4"="Trans", "V1"="Chr",
                                               "V6"="hetSites", "V7"="sex", "V8"="species"))

# Merge heterozygosity data with ZF exon data (output from section 1.1)
exonList.data.8.hetSites <- merge(exonList.data.8, hetSites.data, by=c("Chr", "Trans", "ExonIntron.start", "ExonIntron.end", "ExonIntron.type", "sex", "species"), all = TRUE)
length(exonList.data.8.hetSites$Trans)/(nr_species*2) # Nr of exons: 313334

# ==================================== # 
# 1.4 Merge with loss-of-function data #
# ==================================== # 

# Rename columns
LOFSites.data <- plyr::rename(LOFSites.data, c("V5"="ExonIntron.type","V2"="ExonIntron.start", "V3"="ExonIntron.end", "V4"="Trans", "V1"="Chr",
                                               "V6"="LOFSites", "V7"="sex", "V8"="species"))

# Merge loss-of-function mutation data with ZF+heterozygosity exon data (output from section 1.2)
exonList.data.8.hetSites.LOFSites <- merge(exonList.data.8.hetSites, LOFSites.data, by=c("Chr", "Trans", "ExonIntron.start", "ExonIntron.end", "ExonIntron.type", "sex", "species"), all = TRUE)
length(exonList.data.8.hetSites.LOFSites$Trans)/(nr_species*2) # Length: 313334

# =================================== # 
# 1.4 Merge with genome coverage data #
# =================================== # 

# Rename columns
gencov.data <- plyr::rename(gencov.data, c("V5"="ExonIntron.length","V3"="ExonIntron.start", "V4"="ExonIntron.end", "V1"="Trans", "V2"="Chr",
                                           "V6"="sex", "V7"="minCov", "V8"="maxCov", "V9"="avg_cov", "V10"="median_cov", "V11"="nocoveragebp", "V12"="percentcov", "V13"="species"))

# Merge genome coverage data with ZF+heterozygosity+LOF exon data (output from section 1.3)
exonList.data.8.hetSites.LOFSites.genCov <- merge(exonList.data.8.hetSites.LOFSites, gencov.data, by=c("Chr", "Trans", "ExonIntron.start", "ExonIntron.end", "sex", "species"), all = TRUE)
length(exonList.data.8.hetSites.LOFSites.genCov$Trans)/(nr_species*2) # Length: 313334


# ================================== # 
# 1.5 Merge with callable sites data #
# ================================== # 

callable_sites <- plyr::rename(callable_sites, c("V2"="ExonIntron.start", "V3"="ExonIntron.end", "V4"="Trans", "V1"="Chr",
                                                 "V5"="callableSites", "V7"="sex", "V6"="species"))

callable_sites.df <- merge(exonList.data.8.hetSites.LOFSites.genCov, callable_sites, by=c("Chr", "Trans", "ExonIntron.start", "ExonIntron.end", "sex", "species"), all = TRUE)
length(callable_sites.df$Trans)/(nr_species*2) # Length: 313334

exonList.data.8.hetSites.LOFSites.genCov <- callable_sites.df

length(exonList.data.8.hetSites.LOFSites.genCov$Trans)/(nr_species*2) # Length: 313334

hist(callable_sites.df$callableSites, breaks = 10000, xlim=c(0,1000))

hist(exonList.data.8.hetSites.LOFSites.genCov$ExonIntron.length, breaks = 10000, xlim=c(0,1000))


# ====================================== # 
# 1.6 Long-to-wide format transformation #
# == (female and male on the same row) = #
# ====================================== # 

data.female <- subset(exonList.data.8.hetSites.LOFSites.genCov, exonList.data.8.hetSites.LOFSites.genCov$sex=="female")
data.male <- subset(exonList.data.8.hetSites.LOFSites.genCov, exonList.data.8.hetSites.LOFSites.genCov$sex=="male")
data.female <- plyr::rename(data.female, c("hetSites"="female.hetSites", "LOFSites"="female.LOFSites", "minCov"="female.mincov","maxCov"="female.maxcov", 
                                           "avg_cov"="female.avg_cov", "median_cov"="female.median_cov", "nocoveragebp"="female.nocoveragebp", "percentcov"="female.percentcov", "callableSites"="female.callableSites"))
data.male <- plyr::rename(data.male, c("hetSites"="male.hetSites", "LOFSites"="male.LOFSites", "minCov"="male.mincov","maxCov"="male.maxcov", 
                                       "avg_cov"="male.avg_cov", "median_cov"="male.median_cov", "nocoveragebp"="male.nocoveragebp", "percentcov"="male.percentcov", "callableSites"="male.callableSites"))

data.female <- subset(data.female, select = -c(sex) )
data.male <- subset(data.male, select = -c(sex) )

allData.wide <- merge(data.female, data.male, by=c("Trans", "ExonIntron.type", "species", "Chr", "ExonIntron.start", "ExonIntron.end", "Gene.start", "Gene.end", "ExonIntron.length", "ExonIntron.nr"))

length(allData.wide$Trans)/nr_species # Length: 313334


ggplot(allData.wide, aes(x = ExonIntron.length, y = male.callableSites, fill = species)) +
  geom_point(alpha = 0.1) + scale_fill_jco() +
  text_size_colour + facet_wrap(~species)

long_ones <- subset(allData.wide, allData.wide$ExonIntron.length> 50000)

# ========================================== # 
# ========================================== # 
# SECTION 2 - FILTER AND SUMMARIZE EXON DATA #
# ========================================== # 
# ========================================== # 

# ======================================================== # 
# == 2.1 Normalize genome coverage values between sexes == #
# median values of exons that are autosomal in all species #
# ======================================================== # 

#Calculate the average coverage per sample
allData.wide.test <- allData.wide
fun <- function(x){
  c(mean=mean(x, na.rm=TRUE))} # Mean and standard deviation function

cov_per_sample_autosomes <- summaryBy(female.avg_cov + male.avg_cov ~ species, 
                                      data=subset(allData.wide.test, c(allData.wide.test$Chr!="Z" & 
                                                                       allData.wide.test$Chr!="Z_random" & 
                                                                       allData.wide.test$Chr!="3" & 
                                                                       allData.wide.test$Chr!="4A" & 
                                                                       allData.wide.test$Chr!="4" & 
                                                                       allData.wide.test$Chr!="5" & 
                                                                       allData.wide.test$Chr!="8")), keep.names=TRUE, FUN=mean)

cov_per_sample_autosomes <- plyr::rename(cov_per_sample_autosomes, c("female.avg_cov"="female.genomeWide.avg_cov","male.avg_cov"="male.genomeWide.avg_cov"))
allData.wide.test <- merge(allData.wide.test, cov_per_sample_autosomes, by=c("species"))


# First, replace outlier genome coverage values with NA
allData.wide.test$female.avg_cov.filter <- allData.wide.test$female.avg_cov
allData.wide.test$female.avg_cov.filter[which(allData.wide.test$female.avg_cov.filter > (allData.wide.test$female.genomeWide.avg_cov*2))] <- NA
allData.wide.test$female.avg_cov.filter[which(allData.wide.test$female.avg_cov.filter < 5)] <- NA
allData.wide.test$female.avg_cov.filter[which(allData.wide.test$female.callableSites < 50)] <- NA
allData.wide.test$male.avg_cov.filter <- allData.wide.test$male.avg_cov
allData.wide.test$male.avg_cov.filter[which(allData.wide.test$male.avg_cov.filter > (allData.wide.test$male.genomeWide.avg_cov*2))] <- NA
allData.wide.test$male.avg_cov.filter[which(allData.wide.test$male.avg_cov.filter < 5)] <- NA
allData.wide.test$male.avg_cov.filter[which(allData.wide.test$male.callableSites < 50)] <- NA

allData.wide <- allData.wide.test

fun <- function(x){
  c(median=median(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE))} # Median and standard deviation function

# Second, calculate average female and male genome coverage values for each species
autosome.chr.cov_ratio <- summaryBy(female.avg_cov.filter + male.avg_cov.filter ~ species, data=subset(allData.wide, 
                                                                                                       c(allData.wide$Chr!="Z" & allData.wide$Chr!="Z_random" & allData.wide$Chr!="3" & 
                                                                                                           allData.wide$Chr!="4A" & allData.wide$Chr!="4" & allData.wide$Chr!="5" & 
                                                                                                           allData.wide$Chr!="8")), keep.names=TRUE, FUN=c(fun))

# Third, calculate "genome-wide female-to-male ratio"
autosome.chr.cov_ratio$cov_ratio <- autosome.chr.cov_ratio$female.avg_cov.filter.median / autosome.chr.cov_ratio$male.avg_cov.filter.median
allData.wide.scaling <- merge(allData.wide, autosome.chr.cov_ratio, by=c("species"))

# Fourth, scale the female-to-male ratio value for each exon using the "genome-wide female-to-male ratio"
allData.wide.scaling$avg_cov_ratio_scaled <- (allData.wide.scaling$female.avg_cov.filter / allData.wide.scaling$male.avg_cov.filter) / allData.wide.scaling$cov_ratio

# Fifth, replace outlier values (F:M ratio > 2 OR F:M ratio <0.1) with NA
allData.wide.scaling$avg_cov_ratio_scaled[which(allData.wide.scaling$avg_cov_ratio_scaled > 2)] <- NA
allData.wide.scaling$avg_cov_ratio_scaled[which(allData.wide.scaling$avg_cov_ratio_scaled < 0.1)] <- NA

# Density plots for the non-normalized genome coverage values per species and sex 
pdf("plots/nonNorm_genome_coverage_exonIntron.pdf", height = 10, width = 15)
ggplot(allData.wide.scaling, aes(x = female.avg_cov.filter, fill = "female")) +
  geom_density(alpha = 0.5) + scale_fill_jco() + geom_density(aes(x = male.avg_cov.filter, fill = "male"), alpha = 0.5) +
  text_size_colour + xlab("genome coverage") + facet_wrap(~species)
dev.off()

# Density plot show that normalization worked well (coverage values normal distribution around 1 in all species)
pdf("plots/norm_genome_coverage_exonIntron.pdf", height = 5, width = 10)
ggplot(allData.wide.scaling, aes(x = avg_cov_ratio_scaled, fill = species)) +
  geom_density(alpha = 0.1) + scale_fill_jco() +
  text_size_colour + xlab("normalized female/male genome coverage")
dev.off()


# Summary of female and male coverage per species
summaryBy(female.avg_cov.filter + male.avg_cov.filter ~ species, data=subset(allData.wide), keep.names=TRUE, FUN=c(fun))

length(allData.wide.scaling$Trans) # Length: 2506672

# ======================================================= # 
# 2.2 Calculate proportion of heterozygous sites per exon #
# === (divide heterozygous sites with length of exon) === #
# ======================================================= # 

allData.wide.scaling$female.hetSites.prop <- allData.wide.scaling$female.hetSites / allData.wide.scaling$female.callableSites
allData.wide.scaling$male.hetSites.prop <- allData.wide.scaling$male.hetSites / allData.wide.scaling$male.callableSites
allData.wide.scaling$diff.hetSites.prop <- allData.wide.scaling$female.hetSites.prop - allData.wide.scaling$male.hetSites.prop

length(allData.wide.scaling$Trans) # Length: 2506672

# ========================================================= # 
# 2.3 Summarize data per gene (mean and standard deviation) #
# ========================================================= # 

fun <- function(x){
  c(mean=mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE))}

perGeneData <- summaryBy(avg_cov_ratio_scaled + diff.hetSites.prop ~ species + Trans + Chr , data=allData.wide.scaling, keep.names=TRUE, FUN=c(fun))

perGeneData <- plyr::rename(perGeneData, c("avg_cov_ratio_scaled.mean"="Gene.cov_ratio_scaled", "avg_cov_ratio_scaled.sd"="Gene.cov_ratio_scaled.sd", 
                                           "diff.hetSites.prop.mean"="Gene.diff.hetSites.prop.mean", "diff.hetSites.prop.sd"="Gene.diff.hetSites.prop.sd"))

length(perGeneData$Trans) # Length: 148944

# ============================================= # 
# 2.4 Add loss-of-function mutation information #
# ============================================= # 

# Add number of loss-of-function mutations per gene
perGeneDataLOF <- summaryBy(female.LOFSites + male.LOFSites + ExonIntron.length ~ species + Trans + Chr + Gene.start + Gene.end, data=subset(allData.wide.scaling, allData.wide.scaling$ExonIntron.type=="exon"), keep.names=TRUE, FUN=c(sum))
length(perGeneDataLOF$Trans) # Length: 148944

perGeneData <- merge(perGeneData, perGeneDataLOF, by=(c("Trans", "species", "Chr")), all = TRUE)
length(perGeneData$Trans) # Length: 148944

length(perGeneData$Trans)/ nr_species # 18618 - all genes still present

# Create data table listing all exons with loss-of-function mutations for each gene (females separately)
allData.wide.scaling.femaleLOF <- subset(allData.wide.scaling, c(allData.wide.scaling$female.LOFSites > 0 & allData.wide.scaling$male.LOFSites == 0 & allData.wide.scaling$ExonIntron.type=="exon"))
allData.wide.scaling.femaleLOF <- as.data.table(allData.wide.scaling.femaleLOF)[, toString(ExonIntron.nr), by = list(species, Trans)]
allData.wide.scaling.femaleLOF <- plyr::rename(allData.wide.scaling.femaleLOF, c("V1"="female.LOFSites.ExonNr"))

# Create data table listing all exons with loss-of-function mutations for each gene (males separately)
allData.wide.scaling.maleLOF <- subset(allData.wide.scaling, c(allData.wide.scaling$male.LOFSites > 0 & allData.wide.scaling$female.LOFSites == 0 & allData.wide.scaling$ExonIntron.type=="exon"))
allData.wide.scaling.maleLOF <- as.data.table(allData.wide.scaling.maleLOF)[, toString(ExonIntron.nr), by = list(species, Trans)]
allData.wide.scaling.maleLOF <- plyr::rename(allData.wide.scaling.maleLOF, c("V1"="male.LOFSites.ExonNr"))

# Create data table listing all exons with loss-of-function mutations for each gene (females and males together)
allData.wide.scaling.bothLOF <- subset(allData.wide.scaling, c(allData.wide.scaling$female.LOFSites > 0 & allData.wide.scaling$male.LOFSites > 0  & allData.wide.scaling$ExonIntron.type=="exon"))
allData.wide.scaling.bothLOF <- as.data.table(allData.wide.scaling.bothLOF)[, toString(ExonIntron.nr), by = list(species, Trans)]
allData.wide.scaling.bothLOF <- plyr::rename(allData.wide.scaling.bothLOF, c("V1"="bothSexes.LOFSites.ExonNr"))

# Merge the three data tables above into one
LOF.exonNr <- merge(allData.wide.scaling.femaleLOF, allData.wide.scaling.maleLOF, by=(c("Trans", "species")), all = TRUE)
LOF.exonNr <- merge(LOF.exonNr, allData.wide.scaling.bothLOF, by=(c("Trans", "species")), all = TRUE)
LOF.exonNr$female.LOFSites.ExonNr[which(c(!is.na(LOF.exonNr$bothSexes.LOFSites.ExonNr)))] <- LOF.exonNr$bothSexes.LOFSites.ExonNr[which(c(!is.na(LOF.exonNr$bothSexes.LOFSites.ExonNr)))]
LOF.exonNr$male.LOFSites.ExonNr[which(c(!is.na(LOF.exonNr$bothSexes.LOFSites.ExonNr)))] <- LOF.exonNr$bothSexes.LOFSites.ExonNr[which(c(!is.na(LOF.exonNr$bothSexes.LOFSites.ExonNr)))]

# Merge with per-gene data table (output from section 2.3)
perGeneData <- merge(perGeneData, LOF.exonNr, by=(c("Trans", "species")), all = TRUE)
length(perGeneData$Trans) # Length: 148944

# Categorize genes into having (a) female-only LOF mutations, (b) male-only LOF mutations or (c) both sexes LOF mutations
perGeneData$LOFstatus <- NA
perGeneData$LOFstatus[which(c(perGeneData$female.LOFSites > 0 & perGeneData$male.LOFSites == 0 ))] <- "female.LOF"
perGeneData$LOFstatus[which(c(perGeneData$female.LOFSites == 0 & perGeneData$male.LOFSites > 0 ))] <- "male.LOF"
perGeneData$LOFstatus[which(c(perGeneData$female.LOFSites > 0 & perGeneData$male.LOFSites > 0 ))] <- "bothSexes.LOF"

# Create data table of exons with loss-of-function mutations
LOFexons <- subset(allData.wide.scaling, c((allData.wide.scaling$female.LOFSites > 0 | allData.wide.scaling$male.LOFSites > 0) & allData.wide.scaling$ExonIntron.type=="exon"))
LOFexons <- as.data.table(LOFexons)[, toString(species), by = list(ExonIntron.nr, Trans)]

# Add information to per-exon data table (output section 2.1)
LOFexons.allData <- merge(allData.wide.scaling, LOFexons, by=(c("Trans", "ExonIntron.nr")), all = TRUE)
length(LOFexons.allData$Trans) # Length: 2506672

# If several species have loss-of-function mutations in the same exon, list them in column "V1", and repeat species name in column "Z"
LOFexons.allData$Z <- mapply(function(x, y) {
  temp <- intersect(x, y)
  if(length(temp)) toString(temp) else ""
}, strsplit(LOFexons.allData$V1, ", "), LOFexons.allData$species)

length(LOFexons.allData$Trans) # Length: 2506672

# Introns (with or without loss-of-function mutations) and exons without loss-of-function mutations should be "NA" in column "V1"
LOFexons.allData$V1[which(LOFexons.allData$species != LOFexons.allData$Z)] <- NA

# Subset data table to exons with loss-of-function mutations
LOFexons.allData.LOF <- subset(LOFexons.allData, c(LOFexons.allData$species == LOFexons.allData$Z)) 
LOFexons.allData.LOF <- subset(LOFexons.allData.LOF, select = c(Trans, Chr, species, ExonIntron.nr, V1, Z) )

# Subset data table to exons where there is a loss-of-function mutation in more than one species
LOFexons.allData.LOF.shared <- subset(LOFexons.allData.LOF, c(LOFexons.allData.LOF$V1 != LOFexons.allData.LOF$Z)) 
LOF.shared <- unique(subset(LOFexons.allData.LOF.shared, select = c(Trans, species, V1) ))
LOF.shared.Gene <- unique(subset(LOFexons.allData.LOF.shared, select = c(Trans, species) ))
LOF.shared.Gene$sharedExons <- "sharedExons"

# Merge with per-gene data table (genes with shared loss-of-function exons are marked "sharedExons")
perGeneData <- merge(perGeneData, LOF.shared.Gene, by=(c("Trans", "species")), all = TRUE)

# Mark genes without shared loss-of-function exons as "notSharedExons")
perGeneData$sharedExons[which(c(! is.na(perGeneData$LOFstatus) & is.na(perGeneData$sharedExons)))] <- "notSharedExons"

length(perGeneData$Trans) # Length: 148944

# ================================================ # 
# 2.5 Add information on sex-linked genome regions #
# ================================================ # 

# Label chromosomes that are autosomal in all species as "autosomes", and the rest by their chromosome ID
perGeneData$group <- "autosomes"
perGeneData$group[which(perGeneData$Chr=="Z")] <- "Z"
perGeneData$group[which(perGeneData$Chr=="Z_random")] <- "Z_random"
perGeneData$group[which(perGeneData$Chr=="4A")] <- "4A"
perGeneData$group[which(perGeneData$Chr=="4")] <- "4"
perGeneData$group[which(perGeneData$Chr=="3")] <- "3"
perGeneData$group[which(perGeneData$Chr=="5")] <- "5"
perGeneData$group[which(perGeneData$Chr=="8")] <- "8"

# Subdivide into different sex-linked regions
perGeneData$sc.region <- "Autosomal"
perGeneData$sc.region[which(c(perGeneData$Chr == "3" & ((perGeneData$Gene.start > 5800000 & perGeneData$Gene.start < 24100000) | (perGeneData$Gene.start > 29800000 & 
                                                                                                                                    perGeneData$Gene.start < 88000000))))] <- "chr3_c"
perGeneData$sc.region[which(c(perGeneData$Chr == "3" & ((perGeneData$Gene.start > 8400000 & perGeneData$Gene.start < 10400000) | (perGeneData$Gene.start > 18100000 & 
                                                                                                                                    perGeneData$Gene.start < 24100000))))] <- "chr3_a"
perGeneData$sc.region[which(c(perGeneData$Chr == "3" & (perGeneData$Gene.start > 10400000 & perGeneData$Gene.start < 14000000)))] <- "chr3_b"
perGeneData$sc.region[which(c(perGeneData$Chr == "5" & perGeneData$Gene.start > 9100000 & perGeneData$Gene.start < 45400000))] <- "chr5"
perGeneData$sc.region[which(c(perGeneData$Chr == "4" & perGeneData$Gene.start > 13800000  & perGeneData$Gene.start < 49100000))] <- "chr4"
perGeneData$sc.region[which(c(perGeneData$Chr == "4A" & perGeneData$Gene.start > 5300000  & perGeneData$Gene.start < 9600000))] <- "chr4A_a"
perGeneData$sc.region[which(c(perGeneData$Chr == "4A" & perGeneData$Gene.start < 5300000))] <- "chr4A_b"
perGeneData$sc.region[which(c(perGeneData$Chr == "8" & perGeneData$Gene.start > 7300000  & perGeneData$Gene.start < 21800000))] <- "chr8"
perGeneData$sc.region[which(perGeneData$Chr == "Z")] <- "chrZ"

# Label each gene by if they are sex-linked or not in each species
perGeneData$sex.linkage <- "autosomal"
perGeneData$sex.linkage[which(perGeneData$sc.region == "chrZ")] <- "sex-linked"
perGeneData$sex.linkage[which(c(perGeneData$sc.region == "chr4A_a"))] <- "sex-linked"
perGeneData$sex.linkage[which(c(perGeneData$sc.region == "chr4" & perGeneData$species == "Cisticola_juncidis"))] <- "sex-linked"
perGeneData$sex.linkage[which(c(perGeneData$sc.region == "chr4" & perGeneData$species == "Camaroptera_brevicaudata"))] <- "sex-linked"
perGeneData$sex.linkage[which(c(perGeneData$sc.region == "chr8" & perGeneData$species == "Sylvietta_brachyura"))] <- "sex-linked"
perGeneData$sex.linkage[which(c(perGeneData$sc.region == "chr5" & perGeneData$species == "Alauda_arvensis"))] <- "sex-linked"
perGeneData$sex.linkage[which(c(perGeneData$sc.region == "chr5" & perGeneData$species == "Alauda_razae"))] <- "sex-linked"
perGeneData$sex.linkage[which(c(perGeneData$sc.region == "chr3_c" & perGeneData$species == "Alauda_arvensis"))] <- "sex-linked"
perGeneData$sex.linkage[which(c(perGeneData$sc.region == "chr3_c" & perGeneData$species == "Alauda_razae"))] <- "sex-linked"
perGeneData$sex.linkage[which(c(perGeneData$sc.region == "chr3_a" & perGeneData$species == "Alauda_arvensis"))] <- "sex-linked"
perGeneData$sex.linkage[which(c(perGeneData$sc.region == "chr3_a" & perGeneData$species == "Alauda_razae"))] <- "sex-linked"
perGeneData$sex.linkage[which(c(perGeneData$sc.region == "chr3_a" & perGeneData$species == "Panurus_biarmicus"))] <- "sex-linked"
perGeneData$sex.linkage[which(c(perGeneData$sc.region == "chr3_a" & perGeneData$species == "Eremophila_alpestris"))] <- "sex-linked"
perGeneData$sex.linkage[which(c(perGeneData$sc.region == "chr3_a" & perGeneData$species == "Calandrella_cinerea"))] <- "sex-linked"
perGeneData$sex.linkage[which(c(perGeneData$sc.region == "chr3_b" & perGeneData$species == "Eremophilis_alpestris"))] <- "sex-linked"
perGeneData$sex.linkage[which(c(perGeneData$sc.region == "chr3_b" & perGeneData$species == "Calandrella_cinerea"))] <- "sex-linked"
perGeneData$sex.linkage[which(c(perGeneData$sc.region == "chr3_b" & perGeneData$species == "Alauda_arvensis"))] <- "sex-linked"
perGeneData$sex.linkage[which(c(perGeneData$sc.region == "chr3_b" & perGeneData$species == "Alauda_razae"))] <- "sex-linked"
perGeneData$sex.linkage[which(c(perGeneData$sc.region == "chr4A_b" & perGeneData$species != "Sylvietta_brachyura"))] <- "sex-linked"

# Remove genes belonging to chr "Un" or those ending with "_random
perGeneData.noRandom <- subset(perGeneData, perGeneData$Chr!="Un")
random <- unique(perGeneData.noRandom$Chr[grep("random", perGeneData.noRandom$Chr)])
perGeneData.noRandom <- perGeneData.noRandom[ ! perGeneData.noRandom$Chr %in% random, ]
length(perGeneData.noRandom$Trans)/8 # 14320 genes left

# =============================================== # 
# 2.6 Calculate outlier values                    #
# (used to categorize genes as sex-linked or not) #
# =============================================== # 

# Calculating IQR from heterozygosity data
IQR_mean <- describeBy(perGeneData.noRandom$Gene.diff.hetSites.prop.mean, perGeneData.noRandom$species, mat = TRUE, IQR=TRUE, quant=c(.25,.75)) 
IQR_mean <- subset(IQR_mean, select=c("group1", "IQR", "Q0.25", "Q0.75"))
IQR_mean$IQR.1.5 <- IQR_mean$IQR*1.5
IQR_mean$lower <- IQR_mean$Q0.25-IQR_mean$IQR.1.5
IQR_mean$upper <- IQR_mean$Q0.75+IQR_mean$IQR.1.5
IQR_mean <- subset(IQR_mean, select=c("group1", "lower", "upper"))
IQR_mean <- plyr::rename(IQR_mean, c("group1"="species", "lower"="hetSites_lower_CI", "upper"="hetSites_upper_CI"))
perGeneData.IQR <- merge(perGeneData.noRandom, IQR_mean, by=c("species")) 
IQR_mean_het <- IQR_mean

# Calculating IQR from genome coverage data (not used)
IQR_mean <- describeBy(perGeneData.noRandom$Gene.cov_ratio_scaled, perGeneData.noRandom$species, mat = TRUE, IQR=TRUE, quant=c(.25,.75)) 
IQR_mean <- subset(IQR_mean, select=c("group1", "IQR", "Q0.25", "Q0.75"))
IQR_mean$IQR.1.5 <- IQR_mean$IQR*1.5
IQR_mean$lower <- IQR_mean$Q0.25-IQR_mean$IQR.1.5
IQR_mean$upper <- IQR_mean$Q0.75+IQR_mean$IQR.1.5
IQR_mean <- subset(IQR_mean, select=c("group1", "lower", "upper"))
IQR_mean <- plyr::rename(IQR_mean, c("group1"="species", "lower"="cov_ratio_lower_CI", "upper"="cov_ratio_upper_CI"))
perGeneData.IQR <- merge(perGeneData.IQR, IQR_mean, by=c("species")) 
IQR_mean_cov <- IQR_mean

# Label heterozygosity outliers
perGeneData.IQR$hetSites_outlier <- "no"
perGeneData.IQR$hetSites_outlier[which(c(perGeneData.IQR$Gene.diff.hetSites.prop.mean > perGeneData.IQR$hetSites_upper_CI))] <- "yes"
perGeneData.IQR$hetSites_outlier[which(c(perGeneData.IQR$Gene.diff.hetSites.prop.mean < perGeneData.IQR$hetSites_lower_CI))] <- "yes"

# Label genome coverage outliers
perGeneData.IQR$cov_ratio_outlier <- "no"
perGeneData.IQR$cov_ratio_outlier[which(c(perGeneData.IQR$Gene.cov_ratio_scaled > perGeneData.IQR$cov_ratio_upper_CI))] <- "yes"
perGeneData.IQR$cov_ratio_outlier[which(c(perGeneData.IQR$Gene.cov_ratio_scaled < perGeneData.IQR$cov_ratio_lower_CI))] <- "yes"

length(perGeneData.IQR$Trans)/ nr_species # 14320 genes left

perGeneData.IQR <- merge(zfPep, perGeneData.IQR, by=c("Trans", "Chr")) # Add ZF gene coordinates

# Remove genes where genome coverage were mas masked (<5x or >75x) 
length(which(is.na(perGeneData.IQR$Gene.cov_ratio_scaled))) # 11949
length(unique(perGeneData.IQR$Trans[which(is.na(perGeneData.IQR$Gene.cov_ratio_scaled))])) # 554
perGeneData.IQR.rm.NAcov <- perGeneData.IQR[!is.na(perGeneData.IQR$Gene.cov_ratio_scaled),]
length(which(is.na(perGeneData.IQR.rm.NAcov$Gene.cov_ratio_scaled))) # 0

# Remove genes where both sexes have loss-of-function mutations
length(which(perGeneData.IQR.rm.NAcov$LOFstatus=="bothSexes.LOF"))
perGeneData.IQR.rm.NAcov[perGeneData.IQR.rm.NAcov$LOFstatus=="bothSexes.LOF",] %>% count(species, sex.linkage)
length(unique(perGeneData.IQR.rm.NAcov$Trans[which(perGeneData.IQR.rm.NAcov$LOFstatus=="bothSexes.LOF")])) # 164
perGeneData.IQR.rm.NAcov[perGeneData.IQR.rm.NAcov$LOFstatus=="bothSexes.LOF",] %>% count(species, sex.linkage)

perGeneData.IQR.rm.NAcov.noLOFbothsexes <- perGeneData.IQR.rm.NAcov %>% filter(is.na(LOFstatus) | LOFstatus != "bothSexes.LOF")
perGeneData.IQR.rm.NAcov.noLOFbothsexes %>% count(species, sex.linkage, LOFstatus)

# Subset for genes longer than 700bp
perGeneData.IQR.noShort <- subset(perGeneData.IQR.rm.NAcov.noLOFbothsexes, c(perGeneData.IQR.rm.NAcov.noLOFbothsexes$ExonIntron.length > 700))
length(unique(perGeneData.IQR.noShort$Trans)) # 8452 genes remaining <- now 8440

#testZ <- subset(perGeneData.IQR.rm.NAcov, c(perGeneData.IQR.rm.NAcov$Chr =="Z" & perGeneData.IQR.rm.NAcov$ExonIntron.length < 1000))


# ============================================ # 
# 2.7 Categorize genes as lost W or survived W #
# ============================================ # 

# "diploid": Genes with F:M coverage ratio >0.75 
diploid_F <- subset(perGeneData.IQR.noShort, c(perGeneData.IQR.noShort$Gene.cov_ratio_scaled>0.75 ))    
diploid_F$lost_status <- "diploid_female"
diploid_F <- subset(diploid_F, select=c(Trans, species, lost_status))

# "lost W" (sex-linked genes with no W gene copy in females): Genes with F:M coverage ratio <0.7 and have little difference in heterozygosity between sexes
haploid_F <- subset(perGeneData.IQR.noShort, c(perGeneData.IQR.noShort$Gene.cov_ratio_scaled<0.7 & perGeneData.IQR.noShort$Gene.diff.hetSites.prop.mean < 0.001))
haploid_F$lost_status <- "haploid_female"
haploid_F <- subset(haploid_F, select=c(Trans, species, lost_status))

# Concatenate data tables
lost_W_status <- rbind(diploid_F, haploid_F )

# Check that no genes belong to more than one category
merge(diploid_F, haploid_F, by=(c("Trans", "species"))) # 0 rows - no overlap

# Merge gene status data table (lost_W_status) with main data table (perGeneData.IQR.noShort)
perGeneData.IQR.noShort <- merge(perGeneData.IQR.noShort, lost_W_status, by=(c("Trans", "species")), all = TRUE)
perGeneData.IQR.noShort$lost_status[which(c(perGeneData.IQR.noShort$LOFstatus == "female.LOF"))] <- "haploid_female"

# Label genes without a gene status category as "unsure"
perGeneData.IQR.noShort$lost_status[which(c(is.na(perGeneData.IQR.noShort$lost_status)))] <- "unsure"

# ======================= # 
# 2.8 Write output table #
# ======================= # 

# Per-gene data. This is the table that is analysed further
outname <- sprintf("perGeneData.IQR.noShort.tsv")
write.table(perGeneData.IQR.noShort, file = outname, sep = "\t", quote = FALSE, row.names = F)

