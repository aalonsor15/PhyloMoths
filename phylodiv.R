# ----------------------------------------------------------------------
# El Verde moths project: Calculating phylogenetic diversity metrics
# 26 May 2022
# Aura M. Alonso-Rodríguez
# ----------------------------------------------------------------------
#

# loading libraries

# install.packages("picante", dependencies = TRUE)
# install.packages("ape")
# install.packages("phytools")

library(picante)
library(tidyverse)
library(vegan)
library(ape)
library(phytools)
library(ggplot2)
library(ggpubr)

# picante package resources #
# https://picante.r-forge.r-project.org/picante-intro.pdf
# http://www2.uaem.mx/r-mirror/web/packages/picante/vignettes/picante-intro.pdf
# https://cran.r-project.org/web/packages/picante/picante.pdf


# practicing picante ------------------------------------------------


# data(phylocom)
# names(phylocom)
# 
# phy <- phylocom$phylo
# comm <- phylocom$sample
# traits <- phylocom$traits
# 
# phy 
# 
# prunedphy <- prune.sample(comm, phy)
# 
# rownames(comm)
# 
# pd.result <- pd(comm, phy, include.root=TRUE)
# pd.result


#SUPER IMPORTANT NOTE ----------------------------------------------------------
# The names on the tree newick file need to be the same as in the matrix.
# I need to remove any - or > or any other symbol from the name, the only symbol
# I should have is the underscore _.
# When I run the tree, I need to make sure the name does not contain any symbols
# except the _. And the names in the matrix should be written exactly the same.
# I can edit the names in all files (newick, fafsa, excel) in Notepad using 
# Find/Replace


# Loading my tree -----------------------------------------------------------

#  How to upload tree data -> If you have a phylogeny in Newick or Nexus 
# format it can be imported into R with the read.tree or read.nexus functions

tree <- read.tree("tree_202210.newick")
tree
class(tree)


# Loading my community data --------------------------------------------------

alldata <- read.csv("matrix_202210.csv")
moths <- select(alldata, PRMAR661_21_Nystalea_superciliosa_A001c_Notodontidae:
                  PRMAR716_21_Glaphyria_dolatalis_M110a_Crambidae)
codes <- select(alldata, Code)
matrix <- cbind(codes, moths)
matrix <- data.frame(matrix[,-1], row.names=matrix[,1])

class(matrix)
matrix <- as.matrix(matrix)
colnames(matrix)
rownames(matrix)


# Calculating phylogenetic diversity metrics ------------------------------------

# We can prune the phylogeny to include only the species that actually 
# occurred in some community.
prunedtree <- prune.sample(matrix, tree)
prunedtree

#### MPD and MNTD metrics ####

# First we need a phylogenetic distance matrix for the input
#A phylo object can be converted to a interspecific phylogenetic distance matrix 
# using the cophenetic function
phydist <- cophenetic(tree) 


# Calculating MPD index
# MPD = Mean Pairwise Distance

mpd <- mpd(matrix, phydist, abundance.weighted = FALSE) 
mpd #this is the observed MPD index

ses.mpd.result <- ses.mpd(matrix, phydist, null.model = "taxa.labels",
    abundance.weighted = FALSE, runs = 99) # eventually, change to 
                                          # runs = 999, iterations = 1000
ses.mpd.result #this is the standardized effect size of the MPD index (SESmpd),  
# which compares the observed MPD with a NULL (randomized) community


# Calculating MNTD index
# MNTD = Mean Nearest Taxon Distance

mntd <- mntd(matrix, phydist, abundance.weighted = FALSE)
mntd #this is the observed MNTD index

ses.mntd.result <- ses.mntd(matrix, phydist, null.model = "taxa.labels",
    abundance.weighted = FALSE, runs = 99)  # eventually, change to 
                                            # runs = 999, iterations = 1000
ses.mntd.result #this is the standardized effect size of the MNTD index (SESmntd),  
# which compares the observed MNTD with a NULL (randomized) community

# The output includes the following columns:
# • ntaxa Number of taxa in community
# • mpd.obs Observed mpd in community
# • mpd.rand.mean Mean mpd in null communities
# • mpd.rand.sd Standard deviation of mpd in null communities
# • mpd.obs.rank Rank of observed mpd vs. null communities
# • mpd.obs.z Standardized effect size of mpd vs. null communities 
#   (equivalent to -NRI)
# • mpd.obs.p P-value (quantile) of observed mpd vs. null communities 
#   (= mpd.obs.rank / runs + 1)
# • runs Number of randomizations

# Positive SES values (mpd.obs.z > 0) and high quantiles (mpd.obs.p > 0.95) indicate 
# phylogenetic evenness, or a greater phylogenetic distance among co-occurring
# species than expected. Negative SES values and low quantiles (mpd.obs.p < 0.05)
# indicate phylogenetic clustering, or small phylogenetic distances among co-occurring
# species than expected. MPD is generally thought to be more sensitive to tree-wide
# patterns of phylogenetic clustering and eveness, while MNTD is more sensitive to
# patterns of evenness and clustering closer to the tips of the phylogeny.
# For example,
# community ’clump4’ contains species that are spread randomly across the entire tree
# (SESMP D close to zero) but phylogenetically clustered towards the tips (negative
# SESMNT D and mntd.obs.p in the low quantiles of the null distribution).
# All of these measures can incorporate abundance information when available using
# the abundance.weighted argument. This will change the interpretation of these
# metrics from the mean phylogenetic distances among species, to the mean phylogenetic
# distances among individuals.


#### Faith's phylogenetic diversity (PD) - also calculates species richness (SP) ####

pd.result <- pd(matrix, tree, include.root=TRUE)
pd.result  #this is the observed PD index

ses.pd.result <- ses.pd(matrix, tree, null.model = c("taxa.labels"),
                        runs = 99, include.root=TRUE)
ses.pd.result  #this is the standardized effect size of the PD index (SESpd),  
# which compares the observed PD with a NULL (randomized) community


# Phylogenetic beta diversity ---------------------------------------------

# comdist = Calculates MPD (mean pairwise distance) separating taxa in two 
# communities, a measure of phylogenetic beta diversity
comdist.result <- comdist(matrix, phydist, abundance.weighted = FALSE)
comdist.result

library(cluster)
comdist.clusters <- hclust(comdist.result)
plot(comdist.clusters)

# comdistnt = Calculates MNTD (mean nearest taxon distance) separating taxa in 
# two communities, a measure of phylogenetic beta diversity
comdistnt.result <- comdistnt(matrix, phydist, abundance.weighted = FALSE)
comdistnt.result

comdistnt.clusters <- hclust(comdistnt.result)
plot(comdistnt.clusters)


# Phylogenetic index of beta-diversity PhyloSor (phylosor)

phylosor(matrix, prunedtree) # Fraction of branch-length shared between two communities
# there may not be enough difference in branch length in my data to use this metric


# Species' evolutionary distinctiveness -----------------------------------

evol.distinct(prunedtree, type = c("equal.splits", "fair.proportion"),
              scale = FALSE, use.branch.lengths = TRUE)

# Calculates evolutionary distinctiveness measures for a suite of species by: 
# a) equal splits or b) fair proportions. Returns a dataframe with species identifiers
# and species scores.

# Probably not useful for my purposes. Would serve to identify species of high 
# conservation value due to their distinctiveness.


# Phylogenetic Species Diversity Metrics (PSD variants) -------------------

# Calculate the bounded phylogenetic biodiversity metrics: phylogenetic species 
# variability, richness, evenness and clustering for one or multiple samples.
# more info here: https://rdrr.io/cran/picante/man/psd.html

##### every time I try to run this, my Rstudio freezes :( 
# psv(matrix,tree,compute.var=TRUE,scale.vcv=TRUE)
# psr(matrix,tree,compute.var=TRUE,scale.vcv=TRUE)
# pse(matrix,tree,scale.vcv=TRUE)
# psc(matrix,tree,scale.vcv=TRUE)
# psd(matrix,tree,compute.var=TRUE,scale.vcv=TRUE)
# psv.spp(matrix,tree)



# Null model for randomizing the community data matrix --------------------

null <- randomizeMatrix(matrix, null.model = c("frequency", "richness",
                  "independentswap", "trialswap"), iterations = 99)



# Rao's quadratic entropy (Rao's Q) ---------------------------------------

# Calculates Rao’s quadratic entropy, a measure of within- and among-community 
# diversity taking species dissimilarities into account

utree <- chronos(tree, lambda = 0) #converts tree to ultrametric, needed for Rao function
is.ultrametric(utree)

RaoQ <- raoD(matrix, utree)
RaoQ

# What do the results mean? #
# Dkk Within-community diversity
# Dkl Among-community diversity
# H Among-community diversities excluding within-community diversity
# total Total diversity across all samples
# alpha Alpha diversity - average within-community diversity
# beta Beta diversity - average among-community diversity
# Fst Beta diversity / total diversity


#---------------------------------------------------------------------------#

# Statistical analyses and graphs for SACNAS 2022 poster ---------------------

### Using the full matrix (by night) ###

# Extract indexes

mpd.values <- as.data.frame(mpd)

mntd.values <- as.data.frame(mntd)

pd.values <- as.data.frame(pd.result)

habitat <- select(alldata, Habitat)
period <- select(alldata,Period)


#create table with output of all diversity indexes
phylo.indexes.night <- cbind(habitat, period, mpd.values, mntd.values, pd.values)


#boxplots

mpd.plot <- ggboxplot(phylo.indexes.night, x = "Period", y = "mpd", width = 0.7,
                 palette = c("#00AFBB","#FC4E07"),size =0.3, fill= "Habitat", 
                 order = c("Pre-Hurricane", "Post-Hurricane"),legend = "top",
                 font.x = c(8, "black"),font.y = c(8, "black"),
                 font.ytickslab= c(8, "black"), font.xtickslab= c(8, "black")) + 
                 xlab(NULL) + ylab("Mean Pairwise Distance") + 
                 theme(legend.text = element_text(size=8), 
                       legend.title = element_text(size=10)) #+
                #ylim(0,0.5)
mpd.plot

mntd.plot <- ggboxplot(phylo.indexes.night, x = "Period", y = "mntd", width = 0.7,
                  palette = c("#00AFBB","#FC4E07"),size =0.3, fill= "Habitat", 
                  order = c("Pre-Hurricane", "Post-Hurricane"),legend = "top",
                  font.x = c(8, "black"),font.y = c(8, "black"),
                  font.ytickslab= c(8, "black"), font.xtickslab= c(8, "black")) + 
                  xlab(NULL) + ylab("Mean Nearest\nTaxon Distance") + 
                  theme(legend.text = element_text(size=8), 
                        legend.title = element_text(size=10)) +
                  scale_y_continuous(
                        labels = scales::number_format(accuracy = 0.01))
mntd.plot

pd.plot <- ggboxplot(phylo.indexes.night, x = "Period", y = "PD", width = 0.7,
                    palette = c("#00AFBB","#FC4E07"),size =0.3, fill= "Habitat", 
                    order = c("Pre-Hurricane", "Post-Hurricane"),legend = "top",
                    font.x = c(8, "black"),font.y = c(8, "black"),
                    font.ytickslab= c(8, "black"), font.xtickslab= c(8, "black")) + 
                    xlab("Period") + ylab("Faith's Phylogenetic Diversity") + 
                    theme(legend.text = element_text(size=8), 
                         legend.title = element_text(size=10)) +
                    scale_y_continuous(
                          labels = scales::number_format(accuracy = 0.01))
pd.plot


all.plots <- ggarrange(mpd.plot + rremove("x.text") , mntd.plot+ rremove("x.text"), 
                       pd.plot, align = "hv",
                #labels = c("A", "B", "C"),font.label = list(size = 8, color = "black"),
                ncol = 1, nrow = 3,
                common.legend = TRUE, legend = "top")
all.plots

# Two Way Anova with full data (by night)

aov2.mpd <- aov(mpd ~ Habitat * Period, data = phylo.indexes.night)
summary(aov2.mpd)
TukeyHSD(aov2.mpd)

aov2.mntd <- aov(mntd ~ Habitat * Period, data = phylo.indexes.night)
summary(aov2.mntd)
TukeyHSD(aov2.mntd)

aov2.pd <- aov(PD ~ Habitat * Period, data = phylo.indexes.night)
summary(aov2.pd)
TukeyHSD(aov2.pd)

aov2.out <- group_by(phylo.indexes.night, Habitat, Period) %>%
  summarise(
    count = n(),
    mpd_mean = mean(mpd, na.rm = TRUE),
    mpd_sd = sd(mpd, na.rm = TRUE),
    mntd_mean = mean(mntd, na.rm = TRUE),
    mntd_sd = sd(mntd, na.rm = TRUE),
    pd_mean = mean(PD, na.rm = TRUE),
    pd_sd = sd(PD, na.rm = TRUE),
  )

aov2.out


### Using the summary matrix (by site and period, after removing January 2018
# data to avoid having 6 months in the post period vs 5 months in the pre) ###


# Loading the summary data matrix

sum.data <- read.csv("matrix_202210_summed.csv")
moths <- select(sum.data, PRMAR661_21_Nystalea_superciliosa_A001c_Notodontidae:
                  PRMAR716_21_Glaphyria_dolatalis_M110a_Crambidae)
codes <- select(sum.data, Code)
sum.matrix <- cbind(codes, moths)
sum.matrix <- data.frame(sum.matrix[,-1], row.names=sum.matrix[,1])

class(sum.matrix)
sum.matrix <- as.matrix(sum.matrix)
colnames(sum.matrix)
rownames(sum.matrix)


# Calculating phylogenetic diversity metrics 

# We can prune the phylogeny to include only the species that actually 
# occurred in some community.
prunedtree <- prune.sample(sum.matrix, tree)
prunedtree

#### MPD and MNTD metrics ####

# First we need a phylogenetic distance matrix for the input
#A phylo object can be converted to a interspecific phylogenetic distance matrix 
# using the cophenetic function
phydist <- cophenetic(tree) 


# Calculating MPD index
# MPD = Mean Pairwise Distance

mpd.sum <- mpd(sum.matrix, phydist, abundance.weighted = FALSE) 
mpd.sum #this is the observed MPD index

# ses.mpd.result <- ses.mpd(matrix, phydist, null.model = "taxa.labels",
#                           abundance.weighted = FALSE, runs = 99) # eventually, change to 
# # runs = 999, iterations = 1000
# ses.mpd.result #this is the standardized effect size of the MPD index (SESmpd),  
# # which compares the observed MPD with a NULL (randomized) community


# Calculating MNTD index
# MNTD = Mean Nearest Taxon Distance

mntd.sum <- mntd(sum.matrix, phydist, abundance.weighted = FALSE)
mntd.sum #this is the observed MNTD index

# ses.mntd.result <- ses.mntd(matrix, phydist, null.model = "taxa.labels",
#                             abundance.weighted = FALSE, runs = 99)  # eventually, change to 
# # runs = 999, iterations = 1000
# ses.mntd.result #this is the standardized effect size of the MNTD index (SESmntd),  
# # which compares the observed MNTD with a NULL (randomized) community


#### Faith's phylogenetic diversity (PD) - also calculates species richness (SP) ####

pd.sum <- pd(sum.matrix, tree, include.root=TRUE)
pd.sum  #this is the observed PD index

# ses.pd.result <- ses.pd(matrix, tree, null.model = c("taxa.labels"),
#                         runs = 99, include.root=TRUE)
# ses.pd.result  #this is the standardized effect size of the PD index (SESpd),  
# # which compares the observed PD with a NULL (randomized) community


### Using the summary matrix (by site*period) ###

# Extract indexes

mpd.values.sum <- as.data.frame(mpd.sum)

mntd.values.sum <- as.data.frame(mntd.sum)

pd.values.sum <- as.data.frame(pd.sum)

habitat <- select(sum.data, Habitat)
period <- select(sum.data,Period)


#create table with output of all diversity indexes
phylo.indexes.sum <- cbind(habitat, period, mpd.values.sum, mntd.values.sum, pd.values.sum)


#boxplots

mpd.plot.sum <- ggboxplot(phylo.indexes.sum, x = "Period", y = "mpd.sum", width = 0.7,
                      palette = c("#00AFBB","#FC4E07"),size =0.3, fill= "Habitat", 
                      order = c("Pre-Hurricane", "Post-Hurricane"),legend = "top",
                      font.x = c(8, "black"),font.y = c(8, "black"),
                      font.ytickslab= c(8, "black"), font.xtickslab= c(8, "black")) + 
                      xlab(NULL) + ylab("Mean Pairwise Distance") + 
                      theme(legend.text = element_text(size=8), 
                          legend.title = element_text(size=10)) #+
                      #ylim(0,0.5)
mpd.plot.sum

mntd.plot.sum <- ggboxplot(phylo.indexes.sum, x = "Period", y = "mntd.sum", width = 0.7,
                       palette = c("#00AFBB","#FC4E07"),size =0.3, fill= "Habitat", 
                       order = c("Pre-Hurricane", "Post-Hurricane"),legend = "top",
                       font.x = c(8, "black"),font.y = c(8, "black"),
                       font.ytickslab= c(8, "black"), font.xtickslab= c(8, "black")) + 
                        xlab(NULL) + ylab("Mean Nearest Taxon Distance") + 
                        theme(legend.text = element_text(size=8), 
                            legend.title = element_text(size=10)) +
                        scale_y_continuous(
                            labels = scales::number_format(accuracy = 0.01))
mntd.plot.sum

pd.plot.sum <- ggboxplot(phylo.indexes.sum, x = "Period", y = "PD", width = 0.7,
                     palette = c("#00AFBB","#FC4E07"),size =0.3, fill= "Habitat", 
                     order = c("Pre-Hurricane", "Post-Hurricane"),legend = "top",
                     font.x = c(8, "black"),font.y = c(8, "black"),
                     font.ytickslab= c(8, "black"), font.xtickslab= c(8, "black")) + 
                  xlab("Period") + ylab("Faith's Phylogenetic Diversity") + 
                  theme(legend.text = element_text(size=8), 
                          legend.title = element_text(size=10)) +
                  scale_y_continuous(
                          labels = scales::number_format(accuracy = 0.01))
pd.plot.sum


all.plots.sum <- ggarrange(mpd.plot.sum + rremove("x.text") , mntd.plot.sum + rremove("x.text"), 
                       pd.plot.sum, align = "hv",
                       #labels = c("A", "B", "C"),font.label = list(size = 8, color = "black"),
                       ncol = 1, nrow = 3,
                       common.legend = TRUE, legend = "top")
all.plots.sum

# Two Way Anova with summed data 

aov2.mpd.sum <- aov(mpd.sum ~ Habitat * Period, data = phylo.indexes.sum)
summary(aov2.mpd.sum)
TukeyHSD(aov2.mpd.sum)

summary(aov2.mpd) #to compare with results of the full matrix

aov2.mntd.sum <- aov(mntd.sum ~ Habitat * Period, data = phylo.indexes.sum)
summary(aov2.mntd.sum)
TukeyHSD(aov2.mntd.sum)

summary(aov2.mntd) #to compare with results of the full matrix

aov2.pd.sum <- aov(PD ~ Habitat * Period, data = phylo.indexes.sum)
summary(aov2.pd.sum)
TukeyHSD(aov2.pd.sum)

summary(aov2.pd) #to compare with results of the full matrix

aov2.out.sum <- group_by(phylo.indexes.sum, Habitat, Period) %>%
  summarise(
    count = n(),
    mpd_mean = mean(mpd.sum, na.rm = TRUE),
    mpd_sd = sd(mpd.sum, na.rm = TRUE),
    mntd_mean = mean(mntd.sum, na.rm = TRUE),
    mntd_sd = sd(mntd.sum, na.rm = TRUE),
    pd_mean = mean(PD, na.rm = TRUE),
    pd_sd = sd(PD, na.rm = TRUE),
  )

aov2.out.sum


