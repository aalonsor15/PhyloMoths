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

tree <- read.tree("TestTree.newick")
tree
class(tree)


# Loading my community data --------------------------------------------------

matrix <- read.csv("test_summed.csv")
moths <- select(matrix, PRMAR660_21_Nystalea_superciliosa_Notodontidae:
                  PRMAR999_21_Coptarthria_dasypyga_Pyralidae)
codes <- select(matrix, Code)
matrix <- cbind(codes, moths)
matrix <- data.frame(matrix[,-1], row.names=matrix[,1])

class(matrix)
matrix <- as.matrix(matrix)
colnames(matrix)
rownames(matrix)


# calculating phylogenetic diversity metrics ------------------------------------

# We first need to prune the phylogeny to include only the species that actually 
# occurred in some community.
prunedtree <- prune.sample(matrix, tree)
prunedtree

#### Faith's phylogenetic diversity (PD) - also calculates species richness (SP) ####
pd.result <- pd(matrix, tree, include.root=TRUE)
pd.result

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
                          abundance.weighted = FALSE, runs = 99)
ses.mpd.result #this is the standardized effect size of the MPD index (SESmpd),  
# which compares the observed MPD with a NULL (randomized) community


# Calculating MNTD index
# MNTD = Mean Nearest Taxon Distance

mntd <- mntd(matrix, phydist, abundance.weighted = FALSE)
mntd #this is the observed MNTD index

ses.mntd.result <- ses.mntd(matrix, phydist, null.model = "taxa.labels",
                            abundance.weighted = FALSE, runs = 99)
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



# Species' evolutionary distinctiveness -----------------------------------

evol.distinct(prunedtree, type = c("equal.splits", "fair.proportion"),
              scale = FALSE, use.branch.lengths = TRUE)

# Calculates evolutionary distinctiveness measures for a suite of species by: 
# a) equal splits or b) fair proportions. Returns a dataframe with species identifiers
# and species scores.



# Phylogenetic Species Diversity Metrics (PSD variants) -------------------

# Calculate the bounded phylogenetic biodiversity metrics: phylogenetic species 
# variability, richness, evenness and clustering for one or multiple samples.
# more info here: https://rdrr.io/cran/picante/man/psd.html

# every time I try to run this, my Rstudio freezes :( 
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


