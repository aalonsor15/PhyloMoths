# ----------------------------------------------------------------------
# El Verde moths project: Calculating phylogenetic diversity metrics
# 26 May 2022
# Aura M. Alonso-Rodríguez
# ----------------------------------------------------------------------
#

# loading libraries

#install.packages("picante", dependencies = TRUE)
library(picante)
library(tidyverse)
library(vegan)


# practicing picante ------------------------------------------------

data(phylocom)
names(phylocom)

phy <- phylocom$phylo
comm <- phylocom$sample
traits <- phylocom$traits

phy 

prunedphy <- prune.sample(comm, phy)

rownames(comm)

pd.result <- pd(comm, phy, include.root=TRUE)
pd.result


#SUPER IMPORTANT NOTE ----------------------------------------------------------
# The names on the tree newick file need to be the same as in the matrix.
# I need to remove any - or > or any other symbol from the name, the only symbol
# I should have is the underscore _.
# When I run the tree, I need to make sure the name does not contain any symbols
# except the _. And the names in the matrix should be written exactly the same.
# I can edit the names in all files (newick, fafsa, excel) in Notepad using 
# Find/Replace

# Loading test tree -----------------------------------------------------------

#  How to upload tree data -> If you have a phylogeny in Newick or Nexus 
# format it can be imported into R with the read.tree or read.nexus functions

tree <- read.tree("TestTree.newick")
tree
class(tree)

# Loading community data --------------------------------------------------

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

#### Faith's phylogenetic diversity (PD) - also calculated species richness (SP) ####
pd.result <- pd(matrix, tree, include.root=TRUE)
pd.result

#### MPD and MNTD metrics ####

# First we need a phylogenetic distance matrix for the input
#A phylo object can be converted to a interspecific phylogenetic distance matrix 
# using the cophenetic function
phydist <- cophenetic(tree) 

# Calculating MPD index and comparing with NULL (fake) community
# MPD = Mean Pairwise Distance
ses.mpd.result <- ses.mpd(matrix, phydist, null.model = "taxa.labels",
                          abundance.weighted = FALSE, runs = 99)
ses.mpd.result

# Calculating MNTD index and comparing with NULL (fake) community
# MNTD = Mean Nearest Taxon Distance

ses.mntd.result <- ses.mntd(matrix, phydist, null.model = "taxa.labels",
                            abundance.weighted = FALSE, runs = 99)
ses.mntd.result

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



