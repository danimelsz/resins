# Testing signal using Mantel tests

#### TOTAL MATRIX: PHYLOGENY VS INTERACTIONS (CONTROL: NONE) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(ape)
library(bipartite)
library(dplyr)
library(igraph)
library(reshape2)
library(vegan)

# INTERACTION DISTANCE MATRIX
# Load edge and node lists
nodes = read.delim("data/net1nodes_.txt", header = T) # node traits
links = read.delim("data/list.3MBSDs_v3_binary", header = T) # list of interactions
# Load network
all_agreggated <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
all_agreggated = get.adjacency(all_agreggated,sparse=FALSE)
all_agreggated = all_agreggated[c(0:68), ] # keep only rows of bees
all_agreggated = all_agreggated[, 69:ncol(all_agreggated)] # keep only columns of plants
data_int <- empty(all_agreggated, count = F) # remove nodes with no interactions
# Create a distance matrix using Jaccard
dist_interactions = vegdist (data_int, diag=TRUE, upper=T, method="jaccard", binary=T)

# PHYLOGENETIC DISTANCE MATRIX
# Load phylogeny
tree = read.nexus("data/tree_edited.nex")
plotTree(tree,node.numbers=T,fsize=0.35)
# Find which species are common in dist_interactions and phylogeny
tree_species = tree$tip.label # list species in the tree
dist_interactions = as.matrix(dist_interactions) # convert dist to matrix (rownames does not work properly to dist objects)
matrix_species = rownames(dist_interactions) # list species in the matrix
common_species = intersect(tree_species, matrix_species) # list common species
# Prune tree
pruned_tree = drop.tip(tree, tree_species[!(tree_species %in% common_species)])
plotTree(pruned_tree,node.numbers=T,fsize=0.35)
# Prune dist_interactions
dist_interactions = dist_interactions[common_species, common_species]
# Calculate phylogenetic distance matrix
dist_phylo = cophenetic.phylo (pruned_tree)
# Rows and columns of dist_interactions are not in the same order of dist_phylo and thus we must correct it
sorted_rows <- dist_phylo[order(rownames(dist_phylo)), ] # Sort rows in alphabetical order
dist_phylo <- sorted_rows[, order(colnames(sorted_rows))] # Sort columns in alphabetical order

# MANTEL TEST
mantel_intVSphy = mantel(dist_interactions, dist_phylo, method="spearman", permutations = 10000)
mantel_intVSphy
densityplot(permustats(mantel_intVSphy))

#### TOTAL MATRIX: PHYLOGENY VS MODULES (CONTROL: NONE) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(ape)
library(bipartite)
library(dplyr)
library(igraph)
library(reshape2)
library(vegan)

# MODULE DISTANCE MATRIX
# Read module dataframe
data_mod = read_tsv("data/partitions_all.txt")
data_mod = data_mod[c(0:68), ] # keep only rows of bees
data_mod = data_mod[, 1:2] # keep only columns of plants
data_mod = dcast(data_mod, nodes ~ module, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_mod) <- data_mod$nodes
data_mod <- data_mod[, -1]
# Create a distance matrix using Jaccard
dist_mod = vegdist (data_mod, diag=TRUE, upper=T, method="jaccard", binary=T)

# PHYLOGENETIC DISTANCE MATRIX
# Load phylogeny
tree = read.nexus("data/tree_edited.nex")
plotTree(tree,node.numbers=T,fsize=0.35)
# Find which species are common in dist_mod and phylogeny
tree_species = tree$tip.label # list species in the tree
dist_mod = as.matrix(dist_mod) # convert dist to matrix (rownames does not work properly to dist objects)
matrix_species = rownames(dist_mod) # list species in the matrix
common_species = intersect(tree_species, matrix_species) # list common species
# Prune tree
pruned_tree = drop.tip(tree, tree_species[!(tree_species %in% common_species)])
plotTree(pruned_tree,node.numbers=T,fsize=0.35)
# Prune dist_mod
dist_mod = dist_mod[common_species, common_species]
# Calculate phylogenetic distance matrix
dist_phylo = cophenetic.phylo (pruned_tree)
# Rows and columns of dist_mod are not in the same order of dist_phylo and thus we must correct it
sorted_rows <- dist_phylo[order(rownames(dist_phylo)), ] # Sort rows in alphabetical order
dist_phylo <- sorted_rows[, order(colnames(sorted_rows))] # Sort columns in alphabetical order

# MANTEL TEST
mantel_modVSphy = mantel(dist_mod, dist_phylo, method="spearman", permutations = 10000)
mantel_modVSphy
densityplot(permustats(mantel_modVSphy))

#### TOTAL MATRIX: PHYLOGENY VS LAYERS (CONTROL: NONE) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(ape)
library(bipartite)
library(dplyr)
library(igraph)
library(reshape2)
library(vegan)

# LAYER DISTANCE MATRIX
data_layer = read_tsv("data/layers_total.txt")
data_layer = dcast(data_layer, Species ~ Layer, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_layer) <- data_layer$Species # turn the first column into rownames
data_layer <- data_layer[, -1] # remove first column
dist_layer = vegdist (data_layer, diag=TRUE, upper=T, method="jaccard", binary=T)

# PHYLOGENETIC DISTANCE MATRIX
# Load phylogeny
tree = read.nexus("data/tree_edited.nex")
plotTree(tree,node.numbers=T,fsize=0.35)
# Find which species are common in dist_layer and phylogeny
tree_species = tree$tip.label # list species in the tree
dist_layer = as.matrix(dist_layer) # convert dist to matrix (rownames does not work properly to dist objects)
matrix_species = rownames(dist_layer) # list species in the matrix
common_species = intersect(tree_species, matrix_species) # list common species
# Prune tree
pruned_tree = drop.tip(tree, tree_species[!(tree_species %in% common_species)])
plotTree(pruned_tree,node.numbers=T,fsize=0.35)
# Prune dist_layer
dist_layer = dist_layer[common_species, common_species]
# Calculate phylogenetic distance matrix
dist_phylo = cophenetic.phylo (pruned_tree)
# Rows and columns of dist_mod are not in the same order of dist_phylo and thus we must correct it
sorted_rows <- dist_phylo[order(rownames(dist_phylo)), ] # Sort rows in alphabetical order
dist_phylo <- sorted_rows[, order(colnames(sorted_rows))] # Sort columns in alphabetical order

# MANTEL TEST
mantel_layVSphy = mantel(dist_layer, dist_phylo, method="spearman", permutations = 10000)
mantel_layVSphy
densityplot(permustats(mantel_layVSphy))


#### TOTAL MATRIX: MODULE VS INTERACTIONS (CONTROL: NONE) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(bipartite)
library(igraph)
library(reshape2)
library(vegan)

# INTERACTION DISTANCE MATRIX
# Load edge and node lists
nodes = read.delim("data/net1nodes_.txt", header = T) # node traits
links = read.delim("data/list.3MBSDs_v3_binary", header = T) # list of interactions
# Load network
all_agreggated <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
all_agreggated = get.adjacency(all_agreggated,sparse=FALSE)
all_agreggated = all_agreggated[c(0:68), ] # keep only rows of bees
all_agreggated = all_agreggated[, 69:ncol(all_agreggated)] # keep only columns of plants
data_int <- empty(all_agreggated, count = F) # remove nodes with no interactions
# Create a distance matrix using Jaccard
dist_interactions = vegdist (data_int, diag=TRUE, upper=T, method="jaccard", binary=T)

# MODULE DISTANCE MATRIX
# Read module dataframe
data_mod = read_tsv("data/partitions_all.txt")
data_mod = data_mod[c(0:68), ] # keep only rows of bees
data_mod = data_mod[, 1:2] # keep only columns of plants
data_mod = dcast(data_mod, nodes ~ module, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_mod) <- data_mod$nodes
data_mod <- data_mod[, -1]
# Create a distance matrix using Jaccard
dist_mod = vegdist (data_mod, diag=TRUE, upper=T, method="jaccard", binary=T)

# MANTEL TEST
mantel_intVSmod = mantel(dist_interactions, dist_mod, method="spearman", permutations = 10000)
mantel_intVSmod
densityplot(permustats(mantel_intVSmod))
            
#### TOTAL MATRIX: MODULES VS LAYERS (CONTROL: NONE) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(bipartite)
library(igraph)
library(reshape2)
library(vegan)

# MODULE DISTANCE MATRIX
# Read module dataframe
data_mod = read_tsv("data/partitions_all.txt")
data_mod = data_mod[c(0:68), ] # keep only rows of bees
data_mod = data_mod[, 1:2] # keep only columns of plants
data_mod = dcast(data_mod, nodes ~ module, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_mod) <- data_mod$nodes
data_mod <- data_mod[, -1]
# Create a distance matrix using Jaccard
dist_mod = vegdist (data_mod, diag=TRUE, upper=T, method="jaccard", binary=T)

# LAYER DISTANCE MATRIX
data_layer = read_tsv("data/layers_total.txt")
data_layer = dcast(data_layer, Species ~ Layer, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_layer) <- data_layer$Species # turn the first column into rownames
data_layer <- data_layer[, -1] # remove first column
dist_layer = vegdist (data_layer, diag=TRUE, upper=T, method="jaccard", binary=T)

# MANTEL TEST
mantel_modVSlay = mantel(dist_mod, dist_layer, method="spearman", permutations = 10000)
mantel_modVSlay
densityplot(permustats(mantel_modVSlay))


#### TOTAL MATRIX: INTERACTIONS VS LAYERS (CONTROL: NONE) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 
#Load the required packages and functions.
library(bipartite)
library(igraph)
library(reshape2)
library(vegan)

# INTERACTION DISTANCE MATRIX
# Load edge and node lists
nodes = read.delim("data/net1nodes_.txt", header = T) # node traits
links = read.delim("data/list.3MBSDs_v3_binary", header = T) # list of interactions
# Load network
all_agreggated <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
all_agreggated = get.adjacency(all_agreggated,sparse=FALSE)
all_agreggated = all_agreggated[c(0:68), ] # keep only rows of bees
all_agreggated = all_agreggated[, 69:ncol(all_agreggated)] # keep only columns of plants
data_int <- empty(all_agreggated, count = F) # remove nodes with no interactions
# Create a distance matrix using Jaccard
dist_interactions = vegdist (data_int, diag=TRUE, upper=T, method="jaccard", binary=T)

# LAYER DISTANCE MATRIX
data_layer = read_tsv("data/layers_total.txt")
data_layer = dcast(data_layer, Species ~ Layer, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_layer) <- data_layer$Species # turn the first column into rownames
data_layer <- data_layer[, -1] # remove first column
dist_layer = vegdist (data_layer, diag=TRUE, upper=T, method="jaccard", binary=T)

# MANTEL TEST
mantel_intVSlay = mantel(dist_interactions, dist_layer, method="spearman", permutations = 10000)
mantel_intVSlay
densityplot(permustats(mantel_intVSlay))



#### TOTAL MATRIX: PHYLOGENY VS INTERACTIONS (CONTROL: LAYERS) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(ape)
library(bipartite)
library(dplyr)
library(igraph)
library(reshape2)
library(vegan)

# INTERACTION DISTANCE MATRIX
# Load edge and node lists
nodes = read.delim("data/net1nodes_.txt", header = T) # node traits
links = read.delim("data/list.3MBSDs_v3_binary", header = T) # list of interactions
# Load network
all_agreggated <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
all_agreggated = get.adjacency(all_agreggated,sparse=FALSE)
all_agreggated = all_agreggated[c(0:68), ] # keep only rows of bees
all_agreggated = all_agreggated[, 69:ncol(all_agreggated)] # keep only columns of plants
data_int <- empty(all_agreggated, count = F) # remove nodes with no interactions
# Create a distance matrix using Jaccard
dist_interactions = vegdist (data_int, diag=TRUE, upper=T, method="jaccard", binary=T)

# PHYLOGENETIC DISTANCE MATRIX
# Load phylogeny
tree = read.nexus("data/tree_edited.nex")
plotTree(tree,node.numbers=T,fsize=0.35)
# Find which species are common in dist_interactions and phylogeny
tree_species = tree$tip.label # list species in the tree
dist_interactions = as.matrix(dist_interactions) # convert dist to matrix (rownames does not work properly to dist objects)
matrix_species = rownames(dist_interactions) # list species in the matrix
common_species = intersect(tree_species, matrix_species) # list common species
# Prune tree
pruned_tree = drop.tip(tree, tree_species[!(tree_species %in% common_species)])
plotTree(pruned_tree,node.numbers=T,fsize=0.35)
# Prune dist_interactions
dist_interactions = dist_interactions[common_species, common_species]
# Calculate phylogenetic distance matrix
dist_phylo = cophenetic.phylo (pruned_tree)
# Rows and columns of dist_interactions are not in the same order of dist_phylo and thus we must correct it
sorted_rows <- dist_phylo[order(rownames(dist_phylo)), ] # Sort rows in alphabetical order
dist_phylo <- sorted_rows[, order(colnames(sorted_rows))] # Sort columns in alphabetical order

# LAYER DISTANCE MATRIX
data_layer = read_tsv("data/layers_total.txt")
data_layer = dcast(data_layer, Species ~ Layer, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_layer) <- data_layer$Species # turn the first column into rownames
data_layer <- data_layer[, -1] # remove first column
# Find which species are common between dist_interactions and data_layer
matrix_species1 = rownames(data_layer) # list species in the matrix
matrix_species2 = rownames(dist_interactions) # list species in the matrix
common_species = intersect(matrix_species1, matrix_species2) # list common species
# Prune data_layer
data_layer = data_layer[matrix_species1 %in% common_species, ]
dist_layer = vegdist (data_layer, diag=TRUE, upper=T, method="jaccard", binary=T)

# PARTIAL MANTEL TEST
mantel_phyVSintCOlay = mantel.partial(dist_phylo, dist_interactions, dist_layer, method = "pearson", permutations = 10000)
mantel_phyVSintCOlay

#### TOTAL MATRIX: PHYLOGENY VS MODULES (CONTROL: LAYERS) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(ape)
library(bipartite)
library(dplyr)
library(igraph)
library(reshape2)
library(vegan)

# MODULE DISTANCE MATRIX
# Read module dataframe
data_mod = read_tsv("data/partitions_all.txt")
data_mod = data_mod[c(0:68), ] # keep only rows of bees
data_mod = data_mod[, 1:2] # keep only columns of plants
data_mod = dcast(data_mod, nodes ~ module, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_mod) <- data_mod$nodes
data_mod <- data_mod[, -1]
# Create a distance matrix using Jaccard
dist_mod = vegdist (data_mod, diag=TRUE, upper=T, method="jaccard", binary=T)

# PHYLOGENETIC DISTANCE MATRIX
# Load phylogeny
tree = read.nexus("data/tree_edited.nex")
plotTree(tree,node.numbers=T,fsize=0.35)
# Find which species are common in dist_mod and phylogeny
tree_species = tree$tip.label # list species in the tree
dist_mod = as.matrix(dist_mod) # convert dist to matrix (rownames does not work properly to dist objects)
matrix_species = rownames(dist_mod) # list species in the matrix
common_species = intersect(tree_species, matrix_species) # list common species
# Prune tree
pruned_tree = drop.tip(tree, tree_species[!(tree_species %in% common_species)])
plotTree(pruned_tree,node.numbers=T,fsize=0.35)
# Prune dist_mod
dist_mod = dist_mod[common_species, common_species]
# Calculate phylogenetic distance matrix
dist_phylo = cophenetic.phylo (pruned_tree)
# Rows and columns of dist_mod are not in the same order of dist_phylo and thus we must correct it
sorted_rows <- dist_phylo[order(rownames(dist_phylo)), ] # Sort rows in alphabetical order
dist_phylo <- sorted_rows[, order(colnames(sorted_rows))] # Sort columns in alphabetical order

# LAYER DISTANCE MATRIX
data_layer = read_tsv("data/layers_total.txt")
data_layer = dcast(data_layer, Species ~ Layer, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_layer) <- data_layer$Species # turn the first column into rownames
data_layer <- data_layer[, -1] # remove first column
# Find which species are common between dist_mod and data_layer
matrix_species1 = rownames(data_layer) # list species in the matrix
matrix_species2 = rownames(dist_mod) # list species in the matrix
common_species = intersect(matrix_species1, matrix_species2) # list common species
# Prune data_layer
data_layer = data_layer[matrix_species1 %in% common_species, ]
dist_layer = vegdist (data_layer, diag=TRUE, upper=T, method="jaccard", binary=T)

# PARTIAL MANTEL TEST
mantel_phyVSmodCOlay = mantel.partial(dist_phylo, dist_mod, dist_layer, method = "pearson", permutations = 10000)
mantel_phyVSmodCOlay


#### TOTAL MATRIX: MODULES VS INTERACTIONS (CONTROL: PHYLOGENY) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(ape)
library(bipartite)
library(dplyr)
library(igraph)
library(reshape2)
library(vegan)

# MODULE DISTANCE MATRIX
# Read module dataframe
data_mod = read_tsv("data/partitions_all.txt")
data_mod = data_mod[c(0:68), ] # keep only rows of bees
data_mod = data_mod[, 1:2] # keep only columns of plants
data_mod = dcast(data_mod, nodes ~ module, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_mod) <- data_mod$nodes
data_mod <- data_mod[, -1]
# Create a distance matrix using Jaccard
dist_mod = vegdist (data_mod, diag=TRUE, upper=T, method="jaccard", binary=T)

# PHYLOGENETIC DISTANCE MATRIX
# Load phylogeny
tree = read.nexus("data/tree_edited.nex")
plotTree(tree,node.numbers=T,fsize=0.35)
# Find which species are common in dist_mod and phylogeny
tree_species = tree$tip.label # list species in the tree
dist_mod = as.matrix(dist_mod) # convert dist to matrix (rownames does not work properly to dist objects)
matrix_species = rownames(dist_mod) # list species in the matrix
common_species = intersect(tree_species, matrix_species) # list common species
# Prune tree
pruned_tree = drop.tip(tree, tree_species[!(tree_species %in% common_species)])
plotTree(pruned_tree,node.numbers=T,fsize=0.35)
# Prune dist_mod
dist_mod = dist_mod[common_species, common_species]
# Calculate phylogenetic distance matrix
dist_phylo = cophenetic.phylo (pruned_tree)
# Rows and columns of dist_mod are not in the same order of dist_phylo and thus we must correct it
sorted_rows <- dist_phylo[order(rownames(dist_phylo)), ] # Sort rows in alphabetical order
dist_phylo <- sorted_rows[, order(colnames(sorted_rows))] # Sort columns in alphabetical order

# INTERACTION DISTANCE MATRIX
# Load edge and node lists
nodes = read.delim("data/net1nodes_.txt", header = T) # node traits
links = read.delim("data/list.3MBSDs_v3_binary", header = T) # list of interactions
# Load network
all_agreggated <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
all_agreggated = get.adjacency(all_agreggated,sparse=FALSE)
all_agreggated = all_agreggated[c(0:68), ] # keep only rows of bees
all_agreggated = all_agreggated[, 69:ncol(all_agreggated)] # keep only columns of plants
data_int <- empty(all_agreggated, count = F) # remove nodes with no interactions
# Find which species are common between dist_interactions and data_layer
matrix_species1 = rownames(dist_mod) # list species in the matrix
matrix_species2 = rownames(data_int) # list species in the matrix
common_species = intersect(matrix_species1, matrix_species2) # list common species
# Prune data_layer
data_int = data_int[matrix_species2 %in% common_species, ]
# Create a distance matrix using Jaccard
dist_int = vegdist (data_int, diag=TRUE, upper=T, method="jaccard", binary=T)

# PARTIAL MANTEL TEST
mantel_modVSintCOphy = mantel.partial(dist_int, dist_mod, dist_phylo, method = "pearson", permutations = 10000)
mantel_modVSintCOphy

#### TOTAL MATRIX: MODULES VS LAYERS (CONTROL: PHYLOGENY) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(ape)
library(bipartite)
library(dplyr)
library(igraph)
library(reshape2)
library(vegan)

# MODULE DISTANCE MATRIX
# Read module dataframe
data_mod = read_tsv("data/partitions_all.txt")
data_mod = data_mod[c(0:68), ] # keep only rows of bees
data_mod = data_mod[, 1:2] # keep only columns of plants
data_mod = dcast(data_mod, nodes ~ module, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_mod) <- data_mod$nodes
data_mod <- data_mod[, -1]
# Create a distance matrix using Jaccard
dist_mod = vegdist (data_mod, diag=TRUE, upper=T, method="jaccard", binary=T)

# PHYLOGENETIC DISTANCE MATRIX
# Load phylogeny
tree = read.nexus("data/tree_edited.nex")
plotTree(tree,node.numbers=T,fsize=0.35)
# Find which species are common in dist_mod and phylogeny
tree_species = tree$tip.label # list species in the tree
dist_mod = as.matrix(dist_mod) # convert dist to matrix (rownames does not work properly to dist objects)
matrix_species = rownames(dist_mod) # list species in the matrix
common_species = intersect(tree_species, matrix_species) # list common species
# Prune tree
pruned_tree = drop.tip(tree, tree_species[!(tree_species %in% common_species)])
plotTree(pruned_tree,node.numbers=T,fsize=0.35)
# Prune dist_mod
dist_mod = dist_mod[common_species, common_species]
# Calculate phylogenetic distance matrix
dist_phylo = cophenetic.phylo (pruned_tree)
# Rows and columns of dist_mod are not in the same order of dist_phylo and thus we must correct it
sorted_rows <- dist_phylo[order(rownames(dist_phylo)), ] # Sort rows in alphabetical order
dist_phylo <- sorted_rows[, order(colnames(sorted_rows))] # Sort columns in alphabetical order

# LAYER DISTANCE MATRIX
data_layer = read_tsv("data/layers_total.txt")
data_layer = dcast(data_layer, Species ~ Layer, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_layer) <- data_layer$Species # turn the first column into rownames
data_layer <- data_layer[, -1] # remove first column
# Find which species are common between dist_mod and data_layer
matrix_species1 = rownames(data_layer) # list species in the matrix
matrix_species2 = rownames(dist_mod) # list species in the matrix
common_species = intersect(matrix_species1, matrix_species2) # list common species
# Prune data_layer
data_layer = data_layer[matrix_species1 %in% common_species, ]
dist_layer = vegdist (data_layer, diag=TRUE, upper=T, method="jaccard", binary=T)

# PARTIAL MANTEL TEST
mantel_modVSlayCOphy = mantel.partial(dist_mod, dist_layer, dist_phylo, method = "pearson", permutations = 10000)
mantel_modVSlayCOphy

#### TOTAL MATRIX: PHYLOGENY VS INTERACTIONS (CONTROL: MODULES) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(ape)
library(bipartite)
library(dplyr)
library(igraph)
library(reshape2)
library(vegan)

# MODULE DISTANCE MATRIX
# Read module dataframe
data_mod = read_tsv("data/partitions_all.txt")
data_mod = data_mod[c(0:68), ] # keep only rows of bees
data_mod = data_mod[, 1:2] # keep only columns of plants
data_mod = dcast(data_mod, nodes ~ module, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_mod) <- data_mod$nodes
data_mod <- data_mod[, -1]
# Create a distance matrix using Jaccard
dist_mod = vegdist (data_mod, diag=TRUE, upper=T, method="jaccard", binary=T)

# PHYLOGENETIC DISTANCE MATRIX
# Load phylogeny
tree = read.nexus("data/tree_edited.nex")
plotTree(tree,node.numbers=T,fsize=0.35)
# Find which species are common in dist_mod and phylogeny
tree_species = tree$tip.label # list species in the tree
dist_mod = as.matrix(dist_mod) # convert dist to matrix (rownames does not work properly to dist objects)
matrix_species = rownames(dist_mod) # list species in the matrix
common_species = intersect(tree_species, matrix_species) # list common species
# Prune tree
pruned_tree = drop.tip(tree, tree_species[!(tree_species %in% common_species)])
plotTree(pruned_tree,node.numbers=T,fsize=0.35)
# Prune dist_mod
dist_mod = dist_mod[common_species, common_species]
# Calculate phylogenetic distance matrix
dist_phylo = cophenetic.phylo (pruned_tree)
# Rows and columns of dist_mod are not in the same order of dist_phylo and thus we must correct it
sorted_rows <- dist_phylo[order(rownames(dist_phylo)), ] # Sort rows in alphabetical order
dist_phylo <- sorted_rows[, order(colnames(sorted_rows))] # Sort columns in alphabetical order

# INTERACTION DISTANCE MATRIX
# Load edge and node lists
nodes = read.delim("data/net1nodes_.txt", header = T) # node traits
links = read.delim("data/list.3MBSDs_v3_binary", header = T) # list of interactions
# Load network
all_agreggated <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
all_agreggated = get.adjacency(all_agreggated,sparse=FALSE)
all_agreggated = all_agreggated[c(0:68), ] # keep only rows of bees
all_agreggated = all_agreggated[, 69:ncol(all_agreggated)] # keep only columns of plants
data_int <- empty(all_agreggated, count = F) # remove nodes with no interactions
# Find which species are common between dist_interactions and data_layer
matrix_species1 = rownames(dist_mod) # list species in the matrix
matrix_species2 = rownames(data_int) # list species in the matrix
common_species = intersect(matrix_species1, matrix_species2) # list common species
# Prune data_layer
data_int = data_int[matrix_species2 %in% common_species, ]
# Create a distance matrix using Jaccard
dist_int = vegdist (data_int, diag=TRUE, upper=T, method="jaccard", binary=T)

# PARTIAL MANTEL TEST
mantel_phyVSintCOmod = mantel.partial(dist_phylo, dist_int, dist_mod, method = "pearson", permutations = 10000)
mantel_phyVSintCOmod

#### TOTAL MATRIX: PHYLOGENY VS LAYERS (CONTROL: MODULES) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(ape)
library(bipartite)
library(dplyr)
library(igraph)
library(reshape2)
library(vegan)

# MODULE DISTANCE MATRIX
# Read module dataframe
data_mod = read_tsv("data/partitions_all.txt")
data_mod = data_mod[c(0:68), ] # keep only rows of bees
data_mod = data_mod[, 1:2] # keep only columns of plants
data_mod = dcast(data_mod, nodes ~ module, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_mod) <- data_mod$nodes
data_mod <- data_mod[, -1]
# Create a distance matrix using Jaccard
dist_mod = vegdist (data_mod, diag=TRUE, upper=T, method="jaccard", binary=T)

# PHYLOGENETIC DISTANCE MATRIX
# Load phylogeny
tree = read.nexus("data/tree_edited.nex")
plotTree(tree,node.numbers=T,fsize=0.35)
# Find which species are common in dist_mod and phylogeny
tree_species = tree$tip.label # list species in the tree
dist_mod = as.matrix(dist_mod) # convert dist to matrix (rownames does not work properly to dist objects)
matrix_species = rownames(dist_mod) # list species in the matrix
common_species = intersect(tree_species, matrix_species) # list common species
# Prune tree
pruned_tree = drop.tip(tree, tree_species[!(tree_species %in% common_species)])
plotTree(pruned_tree,node.numbers=T,fsize=0.35)
# Prune dist_mod
dist_mod = dist_mod[common_species, common_species]
# Calculate phylogenetic distance matrix
dist_phylo = cophenetic.phylo (pruned_tree)
# Rows and columns of dist_mod are not in the same order of dist_phylo and thus we must correct it
sorted_rows <- dist_phylo[order(rownames(dist_phylo)), ] # Sort rows in alphabetical order
dist_phylo <- sorted_rows[, order(colnames(sorted_rows))] # Sort columns in alphabetical order

# LAYER DISTANCE MATRIX
data_layer = read_tsv("data/layers_total.txt")
data_layer = dcast(data_layer, Species ~ Layer, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_layer) <- data_layer$Species # turn the first column into rownames
data_layer <- data_layer[, -1] # remove first column
# Find which species are common between dist_mod and data_layer
matrix_species1 = rownames(data_layer) # list species in the matrix
matrix_species2 = rownames(dist_mod) # list species in the matrix
common_species = intersect(matrix_species1, matrix_species2) # list common species
# Prune data_layer
data_layer = data_layer[matrix_species1 %in% common_species, ]
dist_layer = vegdist (data_layer, diag=TRUE, upper=T, method="jaccard", binary=T)

# PARTIAL MANTEL TEST
mantel_phyVSlayCOmod = mantel.partial(dist_phylo, dist_layer, dist_mod, method = "pearson", permutations = 10000)
mantel_phyVSlayCOmod



#### OK CONSERVATIVE MATRIX: PHYLOGENY VS INTERACTIONS (CONTROL: NONE) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(ape)
library(bipartite)
library(dplyr)
library(igraph)
library(reshape2)
library(vegan)

# INTERACTION DISTANCE MATRIX
# Load edge and node lists
nodes = read.delim("data/net1nodes_CF.txt", header = T) # node traits
links = read.delim("data/list_CF", header = T) # list of interactions
# Load network
all_agreggated <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
all_agreggated = get.adjacency(all_agreggated,sparse=FALSE)
all_agreggated = all_agreggated[c(0:62), ] # keep only rows of bees
all_agreggated = all_agreggated[, 63:ncol(all_agreggated)] # keep only columns of plants
data_int <- empty(all_agreggated, count = F) # remove nodes with no interactions
# Create a distance matrix using Jaccard
dist_interactions = vegdist (data_int, diag=TRUE, upper=T, method="jaccard", binary=T)
dim(dist_interactions) # nrows and ncolumns

# PHYLOGENETIC DISTANCE MATRIX
# Load phylogeny
tree = read.nexus("data/tree_edited.nex")
plotTree(tree,node.numbers=T,fsize=0.35)
# Find which species are common in dist_interactions and phylogeny
tree_species = tree$tip.label # list species in the tree
dist_interactions = as.matrix(dist_interactions) # convert dist to matrix (rownames does not work properly to dist objects)
matrix_species = rownames(dist_interactions) # list species in the matrix
common_species = intersect(tree_species, matrix_species) # list common species
# Prune tree
pruned_tree = drop.tip(tree, tree_species[!(tree_species %in% common_species)])
plotTree(pruned_tree,node.numbers=T,fsize=0.35)
# Prune dist_interactions
dist_interactions = dist_interactions[common_species, common_species]
dim(dist_interactions) # nrows and ncolumns
# Calculate phylogenetic distance matrix
dist_phylo = cophenetic.phylo (pruned_tree)
# Rows and columns of dist_interactions are not in the same order of dist_phylo and thus we must correct it
sorted_rows <- dist_phylo[order(rownames(dist_phylo)), ] # Sort rows in alphabetical order
dist_phylo <- sorted_rows[, order(colnames(sorted_rows))] # Sort columns in alphabetical order
dim(dist_phylo)

# MANTEL TEST
mantel_intVSphy = mantel(dist_interactions, dist_phylo, method="spearman", permutations = 10000)
mantel_intVSphy
densityplot(permustats(mantel_intVSphy))

#### OK CONSERVATIVE MATRIX: PHYLOGENY VS MODULES (CONTROL: NONE) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(ape)
library(bipartite)
library(dplyr)
library(igraph)
library(reshape2)
library(vegan)

# MODULE DISTANCE MATRIX
# Read module dataframe
data_mod = read_tsv("data/partitions_cf.txt")
data_mod = data_mod[c(0:62), ] # keep only rows of bees
data_mod = data_mod[, 1:2] # keep only columns of plants
data_mod = dcast(data_mod, nodes ~ module, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_mod) <- data_mod$nodes
data_mod <- data_mod[, -1]
# Create a distance matrix using Jaccard
dist_mod = vegdist (data_mod, diag=TRUE, upper=T, method="jaccard", binary=T)
dim(dist_mod)

# PHYLOGENETIC DISTANCE MATRIX
# Load phylogeny
tree = read.nexus("data/tree_edited.nex")
plotTree(tree,node.numbers=T,fsize=0.35)
# Find which species are common in dist_mod and phylogeny
tree_species = tree$tip.label # list species in the tree
dist_mod = as.matrix(dist_mod) # convert dist to matrix (rownames does not work properly to dist objects)
matrix_species = rownames(dist_mod) # list species in the matrix
common_species = intersect(tree_species, matrix_species) # list common species
# Prune tree
pruned_tree = drop.tip(tree, tree_species[!(tree_species %in% common_species)])
plotTree(pruned_tree,node.numbers=T,fsize=0.35)
# Prune dist_mod
dist_mod = dist_mod[common_species, common_species]
dim(dist_mod)
# Calculate phylogenetic distance matrix
dist_phylo = cophenetic.phylo (pruned_tree)
# Rows and columns of dist_mod are not in the same order of dist_phylo and thus we must correct it
sorted_rows <- dist_phylo[order(rownames(dist_phylo)), ] # Sort rows in alphabetical order
dist_phylo <- sorted_rows[, order(colnames(sorted_rows))] # Sort columns in alphabetical order
dim(dist_phylo)

# MANTEL TEST
mantel_modVSphy = mantel(dist_mod, dist_phylo, method="spearman", permutations = 10000)
mantel_modVSphy
densityplot(permustats(mantel_modVSphy))

#### OK CONSERVATIVE MATRIX: PHYLOGENY VS LAYERS (CONTROL: NONE) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(ape)
library(bipartite)
library(dplyr)
library(igraph)
library(reshape2)
library(vegan)

# LAYER DISTANCE MATRIX
data_layer = read_tsv("data/layers_CF")
data_layer = dcast(data_layer, Species ~ Layer, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_layer) <- data_layer$Species # turn the first column into rownames
data_layer <- data_layer[, -1] # remove first column
dist_layer = vegdist (data_layer, diag=TRUE, upper=T, method="jaccard", binary=T)

# PHYLOGENETIC DISTANCE MATRIX
# Load phylogeny
tree = read.nexus("data/tree_edited.nex")
plotTree(tree,node.numbers=T,fsize=0.35)
# Find which species are common in dist_layer and phylogeny
tree_species = tree$tip.label # list species in the tree
dist_layer = as.matrix(dist_layer) # convert dist to matrix (rownames does not work properly to dist objects)
matrix_species = rownames(dist_layer) # list species in the matrix
common_species = intersect(tree_species, matrix_species) # list common species
# Prune tree
pruned_tree = drop.tip(tree, tree_species[!(tree_species %in% common_species)])
plotTree(pruned_tree,node.numbers=T,fsize=0.35)
# Prune dist_layer
dist_layer = dist_layer[common_species, common_species]
# Calculate phylogenetic distance matrix
dist_phylo = cophenetic.phylo (pruned_tree)
# Rows and columns of dist_mod are not in the same order of dist_phylo and thus we must correct it
sorted_rows <- dist_phylo[order(rownames(dist_phylo)), ] # Sort rows in alphabetical order
dist_phylo <- sorted_rows[, order(colnames(sorted_rows))] # Sort columns in alphabetical order

# MANTEL TEST
mantel_layVSphy = mantel(dist_layer, dist_phylo, method="spearman", permutations = 10000)
mantel_layVSphy
densityplot(permustats(mantel_layVSphy))


#### OK CONSERVATIVE MATRIX: MODULE VS INTERACTIONS (CONTROL: NONE) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(bipartite)
library(igraph)
library(reshape2)
library(vegan)

# INTERACTION DISTANCE MATRIX
# Load edge and node lists
nodes = read.delim("data/net1nodes_CF.txt", header = T) # node traits
links = read.delim("data/list_CF", header = T) # list of interactions
# Load network
all_agreggated <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
all_agreggated = get.adjacency(all_agreggated,sparse=FALSE)
all_agreggated = all_agreggated[c(0:62), ] # keep only rows of bees
all_agreggated = all_agreggated[, 63:ncol(all_agreggated)] # keep only columns of plants
data_int <- empty(all_agreggated, count = F) # remove nodes with no interactions
# Create a distance matrix using Jaccard
dist_interactions = vegdist (data_int, diag=TRUE, upper=T, method="jaccard", binary=T)
dim(dist_interactions) # nrows and ncolumns

# MODULE DISTANCE MATRIX
# Read module dataframe
data_mod = read_tsv("data/partitions_cf.txt")
data_mod = data_mod[c(0:62), ] # keep only rows of bees
data_mod = data_mod[, 1:2] # keep only columns of plants
data_mod = dcast(data_mod, nodes ~ module, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_mod) <- data_mod$nodes
data_mod <- data_mod[, -1]
# Create a distance matrix using Jaccard
dist_mod = vegdist (data_mod, diag=TRUE, upper=T, method="jaccard", binary=T)
dim(dist_mod)

# MANTEL TEST
mantel_intVSmod = mantel(dist_interactions, dist_mod, method="spearman", permutations = 10000)
mantel_intVSmod
densityplot(permustats(mantel_intVSmod))

#### OK CONSERVATIVE MATRIX: MODULES VS LAYERS (CONTROL: NONE) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(bipartite)
library(igraph)
library(reshape2)
library(vegan)

# MODULE DISTANCE MATRIX
# Read module dataframe
data_mod = read_tsv("data/partitions_cf.txt")
data_mod = data_mod[c(0:62), ] # keep only rows of bees
data_mod = data_mod[, 1:2] # keep only columns of plants
data_mod = dcast(data_mod, nodes ~ module, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_mod) <- data_mod$nodes
data_mod <- data_mod[, -1]
# Create a distance matrix using Jaccard
dist_mod = vegdist (data_mod, diag=TRUE, upper=T, method="jaccard", binary=T)
dim(dist_mod)

# LAYER DISTANCE MATRIX
data_layer = read_tsv("data/layers_CF")
data_layer = dcast(data_layer, Species ~ Layer, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_layer) <- data_layer$Species # turn the first column into rownames
data_layer <- data_layer[, -1] # remove first column
dist_layer = vegdist (data_layer, diag=TRUE, upper=T, method="jaccard", binary=T)
dim(dist_layer)

# MANTEL TEST
mantel_modVSlay = mantel(dist_mod, dist_layer, method="spearman", permutations = 10000)
mantel_modVSlay
densityplot(permustats(mantel_modVSlay))

#### OK CONSERVATIVE MATRIX: INTERACTIONS VS LAYERS (CONTROL: NONE) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 
#Load the required packages and functions.
library(bipartite)
library(igraph)
library(reshape2)
library(vegan)

# INTERACTION DISTANCE MATRIX
# Load edge and node lists
nodes = read.delim("data/net1nodes_CF.txt", header = T) # node traits
links = read.delim("data/list_CF", header = T) # list of interactions
# Load network
all_agreggated <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
all_agreggated = get.adjacency(all_agreggated,sparse=FALSE)
all_agreggated = all_agreggated[c(0:62), ] # keep only rows of bees
all_agreggated = all_agreggated[, 63:ncol(all_agreggated)] # keep only columns of plants
data_int <- empty(all_agreggated, count = F) # remove nodes with no interactions
# Create a distance matrix using Jaccard
dist_interactions = vegdist (data_int, diag=TRUE, upper=T, method="jaccard", binary=T)
dim(dist_interactions) # nrows and ncolumns

# LAYER DISTANCE MATRIX
data_layer = read_tsv("data/layers_CF")
data_layer = dcast(data_layer, Species ~ Layer, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_layer) <- data_layer$Species # turn the first column into rownames
data_layer <- data_layer[, -1] # remove first column
dist_layer = vegdist (data_layer, diag=TRUE, upper=T, method="jaccard", binary=T)
dim(dist_layer)

# MANTEL TEST
mantel_intVSlay = mantel(dist_interactions, dist_layer, method="spearman", permutations = 10000)
mantel_intVSlay
densityplot(permustats(mantel_intVSlay))

#### OK CONSERVATIVE MATRIX: PHYLOGENY VS INTERACTIONS (CONTROL: LAYERS) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(ape)
library(bipartite)
library(dplyr)
library(igraph)
library(reshape2)
library(vegan)

# INTERACTION DISTANCE MATRIX
# Load edge and node lists
nodes = read.delim("data/net1nodes_CF.txt", header = T) # node traits
links = read.delim("data/list_CF", header = T) # list of interactions
# Load network
all_agreggated <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
all_agreggated = get.adjacency(all_agreggated,sparse=FALSE)
all_agreggated = all_agreggated[c(0:62), ] # keep only rows of bees
all_agreggated = all_agreggated[, 63:ncol(all_agreggated)] # keep only columns of plants
data_int <- empty(all_agreggated, count = F) # remove nodes with no interactions
# Create a distance matrix using Jaccard
dist_interactions = vegdist (data_int, diag=TRUE, upper=T, method="jaccard", binary=T)
dim(dist_interactions) # nrows and ncolumns

# PHYLOGENETIC DISTANCE MATRIX
# Load phylogeny
tree = read.nexus("data/tree_edited.nex")
plotTree(tree,node.numbers=T,fsize=0.35)
# Find which species are common in dist_interactions and phylogeny
tree_species = tree$tip.label # list species in the tree
dist_interactions = as.matrix(dist_interactions) # convert dist to matrix (rownames does not work properly to dist objects)
matrix_species = rownames(dist_interactions) # list species in the matrix
common_species = intersect(tree_species, matrix_species) # list common species
# Prune tree
pruned_tree = drop.tip(tree, tree_species[!(tree_species %in% common_species)])
plotTree(pruned_tree,node.numbers=T,fsize=0.35)
# Prune dist_interactions
dist_interactions = dist_interactions[common_species, common_species]
dim(dist_interactions) # nrows and ncolumns
# Calculate phylogenetic distance matrix
dist_phylo = cophenetic.phylo (pruned_tree)
# Rows and columns of dist_interactions are not in the same order of dist_phylo and thus we must correct it
sorted_rows <- dist_phylo[order(rownames(dist_phylo)), ] # Sort rows in alphabetical order
dist_phylo <- sorted_rows[, order(colnames(sorted_rows))] # Sort columns in alphabetical order
dim(dist_phylo) # nrows and ncolumns

# LAYER DISTANCE MATRIX
data_layer = read_tsv("data/layers_CF")
data_layer = dcast(data_layer, Species ~ Layer, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_layer) <- data_layer$Species # turn the first column into rownames
data_layer <- data_layer[, -1] # remove first column
# Find which species are common between dist_interactions and data_layer
matrix_species1 = rownames(data_layer) # list species in the matrix
matrix_species2 = rownames(dist_interactions) # list species in the matrix
common_species = intersect(matrix_species1, matrix_species2) # list common species
# Prune data_layer
data_layer = data_layer[matrix_species1 %in% common_species, ]
dist_layer = vegdist (data_layer, diag=TRUE, upper=T, method="jaccard", binary=T)
dim(dist_layer)

# PARTIAL MANTEL TEST
mantel_phyVSintCOlay = mantel.partial(dist_phylo, dist_interactions, dist_layer, method = "pearson", permutations = 10000)
mantel_phyVSintCOlay

#### OK CONSERVATIVE MATRIX: PHYLOGENY VS MODULES (CONTROL: LAYERS) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(ape)
library(bipartite)
library(dplyr)
library(igraph)
library(reshape2)
library(vegan)

# MODULE DISTANCE MATRIX
# Read module dataframe
data_mod = read_tsv("data/partitions_cf.txt")
data_mod = data_mod[c(0:62), ] # keep only rows of bees
data_mod = data_mod[, 1:2] # keep only columns of plants
data_mod = dcast(data_mod, nodes ~ module, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_mod) <- data_mod$nodes
data_mod <- data_mod[, -1]
# Create a distance matrix using Jaccard
dist_mod = vegdist (data_mod, diag=TRUE, upper=T, method="jaccard", binary=T)
dim(dist_mod)

# PHYLOGENETIC DISTANCE MATRIX
# Load phylogeny
tree = read.nexus("data/tree_edited.nex")
plotTree(tree,node.numbers=T,fsize=0.35)
# Find which species are common in dist_mod and phylogeny
tree_species = tree$tip.label # list species in the tree
dist_mod = as.matrix(dist_mod) # convert dist to matrix (rownames does not work properly to dist objects)
matrix_species = rownames(dist_mod) # list species in the matrix
common_species = intersect(tree_species, matrix_species) # list common species
# Prune tree
pruned_tree = drop.tip(tree, tree_species[!(tree_species %in% common_species)])
plotTree(pruned_tree,node.numbers=T,fsize=0.35)
# Prune dist_mod
dist_mod = dist_mod[common_species, common_species]
dim(dist_mod)
# Calculate phylogenetic distance matrix
dist_phylo = cophenetic.phylo (pruned_tree)
# Rows and columns of dist_mod are not in the same order of dist_phylo and thus we must correct it
sorted_rows <- dist_phylo[order(rownames(dist_phylo)), ] # Sort rows in alphabetical order
dist_phylo <- sorted_rows[, order(colnames(sorted_rows))] # Sort columns in alphabetical order
dim(dist_phylo)

# LAYER DISTANCE MATRIX
data_layer = read_tsv("data/layers_CF")
data_layer = dcast(data_layer, Species ~ Layer, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_layer) <- data_layer$Species # turn the first column into rownames
data_layer <- data_layer[, -1] # remove first column
# Find which species are common between dist_interactions and data_layer
matrix_species1 = rownames(data_layer) # list species in the matrix
matrix_species2 = rownames(dist_mod) # list species in the matrix
common_species = intersect(matrix_species1, matrix_species2) # list common species
# Prune data_layer
data_layer = data_layer[matrix_species1 %in% common_species, ]
dist_layer = vegdist (data_layer, diag=TRUE, upper=T, method="jaccard", binary=T)
dim(dist_layer)

# PARTIAL MANTEL TEST
mantel_phyVSmodCOlay = mantel.partial(dist_phylo, dist_mod, dist_layer, method = "pearson", permutations = 10000)
mantel_phyVSmodCOlay


#### OK CONSERVATIVE MATRIX: MODULES VS INTERACTIONS (CONTROL: PHYLOGENY) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(ape)
library(bipartite)
library(dplyr)
library(igraph)
library(reshape2)
library(vegan)

# MODULE DISTANCE MATRIX
# Read module dataframe
data_mod = read_tsv("data/partitions_cf.txt")
data_mod = data_mod[c(0:62), ] # keep only rows of bees
data_mod = data_mod[, 1:2] # keep only columns of plants
data_mod = dcast(data_mod, nodes ~ module, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_mod) <- data_mod$nodes
data_mod <- data_mod[, -1]
# Create a distance matrix using Jaccard
dist_mod = vegdist (data_mod, diag=TRUE, upper=T, method="jaccard", binary=T)
dim(dist_mod)

# PHYLOGENETIC DISTANCE MATRIX
# Load phylogeny
tree = read.nexus("data/tree_edited.nex")
plotTree(tree,node.numbers=T,fsize=0.35)
# Find which species are common in dist_mod and phylogeny
tree_species = tree$tip.label # list species in the tree
dist_mod = as.matrix(dist_mod) # convert dist to matrix (rownames does not work properly to dist objects)
matrix_species = rownames(dist_mod) # list species in the matrix
common_species = intersect(tree_species, matrix_species) # list common species
# Prune tree
pruned_tree = drop.tip(tree, tree_species[!(tree_species %in% common_species)])
plotTree(pruned_tree,node.numbers=T,fsize=0.35)
# Prune dist_mod
dist_mod = dist_mod[common_species, common_species]
dim(dist_mod)
# Calculate phylogenetic distance matrix
dist_phylo = cophenetic.phylo (pruned_tree)
# Rows and columns of dist_mod are not in the same order of dist_phylo and thus we must correct it
sorted_rows <- dist_phylo[order(rownames(dist_phylo)), ] # Sort rows in alphabetical order
dist_phylo <- sorted_rows[, order(colnames(sorted_rows))] # Sort columns in alphabetical order
dim(dist_phylo)

# INTERACTION DISTANCE MATRIX
# Load edge and node lists
nodes = read.delim("data/net1nodes_CF.txt", header = T) # node traits
links = read.delim("data/list_CF", header = T) # list of interactions
# Load network
all_agreggated <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
all_agreggated = get.adjacency(all_agreggated,sparse=FALSE)
all_agreggated = all_agreggated[c(0:62), ] # keep only rows of bees
all_agreggated = all_agreggated[, 63:ncol(all_agreggated)] # keep only columns of plants
data_int <- empty(all_agreggated, count = F) # remove nodes with no interactions
# Find which species are common between dist_interactions and data_layer
matrix_species1 = rownames(dist_mod) # list species in the matrix
matrix_species2 = rownames(data_int) # list species in the matrix
common_species = intersect(matrix_species1, matrix_species2) # list common species
# Prune data_layer
data_int = data_int[matrix_species2 %in% common_species, ]
# Create a distance matrix using Jaccard
dist_int = vegdist (data_int, diag=TRUE, upper=T, method="jaccard", binary=T)
dim(dist_int)

# PARTIAL MANTEL TEST
mantel_modVSintCOphy = mantel.partial(dist_int, dist_mod, dist_phylo, method = "pearson", permutations = 10000)
mantel_modVSintCOphy

#### OK CONSERVATIVE MATRIX: MODULES VS LAYERS (CONTROL: PHYLOGENY) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(ape)
library(bipartite)
library(dplyr)
library(igraph)
library(reshape2)
library(vegan)

# MODULE DISTANCE MATRIX
# Read module dataframe
data_mod = read_tsv("data/partitions_cf.txt")
data_mod = data_mod[c(0:62), ] # keep only rows of bees
data_mod = data_mod[, 1:2] # keep only columns of plants
data_mod = dcast(data_mod, nodes ~ module, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_mod) <- data_mod$nodes
data_mod <- data_mod[, -1]
# Create a distance matrix using Jaccard
dist_mod = vegdist (data_mod, diag=TRUE, upper=T, method="jaccard", binary=T)
dim(dist_mod)

# PHYLOGENETIC DISTANCE MATRIX
# Load phylogeny
tree = read.nexus("data/tree_edited.nex")
plotTree(tree,node.numbers=T,fsize=0.35)
# Find which species are common in dist_mod and phylogeny
tree_species = tree$tip.label # list species in the tree
dist_mod = as.matrix(dist_mod) # convert dist to matrix (rownames does not work properly to dist objects)
matrix_species = rownames(dist_mod) # list species in the matrix
common_species = intersect(tree_species, matrix_species) # list common species
# Prune tree
pruned_tree = drop.tip(tree, tree_species[!(tree_species %in% common_species)])
plotTree(pruned_tree,node.numbers=T,fsize=0.35)
# Prune dist_mod
dist_mod = dist_mod[common_species, common_species]
dim(dist_mod)
# Calculate phylogenetic distance matrix
dist_phylo = cophenetic.phylo (pruned_tree)
# Rows and columns of dist_mod are not in the same order of dist_phylo and thus we must correct it
sorted_rows <- dist_phylo[order(rownames(dist_phylo)), ] # Sort rows in alphabetical order
dist_phylo <- sorted_rows[, order(colnames(sorted_rows))] # Sort columns in alphabetical order
dim(dist_phylo)

# LAYER DISTANCE MATRIX
data_layer = read_tsv("data/layers_CF")
data_layer = dcast(data_layer, Species ~ Layer, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_layer) <- data_layer$Species # turn the first column into rownames
data_layer <- data_layer[, -1] # remove first column
# Find which species are common between dist_mod and data_layer
matrix_species1 = rownames(data_layer) # list species in the matrix
matrix_species2 = rownames(dist_mod) # list species in the matrix
common_species = intersect(matrix_species1, matrix_species2) # list common species
# Prune data_layer
data_layer = data_layer[matrix_species1 %in% common_species, ]
dist_layer = vegdist (data_layer, diag=TRUE, upper=T, method="jaccard", binary=T)
dim(dist_layer)

# PARTIAL MANTEL TEST
mantel_modVSlayCOphy = mantel.partial(dist_mod, dist_layer, dist_phylo, method = "pearson", permutations = 10000)
mantel_modVSlayCOphy

#### OK CONSERVATIVE MATRIX: PHYLOGENY VS INTERACTIONS (CONTROL: MODULES) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(ape)
library(bipartite)
library(dplyr)
library(igraph)
library(reshape2)
library(vegan)

# MODULE DISTANCE MATRIX
# Read module dataframe
data_mod = read_tsv("data/partitions_cf.txt")
data_mod = data_mod[c(0:62), ] # keep only rows of bees
data_mod = data_mod[, 1:2] # keep only columns of plants
data_mod = dcast(data_mod, nodes ~ module, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_mod) <- data_mod$nodes
data_mod <- data_mod[, -1]
# Create a distance matrix using Jaccard
dist_mod = vegdist (data_mod, diag=TRUE, upper=T, method="jaccard", binary=T)
dim(dist_mod)

# PHYLOGENETIC DISTANCE MATRIX
# Load phylogeny
tree = read.nexus("data/tree_edited.nex")
plotTree(tree,node.numbers=T,fsize=0.35)
# Find which species are common in dist_mod and phylogeny
tree_species = tree$tip.label # list species in the tree
dist_mod = as.matrix(dist_mod) # convert dist to matrix (rownames does not work properly to dist objects)
matrix_species = rownames(dist_mod) # list species in the matrix
common_species = intersect(tree_species, matrix_species) # list common species
# Prune tree
pruned_tree = drop.tip(tree, tree_species[!(tree_species %in% common_species)])
plotTree(pruned_tree,node.numbers=T,fsize=0.35)
# Prune dist_mod
dist_mod = dist_mod[common_species, common_species]
dim(dist_mod)
# Calculate phylogenetic distance matrix
dist_phylo = cophenetic.phylo (pruned_tree)
# Rows and columns of dist_mod are not in the same order of dist_phylo and thus we must correct it
sorted_rows <- dist_phylo[order(rownames(dist_phylo)), ] # Sort rows in alphabetical order
dist_phylo <- sorted_rows[, order(colnames(sorted_rows))] # Sort columns in alphabetical order
dim(dist_phylo)

# INTERACTION DISTANCE MATRIX
# Load edge and node lists
nodes = read.delim("data/net1nodes_CF.txt", header = T) # node traits
links = read.delim("data/list_CF", header = T) # list of interactions
# Load network
all_agreggated <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
all_agreggated = get.adjacency(all_agreggated,sparse=FALSE)
all_agreggated = all_agreggated[c(0:62), ] # keep only rows of bees
all_agreggated = all_agreggated[, 63:ncol(all_agreggated)] # keep only columns of plants
data_int <- empty(all_agreggated, count = F) # remove nodes with no interactions
# Find which species are common between dist_interactions and data_layer
matrix_species1 = rownames(dist_mod) # list species in the matrix
matrix_species2 = rownames(data_int) # list species in the matrix
common_species = intersect(matrix_species1, matrix_species2) # list common species
# Prune data_layer
data_int = data_int[matrix_species2 %in% common_species, ]
# Create a distance matrix using Jaccard
dist_int = vegdist (data_int, diag=TRUE, upper=T, method="jaccard", binary=T)

# PARTIAL MANTEL TEST
mantel_phyVSintCOmod = mantel.partial(dist_phylo, dist_int, dist_mod, method = "pearson", permutations = 10000)
mantel_phyVSintCOmod

#### CONSERVATIVE MATRIX: PHYLOGENY VS LAYERS (CONTROL: MODULES) ####
#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(ape)
library(bipartite)
library(dplyr)
library(igraph)
library(reshape2)
library(vegan)

# MODULE DISTANCE MATRIX
# Read module dataframe
data_mod = read_tsv("data/partitions_cf.txt")
data_mod = data_mod[c(0:62), ] # keep only rows of bees
data_mod = data_mod[, 1:2] # keep only columns of plants
data_mod = dcast(data_mod, nodes ~ module, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_mod) <- data_mod$nodes
data_mod <- data_mod[, -1]
# Create a distance matrix using Jaccard
dist_mod = vegdist (data_mod, diag=TRUE, upper=T, method="jaccard", binary=T)
dim(dist_mod)

# PHYLOGENETIC DISTANCE MATRIX
# Load phylogeny
tree = read.nexus("data/tree_edited.nex")
plotTree(tree,node.numbers=T,fsize=0.35)
# Find which species are common in dist_mod and phylogeny
tree_species = tree$tip.label # list species in the tree
dist_mod = as.matrix(dist_mod) # convert dist to matrix (rownames does not work properly to dist objects)
matrix_species = rownames(dist_mod) # list species in the matrix
common_species = intersect(tree_species, matrix_species) # list common species
# Prune tree
pruned_tree = drop.tip(tree, tree_species[!(tree_species %in% common_species)])
plotTree(pruned_tree,node.numbers=T,fsize=0.35)
# Prune dist_mod
dist_mod = dist_mod[common_species, common_species]
dim(dist_mod)
# Calculate phylogenetic distance matrix
dist_phylo = cophenetic.phylo (pruned_tree)
# Rows and columns of dist_mod are not in the same order of dist_phylo and thus we must correct it
sorted_rows <- dist_phylo[order(rownames(dist_phylo)), ] # Sort rows in alphabetical order
dist_phylo <- sorted_rows[, order(colnames(sorted_rows))] # Sort columns in alphabetical order
dim(dist_phylo)

# LAYER DISTANCE MATRIX
data_layer = read_tsv("data/layers_CF")
data_layer = dcast(data_layer, Species ~ Layer, fun.aggregate = length) # create a binary adjacency matrix of modules
rownames(data_layer) <- data_layer$Species # turn the first column into rownames
data_layer <- data_layer[, -1] # remove first column
# Find which species are common between dist_mod and data_layer
matrix_species1 = rownames(data_layer) # list species in the matrix
matrix_species2 = rownames(dist_mod) # list species in the matrix
common_species = intersect(matrix_species1, matrix_species2) # list common species
# Prune data_layer
data_layer = data_layer[matrix_species1 %in% common_species, ]
dist_layer = vegdist (data_layer, diag=TRUE, upper=T, method="jaccard", binary=T)
dim(dist_layer)

# PARTIAL MANTEL TEST
mantel_phyVSlayCOmod = mantel.partial(dist_phylo, dist_layer, dist_mod, method = "pearson", permutations = 10000)
mantel_phyVSlayCOmod
