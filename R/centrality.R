# Centrality ~ Body size

#### 1. BIPARTITE/IGRAPH ESTIMATION ####
#### 1.1 TOTAL NETWORK ####
# Packages
library(bipartite)
library(igraph)

# Load edge and node lists (TOTAL: AGREGGATED)
nodes = read.delim("data/net1nodes_.txt", header = T) # node traits
links = read.delim("data/list.3MBSDs_v3_binary", header = T) # list of interactions
# Load network
all_agreggated <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
all_agreggated = get.adjacency(all_agreggated,sparse=FALSE)
all_agreggated = all_agreggated[c(0:68), ] # keep only rows of bees
all_agreggated = all_agreggated[, 69:ncol(all_agreggated)] # keep only columns of plants
data_tol_ag <- empty(all_agreggated, count = F) # remove nodes with no interactions

# Calculate degree, specialization, betweenness, closeness
tol_ag_degree = bipartite::specieslevel(data_tol_ag, index="degree") # degree
tol_ag_ndegree = bipartite::specieslevel(data_tol_ag, index="normalised degree") # normalized degree
tol_ag_specialization = bipartite::specieslevel(data_tol_ag, index="d") # Bluthgen specialization d
tol_ag_betweenness = bipartite::specieslevel(data_tol_ag, index = "betweenness") # Betweenness
tol_ag_closeness = bipartite::specieslevel(data_tol_ag, index = "closeness") # Closeness
# Calculate within-module degree
Mod <- bipartite::computeModules(data_tol_ag)
tol_ag_withinModuleDegree = bipartite::czvalues(Mod, weighted = F, level = "lower")
# Calculate eigenvector
nodes = read.delim("data/net1nodes_.txt", header = T) # node traits
links = read.delim("data/list.2MBSDs_v3_binary", header = T) # list of interactions
tol_ag_igraph <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
tol_ag_eigen <- igraph::eigen_centrality(tol_ag_igraph, directed = F, scale = T)$vector
tol_ag_eigen = tol_ag_eigen[1:68]

# Merge all bipartite centrality values into a dataframe
df_bipartite_tol = cbind(tol_ag_degree$`lower level`, 
                     tol_ag_ndegree$`lower level`, 
                     tol_ag_specialization$`lower level`, 
                     tol_ag_betweenness$`lower level`, 
                     tol_ag_closeness$`lower level`, 
                     tol_ag_withinModuleDegree$z, 
                     tol_ag_eigen)
# Remove disconnected subnetworks: A.ferruginea, F.silvestri, M.beechei, M. orbynigi, T.williana
df_bipartite_tol = df_bipartite_tol[-c(1, 3,20,27,67), ] 
# Rename columns
colnames(df_bipartite_tol) = c("degr", "ndegr", "d", "betw", "wbetw", "clos", "wclos", "withinMdegr", "eigen")
View(df_bipartite_tol)

#### 1.2 CONSERVATIVE NETWORK ####
# Packages
library(bipartite)
library(igraph)

# Load edge and node lists (CONSERVATIVE: AGGREGATED)
nodes = read.delim("data/net1nodes_CF.txt", header = T) # node traits
links = read.delim("data/list_CF", header = T) # list of interactions
# Load network
all_agreggated <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
all_agreggated = get.adjacency(all_agreggated,sparse=FALSE)
all_agreggated = all_agreggated[c(0:62), ] # keep only rows of bees
all_agreggated = all_agreggated[, 63:ncol(all_agreggated)] # keep only columns of plants
data_con_ag <- empty(all_agreggated, count = F) # remove nodes with no interactions
# Calculate degree, specialization, betweenness, closeness
con_ag_degree = bipartite::specieslevel(data_con_ag, index="degree") # degree
con_ag_ndegree = bipartite::specieslevel(data_con_ag, index="normalised degree") # normalized degree
con_ag_specialization = bipartite::specieslevel(data_con_ag, index="d") # Bluthgen specialization d
con_ag_betweenness = bipartite::specieslevel(data_con_ag, index = "betweenness") # Betweenness
con_ag_closeness = bipartite::specieslevel(data_con_ag, index = "closeness") # Closeness
# Calculate within-module degree
Mod <- bipartite::computeModules(data_con_ag)
con_ag_withinModuleDegree = bipartite::czvalues(Mod, weighted = F, level = "lower")
# Calculate eigenvector
nodes = read.delim("data/net1nodes_CF.txt", header = T) # node traits
links = read.delim("data/list_CF", header = T) # list of interactions
con_ag_igraph <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
con_ag_eigen <- igraph::eigen_centrality(con_ag_igraph, directed = F, scale = T)$vector
con_ag_eigen = con_ag_eigen[1:62]

# Merge all bipartite centrality values into a dataframe
df_bipartite_con = cbind(con_ag_degree$`lower level`, 
                     con_ag_ndegree$`lower level`, 
                     con_ag_specialization$`lower level`, 
                     con_ag_betweenness$`lower level`, 
                     con_ag_closeness$`lower level`, 
                     con_ag_withinModuleDegree$z, 
                     con_ag_eigen)
# Remove disconnected subnetworks: A.ferruginea, F.silvestri, M.beechei, M. orbynigi, T.williana
df_bipartite_con = df_bipartite_con[-c(1, 3,19,24,61), ] 
# Rename columns
colnames(df_bipartite_con) = c("degr", "ndegr", "d", "betw", "wbetw", "clos", "wclos", "withinMdegr", "eigen")
View(df_bipartite_con)

#### 2. EMLN ESTIMATION ####
#### 2.1.1 TOTAL MATRIX: DATA INPUT ####

#Load the required packages and functions.
library(bipartite)
library(igraph)
library(emln)
library(tidyverse)
library(infomapecology)
library(magrittr)

# Load edge and node lists (ASIA)
nodes = read.delim("data/net1nodes_all_asia.txt", header = T) # node traits
links = read.delim("data/list.2MBSDs_asia_binary", header = T) # list of interactions
# Load network
all_asia <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
all_asia = get.adjacency(all_asia,sparse=FALSE)
all_asia = all_asia[c(0:31), ] # keep only rows of bees
all_asia = all_asia[, 32:ncol(all_asia)] # keep only columns of plants
asia <- empty(all_asia, count = F) # remove nodes with no interactions

# Load edge and node lists (NEOTROPICAL)
nodes = read.delim("data/net1nodes_all_neotropical.txt", header = T) # node traits
links = read.delim("data/list.2MBSDs_neotropical_binary", header = T) # list of interactions
# Load network
all_neotrop <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
all_neotrop = get.adjacency(all_neotrop,sparse=FALSE)
all_neotrop = all_neotrop[c(0:36), ] # keep only rows of bees
all_neotrop = all_neotrop[, 37:ncol(all_neotrop)] # keep only columns of plants
neot <- empty(all_neotrop, count = F) # remove nodes with no interactions

# Specify the layers
layer_attributes = tibble(layer_id=1:2, layer_name=c('Indo-Malayan-Australasian',
                                                     'Neotropics'))
# Create interlinks (plants present in more than one layer)
interlinks = read_tsv("data/interlinks_CF.txt")
interlayer <- tibble(layer_from=interlinks[3],
                     node_from=interlinks[1],
                     layer_to=interlinks[4], 
                     node_to=interlinks[2],
                     weight=1)
#interlayer <- tibble(layer_from=c('Neotropics','Neotropics','Neotropics','Neotropics','Neotropics','Neotropics','Neotropics','Neotropics'),
#                     node_from=c('Anacardiaceae', 'Araucariaceae','Clusiaceae','Euphorbiaceae','Fabaceae','Moraceae','Myrtaceae','Pinaceae'),
#                     layer_to=c('Indo-Malayan-Australasia','Indo-Malayan-Australasia','Indo-Malayan-Australasia','Indo-Malayan-Australasia','Indo-Malayan-Australasia','Indo-Malayan-Australasia','Indo-Malayan-Australasia','Indo-Malayan-Australasia'), 
#                     node_to=c('Anacardiaceae', 'Araucariaceae','Clusiaceae','Euphorbiaceae','Fabaceae','Moraceae','Myrtaceae','Pinaceae'),
#                     weight=1)

# Merge layers into a multilayer network
mult_resin = create_multilayer_network(list_of_layers = list(asia,neot),
                                             interlayer_links = NULL,
                                             layer_attributes = layer_attributes,
                                             bipartite = T,
                                             directed = F)

#### 2.1.2 TOTAL MATRIX: DEGREE ####

# MONOLAYER
g_layer <- get_igraph(mult_resin, bipartite = T, directed = F) # get a list of igraph objects as layers igraph 
# Estimate centrality for each layer separately
degr_igraph = c()
for (layer in g_layer$layers_igraph){
  degr_igraph <- append(degr_igraph, igraph::degree(layer, normalized = F))
}
length(degr_igraph) # get the number of nodes

# MULTILAYER
# Get the SAM
sam_multilayer <- get_sam(multilayer = mult_resin, bipartite = T, directed = F, sparse = F, remove_zero_rows_cols = F)
# Converting SAM to igraph object
graph_sam <- graph_from_adjacency_matrix(sam_multilayer$M,mode = "undirected",weighted = NULL,diag = TRUE,add.colnames = NULL,add.rownames = NA)
# Get centrality considering multilayer structure
degr_sam <- igraph::degree(graph_sam, normalized = F)
length(degr_sam) # get the number of nodes

# DATASET
ec_comparison <- data.frame(degr_sam, degr_igraph) # create a df
colnames(ec_comparison) <- c("multilayer","monolayer") # change column names from df
row_names = names(degr_igraph) # extract the names of species
append_continent <- function(x) { # create a function to avoid duplicate names in row_names
  duplicated_names <- duplicated(x)
  if (any(duplicated_names)) {
    x[duplicated_names] <- paste0(x[duplicated_names], "_neot")
    x[!duplicated_names] <- paste0(x[!duplicated_names], "_asia")
  } else {
    x <- paste0(x, "_asia")
  }
  return(x)
}
species_names_with_continent <- append_continent(row_names) # append asia in the first occurrence of a string, neot in the second
rownames(ec_comparison) = species_names_with_continent # change the row names from df
degr_data = ec_comparison

# PLOT CORRELATION BETWEEN MONOLAYER AND MULTILAYER CENTRALITY VALUES
summary(lm(multilayer~monolayer, data=degr_data)) # highly correlated because there are few interlinks
ggplot(data = ec_comparison, aes(monolayer, multilayer)) + 
  geom_point(color = "blue", size = 0.75) + 
  labs(title = "Degree Centrality (EC) comparison", x = "EC (monolayer)", y = "EC (multilayer)")+
  geom_abline()+
  #scale_x_continuous(limits = c(0,1))+
  #scale_y_continuous(limits = c(0,1))+
  coord_fixed() +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=15),
        axis.text = element_text(color='black',size = 10),
        legend.text =  element_text(size=15),
        legend.title = element_text(size=20))

#### 2.1.3 TOTAL MATRIX: BETWEENESS ####

# MONOLAYER
g_layer <- get_igraph(mult_resin, bipartite = T, directed = F) # get a list of igraph objects as layers igraph 
# Estimate centrality for each layer separately
betw_igraph = c()
for (layer in g_layer$layers_igraph){
  betw_igraph <- append(betw_igraph, igraph::betweenness(layer, directed=F))
}
length(betw_igraph) # get the number of nodes

# MULTILAYER
# Get the SAM
sam_multilayer <- get_sam(multilayer = mult_resin, bipartite = T, directed = F, sparse = F, remove_zero_rows_cols = F)
# Converting SAM to igraph object
graph_sam <- graph_from_adjacency_matrix(sam_multilayer$M,mode = "undirected",weighted = NULL,diag = TRUE,add.colnames = NULL,add.rownames = NA)
# Get centrality considering multilayer structure
betw_sam <- igraph::betweenness(graph_sam, directed = F)
length(betw_sam) # get the number of nodes

# DATASET
ec_comparison <- data.frame(betw_sam, betw_igraph) # create a df
colnames(ec_comparison) <- c("multilayer","monolayer") # change column names from df
row_names = names(betw_igraph) # extract the names of species
append_continent <- function(x) { # create a function to avoid duplicate names in row_names
  duplicated_names <- duplicated(x)
  if (any(duplicated_names)) {
    x[duplicated_names] <- paste0(x[duplicated_names], "_neot")
    x[!duplicated_names] <- paste0(x[!duplicated_names], "_asia")
  } else {
    x <- paste0(x, "_asia")
  }
  return(x)
}
species_names_with_continent <- append_continent(row_names) # append asia in the first occurrence of a string, neot in the second
rownames(ec_comparison) = species_names_with_continent # change the row names from df
betw_data = ec_comparison

# PLOT CORRELATION BETWEEN MONOLAYER AND MULTILAYER CENTRALITY VALUES
summary(lm(multilayer~monolayer, data=betw_data)) # highly correlated because there are few interlinks
ggplot(data = ec_comparison, aes(monolayer, multilayer)) + 
  geom_point(color = "blue", size = 0.75) + 
  labs(title = "Betweeness Centrality (EC) comparison", x = "EC (monolayer)", y = "EC (multilayer)")+
  geom_abline()+
  #scale_x_continuous(limits = c(0,1))+
  #scale_y_continuous(limits = c(0,1))+
  coord_fixed() +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=15),
        axis.text = element_text(color='black',size = 10),
        legend.text =  element_text(size=15),
        legend.title = element_text(size=20))

#### 2.1.4 TOTAL MATRIX: CLOSENESS ####

# MONOLAYER
g_layer <- get_igraph(mult_resin, bipartite = T, directed = F) # get a list of igraph objects as layers igraph 
# Estimate centrality for each layer separately
clos_igraph = c()
for (layer in g_layer$layers_igraph){
  clos_igraph <- append(clos_igraph, igraph::closeness(layer))
}
length(clos_igraph) # get the number of nodes

# MULTILAYER
# Get the SAM
sam_multilayer <- get_sam(multilayer = mult_resin, bipartite = T, directed = F, sparse = F, remove_zero_rows_cols = F)
# Converting SAM to igraph object
graph_sam <- graph_from_adjacency_matrix(sam_multilayer$M,mode = "undirected",weighted = NULL,diag = TRUE,add.colnames = NULL,add.rownames = NA)
# Get centrality considering multilayer structure
clos_sam <- igraph::closeness(layer)
length(clos_sam) # get the number of nodes

# DATASET
ec_comparison <- data.frame(clos_sam, clos_igraph) # create a df
colnames(ec_comparison) <- c("multilayer","monolayer") # change column names from df
row_names = names(betw_igraph) # extract the names of species
append_continent <- function(x) { # create a function to avoid duplicate names in row_names
  duplicated_names <- duplicated(x)
  if (any(duplicated_names)) {
    x[duplicated_names] <- paste0(x[duplicated_names], "_neot")
    x[!duplicated_names] <- paste0(x[!duplicated_names], "_asia")
  } else {
    x <- paste0(x, "_asia")
  }
  return(x)
}
species_names_with_continent <- append_continent(row_names) # append asia in the first occurrence of a string, neot in the second
rownames(ec_comparison) = species_names_with_continent # change the row names from df
clos_data = ec_comparison

# PLOT CORRELATION BETWEEN MONOLAYER AND MULTILAYER CENTRALITY VALUES
summary(lm(multilayer~monolayer, data=clos_data)) # highly correlated because there are few interlinks
ggplot(data = clos_data, aes(monolayer, multilayer)) + 
  geom_point(color = "blue", size = 0.75) + 
  labs(title = "Closeness Centrality (EC) comparison", x = "EC (monolayer)", y = "EC (multilayer)")+
  geom_abline()+
  #scale_x_continuous(limits = c(0,1))+
  #scale_y_continuous(limits = c(0,1))+
  coord_fixed() +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=15),
        axis.text = element_text(color='black',size = 10),
        legend.text =  element_text(size=15),
        legend.title = element_text(size=20))

#### 2.1.5 TOTAL MATRIX: EIGENVECTOR ####

# MONOLAYER
g_layer <- get_igraph(mult_resin, bipartite = T, directed = F) # get a list of igraph objects as layers igraph 
# Estimate centrality for each layer separately
eigen_igraph = c()
for (layer in g_layer$layers_igraph){
  eigen_igraph <- append(eigen_igraph, igraph::eigen_centrality(layer, directed = F, scale = T)$vector)
}
length(eigen_igraph) # get the number of nodes

# MULTILAYER
# Get the SAM
sam_multilayer <- get_sam(multilayer = mult_resin, bipartite = F, directed = T, sparse = F, remove_zero_rows_cols = F)
# Converting SAM to igraph object
graph_sam <- graph_from_adjacency_matrix(sam_multilayer$M,mode = "undirected",weighted = NULL,diag = TRUE,add.colnames = NULL,add.rownames = NA)
# Get Eigenvector scores
eigen_sam <- igraph::eigen_centrality(graph_sam, directed = F, scale = T)$vector
length(eigen_sam) # get the number of nodes

# DATASET
ec_comparison <- data.frame(eigen_sam, eigen_igraph) # create a df
colnames(ec_comparison) <- c("multilayer","monolayer") # change column names from df
row_names = names(eigen_igraph) # extract the names of species
append_continent <- function(x) { # create a function to avoid duplicate names in row_names
  duplicated_names <- duplicated(x)
  if (any(duplicated_names)) {
    x[duplicated_names] <- paste0(x[duplicated_names], "_neot")
    x[!duplicated_names] <- paste0(x[!duplicated_names], "_asia")
  } else {
    x <- paste0(x, "_asia")
  }
  return(x)
}
species_names_with_continent <- append_continent(row_names) # append asia in the first occurrence of a string, neot in the second
rownames(ec_comparison) = species_names_with_continent # change the row names from df
eigen_data = ec_comparison

# PLOT CORRELATION BETWEEN MONOLAYER AND MULTILAYER CENTRALITY VALUES
summary(lm(multilayer~monolayer, data=eigen_data)) # highly correlated because there are few interlinks
ggplot(data = eigen_data, aes(monolayer, multilayer)) + 
  geom_point(color = "blue", size = 0.75) + 
  labs(title = "Eigenvector Centrality (EC) comparison", x = "EC (monolayer)", y = "EC (multilayer)")+
  geom_abline()+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  coord_fixed() +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=15),
        axis.text = element_text(color='black',size = 10),
        legend.text =  element_text(size=15),
        legend.title = element_text(size=20))

#### 2.1.6 TOTAL MATRIX: ADD VALUES: df_emln ####
df_emln = cbind(degr_data, betw_data, clos_data, eigen_data) # bind all EMLN centrality metrics
colnames(df_emln) = c("degr_mono", "degr_mult", "betw_mono", "betw_mult", "clos_mono", "clos_mult", "eigen_mono", "eigen_mult") # rename column names
df_emln <- df_emln[grep("_.*_.*", rownames(df_emln)), ] # select only bees and exclude plants
df_emln <- df_emln[df_emln[, 1] != 0, ] # select only rows with degree > 0 (i.e. exclude imaginary neotropical bees occurring in asia, and vice-versa)
rownames(df_emln) <- sub("_asia$|_neot$", "", rownames(df_emln)) # remove _asia and _neot from all rownames
df_emln <- df_emln[order(rownames(df_emln)), ] # sort rownames in alphabetical order
# remove disconnected subnetworks: F.silvestri, M.beechei, M. orbynigi, T.williana
df_emln_tol = df_emln[-c(2,19,26,66), ] 

#### 2.2.1 CONSERVATIVE MATRIX: DATA INPUT ####

#Load the required packages and functions.
library(bipartite)
library(igraph)
library(emln)
library(tidyverse)
library(infomapecology)
library(magrittr)

# Load edge and node lists (ASIA)
nodes = read.delim("data/net1nodes_cf_asia.txt", header = T) # node traits
links = read.delim("data/list_CF_asia", header = T) # list of interactions
# Load network
all_asia <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
all_asia = get.adjacency(all_asia,sparse=FALSE)
all_asia = all_asia[c(0:31), ] # keep only rows of bees
all_asia = all_asia[, 32:ncol(all_asia)] # keep only columns of plants
asia <- empty(all_asia, count = F) # remove nodes with no interactions

# Load edge and node lists (NEOTROPICAL)
nodes = read.delim("data/net1nodes_CF_neotropical.txt", header = T) # node traits
links = read.delim("data/list_CF_neotropical", header = T) # list of interactions
# Load network
all_neotrop <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
all_neotrop = get.adjacency(all_neotrop,sparse=FALSE)
all_neotrop = all_neotrop[c(0:30), ] # keep only rows of bees
all_neotrop = all_neotrop[, 31:ncol(all_neotrop)] # keep only columns of plants
neot <- empty(all_neotrop, count = F) # remove nodes with no interactions

# Specify the layers
layer_attributes = tibble(layer_id=1:2, layer_name=c('Indo-Malayan-Australasian',
                                                     'Neotropics'))
# Create interlinks (plants present in more than one layer)
#interlinks = read_tsv("data/interlinks_CF.txt")
#interlayer <- tibble(layer_from=interlinks[3],
#                     node_from=interlinks[1],
#                     layer_to=interlinks[4], 
#                     node_to=interlinks[2],
#                     weight=1)
#interlayer <- tibble(layer_from=c('Neotropics','Neotropics','Neotropics','Neotropics','Neotropics','Neotropics','Neotropics','Neotropics'),
#                     node_from=c('Anacardiaceae', 'Araucariaceae','Clusiaceae','Euphorbiaceae','Fabaceae','Moraceae','Myrtaceae','Pinaceae'),
#                     layer_to=c('Indo-Malayan-Australasia','Indo-Malayan-Australasia','Indo-Malayan-Australasia','Indo-Malayan-Australasia','Indo-Malayan-Australasia','Indo-Malayan-Australasia','Indo-Malayan-Australasia','Indo-Malayan-Australasia'), 
#                     node_to=c('Anacardiaceae', 'Araucariaceae','Clusiaceae','Euphorbiaceae','Fabaceae','Moraceae','Myrtaceae','Pinaceae'),
#                     weight=1)

# Merge layers into a multilayer network
mult_resin = create_multilayer_network(list_of_layers = list(asia,neot),
                                       interlayer_links = NULL,
                                       layer_attributes = layer_attributes,
                                       bipartite = T,
                                       directed = F)

#### 2.2.2 CONSERVATIVE MATRIX: DEGREE ####

# MONOLAYER
g_layer <- get_igraph(mult_resin, bipartite = T, directed = F) # get a list of igraph objects as layers igraph 
# Estimate centrality for each layer separately
degr_igraph = c()
for (layer in g_layer$layers_igraph){
  degr_igraph <- append(degr_igraph, igraph::degree(layer, normalized = F))
}
length(degr_igraph) # get the number of nodes

# MULTILAYER
# Get the SAM
sam_multilayer <- get_sam(multilayer = mult_resin, bipartite = T, directed = F, sparse = F, remove_zero_rows_cols = F)
# Converting SAM to igraph object
graph_sam <- graph_from_adjacency_matrix(sam_multilayer$M,mode = "undirected",weighted = NULL,diag = TRUE,add.colnames = NULL,add.rownames = NA)
# Get centrality considering multilayer structure
degr_sam <- igraph::degree(graph_sam, normalized = F)
length(degr_sam) # get the number of nodes

# DATASET
ec_comparison <- data.frame(degr_sam, degr_igraph) # create a df
colnames(ec_comparison) <- c("multilayer","monolayer") # change column names from df
row_names = names(degr_igraph) # extract the names of species
append_continent <- function(x) { # create a function to avoid duplicate names in row_names
  duplicated_names <- duplicated(x)
  if (any(duplicated_names)) {
    x[duplicated_names] <- paste0(x[duplicated_names], "_neot")
    x[!duplicated_names] <- paste0(x[!duplicated_names], "_asia")
  } else {
    x <- paste0(x, "_asia")
  }
  return(x)
}
species_names_with_continent <- append_continent(row_names) # append asia in the first occurrence of a string, neot in the second
rownames(ec_comparison) = species_names_with_continent # change the row names from df
degr_data = ec_comparison

# PLOT CORRELATION BETWEEN MONOLAYER AND MULTILAYER CENTRALITY VALUES
summary(lm(multilayer~monolayer, data=degr_data)) # highly correlated because there are few interlinks
ggplot(data = ec_comparison, aes(monolayer, multilayer)) + 
  geom_point(color = "blue", size = 0.75) + 
  labs(title = "Degree Centrality (EC) comparison", x = "EC (monolayer)", y = "EC (multilayer)")+
  geom_abline()+
  #scale_x_continuous(limits = c(0,1))+
  #scale_y_continuous(limits = c(0,1))+
  coord_fixed() +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=15),
        axis.text = element_text(color='black',size = 10),
        legend.text =  element_text(size=15),
        legend.title = element_text(size=20))

#### 2.22.3 CONSERVATIVE MATRIX: BETWEENESS ####

# MONOLAYER
g_layer <- get_igraph(mult_resin, bipartite = T, directed = F) # get a list of igraph objects as layers igraph 
# Estimate centrality for each layer separately
betw_igraph = c()
for (layer in g_layer$layers_igraph){
  betw_igraph <- append(betw_igraph, igraph::betweenness(layer, directed=F))
}
length(betw_igraph) # get the number of nodes

# MULTILAYER
# Get the SAM
sam_multilayer <- get_sam(multilayer = mult_resin, bipartite = T, directed = F, sparse = F, remove_zero_rows_cols = F)
# Converting SAM to igraph object
graph_sam <- graph_from_adjacency_matrix(sam_multilayer$M,mode = "undirected",weighted = NULL,diag = TRUE,add.colnames = NULL,add.rownames = NA)
# Get centrality considering multilayer structure
betw_sam <- igraph::betweenness(graph_sam, directed = F)
length(betw_sam) # get the number of nodes

# DATASET
ec_comparison <- data.frame(betw_sam, betw_igraph) # create a df
colnames(ec_comparison) <- c("multilayer","monolayer") # change column names from df
row_names = names(betw_igraph) # extract the names of species
append_continent <- function(x) { # create a function to avoid duplicate names in row_names
  duplicated_names <- duplicated(x)
  if (any(duplicated_names)) {
    x[duplicated_names] <- paste0(x[duplicated_names], "_neot")
    x[!duplicated_names] <- paste0(x[!duplicated_names], "_asia")
  } else {
    x <- paste0(x, "_asia")
  }
  return(x)
}
species_names_with_continent <- append_continent(row_names) # append asia in the first occurrence of a string, neot in the second
rownames(ec_comparison) = species_names_with_continent # change the row names from df
betw_data = ec_comparison

# PLOT CORRELATION BETWEEN MONOLAYER AND MULTILAYER CENTRALITY VALUES
summary(lm(multilayer~monolayer, data=betw_data)) # highly correlated because there are few interlinks
ggplot(data = ec_comparison, aes(monolayer, multilayer)) + 
  geom_point(color = "blue", size = 0.75) + 
  labs(title = "Betweeness Centrality (EC) comparison", x = "EC (monolayer)", y = "EC (multilayer)")+
  geom_abline()+
  #scale_x_continuous(limits = c(0,1))+
  #scale_y_continuous(limits = c(0,1))+
  coord_fixed() +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=15),
        axis.text = element_text(color='black',size = 10),
        legend.text =  element_text(size=15),
        legend.title = element_text(size=20))

#### 2.2.4 CONSERVATIVE MATRIX: CLOSENESS ####

# MONOLAYER
g_layer <- get_igraph(mult_resin, bipartite = T, directed = F) # get a list of igraph objects as layers igraph 
# Estimate centrality for each layer separately
clos_igraph = c()
for (layer in g_layer$layers_igraph){
  clos_igraph <- append(clos_igraph, igraph::closeness(layer))
}
length(clos_igraph) # get the number of nodes

# MULTILAYER
# Get the SAM
sam_multilayer <- get_sam(multilayer = mult_resin, bipartite = T, directed = F, sparse = F, remove_zero_rows_cols = F)
# Converting SAM to igraph object
graph_sam <- graph_from_adjacency_matrix(sam_multilayer$M,mode = "undirected",weighted = NULL,diag = TRUE,add.colnames = NULL,add.rownames = NA)
# Get centrality considering multilayer structure
clos_sam <- igraph::closeness(layer)
length(clos_sam) # get the number of nodes

# DATASET
ec_comparison <- data.frame(clos_sam, clos_igraph) # create a df
colnames(ec_comparison) <- c("multilayer","monolayer") # change column names from df
row_names = names(betw_igraph) # extract the names of species
append_continent <- function(x) { # create a function to avoid duplicate names in row_names
  duplicated_names <- duplicated(x)
  if (any(duplicated_names)) {
    x[duplicated_names] <- paste0(x[duplicated_names], "_neot")
    x[!duplicated_names] <- paste0(x[!duplicated_names], "_asia")
  } else {
    x <- paste0(x, "_asia")
  }
  return(x)
}
species_names_with_continent <- append_continent(row_names) # append asia in the first occurrence of a string, neot in the second
rownames(ec_comparison) = species_names_with_continent # change the row names from df
clos_data = ec_comparison

# PLOT CORRELATION BETWEEN MONOLAYER AND MULTILAYER CENTRALITY VALUES
summary(lm(multilayer~monolayer, data=clos_data)) # highly correlated because there are few interlinks
ggplot(data = clos_data, aes(monolayer, multilayer)) + 
  geom_point(color = "blue", size = 0.75) + 
  labs(title = "Closeness Centrality (EC) comparison", x = "EC (monolayer)", y = "EC (multilayer)")+
  geom_abline()+
  #scale_x_continuous(limits = c(0,1))+
  #scale_y_continuous(limits = c(0,1))+
  coord_fixed() +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=15),
        axis.text = element_text(color='black',size = 10),
        legend.text =  element_text(size=15),
        legend.title = element_text(size=20))

#### 2.2.5 CONSERVATIVE MATRIX: EIGENVECTOR ####

# MONOLAYER
g_layer <- get_igraph(mult_resin, bipartite = T, directed = F) # get a list of igraph objects as layers igraph 
# Estimate centrality for each layer separately
eigen_igraph = c()
for (layer in g_layer$layers_igraph){
  eigen_igraph <- append(eigen_igraph, igraph::eigen_centrality(layer, directed = F, scale = T)$vector)
}
length(eigen_igraph) # get the number of nodes

# MULTILAYER
# Get the SAM
sam_multilayer <- get_sam(multilayer = mult_resin, bipartite = F, directed = T, sparse = F, remove_zero_rows_cols = F)
# Converting SAM to igraph object
graph_sam <- graph_from_adjacency_matrix(sam_multilayer$M,mode = "undirected",weighted = NULL,diag = TRUE,add.colnames = NULL,add.rownames = NA)
# Get Eigenvector scores
eigen_sam <- igraph::eigen_centrality(graph_sam, directed = F, scale = T)$vector
length(eigen_sam) # get the number of nodes

# DATASET
ec_comparison <- data.frame(eigen_sam, eigen_igraph) # create a df
colnames(ec_comparison) <- c("multilayer","monolayer") # change column names from df
row_names = names(eigen_igraph) # extract the names of species
append_continent <- function(x) { # create a function to avoid duplicate names in row_names
  duplicated_names <- duplicated(x)
  if (any(duplicated_names)) {
    x[duplicated_names] <- paste0(x[duplicated_names], "_neot")
    x[!duplicated_names] <- paste0(x[!duplicated_names], "_asia")
  } else {
    x <- paste0(x, "_asia")
  }
  return(x)
}
species_names_with_continent <- append_continent(row_names) # append asia in the first occurrence of a string, neot in the second
rownames(ec_comparison) = species_names_with_continent # change the row names from df
eigen_data = ec_comparison

# PLOT CORRELATION BETWEEN MONOLAYER AND MULTILAYER CENTRALITY VALUES
summary(lm(multilayer~monolayer, data=eigen_data)) # highly correlated because there are few interlinks
ggplot(data = eigen_data, aes(monolayer, multilayer)) + 
  geom_point(color = "blue", size = 0.75) + 
  labs(title = "Eigenvector Centrality (EC) comparison", x = "EC (monolayer)", y = "EC (multilayer)")+
  geom_abline()+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  coord_fixed() +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=15),
        axis.text = element_text(color='black',size = 10),
        legend.text =  element_text(size=15),
        legend.title = element_text(size=20))

#### 2.2.6 CONSERVATIVE MATRIX: ADD VALUES: df_emln ####
df_emln = cbind(degr_data, betw_data, clos_data, eigen_data) # bind all EMLN centrality metrics
colnames(df_emln) = c("degr_mono", "degr_mult", "betw_mono", "betw_mult", "clos_mono", "clos_mult", "eigen_mono", "eigen_mult") # rename column names
df_emln <- df_emln[grep("_.*_.*", rownames(df_emln)), ] # select only bees and exclude plants
df_emln <- df_emln[df_emln[, 1] != 0, ] # select only rows with degree > 0 (i.e. exclude imaginary neotropical bees occurring in asia, and vice-versa)
rownames(df_emln) <- sub("_asia$|_neot$", "", rownames(df_emln)) # remove _asia and _neot from all rownames
df_emln <- df_emln[order(rownames(df_emln)), ] # sort rownames in alphabetical order
# remove disconnected subnetworks: F.silvestri, M.beechei, M. orbynigi, T.williana
df_emln_cons = df_emln[-c(2,18,23,60), ] 
