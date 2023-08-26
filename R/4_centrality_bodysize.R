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

#### 2.2.1 [to be done] CONSERVATIVE MATRIX: DATA INPUT

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

#### [to be done] 3. MULTINET ####
#### 3.1 DEGREE ####
#### 3.2 NEIGHBORHOOD ####
#### 4. GLMMs ####
#### 4.1 TOTAL ####
library(ape)
library(ggcorrplot)
library(glmmTMB)
library(lme4)
library(phytools)

# Load body size
traits = as.data.frame(read_tsv("data/traits.tsv"))
rownames(traits) <- traits[, 1] # turn the first column into rownames
traits <- traits[, -1] # Remove the first column
traits$module <- as.factor(traits$module) # convert module column to categorical data
dim(traits)

# Load centrality (TOTAL)
df_total = cbind (df_emln_tol, df_bipartite_tol) # merge bipartite and emln estimates
View(df_total)
# Plot pairwise correlations between centrality metrics
df_removeNA <- df_total[complete.cases(df_total), ]
cor_matrix <- cor(df_removeNA)
ggcorrplot(cor_matrix, type = "lower", hc.order = TRUE, lab = TRUE) # plot 
# Merge centrality and body size
data_GLM = (merge(df_total, traits, by = "row.names", all = TRUE))
data_GLM <- data_GLM[!is.na(data_GLM$eigen_mono), ]
data_GLM_tol <- data_GLM[!is.na(data_GLM$body_size_mean), ]
write.table(data_GLM_tol, file = 'data/data_GLM_tol.tsv', sep = "\t", quote = FALSE, row.names = T)

# Load phylogeny
phy_tree = read.nexus("data/tree_edited.nex")
plotTree(phy_tree,node.numbers=T,fsize=0.35) # plot original tree
matching_species <- intersect(phy_tree$tip.label, rownames(df_emln)) # Find which species are common in dist_mod and phylogeny
pruned_phy_tree <- drop.tip(phy_tree, tip = setdiff(phy_tree$tip.label, matching_species))
plotTree(pruned_phy_tree,node.numbers=T,fsize=0.35) # plot pruned tree
df_pruned <- df_emln %>% # prune dataset
  filter(rownames(df_emln) %in% matching_species)
phy_corr <- ape::vcv.phylo(pruned_phy_tree) # convert phylogeny to covariance matrix
class(phy_corr)

# LM
summary(lm(eigen_mult~(body_size_mean), data=data_GLM_tol)) # p = 0.00961
plot(log(data_GLM_tol$body_size_mean), data_GLM_tol$eigen_mult)

# GLM
# Total
glm1_tol = glm(eigen_mult~body_size_mean, data=data_GLM_tol, family = quasipoisson())
summary(glm1_tol)

# GLMM
glmm1_tol = glmmTMB(eigen_mult~body_size_mean + (1|biogeographical_region), 
                data=data_GLM_tol
                )
summary(glmm1_tol)

# FIGURE 5
data_GLM_tol$predicted <- predict(glm1_tol, newdata = data_GLM_tol, type = "response")
data_GLM_tol$se <- sqrt(diag(vcov(glm1_tol)))
predictions_tol <- data.frame(predict(glm1_tol, type = "response", interval = "confidence"))

ggplot(data_GLM_tol, aes(x = body_size_mean, y = eigen_mono)) +
  #geom_point() +
  geom_smooth(method = "glm", method.args = list(family = 'quasipoisson'), se=F, alpha = 0.95) +
  #geom_line(aes(y = predicted)) +  # Add the tendency line
  #geom_errorbar(aes(ymin = predicted - 1.96 * se, ymax = predicted + 1.96 * se), width = 0.2) +  # Add 95% confidence intervals
  labs(x = "\nBody size (mm)", y = "Centrality\n") +
  theme_classic()

png(filename= "figures/fig5B_TOTAL.png", res= 300, height= 2500, width= 3000)
ggplot(data_GLM_tol, aes(x = body_size_mean, y = eigen_mono)) +
  #geom_point() +
  geom_line(aes(y = predicted)) +  # Add the tendency line
  geom_ribbon(aes(ymin = predicted - 1.96*se, ymax = predicted + 1.96*se), alpha = 0.2) +
  labs(x = "\nBody size (mm)", y = "Eigenvector centrality\n") +
  theme_classic()+
  theme(axis.text = element_text(size = 17),
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17))
  #ylim(-.1, .8)
dev.off()

png(filename= "figures/fig5B_TOTAL_pointsNOse.png", res= 300, height= 2500, width= 3000)
ggplot(data_GLM_tol, aes(x = body_size_mean, y = eigen_mono)) +
  #geom_point() +
  geom_line(aes(y = predicted)) +  # Add the tendency line
  #geom_ribbon(aes(ymin = predicted - 1.96*se, ymax = predicted + 1.96*se), alpha = 0.2) +
  labs(x = "\nBody size (mm)", y = "Eigenvector centrality\n") +
  theme_classic()+
  theme(axis.text = element_text(size = 17),
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17))
dev.off()

#### 4.2 GLMM CONSERVATIVE ####
# Load centrality (CONSERVATIVE)
df_conservative = cbind(df_emln_cons, df_bipartite_con)
data_GLM = (merge(df_conservative, traits, by = "row.names", all = TRUE))
data_GLM <- data_GLM[!is.na(data_GLM$eigen_mono), ]
data_GLM_con <- data_GLM[!is.na(data_GLM$body_size_mean), ]
View(data_GLM_con)
write.table(data_GLM_con, file = 'data/data_GLM_con.tsv', sep = "\t", quote = FALSE, row.names = T)


# LM
summary(lm(data_GLM_con$clos_mult~log(data_GLM_con$body_size_mean))) # p = 0.0442
plot(log(data_GLM_con$body_size_mean), data_GLM_con$clos_mult)
summary(lm(data_GLM_con$eigen_mult~log(data_GLM_con$body_size_mean))) # p = 0.3143
plot(log(data_GLM_con$body_size_mean), data_GLM_con$eigen_mult)

# GLM
glm2_con = glm(clos_mult~body_size_mean, data=data_GLM_con, family=quasipoisson())
summary(glm2_con) # p = 0.0424
glm3_con = glm(eigen_mult~body_size_mean, data=data_GLM_con, family = quasipoisson())
summary(glm3_con) # p = 0.348

# GLMM
glmm2_con = glmmTMB(clos_mult~body_size_mean + (1|biogeographical_region),
                    data=data_GLM_con)
summary(glmm2_con) # p = 0.0315
glmm3_con = glmmTMB(eigen_mult~body_size_mean+(1|biogeographical_region),
                    data=data_GLM_con)
summary(glmm3_con) # p = 0.3190

# Predict
data_GLM_con$predicted <- predict(glm2_con, newdata = data_GLM_con, type = "response")
data_GLM_con$se <- sqrt(diag(vcov(glm2_con)))
predictions_con <- data.frame(predict(glm2_con, type = "response", interval = "confidence"))

ggplot(data_GLM_con, aes(x = body_size_mean, y = clos_mult)) +
  #geom_point() +
  geom_smooth(method = "glm", method.args = list(family = 'quasipoisson'), se=F, alpha = 0.95) +
  #geom_line(aes(y = predicted)) +  # Add the tendency line
  #geom_errorbar(aes(ymin = predicted - 1.96 * se, ymax = predicted + 1.96 * se), width = 0.2) +  # Add 95% confidence intervals
  labs(x = "\nBody size (mm)", y = "Centrality\n") +
  theme_classic()

ggplot(data_GLM_con, aes(x = body_size_mean, y = clos_mult)) +
  geom_point() +
  geom_line(aes(y = predicted)) +  # Add the tendency line
  #geom_ribbon(aes(ymin = predicted - 1.96*se, ymax = predicted + 1.96*se), alpha = 0.2) +
  labs(x = "\nBody size (mm)", y = "Closeness centrality\n") +
  theme_classic()+
  theme(axis.text = element_text(size = 17),
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17))

#### 4.3 COMBINED (FIGURE 5A) ####
# Load total
data_GLM = (merge(df_total, traits, by = "row.names", all = TRUE))
data_GLM <- data_GLM[!is.na(data_GLM$eigen_mono), ]
data_GLM_tol <- data_GLM[!is.na(data_GLM$body_size_mean), ]
data_GLM_tol$Approach <- "Total"
dim(data_GLM_tol)
# Load conservative
df_conservative = cbind(df_emln_cons, df_bipartite_con)
data_GLM = (merge(df_conservative, traits, by = "row.names", all = TRUE))
data_GLM <- data_GLM[!is.na(data_GLM$eigen_mono), ]
data_GLM_con <- data_GLM[!is.na(data_GLM$body_size_mean), ]
data_GLM_con$Approach <- "Conservative"
dim(data_GLM_con)
# Merge vertically
combined_data <- rbind(data_GLM_tol, data_GLM_con)
# GLM
glmm_comb1 = glmmTMB(clos_mult~body_size_mean:Approach + (1|biogeographical_region),
                    data=combined_data)
summary(glmm_comb1)
glmm_comb2 = glmmTMB(eigen_mult~body_size_mean:Approach + (1|biogeographical_region),
                     data=combined_data)
summary(glmm_comb2)
# Plot
png(filename= "figures/fig5A.png", res= 300, height= 2200, width= 3000)
ggplot(combined_data, aes(x = body_size_mean, y = eigen_mult, color=Approach)) +
  geom_point(alpha=.15, size= 3) +
  geom_smooth(method = 'glm', method.args = list(family='negative.binomial'(theta=1)), se=F) +
  #geom_line(aes(y = predicted)) +  # Add the tendency line
  #geom_errorbar(aes(ymin = predicted - 1.96 * se, ymax = predicted + 1.96 * se), width = 0.2) +  # Add 95% confidence intervals
  labs(x = "\nIntertegular distance (mm)", y = "Eigenvector centrality\n") +
  theme_classic()+
  scale_color_manual(values=c("yellow4", "black"))+
  theme(axis.text = element_text(size = 15), # values size
        axis.title.x = element_text(size = 20), # x title size
        axis.title.y = element_text(size = 20), # y title size
        legend.text = element_text(size = 14), # legend size
        legend.title = element_text(size = 14)) # legend title size
dev.off()

#### 5. Spider charts ####
library(fmsb)

#### 5.1 TOTAL: TOP 1:5 ####
spider_data_total = data_GLM_tol[c("Row.names", "degr_mult", "betw_mult", "clos_mult", "eigen_mult", "d", "withinMdegr")]
rownames(spider_data_total) <- spider_data_total[, 1]
spider_data_total <- spider_data_total[, -1]
spider_data_total = spider_data_total[order(-spider_data_total$eigen_mult), ]
spider_data_total <- spider_data_total[    #Criamos uma cópia dos resultados
  complete.cases(spider_data_total), ] %>% #Removemos todos os casos com NA
  scale() %>%                              #Escalonamos todos os valores
  as.data.frame() %>%                      #Convertemos para o formato de data.frame
  slice(1:5) %>%                           #Pegamos apenas os 5 primeiros casos
  rbind(apply(., 2, max, na.rm=T),         #Adicionamos uma linha com os valores máximos
        apply(., 2, min, na.rm=T),         #Adicionamos uma linha com os valores mínimos
        .)
row.names(spider_data_total)[1:2] <- c("max","min") #Nomeamos as linhas com valores max e min
new_column_names <- c("Degree", "Betweenness", "Closeness", "Eigenvector", "Specialization", "Module-degree") # Define the new column names
names(spider_data_total) <- new_column_names # Assign the new column names to the dataframe
cores.radar <- rainbow(nrow(spider_data_total))  #Cores a serem usadas no plot

png(filename= "figures/fig5B_total_1_5.png", res= 300, height= 2500, width= 4300)
radarchart(spider_data_total, 
           na.itp = T,
           axistype = 2, 
           maxmin = T,
           pcol = adjustcolor(cores.radar, alpha = 1),
           pfcol = adjustcolor(cores.radar, alpha = 0.2),
           plwd = 1 , 
           plty = 1,
           cglcol = "grey", 
           cglty = 1,
           axislabcol = F, 
           caxislabels = 10, 
           cglwd = 0.5,
           vlcex = 1)
legend(x=1.5, y=1.0, 
       legend = rownames(spider_data_total[-c(1,2),]),
       bty = "n", 
       pch = 20,
       col = adjustcolor(cores.radar, alpha = 0.2),
       text.col = "black",
       cex = 1.2, 
       pt.cex = 3, 
       title = "Species")
dev.off()

#### 5.2 TOTAL: TOP 6:10 ####
spider_data_total = data_GLM_tol[c("Row.names", "degr_mult", "betw_mult", "clos_mult", "eigen_mult", "d", "withinMdegr")]
rownames(spider_data_total) <- spider_data_total[, 1]
spider_data_total <- spider_data_total[, -1]
spider_data_total = spider_data_total[order(-spider_data_total$eigen_mult), ]
spider_data_total <- spider_data_total[    #Criamos uma cópia dos resultados
  complete.cases(spider_data_total), ] %>% #Removemos todos os casos com NA
  scale() %>%                              #Escalonamos todos os valores
  as.data.frame() %>%                      #Convertemos para o formato de data.frame
  slice(6:10) %>%                           #Pegamos apenas os 5 primeiros casos
  rbind(apply(., 2, max, na.rm=T),         #Adicionamos uma linha com os valores máximos
        apply(., 2, min, na.rm=T),         #Adicionamos uma linha com os valores mínimos
        .)
row.names(spider_data_total)[1:2] <- c("max","min") #Nomeamos as linhas com valores max e min
new_column_names <- c("Degree", "Betweenness", "Closeness", "Eigenvector", "Specialization", "Module-degree") # Define the new column names
names(spider_data_total) <- new_column_names # Assign the new column names to the dataframe
cores.radar <- rainbow(nrow(spider_data_total))  #Cores a serem usadas no plot

png(filename= "figures/fig5B_total_6_10.png", res= 300, height= 2500, width= 4300)
radarchart(spider_data_total, 
           na.itp = T,
           axistype = 2, 
           maxmin = T,
           pcol = adjustcolor(cores.radar, alpha = 1),
           pfcol = adjustcolor(cores.radar, alpha = 0.2),
           plwd = 1 , 
           plty = 1,
           cglcol = "grey", 
           cglty = 1,
           axislabcol = F, 
           caxislabels = 10, 
           cglwd = 0.5,
           vlcex = 1)
legend(x=1.5, y=1.0, 
       legend = rownames(spider_data_total[-c(1,2),]),
       bty = "n", 
       pch = 20,
       col = adjustcolor(cores.radar, alpha = 0.2),
       text.col = "black",
       cex = 1.2, 
       pt.cex = 3, 
       title = "Species")
dev.off()
#### 5.3 CONSERVATIVE: TOP 1:5 ####
spider_data_con = data_GLM_con[c("Row.names", "degr_mult", "betw_mult", "clos_mult", "eigen_mult", "d", "withinMdegr")]
rownames(spider_data_con) <- spider_data_con[, 1]
spider_data_con <- spider_data_con[, -1]
spider_data_con = spider_data_con[order(-spider_data_con$eigen_mult), ]
spider_data_con <- spider_data_con[    #Criamos uma cópia dos resultados
  complete.cases(spider_data_con), ] %>% #Removemos todos os casos com NA
  scale() %>%                              #Escalonamos todos os valores
  as.data.frame() %>%                      #Convertemos para o formato de data.frame
  slice(1:5) %>%                           #Pegamos apenas os 5 primeiros casos
  rbind(apply(., 2, max, na.rm=T),         #Adicionamos uma linha com os valores máximos
        apply(., 2, min, na.rm=T),         #Adicionamos uma linha com os valores mínimos
        .)
row.names(spider_data_con)[1:2] <- c("max","min") #Nomeamos as linhas com valores max e min
new_column_names <- c("Degree", "Betweenness", "Closeness", "Eigenvector", "Specialization", "Module-degree") # Define the new column names
names(spider_data_con) <- new_column_names # Assign the new column names to the dataframe
cores.radar <- rainbow(nrow(spider_data_con))  #Cores a serem usadas no plot

png(filename= "figures/fig5B_con_1_5.png", res= 300, height= 2500, width= 4300)
radarchart(spider_data_con, 
           na.itp = T,
           axistype = 2, 
           maxmin = T,
           pcol = adjustcolor(cores.radar, alpha = 1),
           pfcol = adjustcolor(cores.radar, alpha = 0.2),
           plwd = 1 , 
           plty = 1,
           cglcol = "grey", 
           cglty = 1,
           axislabcol = F, 
           caxislabels = 10, 
           cglwd = 0.5,
           vlcex = 1)
legend(x=1.5, y=1.0, 
       legend = rownames(spider_data_con[-c(1,2),]),
       bty = "n", 
       pch = 20,
       col = adjustcolor(cores.radar, alpha = 0.2),
       text.col = "black",
       cex = 1.2, 
       pt.cex = 3, 
       title = "Species")
dev.off()

#### 5.4 CONSERVATIVE TOP 6:10 ####
spider_data_con = data_GLM_con[c("Row.names", "degr_mult", "betw_mult", "clos_mult", "eigen_mult", "d", "withinMdegr")]
rownames(spider_data_con) <- spider_data_con[, 1]
spider_data_con <- spider_data_con[, -1]
spider_data_con = spider_data_con[order(-spider_data_con$eigen_mult), ]
spider_data_con <- spider_data_con[    #Criamos uma cópia dos resultados
  complete.cases(spider_data_con), ] %>% #Removemos todos os casos com NA
  scale() %>%                              #Escalonamos todos os valores
  as.data.frame() %>%                      #Convertemos para o formato de data.frame
  slice(6:10) %>%                           #Pegamos apenas os 5 primeiros casos
  rbind(apply(., 2, max, na.rm=T),         #Adicionamos uma linha com os valores máximos
        apply(., 2, min, na.rm=T),         #Adicionamos uma linha com os valores mínimos
        .)
row.names(spider_data_con)[1:2] <- c("max","min") #Nomeamos as linhas com valores max e min
new_column_names <- c("Degree", "Betweenness", "Closeness", "Eigenvector", "Specialization", "Module-degree") # Define the new column names
names(spider_data_con) <- new_column_names # Assign the new column names to the dataframe
cores.radar <- rainbow(nrow(spider_data_con))  #Cores a serem usadas no plot

png(filename= "figures/fig5B_con_6_10.png", res= 300, height= 2500, width= 4300)
radarchart(spider_data_con, 
           na.itp = T,
           axistype = 2, 
           maxmin = T,
           pcol = adjustcolor(cores.radar, alpha = 1),
           pfcol = adjustcolor(cores.radar, alpha = 0.2),
           plwd = 1 , 
           plty = 1,
           cglcol = "grey", 
           cglty = 1,
           axislabcol = F, 
           caxislabels = 10, 
           cglwd = 0.5,
           vlcex = 1)
legend(x=1.5, y=1.0, 
       legend = rownames(spider_data_con[-c(1,2),]),
       bty = "n", 
       pch = 20,
       col = adjustcolor(cores.radar, alpha = 0.2),
       text.col = "black",
       cex = 1.2, 
       pt.cex = 3, 
       title = "Species")
dev.off()