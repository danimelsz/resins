###########
# 1 SETUP #
###########

# 1.1 Setting directories

# Set directory where packages should be installed
#.libPaths('/home/danimelsz/Downloads/R/x86_64-pc-linux-gnu-library/4.1')

# Set directory
# setwd("~/Desktop/Stingless_Bees/pub_abelha/dados")
setwd("~/Downloads/meliponini")

# 1.2 Load Packages

# The packages 'car' and 'tidyverse' require the following libraries in Ubuntu: 
# $ sudo apt-get install libcurl4-openssl-dev
# $ sudo apt-get install -y libssl-dev
# $ sudo apt-get install libxml2-dev

library(ape)
library(bipartite)
library(brainGraph)
library(brms)
library(caper)
library(car)
library(corrgram)
library(devtools)
library(factoextra)
library(FactoMineR)
library(fmsb)
library(geiger)
library(ggplot2)
library(igraph)
library(iNEXT)
library(intergraph)
library(MCMCglmm)
library(phytools)
library(picante)
library(reshape2)
library(RInSp)
library(tidyverse)
library(gdata)
library("Matrix")
library("paco")
#library("tcltk2")
library("glm2")
library("lme4")
library("gridExtra")
library("dplyr")


###########################
# 2 SINGLELAYER ANALYSIS #
###########################

# A singlelayer analysis is required to estimate module membership of each node.

##################################
### 2.1 Singlelayer: All BMSDs ###

# Load edge and node lists
nodes = read.delim("data/net1nodes_.txt", header = T) # node traits
links = read.delim("data/list.3MBSDs_v3_binary", header = T) # list of interactions
# Inspect objects
head(nodes)
class(nodes)
head(links)
class(links)
unique(links$layer)
# Load network
singlelayer <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Inspect network
class(singlelayer)
singlelayer
V(singlelayer)$name
V(singlelayer)$taxon
# Create an adjacent matrix for bipartite
singlelayer = get.adjacency(singlelayer,sparse=FALSE)
singlelayer = singlelayer[-c(0:68), ]
singlelayer <- empty(singlelayer, count = F) # remove nodes with no interactions
# NODF
observado_NODF_g = networklevel(singlelayer, index = "NODF") # empirical NODF = 33.88934 
# Modules
modules = computeModules(singlelayer, method = "Beckett")
modules@likelihood # M = 0.3835738
# Compute module membership
modules.list = listModuleInformation(modules) # create empty list of modules
source("R/Extracting_modules_from_bipartite.R") # fuction for extracting modules
modules.dataframe = modules_from_bipartite(modules.list)
write.table(modules.dataframe$Rows_modules, "data/modulemembership.rows.plants.txt", sep="\t", quote=F) # export rows (plants)
write.table(modules.dataframe$Cols_modules, "data/modulemembership.cols.bees.txt", sep="\t", quote=F) # export columns (bees)

#################################
### 2.2 SINGLELAYER: CHEMICAL ###

# Load edge and node lists
nodes_c = read.delim("data/net1nodes_c.txt", header = T) # node traits
links_c = read.delim("data/list_chemical", header = T) # list of interactions
# Inspect objects
head(nodes_c)
class(nodes_c)
head(links_c)
class(links_c)
unique(links_c$layer)
# Load network
singlelayer_c <- graph_from_data_frame(d=links_c, vertices=nodes_c, directed=F)
# Inspect network
class(singlelayer_c)
singlelayer_c
V(singlelayer_c)$name
V(singlelayer_c)$taxon
# Create an adjacent matrix for bipartite
singlelayer_c = get.adjacency(singlelayer_c,sparse=FALSE)
singlelayer_c = singlelayer_c[-c(0:48), ] # remove rows
singlelayer_c = singlelayer_c[,-49:-77] # remove columns
singlelayer_c = empty(singlelayer_c, count = F) # remove nodes with no interactions
# NODF
observado_NODF_c = networklevel(singlelayer_c, index = "NODF") # empirical NODF = 16.34849
# Modules
modules_c = computeModules(singlelayer_c, method = "Beckett")
modules_c@likelihood # M = 0.6148095
# Compute module membership
modules.list = listModuleInformation(modules_c) # create empty list of modules
source("R/Extracting_modules_from_bipartite.R") # fuction for extracting modules
modules.dataframe = modules_from_bipartite(modules.list)
write.table(modules.dataframe$Rows_modules, "data/modulemembership.rows.plants_c.txt", sep="\t", quote=F) # export rows (plants)
write.table(modules.dataframe$Cols_modules, "data/modulemembership.cols.bees_c.txt", sep="\t", quote=F) # export columns (bees)

##################################
### 2.3 SINGLELAYER: FIELDWORK ###

# Load edge and node lists
nodes_f = read.delim("data/net1nodes_f.txt", header = T) # node traits
links_f = read.delim("data/list_fieldwork", header = T) # list of interactions
# Inspect objects
head(nodes_f)
class(nodes_f)
head(links_f)
class(links_f)
unique(links_f$layer)
# Load network
singlelayer_f <- graph_from_data_frame(d=links_f, vertices=nodes_f, directed=F)
# Inspect network
class(singlelayer_f)
singlelayer_f
V(singlelayer_f)$name
V(singlelayer_f)$taxon
# Create an adjacent matrix for bipartite
singlelayer_f = get.adjacency(singlelayer_f,sparse=FALSE)
singlelayer_f = singlelayer_f[-c(0:38), ] # remove rows
singlelayer_f = (singlelayer_f[,-39:-55]) # remove columns
singlelayer_f = empty(singlelayer_f, count = F) # remove nodes with no interactions
# NODF
observado_NODF_f = networklevel(singlelayer_f, index = "NODF") # empirical NODF = 22.45398
# Modules
modules_f = computeModules(singlelayer_f, method = "Beckett")
modules_f@likelihood # M = 0.5894409
# Compute module membership
modules.list = listModuleInformation(modules_f) # create empty list of modules
source("R/Extracting_modules_from_bipartite.R") # fuction for extracting modules
modules.dataframe = modules_from_bipartite(modules.list)
write.table(modules.dataframe$Rows_modules, "data/modulemembership.rows.plants_f.txt", sep="\t", quote=F) # export rows (plants)
write.table(modules.dataframe$Cols_modules, "data/modulemembership.cols.bees_f.txt", sep="\t", quote=F) # export columns (bees)


######################################
### 2.4 SINGLELAYER: PALYNOLOGICAL ###

# Load edge and node lists
nodes_p = read.delim("data/net1nodes_p.txt", header = T) # node traits
links_p = read.delim("data/list_palynological", header = T) # list of interactions
# Inspect objects
head(nodes_p)
class(nodes_p)
head(links_p)
class(links_p)
unique(links_p$layer)
# Load network
singlelayer_p <- graph_from_data_frame(d=links_p, vertices=nodes_p, directed=F)
# Inspect network
class(singlelayer_p)
singlelayer_p
V(singlelayer_p)$name
V(singlelayer_p)$taxon
# Create an adjacent matrix for bipartite
singlelayer_p = get.adjacency(singlelayer_p,sparse=FALSE)
singlelayer_p = singlelayer_p[-c(0:14), ] # remove rows
singlelayer_p = singlelayer_p[,-15:-100] # remove columns
singlelayer_p = empty(singlelayer_p, count = F) # remove nodes with no interactions
# NODF
observado_NODF_p = networklevel(singlelayer_p, index = "NODF") # empirical NODF = 47.66621 
# Modules
modules_p = computeModules(singlelayer_p, method = "Beckett")
modules_p@likelihood # M = 0.2545233
# Compute module membership
modules.list = listModuleInformation(modules_p) # create empty list of modules
source("R/Extracting_modules_from_bipartite.R") # fuction for extracting modules
modules.dataframe = modules_from_bipartite(modules.list)
write.table(modules.dataframe$Rows_modules, "data/modulemembership.rows.plants_p.txt", sep="\t", quote=F) # export rows (plants)
write.table(modules.dataframe$Cols_modules, "data/modulemembership.cols.bees_p.txt", sep="\t", quote=F) # export columns (bees)


######################################
### 2.5 SINGLELAYER: CHEMICAL + FIELDWORK ###

# Load edge and node lists
nodes_cf = read.delim("data/net1nodes_CF.txt", header = T) # node traits
links_cf = read.delim("data/list_CF", header = T) # list of interactions
# Inspect objects
head(nodes_cf)
class(nodes_cf)
head(links_cf)
class(links_cf)
unique(links_cf$layer)
# Load network
singlelayer_cf <- graph_from_data_frame(d=links_cf, vertices=nodes_cf, directed=F)
# Inspect network
class(singlelayer_cf)
singlelayer_cf
V(singlelayer_cf)$name
V(singlelayer_cf)$taxon
# Create an adjacent matrix for bipartite
singlelayer_cf = get.adjacency(singlelayer_cf,sparse=FALSE)
singlelayer_cf = singlelayer_cf[c(0:62), ] # keep only rows of bees
singlelayer_cf = singlelayer_cf[, 63:ncol(singlelayer_cf)] # keep only columns of plants
singlelayer_cf = empty(singlelayer_cf, count = F) # remove nodes with no interactions
dim(singlelayer_cf)
# NODF
observado_NODF_cf = networklevel(singlelayer_cf, index = "NODF") # empirical NODF = 47.66621 
# Modules
modules_cf = computeModules(singlelayer_cf, method = "Beckett")
modules_cf@likelihood # M = 0.56
# Compute module membership
modules.list = listModuleInformation(modules_cf) # create empty list of modules
source("R/Extracting_modules_from_bipartite.R") # fuction for extracting modules
modules.dataframe = modules_from_bipartite(modules.list)
write.table(modules.dataframe$Rows_modules, "data/modulemembership.rows.plants_cf.txt", sep="\t", quote=F) # export rows (plants)
write.table(modules.dataframe$Cols_modules, "data/modulemembership.cols.bees_cf.txt", sep="\t", quote=F) # export columns (bees)

#########################
# 3 MULTILAYER ANALYSIS #
#########################

#################################
### 3.1 MULTILAYER: ALL MBSDs ###

# Create network
net1_multi = graph_from_data_frame(d=links, vertices=nodes, directed=F) 

#Check the multilayer network
net1_multi
attributes(V(net1_multi))
attributes(E(net1_multi))
V(net1_multi)$name # node names
V(net1_multi)$taxon # plant or bee
V(net1_multi)$type # 1 = bee, 2 = plant
E(net1_multi)$layer # Neotropics or Indo-Malayan-Australasia or Afrotropics
E(net1_multi)$weight # 1 for all

# Import module membership
modules=read.table("data/partitions_all.txt", h=T)
head(modules)
unique(modules$module) # check the number of modules (= 9)

# Draw the network
# Set the same layout for all graphs
set.seed(16381)
l <- layout_nicely(net1_multi)
# Set NODE shape for bipartite structure (bee = square; plant = circle)
V(net1_multi)$shape = V(net1_multi)$taxon
V(net1_multi)$shape = gsub("bee","square",V(net1_multi)$shape)
V(net1_multi)$shape = gsub("plant","circle",V(net1_multi)$shape)
# Set NODE colors for the bipartite structure (bee = brown; plant = green) 
#V(net1_multi)$color = V(net1_multi)$taxon
#V(net1_multi)$color = gsub("bee","#D3802B",V(net1_multi)$color)
#V(net1_multi)$color = gsub("plant","#2C8437",V(net1_multi)$color)
# Set EDGE colors for the layers (biogeographic regions)
E(net1_multi)$color = E(net1_multi)$layer
E(net1_multi)$color = gsub("Neotropics","#B3DCF9",E(net1_multi)$color)
E(net1_multi)$color = gsub("Indo-Malayan-Australasia","#FFD26B",E(net1_multi)$color)
E(net1_multi)$color = gsub("Afrotropics","#A151A1",E(net1_multi)$color)
# Set polarity
E(net1_multi)$arrow.mode = 0
# Set EDGE width
E(net1_multi)$width = E(net1_multi)$layer
E(net1_multi)$width = gsub("Neotropics", 6, E(net1_multi)$width)
E(net1_multi)$width = gsub("Indo-Malayan-Australasia", 1, E(net1_multi)$width)
E(net1_multi)$width = gsub("Afrotropics", 1, E(net1_multi)$width)
# Set NODE colors for modules (= 10 modules)
vertex_name_order = data.frame(nodes=V(net1_multi)$name)
modules_order = merge(vertex_name_order, modules, by="nodes", sort = F)
# colrs = rainbow(9, alpha=1) #rainbow palette with 9 colors
colrs = c("#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#B41572", "#FF7F00", "#CAB2D6", "#FFFF99", "#B15928") #colorblind-friendly palette with 12 colors.
V(net1_multi)$color <- colrs[modules_order$module]
# Plot
png(filename= "figures/multilayer_all.mod.png", res= 300, height= 3000, width= 4900)
par(mfrow=c(1,1),mar=c(1,1,3,17))
plot(net1_multi, 
     vertex.color = V(net1_multi)$color, 
     vertex.frame.color= V(net1_multi)$color, 
     vertex.shape = V(net1_multi)$shape, 
     vertex.size=3,
     vertex.label=V(net1_multi)$name,
     vertex.label.color="white", # white
     vertex.label.cex=.02, # .02
     edge.color = E(net1_multi)$color, 
     edge.curved=0.3, 
     layout=l)
title(main = "Plant-bee resin foraging multilayer network (All MBSDs)", cex.main=2)
legend(x = 1.3,y = 0.9, title="Taxon",
       legend = c("Bees", "Plants"), pch = c(15,16),   
       text.col = "gray20", title.col = "black", 
       box.lwd = 0, cex = 1.5, col=c("grey", "grey"))
legend(x = 1.3,y = 0.1, title="Layers (biogeographic regions)",
       legend = c("Afrotropics", "Indo-Malayan-Australasia", "Neotropics"),
       fill = c("#A151A1", "#FFD26B", "#B3DCF9"), border = "white", 
       text.col = "gray20", title.col = "black", box.lwd = 0, cex = 1.5)
legend(x = 1.3,y = -0.6, title="Modules",
       legend = c("node colors = modules"),
       fill = c("grey"), border = "white", 
       text.col = "gray20", title.col = "black", box.lwd = 0, cex = 1.5)
par(mfrow=c(1,1))
dev.off()

################################
### 3.2 MULTILAYER: CHEMICAL ###

# Create chemical network
net_multi_c = graph_from_data_frame(d=links_c, vertices=nodes_c, directed=F) 

#Check the multilayer network
net_multi_c
attributes(V(net_multi_c))
attributes(E(net_multi_c))
V(net_multi_c)$name # node names
V(net_multi_c)$taxon # plant or bee
V(net_multi_c)$type # 1 = bee, 2 = plant
E(net_multi_c)$layer # Neotropics or Indo-Malayan-Australasia or Afrotropics
E(net_multi_c)$weight # 1 for all

# Import module membership
modules_c=read.table("data/partitions_c.txt", h=T)
head(modules_c)
unique(modules_c$module) # check the number of modules (= 8)

# Draw the network
# Set the same layout for all graphs
set.seed(16381)
l <- layout_nicely(net_multi_c)
# Set NODE shape for bipartite structure (bee = square; plant = circle)
V(net_multi_c)$shape = V(net_multi_c)$taxon
V(net_multi_c)$shape = gsub("bee","square",V(net_multi_c)$shape)
V(net_multi_c)$shape = gsub("plant","circle",V(net_multi_c)$shape)
# Set EDGE colors for the layers (biogeographic regions)
E(net_multi_c)$color = E(net_multi_c)$layer
E(net_multi_c)$color = gsub("Neotropics","#B3DCF9",E(net_multi_c)$color)
E(net_multi_c)$color = gsub("Indo-Malayan-Australasia","#FFD26B",E(net_multi_c)$color)
E(net_multi_c)$color = gsub("Afrotropics","#A151A1",E(net_multi_c)$color)
# Set polarity
E(net_multi_c)$arrow.mode = 0
# Set EDGE width
E(net_multi_c)$width = E(net_multi_c)$layer
E(net_multi_c)$width = gsub("Neotropics", 1.5, E(net_multi_c)$width)
E(net_multi_c)$width = gsub("Indo-Malayan-Australasia", 1.5, E(net_multi_c)$width)
E(net_multi_c)$width = gsub("Afrotropics", 1.5, E(net_multi_c)$width)
# Set NODE colors for modules (= 10 modules)
vertex_name_order = data.frame(nodes=V(net_multi_c)$name)
modules_order = merge(vertex_name_order, modules_c, by="nodes", sort = F)
# colrs = rainbow(8, alpha=1) #rainbow palette with 9 colors
colrs = c("#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#FF7F00", "#CAB2D6", "#FFFF99", "#B15928") #colorblind-friendly palette with 12 colors.
V(net_multi_c)$color <- colrs[modules_order$module]
# Plot
png(filename= "figures/multilayer_c.mod.png", res= 300, height= 3000, width= 4900)
par(mfrow=c(1,1),mar=c(1,1,3,17))
plot(net_multi_c, 
     vertex.color = V(net_multi_c)$color, 
     vertex.frame.color= V(net_multi_c)$color, 
     vertex.shape = V(net_multi_c)$shape, 
     vertex.size=3,
     vertex.label=V(net_multi_c)$name,
     vertex.label.color="white", # white or black
     vertex.label.cex=.02, # .02 or .5
     edge.color = E(net_multi_c)$color, 
     edge.curved=0.3, 
     layout=l)
title(main = "Plant-bee resin foraging multilayer network (MBSD: Chemical)", cex.main=2)
legend(x = 1.3,y = 0.9, title="Taxon",
       legend = c("Bees", "Plants"), pch = c(15,16),   
       text.col = "gray20", title.col = "black", 
       box.lwd = 0, cex = 1.5, col=c("grey", "grey"))
legend(x = 1.3,y = 0.1, title="Layers (biogeographic regions)",
       legend = c("Afrotropics", "Indo-Malayan-Australasia", "Neotropics"),
       fill = c("#A151A1", "#FFD26B", "#B3DCF9"), border = "white", 
       text.col = "gray20", title.col = "black", box.lwd = 0, cex = 1.5)
legend(x = 1.3,y = -0.6, title="Modules",
       legend = c("node colors = modules"),
       fill = c("grey"), border = "white", 
       text.col = "gray20", title.col = "black", box.lwd = 0, cex = 1.5)
par(mfrow=c(1,1))
dev.off()

#################################
### 3.3 MULTILAYER: FIELDWORK ###

# Create chemical network
net_multi_f = graph_from_data_frame(d=links_f, vertices=nodes_f, directed=F) 

#Check the multilayer network
net_multi_f
attributes(V(net_multi_f))
attributes(E(net_multi_f))
V(net_multi_f)$name # node names
V(net_multi_f)$taxon # plant or bee
V(net_multi_f)$type # 1 = bee, 2 = plant
E(net_multi_f)$layer # Neotropics or Indo-Malayan-Australasia or Afrotropics
E(net_multi_f)$weight # 1 for all

# Import module membership
modules_f=read.table("data/partitions_f.txt", h=T)
head(modules_f)
unique(modules_f$module) # check the number of modules (= 6)

# Draw the network
# Set the same layout for all graphs
set.seed(16381)
l <- layout_nicely(net_multi_f)
# Set NODE shape for bipartite structure (bee = square; plant = circle)
V(net_multi_f)$shape = V(net_multi_f)$taxon
V(net_multi_f)$shape = gsub("bee","square",V(net_multi_f)$shape)
V(net_multi_f)$shape = gsub("plant","circle",V(net_multi_f)$shape)
# Set EDGE colors for the layers (biogeographic regions)
E(net_multi_f)$color = E(net_multi_f)$layer
E(net_multi_f)$color = gsub("Neotropics","#B3DCF9",E(net_multi_f)$color)
E(net_multi_f)$color = gsub("Indo-Malayan-Australasia","#FFD26B",E(net_multi_f)$color)
E(net_multi_f)$color = gsub("Afrotropics","#A151A1",E(net_multi_f)$color)
# Set polarity
E(net_multi_f)$arrow.mode = 0
# Set EDGE width
E(net_multi_f)$width = E(net_multi_f)$layer
E(net_multi_f)$width = gsub("Neotropics", 1.5, E(net_multi_f)$width)
E(net_multi_f)$width = gsub("Indo-Malayan-Australasia", 1.5, E(net_multi_f)$width)
E(net_multi_f)$width = gsub("Afrotropics", 1.5, E(net_multi_f)$width)
# Set NODE colors for modules (= 10 modules)
vertex_name_order = data.frame(nodes=V(net_multi_f)$name)
modules_order = merge(vertex_name_order, modules_f, by="nodes", sort = F)
# colrs = rainbow(6, alpha=1) #rainbow palette with 6 colors
colrs = c("#1F78B4", "#B2DF8A", "#FB9A99", "#CAB2D6", "#FFFF99", "#B15928") #colorblind-friendly palette with 12 colors.
V(net_multi_f)$color <- colrs[modules_order$module]
# Plot
png(filename= "figures/multilayer_f.mod.png", res= 300, height= 3000, width= 4900)
par(mfrow=c(1,1),mar=c(1,1,3,17))
plot(net_multi_f, 
     vertex.color = V(net_multi_f)$color, 
     vertex.frame.color= V(net_multi_f)$color, 
     vertex.shape = V(net_multi_f)$shape, 
     vertex.size=3,
     vertex.label=V(net_multi_f)$name,
     vertex.label.color="white", # white or black
     vertex.label.cex=.02, # .02 or .5
     edge.color = E(net_multi_f)$color, 
     edge.curved=0.3, 
     layout=l)
title(main = "Plant-bee resin foraging multilayer network (MBSD: Fieldwork)", cex.main=2)
legend(x = 1.3,y = 0.9, title="Taxon",
       legend = c("Bees", "Plants"), pch = c(15,16),   
       text.col = "gray20", title.col = "black", 
       box.lwd = 0, cex = 1.5, col=c("grey", "grey"))
legend(x = 1.3,y = 0.1, title="Layers (biogeographic regions)",
       legend = c("Indo-Malayan-Australasia", "Neotropics"),
       fill = c("#FFD26B", "#B3DCF9"), border = "white", 
       text.col = "gray20", title.col = "black", box.lwd = 0, cex = 1.5)
legend(x = 1.3,y = -0.6, title="Modules",
       legend = c("node colors = modules"),
       fill = c("grey"), border = "white", 
       text.col = "gray20", title.col = "black", box.lwd = 0, cex = 1.5)
par(mfrow=c(1,1))
dev.off()

#####################################
### 3.4 MULTILAYER: PALYNOLOGICAL ###

# Create chemical network
net_multi_p = graph_from_data_frame(d=links_p, vertices=nodes_p, directed=F) 

#Check the multilayer network
net_multi_p
attributes(V(net_multi_p))
attributes(E(net_multi_p))
V(net_multi_p)$name # node names
V(net_multi_p)$taxon # plant or bee
V(net_multi_p)$type # 1 = bee, 2 = plant
E(net_multi_p)$layer # Neotropics or Indo-Malayan-Australasia or Afrotropics
E(net_multi_p)$weight # 1 for all

# Import module membership
modules_p=read.table("data/partitions_p.txt", h=T)
head(modules_p)
unique(modules_p$module) # check the number of modules (= 5)

# Draw the network
# Set the same layout for all graphs
set.seed(16381)
l <- layout_nicely(net_multi_p)
# Set NODE shape for bipartite structure (bee = square; plant = circle)
V(net_multi_p)$shape = V(net_multi_p)$taxon
V(net_multi_p)$shape = gsub("bee","square",V(net_multi_p)$shape)
V(net_multi_p)$shape = gsub("plant","circle",V(net_multi_p)$shape)
# Set EDGE colors for the layers (biogeographic regions)
E(net_multi_p)$color = E(net_multi_p)$layer
E(net_multi_p)$color = gsub("Neotropics","#B3DCF9",E(net_multi_p)$color)
E(net_multi_p)$color = gsub("Indo-Malayan-Australasia","#FFD26B",E(net_multi_p)$color)
E(net_multi_p)$color = gsub("Afrotropics","#A151A1",E(net_multi_p)$color)
# Set polarity
E(net_multi_p)$arrow.mode = 0
# Set EDGE width
E(net_multi_p)$width = E(net_multi_p)$layer
E(net_multi_p)$width = gsub("Neotropics", 1.5, E(net_multi_p)$width)
E(net_multi_p)$width = gsub("Indo-Malayan-Australasia", 1.5, E(net_multi_p)$width)
E(net_multi_p)$width = gsub("Afrotropics", 1.5, E(net_multi_p)$width)
# Set NODE colors for modules (= 5 modules)
vertex_name_order = data.frame(nodes=V(net_multi_p)$name)
modules_order = merge(vertex_name_order, modules_p, by="nodes", sort = F)
# colrs = rainbow(6, alpha=1) #rainbow palette with 6 colors
colrs = c("#1F78B4", "#B2DF8A", "#FB9A99", "#CAB2D6", "#FFFF99") #colorblind-friendly palette with 12 colors.
V(net_multi_p)$color <- colrs[modules_order$module]
# Plot
png(filename= "figures/multilayer_p_.mod.png", res= 300, height= 3000, width= 4900)
par(mfrow=c(1,1),mar=c(1,1,3,17))
plot(net_multi_p, 
     vertex.color = V(net_multi_p)$color, 
     vertex.frame.color= V(net_multi_p)$color, 
     vertex.shape = V(net_multi_p)$shape, 
     vertex.size=3,
     vertex.label=V(net_multi_p)$name,
     vertex.label.color="black", # white or black
     vertex.label.cex=.5, # .02 or .5
     edge.color = E(net_multi_p)$color, 
     edge.curved=0.3, 
     layout=l)
title(main = "Plant-bee resin foraging multilayer network (MBSD: Palynological)", cex.main=2)
legend(x = 1.3,y = 0.9, title="Taxon",
       legend = c("Bees", "Plants"), pch = c(15,16),   
       text.col = "gray20", title.col = "black", 
       box.lwd = 0, cex = 1.5, col=c("grey", "grey"))
legend(x = 1.3,y = 0.1, title="Layers (biogeographic regions)",
       legend = c("Indo-Malayan-Australasia", "Neotropics"),
       fill = c("#FFD26B", "#B3DCF9"), border = "white", 
       text.col = "gray20", title.col = "black", box.lwd = 0, cex = 1.5)
legend(x = 1.3,y = -0.6, title="Modules",
       legend = c("node colors = modules"),
       fill = c("grey"), border = "white", 
       text.col = "gray20", title.col = "black", box.lwd = 0, cex = 1.5)
par(mfrow=c(1,1))
dev.off()


#################################
### 3.5 MULTILAYER: CHEMICAL + FIELDWORK ###

# Create network
net1_multi = graph_from_data_frame(d=links_cf, vertices=nodes_cf, directed=F) 

#Check the multilayer network
net1_multi
attributes(V(net1_multi))
attributes(E(net1_multi))
V(net1_multi)$name # node names
V(net1_multi)$taxon # plant or bee
V(net1_multi)$type # 1 = bee, 2 = plant
E(net1_multi)$layer # Neotropics or Indo-Malayan-Australasia or Afrotropics
E(net1_multi)$weight # 1 for all

# Import module membership
modules=read.table("data/partitions_cf.txt", h=T)
head(modules)
unique(modules$module) # check the number of modules (= 9)

# Draw the network
# Set the same layout for all graphs
set.seed(16381)
l <- layout_nicely(net1_multi)
# Set NODE shape for bipartite structure (bee = square; plant = circle)
V(net1_multi)$shape = V(net1_multi)$taxon
V(net1_multi)$shape = gsub("bee","square",V(net1_multi)$shape)
V(net1_multi)$shape = gsub("plant","circle",V(net1_multi)$shape)
# Set NODE colors for the bipartite structure (bee = brown; plant = green) 
#V(net1_multi)$color = V(net1_multi)$taxon
#V(net1_multi)$color = gsub("bee","#D3802B",V(net1_multi)$color)
#V(net1_multi)$color = gsub("plant","#2C8437",V(net1_multi)$color)
# Set EDGE colors for the layers (biogeographic regions)
E(net1_multi)$color = E(net1_multi)$layer
E(net1_multi)$color = gsub("Neotropics","#B3DCF9",E(net1_multi)$color)
E(net1_multi)$color = gsub("Indo-Malayan-Australasia","#FFD26B",E(net1_multi)$color)
E(net1_multi)$color = gsub("Afrotropics","#A151A1",E(net1_multi)$color)
# Set polarity
E(net1_multi)$arrow.mode = 0
# Set EDGE width
E(net1_multi)$width = E(net1_multi)$layer
E(net1_multi)$width = gsub("Neotropics", 1, E(net1_multi)$width)
E(net1_multi)$width = gsub("Indo-Malayan-Australasia", 1, E(net1_multi)$width)
E(net1_multi)$width = gsub("Afrotropics", 1, E(net1_multi)$width)
# Set NODE colors for modules (= 10 modules)
vertex_name_order = data.frame(nodes=V(net1_multi)$name)
modules_order = merge(vertex_name_order, modules, by="nodes", sort = F)
#colrs = rainbow(12, alpha=1) #rainbow palette with 9 colors
colrs = c("#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#B41572", "#FF7F00", "#CAB2D6", "#FFFF99", "#B15928", "#AB2737", "#CCA737", "#1A1CC7") #colorblind-friendly palette with 12 colors.
V(net1_multi)$color <- colrs[modules_order$module]
# Plot
png(filename= "figures/multilayer_cf_.mod.png", res= 300, height= 3000, width= 4900)
par(mfrow=c(1,1),mar=c(1,1,3,17))
plot(net1_multi, 
     vertex.color = V(net1_multi)$color, 
     vertex.frame.color= V(net1_multi)$color, 
     vertex.shape = V(net1_multi)$shape, 
     vertex.size=3,
     vertex.label=V(net1_multi)$name,
     vertex.label.color="black", # white
     vertex.label.cex=.7, # .02
     edge.color = E(net1_multi)$color, 
     edge.curved=0.3, 
     layout=l)
title(main = "Plant-bee resin foraging multilayer network (Chemical + Fieldwork)", cex.main=2)
legend(x = 1.3,y = 0.9, title="Taxon",
       legend = c("Bees", "Plants"), pch = c(15,16),   
       text.col = "gray20", title.col = "black", 
       box.lwd = 0, cex = 1.5, col=c("grey", "grey"))
legend(x = 1.3,y = 0.1, title="Layers (biogeographic regions)",
       legend = c("Afrotropics", "Indo-Malayan-Australasia", "Neotropics"),
       fill = c("#A151A1", "#FFD26B", "#B3DCF9"), border = "white", 
       text.col = "gray20", title.col = "black", box.lwd = 0, cex = 1.5)
legend(x = 1.3,y = -0.6, title="Modules",
       legend = c("node colors = modules"),
       fill = c("grey"), border = "white", 
       text.col = "gray20", title.col = "black", box.lwd = 0, cex = 1.5)
par(mfrow=c(1,1))
dev.off()
