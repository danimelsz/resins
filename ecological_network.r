# Load Packages
library(ape)
library(bipartite)
library(brainGraph)
library(corrgram)
library(fmsb)
library(geiger)
library(ggplot2)
library(igraph)
library(phytools)
library(picante)
library(reshape2)
library(RInSp)
library(tidyverse)

################################################
################ 1. Data input #################
################################################

# 1.1 Input for bipartite analysis

# Load resin matrix with plant families
resin_family.mat = read.delim("matrix_family.txt", row.names=1)%>% as.matrix()

# Load resin matrix with plant species
resin_species.mat <- read.delim("matrix_sp.txt", row.names=1)%>% as.matrix()

# Remove nodes with no interactions
resin_family2.mat <- empty(resin_family.mat, count = F)
resin_species2.mat <- empty(resin_species.mat, count = F)

# Check the object class (is it "matrix array"?)
class(resin_family2.mat)
class(resin_species2.mat)

# 1.2 Input for igraph analyses

#Convert resin_family matrix from bipartite to igraph
resin_family.igr <- graph_from_incidence_matrix(resin_family.mat, 
                                              directed = F, 
                                              weighted = NULL)
attributes(V(resin_family.igr)) #Check nodes
attributes(E(resin_family.igr)) #check edges

#Convert resin_species matrix from bipartite to igraph
resin_species.igr <- graph_from_incidence_matrix(resin_species.mat, 
                                              directed = F, 
                                              weighted = NULL)
attributes(V(resin_species.igr)) #Check nodes
attributes(E(resin_species.igr)) #check edges

# Remove nodes with no interactions
resin_family2.igr <- delete.vertices(resin_family.igr, degree(resin_family.igr)==0)
resin_species2.igr <- delete.vertices(resin_species.igr, degree(resin_species.igr)==0)

# 1.3 Phylogeny

# Load tree
tree <- read.nexus(file = "tree_edited.nex")
plotTree(tree,node.numbers=T,fsize=0.35)

########################################################
################# 2. Network structure #################
########################################################

# 2.1 Exploratory Visualization

# Interaction matrix
visweb(resin_family2.mat,type="nested") #family-level
visweb(resin_species2.mat, type="nested") #species-level

# Bipartite plotting
plotweb(resin_family2.mat,
        method = "cca", 
        text.rot = 90, 
        bor.col.interaction = F, 
        col.interaction = "grey90", 
        col.high = "#D4B36A",      
        col.low = "#87468B")
plotweb(resin_species2.mat,
        method = "cca", 
        text.rot = 90, 
        bor.col.interaction = F, 
        col.interaction = "grey90", 
        col.high = "#D4B36A",      
        col.low = "#87468B",
        labsize=0.3)  

#Setting colors for plants and bees
V(resin_family2.igr)$color <- V(resin_family2.igr)$type
V(resin_family2.igr)$color <- gsub(FALSE,"purple",V(resin_family2.igr)$color)
V(resin_family2.igr)$color <- gsub(TRUE,"gold2",V(resin_family2.igr)$color)

V(resin_species2.igr)$color <- V(resin_species2.igr)$type
V(resin_species2.igr)$color <- gsub(FALSE,"purple",V(resin_species2.igr)$color)
V(resin_species2.igr)$color <- gsub(TRUE,"gold2",V(resin_species2.igr)$color)

# Fruchterman Reingold
algoritmo_family = layout.fruchterman.reingold(resin_family2.igr)
plot(resin_family2.igr, # note que há apenas 1 componente
     layout = algoritmo_family,
     vertex.color = V(resin_species2.igr)$color, 
     vertex.size = 6, 
     vertex.label.cex = 0.3,
     vertex.label.color = "white",
     vertex.frame.color = F,
     edge.color = "grey90",
     edge.width = 1)

algoritmo_species = layout.fruchterman.reingold(resin_species2.igr)
plot(resin_species2.igr, # note que há 7 componentes
     layout = algoritmo_species,
     vertex.color = V(resin_species2.igr)$color, 
     vertex.size = 4, 
     vertex.label.cex = 0.3,
     vertex.label.color = "white",
     vertex.frame.color = F,
     edge.color = "grey90",
     edge.width = 3)
     
# 2.2 Network-level metrics

# C-score
grouplevel(resin_family2.mat, index = "C score")
grouplevel(resin_species2.mat, index = "C score")

# NODF
networklevel(resin_family2.mat, index = "NODF")
networklevel(resin_species2.mat, index = "NODF")

# Modularity for family-level botanical identification
modulos_family = computeModules(resin_family2.mat, method = "Beckett") # calculate with bipartite
modulos_family@likelihood 
plotModuleWeb(modulos_family, labsize=0.2) # plot module composition
modulos_family.lou = cluster_louvain(resin_family2.igr) # calculate with igraph
modulos_family_verificar.lis = listModuleInformation(modulos_family) # list bee and plant species in each module
modulos_family.vet = module2constraints(modulos_family) # save a vector with the module within which each node is part of
modulos_family.df = data.frame(c(rownames(resin_family2.mat),  
                                 colnames(resin_family2.mat)), 
                               modulos_family.vet)             # create a dataframe with nodes and modules
colnames(modulos_family.df) = c("vertices", "modulos") # rename columns from the previous dataframe
modulos_family.lis = split(modulos_family.df$vertices, 
                           modulos_family.df$modulos)  # convert dataframe into a list
cores_family = rainbow(length(modulos_family.lis), 
                       alpha = 0.7,                
                       s = 1, v = 0.8)                 # create a vector with colors for nodes
V(resin_family2.igr)$color = cores_family[modulos_family.df$modulos]
nuvens_family = cores_family
aresta_family.inicio = ends(resin_family2.igr,       
                            es=E(resin_family2.igr), 
                            names=F)[,1]             
aresta_family.cor = V(resin_family2.igr)$color[aresta_family.inicio]
V(resin_family2.igr)$shape <- V(resin_family2.igr)$type
V(resin_family2.igr)$shape <- gsub(FALSE,"square",V(resin_family2.igr)$shape) #Plantas
V(resin_family2.igr)$shape <- gsub(TRUE,"circle",V(resin_family2.igr)$shape)  #Abelhas
plot (resin_family2.igr,
     layout = algoritmo_family, 
     vertex.color = V(resin_family2.igr)$color,         #nós com as cores dos módulos
     vertex.size = 4, 
     vertex.label.cex = 0.2,
     vertex.label.color = "white",
     vertex.frame.color = V(resin_family2.igr)$color,   #Borda com a mesma cor do nó
     vertex.shape = V(resin_family2.igr)$shape,          #Forma dos 2 tipos de nós
     edge.color = adjustcolor(aresta_family.cor, alpha = 0.3), #Arestas c/cores de módulos
     mark.groups = modulos_family.lis,                  #Nuvens baseadas nos módulos
     mark.border = "white",                             #Nuvens com bordas brancas
     mark.col = adjustcolor(nuvens_family, alpha = 0.2) #Nuvens com as cores dos módulos
)

# 2.3 p-values for network-level metrics

# Create null models
set.seed(14) # create a seed
permutacoes = 100 # change this for 1000 to improve the null models
modelo = "vaznull" # null model method
aleatorizadas_family = nullmodel(resin_family2.mat,
                                      N = permutacoes,
                                      method = modelo)
aleatorizadas_species = nullmodel(resin_species2.mat,
                                      N = permutacoes,
                                      method = modelo)
                     
# p-value for NODF                    
metrica_family_NODF = c("NODF")
observado_family_NODF = networklevel(resin_family2.mat, index = metrica_family_NODF)
nulos_family_NODF = unlist(sapply(aleatorizadas_family,
                                  networklevel, 
                                  index=metrica_family_NODF))
par(mar = c(4,4,5,4))
plot(density(nulos_family_NODF),
     main = "Observed (red) vs random values (black)",
     xlim = c(min((observado_family_NODF), min(nulos_family_NODF)),
              max((observado_family_NODF), max(nulos_family_NODF))))
abline(v=observado_family_NODF, col="red", lwd=2,xlab="")
sum(nulos_family_NODF>=(observado_family_NODF)) / length(nulos_family_NODF)
sum(nulos_family_NODF<=(observado_family_NODF)) / length(nulos_family_NODF)

# p-value for modularity                  
algoritmo_family_modulo = c("Beckett")
mod_family = computeModules(resin_family2.mat, method = "Beckett")
mod_family@likelihood
part_family_modulo = bipartite::module2constraints(mod_family)
row_family.part = part_family_modulo[1:nrow(resin_family2.mat)]
col_family.part = part_family_modulo[(nrow(resin_family2.mat)+1):(nrow(resin_family2.mat)+ncol(resin_family2.mat))]
length(unique(part_family_modulo))
nullmod_family = sapply(aleatorizadas_family, computeModules, method = algoritmo_family_modulo)
modnull_family = sapply (nullmod_family, function(x) x@likelihood)
par(mar = c(4,4,5,4))
plot(density(modnull_family),
     main = "Observed vs randomized",
     xlim = c(min((mod_family@likelihood), min(modnull_family)),
              max((mod_family@likelihood), max(modnull_family))))
abline(v=mod_family@likelihood, col="red", lwd=2, xlab="")
mod_family@likelihood 
mean(modnull_family)  
sd(modnull_family)    
(mod_family@likelihood - mean(modnull_family)) / sd(modnull_family)   #valor Z
sum(modnull_family>=(mod_family@likelihood)) / length(modnull_family)
sum(modnull_family<=(mod_family@likelihood)) / length(modnull_family) 

#########################################################
################# 3. Node-level metrics #################
#########################################################

# Specialization d'
specieslevel(resin_family2.mat, index="d")
specieslevel(resin_species2.mat, index="d")

# Closeness centrality
proximidade.bip.family = specieslevel(resin_family2.mat, index = "closeness") # bipartite
igraph::closeness(resin_family2.igr, weights = NULL) # igraph

# Betweeness centrality
intermedio.bip.family = specieslevel(resin_family2.mat, index = "betweenness")
