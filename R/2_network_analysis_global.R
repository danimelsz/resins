###########
# 1 SETUP #
###########

# 1.1 Setting directories

# Set directory where packages should be installed
.libPaths('/home/danimelsz/Downloads/R/x86_64-pc-linux-gnu-library/4.1')

# Set directory
setwd("~/Desktop/Stingless_Bees/pub_abelha/dados")

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
library(MCMCglmm)
library(phytools)
library(picante)
library(reshape2)
library(RInSp)
library(tidyverse)

################
# 2 DATA INPUT #
################

# 2.1 Global

# 2.1.1 Bipartite
resin_family.mat = read.delim("data/matrix_family_full.txt", row.names=1)%>% as.matrix()
# Remove nodes with no interactions
resin_family2.mat <- empty(resin_family.mat, count = F)
# Remove nodes from disconnected, isolated subnetworks
resin_family2.mat = resin_family2.mat[, !colnames(resin_family2.mat) %in% c("Friseomelita_silvestrii", 'Lepidotrigona_terminata', "Melipona_orbignyi")]
problematic.plants.f = c("Podocarpaceae","Balanophoraceae","Orchidaceae")
resin_family2.mat = resin_family2.mat[!rownames(resin_family2.mat) %in% problematic.plants.f, ]

# 2.1.2 igraph
resin_family.igr <- graph_from_incidence_matrix(resin_family2.mat, directed = F, weighted = NULL)
attributes(V(resin_family.igr)) #Check nodes
attributes(E(resin_family.igr)) #check edges
# Remove nodes with no interactions
resin_family2.igr <- delete.vertices(resin_family.igr, degree(resin_family.igr)==0)

# 2.2 Chemical
global_family_chem.mat = read.delim("data/matrix_family_chem.txt", row.names=1)%>% as.matrix()
# Remove nodes with no interactions
global_family_chem2.mat <- empty(global_family_chem.mat, count = F)
# Remove problematic bee species from disconnected subnetworks
global_family_chem2.mat = global_family_chem2.mat[, !colnames(global_family_chem2.mat) %in% c("Friseomelita_silvestrii", "Melipona_orbignyi")]
# Remove problematic plants from disconnected subnetworks
problematic.plants.f = c("Podocarpaceae","Balanophoraceae","Orchidaceae")
global_family_chem2.mat = global_family_chem2.mat[!rownames(global_family_chem2.mat) %in% problematic.plants.f, ]
# Remove species that are not connected with the main network
global_family_chem2.mat = global_family_chem2.mat[, !colnames(global_family_chem2.mat) %in% c("Melipona_beecheii")]

# 2.3 Fieldwork
global_family_field.mat = read.delim("data/matrix_family_fieldwork.txt", row.names=1) %>% as.matrix()
# Remove nodes with no interactions
global_family_field2.mat = empty(global_family_field.mat, count = F)
# Remove nodes with no interactions
global_family_field2.igr <- delete.vertices(global_family_field.igr, degree(global_family_field.igr)==0)

# 2.4 Palynological
global_family_pal.mat = read.delim("data/matrix_family_pal.txt", row.names=1)%>% as.matrix()
# Remove nodes with no interactions
global_family_pal2.mat = empty(global_family_pal.mat, count = F)
# Remove problematic bee species with a single interaction reported
global_family_pal2.mat = global_family_pal2.mat[, !colnames(global_family_pal2.mat) %in% c("Friseomelita_silvestrii", "Melipona_orbignyi")]
# Remove problematic plants with a single interaction reported
problematic.plants.f = c("Podocarpaceae","Balanophoraceae","Orchidaceae")
global_family_pal2.mat = global_family_pal2.mat[!rownames(global_family_pal2.mat) %in% problematic.plants.f, ]

################
# 3 NESTEDNESS #
################

# 3.1 Global
# Empirical NODF
observado_family_NODF_g = networklevel(resin_family2.mat, index = "NODF")
# NODF from null models
set.seed(14)
permutacoes = 10 # change to 1000
modelo = "vaznull"
aleatorizadas_family_g = nullmodel(resin_family2.mat, N = permutacoes, method = modelo)
nulos_family_NODF_g = unlist(sapply(aleatorizadas_family_g,
                                  networklevel,
                                  index=metrica_family_NODF))
mean(nulos_family_NODF_g)
# Z value
(observado_family_NODF_g-mean(nulos_family_NODF_g)) / sd(nulos_family_NODF_g)
# p-value left
sum(nulos_family_NODF_g>=(observado_family_NODF_g)) / length(nulos_family_NODF_g)
# p-value right
sum(nulos_family_NODF_g<=(observado_family_NODF_g)) / length(nulos_family_NODF_g)
# Plot NODF values from null models vs empirical NODF
pdf(file = "figures/NODF_nullmodels_global.pdf",   # The directory you want to save the file in
    width = 20, # The width of the plot in inches
    height = 10) # The height of the plot in inches
par(mar = c(4,4,5,4))
plot(density(nulos_family_NODF_g),
     main = "Observed (red) vs random values (black)",
     xlim = c(min((observado_family_NODF_g), min(nulos_family_NODF_g)),
              max((observado_family_NODF_g), max(nulos_family_NODF_g))))
abline(v=observado_family_NODF_g, col="red", lwd=2,xlab="")
dev.off()

# 3.2 Chemical
# Empirical NODF
observado_family_NODF_c = networklevel(global_family_chem2.mat, index = "NODF")
# NODF from null models
set.seed(14)
permutacoes = 10 # change to 1000
modelo = "vaznull"
aleatorizadas_family_c = nullmodel(global_family_chem2.mat, N= permutacoes,method = modelo)
nulos_family_NODF_c = unlist(sapply(aleatorizadas_family_c,
                                  networklevel, 
                                  index=metrica_family_NODF))
mean(nulos_family_NODF_c)
# Z value
(observado_family_NODF_c-mean(nulos_family_NODF_c)) / sd(nulos_family_NODF_c)
# p-value left
sum(nulos_family_NODF_c>=(observado_family_NODF_c)) / length(nulos_family_NODF_c)
# p-value right
sum(nulos_family_NODF_c<=(observado_family_NODF_c)) / length(nulos_family_NODF_c)
# Plot NODF values from null models vs empirical NODF
pdf(file = "figures/NODF_nullmodels_chem.pdf",   # The directory you want to save the file in
    width = 20, # The width of the plot in inches
    height = 10) # The height of the plot in inches
par(mar = c(4,4,5,4))
plot(density(nulos_family_NODF_c),
     main = "Observed (red) vs random values (black)",
     xlim = c(min((observado_family_NODF_c), min(nulos_family_NODF_c)),
              max((observado_family_NODF_c), max(nulos_family_NODF_c))))
abline(v=observado_family_NODF_c, col="red", lwd=2,xlab="")
dev.off()

# 3.3 Fieldwork
# Empirical NODF
observado_family_NODF_f = networklevel(global_family_field2.mat, index = "NODF")
# NODF from null models
set.seed(14)
permutacoes = 10 # change to 1000
modelo = "vaznull"
aleatorizadas_family_f = nullmodel(global_family_field2.mat, N= permutacoes,method = modelo)
nulos_family_NODF_f = unlist(sapply(aleatorizadas_family_f,
                                  networklevel, 
                                  index=metrica_family_NODF))
mean(nulos_family_NODF_f)
# Z value
(observado_family_NODF_f-mean(nulos_family_NODF_f)) / sd(nulos_family_NODF_f)
# p-value left
sum(nulos_family_NODF_f>=(observado_family_NODF_f)) / length(nulos_family_NODF_f)
# p-value right
sum(nulos_family_NODF_f<=(observado_family_NODF_f)) / length(nulos_family_NODF_f)
# Plot NODF values from null models vs empirical NODF
pdf(file = "figures/NODF_nullmodels_field.pdf",   # The directory you want to save the file in
    width = 20, # The width of the plot in inches
    height = 10) # The height of the plot in inches
par(mar = c(4,4,5,4))
plot(density(nulos_family_NODF_f),
     main = "Observed (red) vs random values (black)",
     xlim = c(min((observado_family_NODF_f), min(nulos_family_NODF_f)),
              max((observado_family_NODF_f), max(nulos_family_NODF_f))))
abline(v=observado_family_NODF_f, col="red", lwd=2,xlab="")
dev.off()

# 3.4 Palynological
# Empirical NODF
observado_family_NODF_p = networklevel(global_family_pal2.mat, index = "NODF")
# NODF from null models
set.seed(14)
permutacoes = 10 # change to 1000
modelo = "vaznull"
aleatorizadas_family_p = nullmodel(global_family_pal2.mat, N= permutacoes,method = modelo)
nulos_family_NODF_p = unlist(sapply(aleatorizadas_family_p,
                                  networklevel, 
                                  index=metrica_family_NODF))
mean(nulos_family_NODF_p)
# Z value
(observado_family_NODF_p-mean(nulos_family_NODF_p)) / sd(nulos_family_NODF_p)
# p-value left
sum(nulos_family_NODF_p>=(observado_family_NODF_p)) / length(nulos_family_NODF_p)
# p-value right
sum(nulos_family_NODF_p<=(observado_family_NODF_p)) / length(nulos_family_NODF_p)
# Plot NODF values from null models vs empirical NODF
pdf(file = "figures/NODF_nullmodels_field.pdf",   # The directory you want to save the file in
    width = 20, # The width of the plot in inches
    height = 10) # The height of the plot in inches
par(mar = c(4,4,5,4))
plot(density(nulos_family_NODF_p),
     main = "Observed (red) vs random values (black)",
     xlim = c(min((observado_family_NODF_p), min(nulos_family_NODF_p)),
              max((observado_family_NODF_p), max(nulos_family_NODF_p))))
abline(v=observado_family_NODF_p, col="red", lwd=2,xlab="")
dev.off()

################
# 4 MODULARITY #
################

# 4.1 Global
# Empirical M
mod_family_g = computeModules(resin_family2.mat, method = "Beckett")
mod_family_g@likelihood
# Extract modules
part_family_modulo_g = bipartite::module2constraints(mod_family_g)
row_family.part_g = part_family_modulo_g[1:nrow(resin_family2.mat)]
col_family.part_g = part_family_modulo_g[(nrow(resin_family2.mat)+1):(nrow(resin_family2.mat>
length(unique(part_family_modulo_g)) # mostra o número de módulos
# M from null models
nullmod_family_g = sapply(aleatorizadas_family_g, computeModules, method = "Beckett")
modnull_family_g = sapply (nullmod_family_g, function(x) x@likelihood)
# Z value
(mod_family_g@likelihood - mean(modnull_family_g)) / sd(modnull_family_g)
# p value left
sum(modnull_family_g>=(mod_family_g@likelihood)) / length(modnull_family_g)
# p value right
sum(modnull_family_g<=(mod_family_g@likelihood)) / length(modnull_family_g)
# Plot M values from null models vs empirical M
pdf(file = "figures/modularity_nullmodels_g.pdf",
    width = 20, # The width of the plot in inches
    height = 10) # The height of the plot in inches
par(mar = c(4,4,5,4))
plot(density(modnull_family_g),
     main = "Observed vs randomized",
     xlim = c(min((mod_family_g@likelihood), min(modnull_family_g)),
              max((mod_family_g@likelihood), max(modnull_family_g))))
abline(v=mod_family_g@likelihood, col="red", lwd=2, xlab="")
dev.off()

# 4.2 Chemical
# Empirical M
mod_family_c = computeModules(global_family_chem2.mat, method = "Beckett")
mod_family_c@likelihood
# Extract modules
part_family_modulo_c = bipartite::module2constraints(mod_family_c)
row_family.part_c = part_family_modulo_c[1:nrow(global_family_chem2.mat)]
col_family.part_c = part_family_modulo_c[(nrow(global_family_chem2.mat)+1):(nrow(global_family_chem2.mat)+ncol(global_family_chem2.mat))]
length(unique(part_family_modulo_c)) # mostra o número de módulos
# M from null models
nullmod_family_c = sapply(aleatorizadas_family_c, computeModules, method = "Beckett")
modnull_family_c = sapply (nullmod_family_c, function(x) x@likelihood)
# Z value
(mod_family_c@likelihood - mean(modnull_family_c)) / sd(modnull_family_c)
# p value left
sum(modnull_family_c>=(mod_family_c@likelihood)) / length(modnull_family_c)
# p value right
sum(modnull_family_c<=(mod_family_c@likelihood)) / length(modnull_family_c)
# Plot M values from null models vs empirical M
pdf(file = "figures/modularity_nullmodels_chem.pdf",
    width = 20, # The width of the plot in inches
    height = 10) # The height of the plot in inches
par(mar = c(4,4,5,4))
plot(density(modnull_family_c),
     main = "Observed vs randomized",
     xlim = c(min((mod_family_c@likelihood), min(modnull_family_c)),
              max((mod_family_c@likelihood), max(modnull_family_c))))
abline(v=mod_family_c@likelihood, col="red", lwd=2, xlab="")
dev.off()

# 4.3 Fieldwork
# Empirical M
mod_family_f = computeModules(global_family_field2.mat, method = "Beckett")
mod_family_f@likelihood
# Extract modules
part_family_modulo_f = bipartite::module2constraints(mod_family_f)
row_family.part_f = part_family_modulo_f[1:nrow(global_family_field2.mat)]
col_family.part_f = part_family_modulo_f[(nrow(global_family_field2.mat)+1):(nrow(global_family_field2.mat)+ncol(global_family_field2.mat))]
length(unique(part_family_modulo_f)) # mostra o número de módulos
# M from null models
nullmod_family_f = sapply(aleatorizadas_family_f, computeModules, method = "Beckett")
modnull_family_f = sapply (nullmod_family_f, function(x) x@likelihood)
# Z value
(mod_family_f@likelihood - mean(modnull_family_f)) / sd(modnull_family_f)
# p value left
sum(modnull_family_f>=(mod_family_f@likelihood)) / length(modnull_family_f)
# p value right
sum(modnull_family_f<=(mod_family_f@likelihood)) / length(modnull_family_f)
# Plot M values from null models vs empirical M
pdf(file = "figures/modularity_nullmodels_field.pdf",
    width = 20, # The width of the plot in inches
    height = 10) # The height of the plot in inches
par(mar = c(4,4,5,4))
plot(density(modnull_family_f),
     main = "Observed vs randomized",
     xlim = c(min((mod_family_f@likelihood), min(modnull_family_f)),
              max((mod_family_f@likelihood), max(modnull_family_f))))
abline(v=mod_family_f@likelihood, col="red", lwd=2, xlab="")
dev.off()

# 4.4 Palynological
# Empirical M
mod_family_p = computeModules(global_family_pal2.mat, method = "Beckett")
mod_family_p@likelihood
# Extract modules
part_family_modulo_p = bipartite::module2constraints(mod_family_p)
row_family.part_p = part_family_modulo_p[1:nrow(global_family_pal2.mat)]
col_family.part_p = part_family_modulo_p[(nrow(global_family_pal2.mat)+1):(nrow(global_family_pal2.mat)+ncol(global_family_pal2.mat))]
length(unique(part_family_modulo_p)) # mostra o número de módulos
# M from null models
nullmod_family_p = sapply(aleatorizadas_family_p, computeModules, method = "Beckett")
modnull_family_p = sapply (nullmod_family_p, function(x) x@likelihood)
# Z value
(mod_family_p@likelihood - mean(modnull_family_p)) / sd(modnull_family_p)
# p value left
sum(modnull_family_p>=(mod_family_p@likelihood)) / length(modnull_family_p)
# p value right
sum(modnull_family_p<=(mod_family_p@likelihood)) / length(modnull_family_p)
# Plot M values from null models vs empirical M
pdf(file = "figures/modularity_nullmodels_pal.pdf",
    width = 20, # The width of the plot in inches
    height = 10) # The height of the plot in inches
par(mar = c(4,4,5,4))
plot(density(modnull_family_p),
     main = "Observed vs randomized",
     xlim = c(min((mod_family_p@likelihood), min(modnull_family_p)),
              max((mod_family_p@likelihood), max(modnull_family_p))))
abline(v=mod_family_p@likelihood, col="red", lwd=2, xlab="")
dev.off()
