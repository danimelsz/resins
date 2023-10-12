# Testing modularity, nestedness and compound topology

############## ANALYSIS 1.0 ALL MBSDs: AGGREGATED ##############
############## 1. PREPARING THE DATA ##############

#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(bipartite)
library(igraph)
source("R/RestNullModel.R")
source("R/PosteriorProb.R")

# Load edge and node lists
nodes = read.delim("data/net1nodes_.txt", header = T) # node traits
links = read.delim("data/list.3MBSDs_v3_binary", header = T) # list of interactions
# Load network
all_agreggated <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
all_agreggated = get.adjacency(all_agreggated,sparse=FALSE)
all_agreggated = all_agreggated[c(0:68), ] # keep only rows of bees
all_agreggated = all_agreggated[, 69:ncol(all_agreggated)] # keep only columns of plants
data <- empty(all_agreggated, count = F) # remove nodes with no interactions

#Visualize the raw matrix
visweb(data)

############### 2. MODULARITY ANALYSIS  ############### 

#Compute modularity
Mod <- bipartite::computeModules(data)
Mod@likelihood

#Recover the partitions
Part <- bipartite::module2constraints(Mod)
row.Part <- Part[1:nrow(data)]
col.Part <- Part[(nrow(data)+1):(nrow(data)+ncol(data))]

#Test for the significance of modularity with a Monte Carlo procedure

#Generate randomized matrices
set.seed(1)
nulls <- nullmodel(data, N=10, method="r2d")

#Calculate the modularity of the randomized matrices
mod.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(mod.nulls, function(x) x@likelihood)

#Calculate the z-score of the randomized distribution
(z <- (Mod@likelihood - mean(like.nulls))/sd(like.nulls))

#Plot the observed modularity value against the distribution of randomized values
png(filename= "figures/figS5_modularity_AllMBSDs_aggregated.png", res= 300, height= 3000, width= 3000)
plot(density(like.nulls), xlim=c(min((Mod@likelihood), min(like.nulls)), max((Mod@likelihood), max(like.nulls))), 
     main="Observed vs. randomized")
abline(v=(Mod@likelihood), col="red", lwd=2)    
dev.off()

#Estimate the P-value
mean(like.nulls) # E (expected modularity)
z # z
sd(like.nulls)
Mod@likelihood
praw <- sum(like.nulls>(Mod@likelihood)) / length(like.nulls)
ifelse(praw > 0.5, 1-praw, praw)

############### 3a. NESTEDNESS ANALYSIS (RESTRICTED) ############### 

#Calculate the desired nestedness metric (here NODF) for the original network.
obs <- unlist(bipartite::nest.smdm(x = data, 
                                   constraints = Part, #Input the modular structured recovered from step 2
                                   weighted = F, #By considering the edge weights, you are choosing WNODA
                                   decreasing = "fill")) # needs to be fill for binary matrices
#Check the scores
obs

#Calculate constrained interaction probabilities considering the network's modular structure
Pij <- PosteriorProb(M = data, 
                     R.partitions = row.Part, C.partitions = col.Part, #Input the modular structured recovered from step 2
                     Prior.Pij = "degreeprob", #Choose the null model
                     Conditional.level = "modules") #Choose the kind of constraints

#Generate randomized networks with the null model (restricted)
set.seed(1)
nulls <- RestNullModel(M = data, 
                       Pij.Prob = Pij, #Recover the probabilities calculated in the previous command
                       Numbernulls = 100, #This step may take long, so start experimenting with low values
                       Print.null = T, allow.degeneration = F, #Choose whether you allow orphan rows and columns to be removed or not
                       return.nonrm.species = F, 
                       connectance = T, byarea = T, 
                       R.partitions = row.Part, C.partitions = col.Part)

#Calculate the same nestedness metric for all randomized networks
null <- sapply(nulls, function(x) bipartite::nest.smdm(x = x, 
                                                       constraints = Part, 
                                                       weighted = F, 
                                                       decreasing = "fill"))
NODF.null <- unlist(null[3,]) # NODF for entire matrix
NODFsm.null <- unlist(null[8,]) # NODF for nodes belonging to the same module
NODFdm.null <- unlist(null[9,]) # NODF for nodes belonging to different modules

#Plot the observed nestedness value against the distribution of randomized values
png(filename= "figures/figS5_nestedness_AllMBSDs_aggregated.png", res= 300, height= 1800, width= 6000)
par(mfrow = c(1,3))
plot(density(NODF.null), xlim=c(min(obs[3], min(NODF.null)), max(obs[3], max(NODF.null))), 
     main="Observed vs. randomized", xlab = "NODF matrix")
abline(v=obs[3], col="red", lwd=2)    
plot(density(NODFsm.null), xlim=c(min(obs[8], min(NODFsm.null)), max(obs[8], max(NODFsm.null))), 
     main="Observed vs. randomized", xlab = "NODFsm matrix")
abline(v=obs[8], col="red", lwd=2)    
plot(density(NODFdm.null), xlim=c(min(obs[9], min(NODFdm.null)), max(obs[9], max(NODFdm.null))), 
     main="Observed vs. randomized", xlab = "NODFdm matrix")
abline(v=obs[9], col="red", lwd=2)    
dev.off()

#Estimate the E, Z-, P-values

#Nestedness in the entire network
mean(NODF.null) # Expected value
obs[3] - mean(NODF.null)/sd(NODF.null) # Z-value
praw.NODF <- sum(NODF.null>obs[3]) / length(NODF.null)
p.NODF <- ifelse(praw.NODF > 0.5, 1- praw.NODF, praw.NODF)    # P-value
p.NODF

#Nestedness within the modules
mean(NODFsm.null) # Expected value
obs[8] - mean(NODFsm.null)/sd(NODFsm.null) # Z-value
praw.NODFsm <- sum(NODFsm.null>obs[8]) / length(NODFsm.null)
p.NODFsm <- ifelse(praw.NODFsm > 0.5, 1- praw.NODFsm, praw.NODFsm)    # P-value
p.NODFsm

#Nestedness between the modules
mean(NODFdm.null) # Expected value
obs[9] - mean(NODFdm.null)/sd(NODFdm.null) # Z-value
praw.NODFdm <- sum(NODFdm.null>obs[9]) / length(NODFdm.null)
p.NODFdm <- ifelse(praw.NODFdm > 0.5, 1- praw.NODFdm, praw.NODFdm)    # P-value
p.NODFdm

############### 3b. NESTEDNESS ANALYSIS (FREE) ############### 

#Calculate the desired nestedness metric (here NODF) for the original network.
obs <- unlist(bipartite::nest.smdm(x = data, 
                                   constraints = Part, #Input the modular structured recovered from step 2
                                   weighted = F, #By considering the edge weights, you are choosing WNODA
                                   decreasing = "fill")) # needs to be fill for binary matrices
#Check the scores
obs

#Generate randomized networks with the null model (free)
set.seed(1)
nulls <- nullmodel(data, N = 100, method = "vaznull")

#Calculate the same nestedness metric for all randomized networks
null <- (sapply(nulls, function(x) bipartite::nest.smdm(x = x,
                                                        constraints = Part,
                                                        weighted = F,
                                                        decreasing = "fill")))
NODF.null <- unlist(null[3,]) # NODF for entire matrix
NODFsm.null <- unlist(null[8,]) # NODF for nodes belonging to the same module
NODFdm.null <- unlist(null[9,]) # NODF for nodes belonging to different modules

#Plot the observed nestedness value against the distribution of randomized values
png(filename= "figures/figS5_nestednessFree_allMBSDs_aggregated.png", res= 300, height= 1800, width= 6000)
par(mfrow = c(1,3))
plot(density(NODF.null), xlim=c(min(obs[3], min(NODF.null)), max(obs[3], max(NODF.null))), 
     main="Observed vs. randomized", xlab = "NODF matrix")
abline(v=obs[3], col="red", lwd=2)    
plot(density(NODFsm.null), xlim=c(min(obs[8], min(NODFsm.null)), max(obs[8], max(NODFsm.null))), 
     main="Observed vs. randomized", xlab = "NODFsm matrix")
abline(v=obs[8], col="red", lwd=2)    
plot(density(NODFdm.null), xlim=c(min(obs[9], min(NODFdm.null)), max(obs[9], max(NODFdm.null))), 
     main="Observed vs. randomized", xlab = "NODFdm matrix")
abline(v=obs[9], col="red", lwd=2)    
dev.off()

#Estimate the E, Z-, P-values

#Nestedness in the entire network
mean(NODF.null) # Expected value
obs[3] - mean(NODF.null)/sd(NODF.null) # Z-value
praw.NODF <- sum(NODF.null>obs[3]) / length(NODF.null)
p.NODF <- ifelse(praw.NODF > 0.5, 1- praw.NODF, praw.NODF)    # P-value
p.NODF

#Nestedness within the modules
mean(NODFsm.null) # Expected value
obs[8] - mean(NODFsm.null)/sd(NODFsm.null) # Z-value
praw.NODFsm <- sum(NODFsm.null>obs[8]) / length(NODFsm.null)
p.NODFsm <- ifelse(praw.NODFsm > 0.5, 1- praw.NODFsm, praw.NODFsm)    # P-value
p.NODFsm

#Nestedness between the modules
mean(NODFdm.null) # Expected value
obs[9] - mean(NODFdm.null)/sd(NODFdm.null) # Z-value
praw.NODFdm <- sum(NODFdm.null>obs[9]) / length(NODFdm.null)
p.NODFdm <- ifelse(praw.NODFdm > 0.5, 1- praw.NODFdm, praw.NODFdm)    # P-value
p.NODFdm


############### 4. PLOTTING THE NETWORK ############### 
par(mfrow = c(1,1))
#Sort the matrix in a way that facilitates visualizing the compound topology
data.comp <- bipartite::sortmatrix(matrix = data, topology = "compound", sort_by = "weights", row_partitions = row.Part, col_partitions = col.Part)

#Assign colors for the modules
modcol <- rainbow((length(unique(Part))), alpha=1)

#Plot the matrix
png(filename= "figures/figSX_compound_all_aggregated.png", res= 300, height= 2000, width= 3700)
plotmatrix(data.comp$matrix, 
           row_partitions = data.comp$row_partitions, 
           col_partitions = data.comp$col_partitions, 
           border = T,
           binary = F,
           modules_colors = modcol,
           within_color = "black", 
           between_color = "lightgrey",
           plot_labels = T,
           cex.axis = 0.05, srt=90)
dev.off()

############## ANALYSIS 1.1 ALL MBSDs: NEOTROPICAL ##############
############## 1. PREPARING THE DATA ##############

#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(bipartite)
library(igraph)
source("R/RestNullModel.R")
source("R/PosteriorProb.R")

# Load edge and node lists
nodes = read.delim("data/net1nodes_all_neotropical.txt", header = T) # node traits
links = read.delim("data/list.2MBSDs_neotropical_binary", header = T) # list of interactions
# Load network
all_neotrop <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
all_neotrop = get.adjacency(all_neotrop,sparse=FALSE)
all_neotrop = all_neotrop[c(0:36), ] # keep only rows of bees
all_neotrop = all_neotrop[, 37:ncol(all_neotrop)] # keep only columns of plants
data <- empty(all_neotrop, count = F) # remove nodes with no interactions

#Visualize the raw matrix
visweb(data)

############### 2. MODULARITY ANALYSIS ############### 

#Compute modularity
Mod <- bipartite::computeModules(data)
Mod@likelihood

#Recover the partitions
Part <- bipartite::module2constraints(Mod)
row.Part <- Part[1:nrow(data)]
col.Part <- Part[(nrow(data)+1):(nrow(data)+ncol(data))]

#Test for the significance of modularity with a Monte Carlo procedure

#Generate randomized matrices
set.seed(1)
nulls <- nullmodel(data, N=10, method="r2d")

#Calculate the modularity of the randomized matrices
mod.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(mod.nulls, function(x) x@likelihood)

#Calculate the z-score of the randomized distribution
(z <- (Mod@likelihood - mean(like.nulls))/sd(like.nulls))

#Plot the observed modularity value against the distribution of randomized values
png(filename= "figures/figS5_modularity_AllMBSDs_neotropical.png", res= 300, height= 3000, width= 3000)
plot(density(like.nulls), xlim=c(min((Mod@likelihood), min(like.nulls)), max((Mod@likelihood), max(like.nulls))), 
     main="Observed vs. randomized")
abline(v=(Mod@likelihood), col="red", lwd=2)    
dev.off()

#Estimate the P-value
mean(like.nulls) # E (expected modularity)
z # z
sd(like.nulls)
Mod@likelihood
praw <- sum(like.nulls>(Mod@likelihood)) / length(like.nulls)
ifelse(praw > 0.5, 1-praw, praw)


############### 3a. NESTEDNESS ANALYSIS (RESTRICTED) ############### 

#Calculate the desired nestedness metric (here NODF) for the original network.
obs <- unlist(bipartite::nest.smdm(x = data, 
                                   constraints = Part, #Input the modular structured recovered from step 2
                                   weighted = F, #By considering the edge weights, you are choosing WNODA
                                   decreasing = "fill")) # needs to be fill for binary matrices
#Check the scores
obs

#Calculate constrained interaction probabilities considering the network's modular structure
Pij <- PosteriorProb(M = data, 
                     R.partitions = row.Part, C.partitions = col.Part, #Input the modular structured recovered from step 2
                     Prior.Pij = "degreeprob", #Choose the null model
                     Conditional.level = "modules") #Choose the kind of constraints

#Generate randomized networks with the null model (restricted)
set.seed(1)
nulls <- RestNullModel(M = data, 
                       Pij.Prob = Pij, #Recover the probabilities calculated in the previous command
                       Numbernulls = 100, #This step may take long, so start experimenting with low values
                       Print.null = T, allow.degeneration = F, #Choose whether you allow orphan rows and columns to be removed or not
                       return.nonrm.species = F, 
                       connectance = T, byarea = T, 
                       R.partitions = row.Part, C.partitions = col.Part)

#Calculate the same nestedness metric for all randomized networks
null <- sapply(nulls, function(x) bipartite::nest.smdm(x = x, 
                                                       constraints = Part, 
                                                       weighted = F, 
                                                       decreasing = "fill"))
NODF.null <- unlist(null[3,]) # NODF for entire matrix
NODFsm.null <- unlist(null[8,]) # NODF for nodes belonging to the same module
NODFdm.null <- unlist(null[9,]) # NODF for nodes belonging to different modules

#Plot the observed nestedness value against the distribution of randomized values
png(filename= "figures/figS5_nestedness_all_neotropical.png", res= 300, height= 1800, width= 6000)
par(mfrow = c(1,3))
plot(density(NODF.null), xlim=c(min(obs[3], min(NODF.null)), max(obs[3], max(NODF.null))), 
     main="Observed vs. randomized", xlab = "NODF matrix")
abline(v=obs[3], col="red", lwd=2)    
plot(density(NODFsm.null), xlim=c(min(obs[8], min(NODFsm.null)), max(obs[8], max(NODFsm.null))), 
     main="Observed vs. randomized", xlab = "NODFsm matrix")
abline(v=obs[8], col="red", lwd=2)    
plot(density(NODFdm.null), xlim=c(min(obs[9], min(NODFdm.null)), max(obs[9], max(NODFdm.null))), 
     main="Observed vs. randomized", xlab = "NODFdm matrix")
abline(v=obs[9], col="red", lwd=2)
dev.off()

#Estimate the E, Z-, P-values

#Nestedness in the entire network
mean(NODF.null) # Expected value
obs[3] - mean(NODF.null)/sd(NODF.null) # Z-value
praw.NODF <- sum(NODF.null>obs[3]) / length(NODF.null)
p.NODF <- ifelse(praw.NODF > 0.5, 1- praw.NODF, praw.NODF)    # P-value
p.NODF

#Nestedness within the modules
mean(NODFsm.null) # Expected value
obs[8] - mean(NODFsm.null)/sd(NODFsm.null) # Z-value
praw.NODFsm <- sum(NODFsm.null>obs[8]) / length(NODFsm.null)
p.NODFsm <- ifelse(praw.NODFsm > 0.5, 1- praw.NODFsm, praw.NODFsm)    # P-value
p.NODFsm

#Nestedness between the modules
mean(NODFdm.null) # Expected value
obs[9] - mean(NODFdm.null)/sd(NODFdm.null) # Z-value
praw.NODFdm <- sum(NODFdm.null>obs[9]) / length(NODFdm.null)
p.NODFdm <- ifelse(praw.NODFdm > 0.5, 1- praw.NODFdm, praw.NODFdm)    # P-value
p.NODFdm

############### 3b. NESTEDNESS ANALYSIS (FREE) ############### 

#Calculate the desired nestedness metric (here NODF) for the original network.
obs <- unlist(bipartite::nest.smdm(x = data, 
                                   constraints = Part, #Input the modular structured recovered from step 2
                                   weighted = F, #By considering the edge weights, you are choosing WNODA
                                   decreasing = "fill")) # needs to be fill for binary matrices
#Check the scores
obs

#Generate randomized networks with the null model (free)
set.seed(1)
nulls <- nullmodel(data, N = 100, method = "vaznull")

#Calculate the same nestedness metric for all randomized networks
null <- (sapply(nulls, function(x) bipartite::nest.smdm(x = x,
                                                        constraints = Part,
                                                        weighted = F,
                                                        decreasing = "fill")))
NODF.null <- unlist(null[3,]) # NODF for entire matrix
NODFsm.null <- unlist(null[8,]) # NODF for nodes belonging to the same module
NODFdm.null <- unlist(null[9,]) # NODF for nodes belonging to different modules

#Plot the observed nestedness value against the distribution of randomized values
png(filename= "figures/figS5_nestednessFree_all_neotropical.png", res= 300, height= 1800, width= 6000)
par(mfrow = c(1,3))
plot(density(NODF.null), xlim=c(min(obs[3], min(NODF.null)), max(obs[3], max(NODF.null))), 
     main="Observed vs. randomized", xlab = "NODF matrix")
abline(v=obs[3], col="red", lwd=2)    
plot(density(NODFsm.null), xlim=c(min(obs[8], min(NODFsm.null)), max(obs[8], max(NODFsm.null))), 
     main="Observed vs. randomized", xlab = "NODFsm matrix")
abline(v=obs[8], col="red", lwd=2)    
plot(density(NODFdm.null), xlim=c(min(obs[9], min(NODFdm.null)), max(obs[9], max(NODFdm.null))), 
     main="Observed vs. randomized", xlab = "NODFdm matrix")
abline(v=obs[9], col="red", lwd=2)    
dev.off()

#Estimate the E, Z-, P-values

#Nestedness in the entire network
mean(NODF.null) # Expected value
obs[3] - mean(NODF.null)/sd(NODF.null) # Z-value
praw.NODF <- sum(NODF.null>obs[3]) / length(NODF.null)
p.NODF <- ifelse(praw.NODF > 0.5, 1- praw.NODF, praw.NODF)    # P-value
p.NODF

#Nestedness within the modules
mean(NODFsm.null) # Expected value
obs[8] - mean(NODFsm.null)/sd(NODFsm.null) # Z-value
praw.NODFsm <- sum(NODFsm.null>obs[8]) / length(NODFsm.null)
p.NODFsm <- ifelse(praw.NODFsm > 0.5, 1- praw.NODFsm, praw.NODFsm)    # P-value
p.NODFsm

#Nestedness between the modules
mean(NODFdm.null) # Expected value
obs[9] - mean(NODFdm.null)/sd(NODFdm.null) # Z-value
praw.NODFdm <- sum(NODFdm.null>obs[9]) / length(NODFdm.null)
p.NODFdm <- ifelse(praw.NODFdm > 0.5, 1- praw.NODFdm, praw.NODFdm)    # P-value
p.NODFdm


############### 4. PLOTTING THE NETWORK ############### 
par(mfrow = c(1,1))
#Sort the matrix in a way that facilitates visualizing the compound topology
data.comp <- bipartite::sortmatrix(matrix = data, topology = "compound", sort_by = "weights", row_partitions = row.Part, col_partitions = col.Part)

#Assign colors for the modules
modcol <- rainbow((length(unique(Part))), alpha=1)

#Plot the matrix
png(filename= "figures/figSX_compound_all_neotropical.png", res= 300, height= 2000, width= 3500)
plotmatrix(data.comp$matrix, 
           row_partitions = data.comp$row_partitions, 
           col_partitions = data.comp$col_partitions, 
           border = T,
           binary = F,
           modules_colors = modcol,
           within_color = "black", 
           between_color = "lightgrey")
dev.off()
############## ANALYSIS 1.2 ALL MBSDs: INDO-MALAYAN-AUSTRALASIA ##############
############## 1. PREPARING THE DATA ##############

#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(bipartite)
library(igraph)
source("R/RestNullModel.R")
source("R/PosteriorProb.R")

# Load edge and node lists
nodes = read.delim("data/net1nodes_all_asia.txt", header = T) # node traits
links = read.delim("data/list.2MBSDs_asia_binary", header = T) # list of interactions
# Load network
all_asia <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
all_asia = get.adjacency(all_asia,sparse=FALSE)
all_asia = all_asia[c(0:31), ] # keep only rows of bees
all_asia = all_asia[, 32:ncol(all_asia)] # keep only columns of plants
data <- empty(all_asia, count = F) # remove nodes with no interactions

#Visualize the raw matrix
visweb(data)

############### 2. MODULARITY ANALYSIS ############### 

#Compute modularity
Mod <- bipartite::computeModules(data)
Mod@likelihood

#Recover the partitions
Part <- bipartite::module2constraints(Mod)
row.Part <- Part[1:nrow(data)]
col.Part <- Part[(nrow(data)+1):(nrow(data)+ncol(data))]

#Test for the significance of modularity with a Monte Carlo procedure

#Generate randomized matrices
set.seed(1)
nulls <- nullmodel(data, N=10, method="r2d")

#Calculate the modularity of the randomized matrices
mod.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(mod.nulls, function(x) x@likelihood)

#Calculate the z-score of the randomized distribution
(z <- (Mod@likelihood - mean(like.nulls))/sd(like.nulls))

#Plot the observed modularity value against the distribution of randomized values
png(filename= "figures/figS5_modularity_AllMBSDs_asia.png", res= 300, height= 3000, width= 3000)
plot(density(like.nulls), xlim=c(min((Mod@likelihood), min(like.nulls)), max((Mod@likelihood), max(like.nulls))), 
     main="Observed vs. randomized")
abline(v=(Mod@likelihood), col="red", lwd=2)    
dev.off()

#Estimate the P-value
mean(like.nulls) # E (expected modularity)
z # z
sd(like.nulls)
Mod@likelihood
praw <- sum(like.nulls>(Mod@likelihood)) / length(like.nulls)
ifelse(praw > 0.5, 1-praw, praw)


############### 3a. NESTEDNESS ANALYSIS (RESTRICTED) ############### 

#Calculate the desired nestedness metric (here NODF) for the original network.
obs <- unlist(bipartite::nest.smdm(x = data, 
                                   constraints = Part, #Input the modular structured recovered from step 2
                                   weighted = F, #By considering the edge weights, you are choosing WNODA
                                   decreasing = "fill")) # needs to be fill for binary matrices
#Check the scores
obs

#Calculate constrained interaction probabilities considering the network's modular structure
Pij <- PosteriorProb(M = data, 
                     R.partitions = row.Part, C.partitions = col.Part, #Input the modular structured recovered from step 2
                     Prior.Pij = "degreeprob", #Choose the null model
                     Conditional.level = "modules") #Choose the kind of constraints

#Generate randomized networks with the null model (restricted)
set.seed(1)
nulls <- RestNullModel(M = data, 
                       Pij.Prob = Pij, #Recover the probabilities calculated in the previous command
                       Numbernulls = 1000, #This step may take long, so start experimenting with low values
                       Print.null = T, allow.degeneration = F, #Choose whether you allow orphan rows and columns to be removed or not
                       return.nonrm.species = F, 
                       connectance = T, byarea = T, 
                       R.partitions = row.Part, C.partitions = col.Part)

#Calculate the same nestedness metric for all randomized networks
null <- sapply(nulls, function(x) bipartite::nest.smdm(x = x, 
                                                       constraints = Part, 
                                                       weighted = F, 
                                                       decreasing = "fill"))
NODF.null <- unlist(null[3,]) # NODF for entire matrix
NODFsm.null <- unlist(null[8,]) # NODF for nodes belonging to the same module
NODFdm.null <- unlist(null[9,]) # NODF for nodes belonging to different modules

#Plot the observed nestedness value against the distribution of randomized values
png(filename= "figures/figS5_nestedness_all_asia.png", res= 300, height= 1800, width= 6000)
par(mfrow = c(1,3))
plot(density(NODF.null), xlim=c(min(obs[3], min(NODF.null)), max(obs[3], max(NODF.null))), 
     main="Observed vs. randomized", xlab = "NODF matrix")
abline(v=obs[3], col="red", lwd=2)    
plot(density(NODFsm.null), xlim=c(min(obs[8], min(NODFsm.null)), max(obs[8], max(NODFsm.null))), 
     main="Observed vs. randomized", xlab = "NODFsm matrix")
abline(v=obs[8], col="red", lwd=2)    
plot(density(NODFdm.null), xlim=c(min(obs[9], min(NODFdm.null)), max(obs[9], max(NODFdm.null))), 
     main="Observed vs. randomized", xlab = "NODFdm matrix")
abline(v=obs[9], col="red", lwd=2)
dev.off()

#Estimate the E, Z-, P-values

#Nestedness in the entire network
mean(NODF.null) # Expected value
obs[3] - mean(NODF.null)/sd(NODF.null) # Z-value
praw.NODF <- sum(NODF.null>obs[3]) / length(NODF.null)
p.NODF <- ifelse(praw.NODF > 0.5, 1- praw.NODF, praw.NODF)    # P-value
p.NODF

#Nestedness within the modules
mean(NODFsm.null) # Expected value
obs[8] - mean(NODFsm.null)/sd(NODFsm.null) # Z-value
praw.NODFsm <- sum(NODFsm.null>obs[8]) / length(NODFsm.null)
p.NODFsm <- ifelse(praw.NODFsm > 0.5, 1- praw.NODFsm, praw.NODFsm)    # P-value
p.NODFsm

#Nestedness between the modules
mean(NODFdm.null) # Expected value
obs[9] - mean(NODFdm.null)/sd(NODFdm.null) # Z-value
praw.NODFdm <- sum(NODFdm.null>obs[9]) / length(NODFdm.null)
p.NODFdm <- ifelse(praw.NODFdm > 0.5, 1- praw.NODFdm, praw.NODFdm)    # P-value
p.NODFdm

############### 3b. NESTEDNESS ANALYSIS (FREE) ############### 

#Calculate the desired nestedness metric (here NODF) for the original network.
obs <- unlist(bipartite::nest.smdm(x = data, 
                                   constraints = Part, #Input the modular structured recovered from step 2
                                   weighted = F, #By considering the edge weights, you are choosing WNODA
                                   decreasing = "fill")) # needs to be fill for binary matrices
#Check the scores
obs

#Generate randomized networks with the null model (free)
set.seed(1)
nulls <- nullmodel(data, N = 100, method = "vaznull")

#Calculate the same nestedness metric for all randomized networks
null <- (sapply(nulls, function(x) bipartite::nest.smdm(x = x,
                                                        constraints = Part,
                                                        weighted = F,
                                                        decreasing = "fill")))
NODF.null <- unlist(null[3,]) # NODF for entire matrix
NODFsm.null <- unlist(null[8,]) # NODF for nodes belonging to the same module
NODFdm.null <- unlist(null[9,]) # NODF for nodes belonging to different modules

#Plot the observed nestedness value against the distribution of randomized values
png(filename= "figures/figS5_nestednessFree_all_asia.png", res= 300, height= 1800, width= 6000)
par(mfrow = c(1,3))
plot(density(NODF.null), xlim=c(min(obs[3], min(NODF.null)), max(obs[3], max(NODF.null))), 
     main="Observed vs. randomized", xlab = "NODF matrix")
abline(v=obs[3], col="red", lwd=2)    
plot(density(NODFsm.null), xlim=c(min(obs[8], min(NODFsm.null)), max(obs[8], max(NODFsm.null))), 
     main="Observed vs. randomized", xlab = "NODFsm matrix")
abline(v=obs[8], col="red", lwd=2)    
plot(density(NODFdm.null), xlim=c(min(obs[9], min(NODFdm.null)), max(obs[9], max(NODFdm.null))), 
     main="Observed vs. randomized", xlab = "NODFdm matrix")
abline(v=obs[9], col="red", lwd=2)    
dev.off()

#Estimate the E, Z-, P-values

#Nestedness in the entire network
mean(NODF.null) # Expected value
obs[3] - mean(NODF.null)/sd(NODF.null) # Z-value
praw.NODF <- sum(NODF.null>obs[3]) / length(NODF.null)
p.NODF <- ifelse(praw.NODF > 0.5, 1- praw.NODF, praw.NODF)    # P-value
p.NODF

#Nestedness within the modules
mean(NODFsm.null) # Expected value
obs[8] - mean(NODFsm.null)/sd(NODFsm.null) # Z-value
praw.NODFsm <- sum(NODFsm.null>obs[8]) / length(NODFsm.null)
p.NODFsm <- ifelse(praw.NODFsm > 0.5, 1- praw.NODFsm, praw.NODFsm)    # P-value
p.NODFsm

#Nestedness between the modules
mean(NODFdm.null) # Expected value
obs[9] - mean(NODFdm.null)/sd(NODFdm.null) # Z-value
praw.NODFdm <- sum(NODFdm.null>obs[9]) / length(NODFdm.null)
p.NODFdm <- ifelse(praw.NODFdm > 0.5, 1- praw.NODFdm, praw.NODFdm)    # P-value
p.NODFdm


############### 4. PLOTTING THE NETWORK ############### 
par(mfrow = c(1,1))
#Sort the matrix in a way that facilitates visualizing the compound topology
data.comp <- bipartite::sortmatrix(matrix = data, topology = "compound", sort_by = "weights", row_partitions = row.Part, col_partitions = col.Part)

#Assign colors for the modules
modcol <- rainbow((length(unique(Part))), alpha=1)

#Plot the matrix
png(filename= "figures/figSX_compound_all_asia.png", res= 300, height= 2500, width= 3500)
plotmatrix(data.comp$matrix, 
           row_partitions = data.comp$row_partitions, 
           col_partitions = data.comp$col_partitions, 
           border = T,
           binary = F,
           modules_colors = modcol,
           within_color = "black", 
           between_color = "lightgrey")
dev.off()

############## ANALYSIS 2.0 CHEMICAL+FIELDWORK: AGGREGATED ##############
############## 1. PREPARING THE DATA ##############

#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(bipartite)
library(igraph)
source("R/RestNullModel.R")
source("R/PosteriorProb.R")

# Load edge and node lists
nodes = read.delim("data/net1nodes_CF.txt", header = T) # node traits
links = read.delim("data/list_CF", header = T) # list of interactions
# Load network
data <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
data = get.adjacency(data,sparse=FALSE)
data = data[c(0:62), ] # keep rows with bees
data = data[, 63:ncol(data)] # keep only columns of plants
data <- empty(data, count = F) # remove nodes with no interactions

#Visualize the raw matrix
visweb(data)

############### 2. MODULARITY ANALYSIS ############### 

#Compute modularity
Mod <- bipartite::computeModules(data)
Mod@likelihood

#Recover the partitions
Part <- bipartite::module2constraints(Mod)
row.Part <- Part[1:nrow(data)]
col.Part <- Part[(nrow(data)+1):(nrow(data)+ncol(data))]

#Test for the significance of modularity with a Monte Carlo procedure

#Generate randomized matrices
set.seed(1)
nulls <- nullmodel(data, N=10, method="r2d")

#Calculate the modularity of the randomized matrices
mod.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(mod.nulls, function(x) x@likelihood)

#Calculate the z-score of the randomized distribution
(z <- (Mod@likelihood - mean(like.nulls))/sd(like.nulls))

#Plot the observed modularity value against the distribution of randomized values
png(filename= "figures/figS5_modularity_CF_aggregated.png", res= 300, height= 3000, width= 3000)
plot(density(like.nulls), xlim=c(min((Mod@likelihood), min(like.nulls)), max((Mod@likelihood), max(like.nulls))), 
     main="Observed vs. randomized")
abline(v=(Mod@likelihood), col="red", lwd=2)
dev.off()

#Estimate the P-value
Mod@likelihood # Obs. value
mean(like.nulls) # E (expected modularity)
z # z
sd(like.nulls)
praw <- sum(like.nulls>(Mod@likelihood)) / length(like.nulls)
ifelse(praw > 0.5, 1-praw, praw)

############### 3a. NESTEDNESS ANALYSIS (RESTRICTED) ############### 

#Calculate the desired nestedness metric (here NODF) for the original network.
obs <- unlist(bipartite::nest.smdm(x = data, 
                                   constraints = Part, #Input the modular structured recovered from step 2
                                   weighted = F, #By considering the edge weights, you are choosing WNODA
                                   decreasing = "fill")) # needs to be fill for binary matrices
#Check the scores
obs

#Calculate constrained interaction probabilities considering the network's modular structure
Pij <- PosteriorProb(M = data, 
                     R.partitions = row.Part, C.partitions = col.Part, #Input the modular structured recovered from step 2
                     Prior.Pij = "degreeprob", #Choose the null model
                     Conditional.level = "modules") #Choose the kind of constraints

#Generate randomized networks with the null model (restricted)
set.seed(1)
nulls <- RestNullModel(M = data, 
                       Pij.Prob = Pij, #Recover the probabilities calculated in the previous command
                       Numbernulls = 100, #This step may take long, so start experimenting with low values
                       Print.null = T, allow.degeneration = F, #Choose whether you allow orphan rows and columns to be removed or not
                       return.nonrm.species = F, 
                       connectance = T, byarea = T, 
                       R.partitions = row.Part, C.partitions = col.Part)

#Calculate the same nestedness metric for all randomized networks
null <- sapply(nulls, function(x) bipartite::nest.smdm(x = x, 
                                                       constraints = Part, 
                                                       weighted = F, 
                                                       decreasing = "fill"))
NODF.null <- unlist(null[3,]) # NODF for entire matrix
NODFsm.null <- unlist(null[8,]) # NODF for nodes belonging to the same module
NODFdm.null <- unlist(null[9,]) # NODF for nodes belonging to different modules

#Plot the observed nestedness value against the distribution of randomized values
png(filename= "figures/figS5_nestedness_CF_aggregated.png", res= 300, height= 1800, width= 6000)
par(mfrow = c(1,3))
plot(density(NODF.null), xlim=c(min(obs[3], min(NODF.null)), max(obs[3], max(NODF.null))), 
     main="Observed vs. randomized", xlab = "NODF matrix")
abline(v=obs[3], col="red", lwd=2)    
plot(density(NODFsm.null), xlim=c(min(obs[8], min(NODFsm.null)), max(obs[8], max(NODFsm.null))), 
     main="Observed vs. randomized", xlab = "NODFsm matrix")
abline(v=obs[8], col="red", lwd=2)    
plot(density(NODFdm.null), xlim=c(min(obs[9], min(NODFdm.null)), max(obs[9], max(NODFdm.null))), 
     main="Observed vs. randomized", xlab = "NODFdm matrix")
abline(v=obs[9], col="red", lwd=2)    
dev.off()

#Estimate the E, Z-, P-values

#Nestedness in the entire network
mean(NODF.null) # Expected value
obs[3] - mean(NODF.null)/sd(NODF.null) # Z-value
praw.NODF <- sum(NODF.null>obs[3]) / length(NODF.null)
p.NODF <- ifelse(praw.NODF > 0.5, 1- praw.NODF, praw.NODF)    # P-value
p.NODF

#Nestedness within the modules
mean(NODFsm.null) # Expected value
obs[8] - mean(NODFsm.null)/sd(NODFsm.null) # Z-value
praw.NODFsm <- sum(NODFsm.null>obs[8]) / length(NODFsm.null)
p.NODFsm <- ifelse(praw.NODFsm > 0.5, 1- praw.NODFsm, praw.NODFsm)    # P-value
p.NODFsm

#Nestedness between the modules
mean(NODFdm.null) # Expected value
obs[9] - mean(NODFdm.null)/sd(NODFdm.null) # Z-value
praw.NODFdm <- sum(NODFdm.null>obs[9]) / length(NODFdm.null)
p.NODFdm <- ifelse(praw.NODFdm > 0.5, 1- praw.NODFdm, praw.NODFdm)    # P-value
p.NODFdm


############### 3b. NESTEDNESS ANALYSIS (FREE) ############### 

#Calculate the desired nestedness metric (here NODF) for the original network.
obs <- unlist(bipartite::nest.smdm(x = data, 
                                   constraints = Part, #Input the modular structured recovered from step 2
                                   weighted = F, #By considering the edge weights, you are choosing WNODA
                                   decreasing = "fill")) # needs to be fill for binary matrices
#Check the scores
obs

#Generate randomized networks with the null model (free)
set.seed(1)
nulls <- nullmodel(data, N = 100, method = "vaznull")

#Calculate the same nestedness metric for all randomized networks
null <- (sapply(nulls, function(x) bipartite::nest.smdm(x = x,
                                                        constraints = Part,
                                                        weighted = F,
                                                        decreasing = "fill")))
NODF.null <- unlist(null[3,]) # NODF for entire matrix
NODFsm.null <- unlist(null[8,]) # NODF for nodes belonging to the same module
NODFdm.null <- unlist(null[9,]) # NODF for nodes belonging to different modules

#Plot the observed nestedness value against the distribution of randomized values
png(filename= "figures/figS5_nestednessFree_CF_aggregated.png", res= 300, height= 1800, width= 6000)
par(mfrow = c(1,3))
plot(density(NODF.null), xlim=c(min(obs[3], min(NODF.null)), max(obs[3], max(NODF.null))), 
     main="Observed vs. randomized", xlab = "NODF matrix")
abline(v=obs[3], col="red", lwd=2)    
plot(density(NODFsm.null), xlim=c(min(obs[8], min(NODFsm.null)), max(obs[8], max(NODFsm.null))), 
     main="Observed vs. randomized", xlab = "NODFsm matrix")
abline(v=obs[8], col="red", lwd=2)    
plot(density(NODFdm.null), xlim=c(min(obs[9], min(NODFdm.null)), max(obs[9], max(NODFdm.null))), 
     main="Observed vs. randomized", xlab = "NODFdm matrix")
abline(v=obs[9], col="red", lwd=2)    
dev.off()

#Estimate the E, Z-, P-values

#Nestedness in the entire network
mean(NODF.null) # Expected value
obs[3] - mean(NODF.null)/sd(NODF.null) # Z-value
praw.NODF <- sum(NODF.null>obs[3]) / length(NODF.null)
p.NODF <- ifelse(praw.NODF > 0.5, 1- praw.NODF, praw.NODF)    # P-value
p.NODF

#Nestedness within the modules
mean(NODFsm.null) # Expected value
obs[8] - mean(NODFsm.null)/sd(NODFsm.null) # Z-value
praw.NODFsm <- sum(NODFsm.null>obs[8]) / length(NODFsm.null)
p.NODFsm <- ifelse(praw.NODFsm > 0.5, 1- praw.NODFsm, praw.NODFsm)    # P-value
p.NODFsm

#Nestedness between the modules
mean(NODFdm.null) # Expected value
obs[9] - mean(NODFdm.null)/sd(NODFdm.null) # Z-value
praw.NODFdm <- sum(NODFdm.null>obs[9]) / length(NODFdm.null)
p.NODFdm <- ifelse(praw.NODFdm > 0.5, 1- praw.NODFdm, praw.NODFdm)    # P-value
p.NODFdm


############### 4. PLOTTING THE NETWORK ############### 
par(mfrow = c(1,1))
#Sort the matrix in a way that facilitates visualizing the compound topology
data.comp <- bipartite::sortmatrix(matrix = data, topology = "compound", sort_by = "weights", row_partitions = row.Part, col_partitions = col.Part)

#Assign colors for the modules
modcol <- rainbow((length(unique(Part))), alpha=1)

#Plot the matrix
png(filename= "figures/figSX_compound_CF.png", res= 300, height= 2000, width= 3700)
plotmatrix(data.comp$matrix, 
           row_partitions = data.comp$row_partitions, 
           col_partitions = data.comp$col_partitions, 
           border = T,
           binary = F,
           modules_colors = modcol,
           within_color = modcol, 
           between_color = "lightgrey")
dev.off()

############## ANALYSIS 2.1 CHEMICAL+FIELDWORK: NEOTROPICAL ##############
############## 1. PREPARING THE DATA ##############

#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(bipartite)
library(igraph)
source("R/RestNullModel.R")
source("R/PosteriorProb.R")

# Load edge and node lists
nodes = read.delim("data/net1nodes_CF_neotropical.txt", header = T) # node traits
links = read.delim("data/list_CF_neotropical", header = T) # list of interactions
# Load network
data <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
data = get.adjacency(data,sparse=FALSE)
data = data[c(0:30), ] # keep rows with bees
data = data[, 31:ncol(data)] # keep only columns of plants
data <- empty(data, count = F) # remove nodes with no interactions

#Visualize the raw matrix
visweb(data)

############### 2. MODULARITY ANALYSIS ############### 

#Compute modularity
Mod <- bipartite::computeModules(data)
Mod@likelihood

#Recover the partitions
Part <- bipartite::module2constraints(Mod)
row.Part <- Part[1:nrow(data)]
col.Part <- Part[(nrow(data)+1):(nrow(data)+ncol(data))]

#Test for the significance of modularity with a Monte Carlo procedure

#Generate randomized matrices
set.seed(1)
nulls <- nullmodel(data, N=10, method="r2d")

#Calculate the modularity of the randomized matrices
mod.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(mod.nulls, function(x) x@likelihood)

#Calculate the z-score of the randomized distribution
(z <- (Mod@likelihood - mean(like.nulls))/sd(like.nulls))

#Plot the observed modularity value against the distribution of randomized values
png(filename= "figures/figS5_modularity_CF_neotropical.png", res= 300, height= 3000, width= 3000)
plot(density(like.nulls), xlim=c(min((Mod@likelihood), min(like.nulls)), max((Mod@likelihood), max(like.nulls))), 
     main="Observed vs. randomized")
abline(v=(Mod@likelihood), col="red", lwd=2)
dev.off()

#Estimate the P-value
Mod@likelihood # Obs. value
mean(like.nulls) # E (expected modularity)
z # z
sd(like.nulls)
praw <- sum(like.nulls>(Mod@likelihood)) / length(like.nulls)
ifelse(praw > 0.5, 1-praw, praw)

############### 3a. NESTEDNESS ANALYSIS (RESTRICTED) ############### 

#Calculate the desired nestedness metric (here NODF) for the original network.
obs <- unlist(bipartite::nest.smdm(x = data, 
                                   constraints = Part, #Input the modular structured recovered from step 2
                                   weighted = F, #By considering the edge weights, you are choosing WNODA
                                   decreasing = "fill")) # needs to be fill for binary matrices
#Check the scores
obs

#Calculate constrained interaction probabilities considering the network's modular structure
Pij <- PosteriorProb(M = data, 
                     R.partitions = row.Part, C.partitions = col.Part, #Input the modular structured recovered from step 2
                     Prior.Pij = "degreeprob", #Choose the null model
                     Conditional.level = "modules") #Choose the kind of constraints

#Generate randomized networks with the null model (restricted)
set.seed(1)
nulls <- RestNullModel(M = data, 
                       Pij.Prob = Pij, #Recover the probabilities calculated in the previous command
                       Numbernulls = 100, #This step may take long, so start experimenting with low values
                       Print.null = T, allow.degeneration = F, #Choose whether you allow orphan rows and columns to be removed or not
                       return.nonrm.species = F, 
                       connectance = T, byarea = T, 
                       R.partitions = row.Part, C.partitions = col.Part)

#Calculate the same nestedness metric for all randomized networks
null <- sapply(nulls, function(x) bipartite::nest.smdm(x = x, 
                                                       constraints = Part, 
                                                       weighted = F, 
                                                       decreasing = "fill"))
NODF.null <- unlist(null[3,]) # NODF for entire matrix
NODFsm.null <- unlist(null[8,]) # NODF for nodes belonging to the same module
NODFdm.null <- unlist(null[9,]) # NODF for nodes belonging to different modules

#Plot the observed nestedness value against the distribution of randomized values
png(filename= "figures/figS5_nestedness_CF_neotropical.png", res= 300, height= 1800, width= 6000)
par(mfrow = c(1,3))
plot(density(NODF.null), xlim=c(min(obs[3], min(NODF.null)), max(obs[3], max(NODF.null))), 
     main="Observed vs. randomized", xlab = "NODF matrix")
abline(v=obs[3], col="red", lwd=2)    
plot(density(NODFsm.null), xlim=c(min(obs[8], min(NODFsm.null)), max(obs[8], max(NODFsm.null))), 
     main="Observed vs. randomized", xlab = "NODFsm matrix")
abline(v=obs[8], col="red", lwd=2)    
plot(density(NODFdm.null), xlim=c(min(obs[9], min(NODFdm.null)), max(obs[9], max(NODFdm.null))), 
     main="Observed vs. randomized", xlab = "NODFdm matrix")
abline(v=obs[9], col="red", lwd=2)    
dev.off()

#Estimate the E, Z-, P-values

#Nestedness in the entire network
mean(NODF.null) # Expected value
obs[3] - mean(NODF.null)/sd(NODF.null) # Z-value
praw.NODF <- sum(NODF.null>obs[3]) / length(NODF.null)
p.NODF <- ifelse(praw.NODF > 0.5, 1- praw.NODF, praw.NODF)    # P-value
p.NODF

#Nestedness within the modules
mean(NODFsm.null) # Expected value
obs[8] - mean(NODFsm.null)/sd(NODFsm.null) # Z-value
praw.NODFsm <- sum(NODFsm.null>obs[8]) / length(NODFsm.null)
p.NODFsm <- ifelse(praw.NODFsm > 0.5, 1- praw.NODFsm, praw.NODFsm)    # P-value
p.NODFsm

#Nestedness between the modules
mean(NODFdm.null) # Expected value
obs[9] - mean(NODFdm.null)/sd(NODFdm.null) # Z-value
praw.NODFdm <- sum(NODFdm.null>obs[9]) / length(NODFdm.null)
p.NODFdm <- ifelse(praw.NODFdm > 0.5, 1- praw.NODFdm, praw.NODFdm)    # P-value
p.NODFdm


############### 3b. NESTEDNESS ANALYSIS (FREE) ############### 

#Calculate the desired nestedness metric (here NODF) for the original network.
obs <- unlist(bipartite::nest.smdm(x = data, 
                                   constraints = Part, #Input the modular structured recovered from step 2
                                   weighted = F, #By considering the edge weights, you are choosing WNODA
                                   decreasing = "fill")) # needs to be fill for binary matrices
#Check the scores
obs

#Generate randomized networks with the null model (free)
set.seed(1)
nulls <- nullmodel(data, N = 100, method = "vaznull")

#Calculate the same nestedness metric for all randomized networks
null <- (sapply(nulls, function(x) bipartite::nest.smdm(x = x,
                                                        constraints = Part,
                                                        weighted = F,
                                                        decreasing = "fill")))
NODF.null <- unlist(null[3,]) # NODF for entire matrix
NODFsm.null <- unlist(null[8,]) # NODF for nodes belonging to the same module
NODFdm.null <- unlist(null[9,]) # NODF for nodes belonging to different modules

#Plot the observed nestedness value against the distribution of randomized values
png(filename= "figures/figS5_nestednessFree_CF_neotropical.png", res= 300, height= 1800, width= 6000)
par(mfrow = c(1,3))
plot(density(NODF.null), xlim=c(min(obs[3], min(NODF.null)), max(obs[3], max(NODF.null))), 
     main="Observed vs. randomized", xlab = "NODF matrix")
abline(v=obs[3], col="red", lwd=2)    
plot(density(NODFsm.null), xlim=c(min(obs[8], min(NODFsm.null)), max(obs[8], max(NODFsm.null))), 
     main="Observed vs. randomized", xlab = "NODFsm matrix")
abline(v=obs[8], col="red", lwd=2)    
plot(density(NODFdm.null), xlim=c(min(obs[9], min(NODFdm.null)), max(obs[9], max(NODFdm.null))), 
     main="Observed vs. randomized", xlab = "NODFdm matrix")
abline(v=obs[9], col="red", lwd=2)    
dev.off()

#Estimate the E, Z-, P-values

#Nestedness in the entire network
mean(NODF.null) # Expected value
obs[3] - mean(NODF.null)/sd(NODF.null) # Z-value
praw.NODF <- sum(NODF.null>obs[3]) / length(NODF.null)
p.NODF <- ifelse(praw.NODF > 0.5, 1- praw.NODF, praw.NODF)    # P-value
p.NODF

#Nestedness within the modules
mean(NODFsm.null) # Expected value
obs[8] - mean(NODFsm.null)/sd(NODFsm.null) # Z-value
praw.NODFsm <- sum(NODFsm.null>obs[8]) / length(NODFsm.null)
p.NODFsm <- ifelse(praw.NODFsm > 0.5, 1- praw.NODFsm, praw.NODFsm)    # P-value
p.NODFsm

#Nestedness between the modules
mean(NODFdm.null) # Expected value
obs[9] - mean(NODFdm.null)/sd(NODFdm.null) # Z-value
praw.NODFdm <- sum(NODFdm.null>obs[9]) / length(NODFdm.null)
p.NODFdm <- ifelse(praw.NODFdm > 0.5, 1- praw.NODFdm, praw.NODFdm)    # P-value
p.NODFdm


############### 4. PLOTTING THE NETWORK ############### 
par(mfrow = c(1,1))
#Sort the matrix in a way that facilitates visualizing the compound topology
data.comp <- bipartite::sortmatrix(matrix = data, topology = "compound", sort_by = "weights", row_partitions = row.Part, col_partitions = col.Part)

#Assign colors for the modules
modcol <- rainbow((length(unique(Part))), alpha=1)

#Plot the matrix
png(filename= "figures/figSX_compound_CF_neotropical.png", res= 300, height= 2500, width= 3500)
plotmatrix(data.comp$matrix, 
           row_partitions = data.comp$row_partitions, 
           col_partitions = data.comp$col_partitions, 
           border = T,
           binary = F,
           modules_colors = modcol,
           within_color = modcol, 
           between_color = "lightgrey")
dev.off()

############## ANALYSIS 2.1 CHEMICAL+FIELDWORK: ASIA ##############
############## 1. PREPARING THE DATA ##############

#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(bipartite)
library(igraph)
source("R/RestNullModel.R")
source("R/PosteriorProb.R")

# Load edge and node lists
nodes = read.delim("data/net1nodes_CF_asia.txt", header = T) # node traits
links = read.delim("data/list_CF_asia", header = T) # list of interactions
# Load network
data <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
data = get.adjacency(data,sparse=FALSE)
data = data[c(0:31), ] # keep rows with bees
data = data[, 32:ncol(data)] # keep only columns of plants
data <- empty(data, count = F) # remove nodes with no interactions

#Visualize the raw matrix
visweb(data)

############### 2. MODULARITY ANALYSIS ############### 

#Compute modularity
Mod <- bipartite::computeModules(data)
Mod@likelihood

#Recover the partitions
Part <- bipartite::module2constraints(Mod)
row.Part <- Part[1:nrow(data)]
col.Part <- Part[(nrow(data)+1):(nrow(data)+ncol(data))]

#Test for the significance of modularity with a Monte Carlo procedure

#Generate randomized matrices
set.seed(1)
nulls <- nullmodel(data, N=10, method="r2d")

#Calculate the modularity of the randomized matrices
mod.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(mod.nulls, function(x) x@likelihood)

#Calculate the z-score of the randomized distribution
(z <- (Mod@likelihood - mean(like.nulls))/sd(like.nulls))

#Plot the observed modularity value against the distribution of randomized values
png(filename= "figures/figS5_modularity_CF_asia.png", res= 300, height= 3000, width= 3000)
plot(density(like.nulls), xlim=c(min((Mod@likelihood), min(like.nulls)), max((Mod@likelihood), max(like.nulls))), 
     main="Observed vs. randomized")
abline(v=(Mod@likelihood), col="red", lwd=2)
dev.off()

#Estimate the P-value
Mod@likelihood # Obs. value
mean(like.nulls) # E (expected modularity)
z # z
sd(like.nulls)
praw <- sum(like.nulls>(Mod@likelihood)) / length(like.nulls)
ifelse(praw > 0.5, 1-praw, praw)

############### 3a. NESTEDNESS ANALYSIS (RESTRICTED) ############### 

#Calculate the desired nestedness metric (here NODF) for the original network.
obs <- unlist(bipartite::nest.smdm(x = data, 
                                   constraints = Part, #Input the modular structured recovered from step 2
                                   weighted = F, #By considering the edge weights, you are choosing WNODA
                                   decreasing = "fill")) # needs to be fill for binary matrices
#Check the scores
obs

#Calculate constrained interaction probabilities considering the network's modular structure
Pij <- PosteriorProb(M = data, 
                     R.partitions = row.Part, C.partitions = col.Part, #Input the modular structured recovered from step 2
                     Prior.Pij = "degreeprob", #Choose the null model
                     Conditional.level = "modules") #Choose the kind of constraints

#Generate randomized networks with the null model (restricted)
set.seed(1)
nulls <- RestNullModel(M = data, 
                       Pij.Prob = Pij, #Recover the probabilities calculated in the previous command
                       Numbernulls = 100, #This step may take long, so start experimenting with low values
                       Print.null = T, allow.degeneration = F, #Choose whether you allow orphan rows and columns to be removed or not
                       return.nonrm.species = F, 
                       connectance = T, byarea = T, 
                       R.partitions = row.Part, C.partitions = col.Part)

#Calculate the same nestedness metric for all randomized networks
null <- sapply(nulls, function(x) bipartite::nest.smdm(x = x, 
                                                       constraints = Part, 
                                                       weighted = F, 
                                                       decreasing = "fill"))
NODF.null <- unlist(null[3,]) # NODF for entire matrix
NODFsm.null <- unlist(null[8,]) # NODF for nodes belonging to the same module
NODFdm.null <- unlist(null[9,]) # NODF for nodes belonging to different modules

#Plot the observed nestedness value against the distribution of randomized values
png(filename= "figures/figS5_nestedness_CF_asia.png", res= 300, height= 1800, width= 6000)
par(mfrow = c(1,3))
plot(density(NODF.null), xlim=c(min(obs[3], min(NODF.null)), max(obs[3], max(NODF.null))), 
     main="Observed vs. randomized", xlab = "NODF matrix")
abline(v=obs[3], col="red", lwd=2)    
plot(density(NODFsm.null), xlim=c(min(obs[8], min(NODFsm.null)), max(obs[8], max(NODFsm.null))), 
     main="Observed vs. randomized", xlab = "NODFsm matrix")
abline(v=obs[8], col="red", lwd=2)    
plot(density(NODFdm.null), xlim=c(min(obs[9], min(NODFdm.null)), max(obs[9], max(NODFdm.null))), 
     main="Observed vs. randomized", xlab = "NODFdm matrix")
abline(v=obs[9], col="red", lwd=2)    
dev.off()

#Estimate the E, Z-, P-values

#Nestedness in the entire network
mean(NODF.null) # Expected value
obs[3] - mean(NODF.null)/sd(NODF.null) # Z-value
praw.NODF <- sum(NODF.null>obs[3]) / length(NODF.null)
p.NODF <- ifelse(praw.NODF > 0.5, 1- praw.NODF, praw.NODF)    # P-value
p.NODF

#Nestedness within the modules
mean(NODFsm.null) # Expected value
obs[8] - mean(NODFsm.null)/sd(NODFsm.null) # Z-value
praw.NODFsm <- sum(NODFsm.null>obs[8]) / length(NODFsm.null)
p.NODFsm <- ifelse(praw.NODFsm > 0.5, 1- praw.NODFsm, praw.NODFsm)    # P-value
p.NODFsm

#Nestedness between the modules
mean(NODFdm.null) # Expected value
obs[9] - mean(NODFdm.null)/sd(NODFdm.null) # Z-value
praw.NODFdm <- sum(NODFdm.null>obs[9]) / length(NODFdm.null)
p.NODFdm <- ifelse(praw.NODFdm > 0.5, 1- praw.NODFdm, praw.NODFdm)    # P-value
p.NODFdm


############### 3b. NESTEDNESS ANALYSIS (FREE) ############### 

#Calculate the desired nestedness metric (here NODF) for the original network.
obs <- unlist(bipartite::nest.smdm(x = data, 
                                   constraints = Part, #Input the modular structured recovered from step 2
                                   weighted = F, #By considering the edge weights, you are choosing WNODA
                                   decreasing = "fill")) # needs to be fill for binary matrices
#Check the scores
obs

#Generate randomized networks with the null model (free)
set.seed(1)
nulls <- nullmodel(data, N = 100, method = "vaznull")

#Calculate the same nestedness metric for all randomized networks
null <- (sapply(nulls, function(x) bipartite::nest.smdm(x = x,
                                                        constraints = Part,
                                                        weighted = F,
                                                        decreasing = "fill")))
NODF.null <- unlist(null[3,]) # NODF for entire matrix
NODFsm.null <- unlist(null[8,]) # NODF for nodes belonging to the same module
NODFdm.null <- unlist(null[9,]) # NODF for nodes belonging to different modules

#Plot the observed nestedness value against the distribution of randomized values
png(filename= "figures/figS5_nestednessFree_CF_asia.png", res= 300, height= 1800, width= 6000)
par(mfrow = c(1,3))
plot(density(NODF.null), xlim=c(min(obs[3], min(NODF.null)), max(obs[3], max(NODF.null))), 
     main="Observed vs. randomized", xlab = "NODF matrix")
abline(v=obs[3], col="red", lwd=2)    
plot(density(NODFsm.null), xlim=c(min(obs[8], min(NODFsm.null)), max(obs[8], max(NODFsm.null))), 
     main="Observed vs. randomized", xlab = "NODFsm matrix")
abline(v=obs[8], col="red", lwd=2)    
plot(density(NODFdm.null), xlim=c(min(obs[9], min(NODFdm.null)), max(obs[9], max(NODFdm.null))), 
     main="Observed vs. randomized", xlab = "NODFdm matrix")
abline(v=obs[9], col="red", lwd=2)    
dev.off()

#Estimate the E, Z-, P-values

#Nestedness in the entire network
mean(NODF.null) # Expected value
obs[3] - mean(NODF.null)/sd(NODF.null) # Z-value
praw.NODF <- sum(NODF.null>obs[3]) / length(NODF.null)
p.NODF <- ifelse(praw.NODF > 0.5, 1- praw.NODF, praw.NODF)    # P-value
p.NODF

#Nestedness within the modules
mean(NODFsm.null) # Expected value
obs[8] - mean(NODFsm.null)/sd(NODFsm.null) # Z-value
praw.NODFsm <- sum(NODFsm.null>obs[8]) / length(NODFsm.null)
p.NODFsm <- ifelse(praw.NODFsm > 0.5, 1- praw.NODFsm, praw.NODFsm)    # P-value
p.NODFsm

#Nestedness between the modules
mean(NODFdm.null) # Expected value
obs[9] - mean(NODFdm.null)/sd(NODFdm.null) # Z-value
praw.NODFdm <- sum(NODFdm.null>obs[9]) / length(NODFdm.null)
p.NODFdm <- ifelse(praw.NODFdm > 0.5, 1- praw.NODFdm, praw.NODFdm)    # P-value
p.NODFdm


############### 4. PLOTTING THE NETWORK ############### 
par(mfrow = c(1,1))
#Sort the matrix in a way that facilitates visualizing the compound topology
data.comp <- bipartite::sortmatrix(matrix = data, topology = "compound", sort_by = "weights", row_partitions = row.Part, col_partitions = col.Part)

#Assign colors for the modules
modcol <- rainbow((length(unique(Part))), alpha=1)

#Plot the matrix
png(filename= "figures/figSX_compound_CF_asia.png", res= 300, height= 2500, width= 3500)
plotmatrix(data.comp$matrix, 
           row_partitions = data.comp$row_partitions, 
           col_partitions = data.comp$col_partitions, 
           border = T,
           binary = F,
           modules_colors = modcol,
           within_color = modcol, 
           between_color = "lightgrey")
dev.off()

############## ANALYSIS 3.0 CHEMICAL: AGGREGATED ##############
############## 1. PREPARING THE DATA ##############

#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(bipartite)
library(igraph)
source("R/RestNullModel.R")
source("R/PosteriorProb.R")

# Load edge and node lists
nodes = read.delim("data/net1nodes_c.txt", header = T) # node traits
links = read.delim("data/list_chemical", header = T) # list of interactions
# Load network
data <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
data = get.adjacency(data,sparse=FALSE)
data = data[-c(0:48), ] # remove rows
data = data[,-49:-77] # remove columns
data <- empty(data, count = F) # remove nodes with no interactions

#Visualize the raw matrix
visweb(data)

############### 2. MODULARITY ANALYSIS ############### 

#Compute modularity
Mod <- bipartite::computeModules(data)
Mod@likelihood

#Recover the partitions
Part <- bipartite::module2constraints(Mod)
row.Part <- Part[1:nrow(data)]
col.Part <- Part[(nrow(data)+1):(nrow(data)+ncol(data))]

#Test for the significance of modularity with a Monte Carlo procedure

#Generate randomized matrices
set.seed(1)
nulls <- nullmodel(data, N=10, method="r2d")

#Calculate the modularity of the randomized matrices
mod.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(mod.nulls, function(x) x@likelihood)

#Calculate the z-score of the randomized distribution
(z <- (Mod@likelihood - mean(like.nulls))/sd(like.nulls))

#Plot the observed modularity value against the distribution of randomized values
png(filename= "figures/figS5_modularity_chemical_aggregated.png", res= 300, height= 3000, width= 3000)
plot(density(like.nulls), xlim=c(min((Mod@likelihood), min(like.nulls)), max((Mod@likelihood), max(like.nulls))), 
     main="Observed vs. randomized")
abline(v=(Mod@likelihood), col="red", lwd=2)    
dev.off()

#Estimate the P-value
mean(like.nulls) # E (expected modularity)
z # z
sd(like.nulls)
Mod@likelihood
praw <- sum(like.nulls>(Mod@likelihood)) / length(like.nulls)
ifelse(praw > 0.5, 1-praw, praw)

############### 3a. NESTEDNESS ANALYSIS (RESTRICTED) ############### 

#Calculate the desired nestedness metric (here NODF) for the original network.
obs <- unlist(bipartite::nest.smdm(x = data, 
                                   constraints = Part, #Input the modular structured recovered from step 2
                                   weighted = F, #By considering the edge weights, you are choosing WNODA
                                   decreasing = "fill")) # needs to be fill for binary matrices
#Check the scores
obs

#Calculate constrained interaction probabilities considering the network's modular structure
Pij <- PosteriorProb(M = data, 
                     R.partitions = row.Part, C.partitions = col.Part, #Input the modular structured recovered from step 2
                     Prior.Pij = "degreeprob", #Choose the null model
                     Conditional.level = "modules") #Choose the kind of constraints

#Generate randomized networks with the null model (restricted)
set.seed(1)
nulls <- RestNullModel(M = data, 
                       Pij.Prob = Pij, #Recover the probabilities calculated in the previous command
                       Numbernulls = 100, #This step may take long, so start experimenting with low values
                       Print.null = T, allow.degeneration = F, #Choose whether you allow orphan rows and columns to be removed or not
                       return.nonrm.species = F, 
                       connectance = T, byarea = T, 
                       R.partitions = row.Part, C.partitions = col.Part)

#Calculate the same nestedness metric for all randomized networks
null <- sapply(nulls, function(x) bipartite::nest.smdm(x = x, 
                                                       constraints = Part, 
                                                       weighted = F, 
                                                       decreasing = "fill"))
NODF.null <- unlist(null[3,]) # NODF for entire matrix
NODFsm.null <- unlist(null[8,]) # NODF for nodes belonging to the same module
NODFdm.null <- unlist(null[9,]) # NODF for nodes belonging to different modules

#Plot the observed nestedness value against the distribution of randomized values
png(filename= "figures/figS5_nestedness_chemical_aggregated.png", res= 300, height= 1800, width= 6000)
par(mfrow = c(1,3))
plot(density(NODF.null), xlim=c(min(obs[3], min(NODF.null)), max(obs[3], max(NODF.null))), 
     main="Observed vs. randomized", xlab = "NODF matrix")
abline(v=obs[3], col="red", lwd=2)    
plot(density(NODFsm.null), xlim=c(min(obs[8], min(NODFsm.null)), max(obs[8], max(NODFsm.null))), 
     main="Observed vs. randomized", xlab = "NODFsm matrix")
abline(v=obs[8], col="red", lwd=2)    
plot(density(NODFdm.null), xlim=c(min(obs[9], min(NODFdm.null)), max(obs[9], max(NODFdm.null))), 
     main="Observed vs. randomized", xlab = "NODFdm matrix")
abline(v=obs[9], col="red", lwd=2)
dev.off()

#Estimate the E, Z-, P-values

#Nestedness in the entire network
mean(NODF.null) # Expected value
obs[3] - mean(NODF.null)/sd(NODF.null) # Z-value
praw.NODF <- sum(NODF.null>obs[3]) / length(NODF.null)
p.NODF <- ifelse(praw.NODF > 0.5, 1- praw.NODF, praw.NODF)    # P-value
p.NODF

#Nestedness within the modules
mean(NODFsm.null) # Expected value
obs[8] - mean(NODFsm.null)/sd(NODFsm.null) # Z-value
praw.NODFsm <- sum(NODFsm.null>obs[8]) / length(NODFsm.null)
p.NODFsm <- ifelse(praw.NODFsm > 0.5, 1- praw.NODFsm, praw.NODFsm)    # P-value
p.NODFsm

#Nestedness between the modules
mean(NODFdm.null) # Expected value
obs[9] - mean(NODFdm.null)/sd(NODFdm.null) # Z-value
praw.NODFdm <- sum(NODFdm.null>obs[9]) / length(NODFdm.null)
p.NODFdm <- ifelse(praw.NODFdm > 0.5, 1- praw.NODFdm, praw.NODFdm)    # P-value
p.NODFdm

############### 3b. NESTEDNESS ANALYSIS (FREE) ############### 

#Calculate the desired nestedness metric (here NODF) for the original network.
obs <- unlist(bipartite::nest.smdm(x = data, 
                                   constraints = Part, #Input the modular structured recovered from step 2
                                   weighted = F, #By considering the edge weights, you are choosing WNODA
                                   decreasing = "fill")) # needs to be fill for binary matrices
#Check the scores
obs

#Generate randomized networks with the null model (free)
set.seed(1)
nulls <- nullmodel(data, N = 100, method = "vaznull")

#Calculate the same nestedness metric for all randomized networks
null <- (sapply(nulls, function(x) bipartite::nest.smdm(x = x,
                                                        constraints = Part,
                                                        weighted = F,
                                                        decreasing = "fill")))
NODF.null <- unlist(null[3,]) # NODF for entire matrix
NODFsm.null <- unlist(null[8,]) # NODF for nodes belonging to the same module
NODFdm.null <- unlist(null[9,]) # NODF for nodes belonging to different modules

#Plot the observed nestedness value against the distribution of randomized values
png(filename= "figures/figS5_nestednessFree_chemical_aggregated.png", res= 300, height= 1800, width= 6000)
par(mfrow = c(1,3))
plot(density(NODF.null), xlim=c(min(obs[3], min(NODF.null)), max(obs[3], max(NODF.null))), 
     main="Observed vs. randomized", xlab = "NODF matrix")
abline(v=obs[3], col="red", lwd=2)    
plot(density(NODFsm.null), xlim=c(min(obs[8], min(NODFsm.null)), max(obs[8], max(NODFsm.null))), 
     main="Observed vs. randomized", xlab = "NODFsm matrix")
abline(v=obs[8], col="red", lwd=2)    
plot(density(NODFdm.null), xlim=c(min(obs[9], min(NODFdm.null)), max(obs[9], max(NODFdm.null))), 
     main="Observed vs. randomized", xlab = "NODFdm matrix")
abline(v=obs[9], col="red", lwd=2)    
dev.off()

#Estimate the E, Z-, P-values

#Nestedness in the entire network
mean(NODF.null) # Expected value
obs[3] - mean(NODF.null)/sd(NODF.null) # Z-value
praw.NODF <- sum(NODF.null>obs[3]) / length(NODF.null)
p.NODF <- ifelse(praw.NODF > 0.5, 1- praw.NODF, praw.NODF)    # P-value
p.NODF

#Nestedness within the modules
mean(NODFsm.null) # Expected value
obs[8] - mean(NODFsm.null)/sd(NODFsm.null) # Z-value
praw.NODFsm <- sum(NODFsm.null>obs[8]) / length(NODFsm.null)
p.NODFsm <- ifelse(praw.NODFsm > 0.5, 1- praw.NODFsm, praw.NODFsm)    # P-value
p.NODFsm

#Nestedness between the modules
mean(NODFdm.null) # Expected value
obs[9] - mean(NODFdm.null)/sd(NODFdm.null) # Z-value
praw.NODFdm <- sum(NODFdm.null>obs[9]) / length(NODFdm.null)
p.NODFdm <- ifelse(praw.NODFdm > 0.5, 1- praw.NODFdm, praw.NODFdm)    # P-value
p.NODFdm


############### 4. PLOTTING THE NETWORK ############### 
par(mfrow = c(1,1))
#Sort the matrix in a way that facilitates visualizing the compound topology
data.comp <- bipartite::sortmatrix(matrix = data, topology = "compound", sort_by = "weights", row_partitions = row.Part, col_partitions = col.Part)

#Assign colors for the modules
modcol <- rainbow((length(unique(Part))), alpha=1)

#Plot the matrix
png(filename= "figures/figSX_compound_chemical.png", res= 300, height= 2500, width= 3500)
plotmatrix(data.comp$matrix, 
           row_partitions = data.comp$row_partitions, 
           col_partitions = data.comp$col_partitions, 
           border = T,
           binary = F,
           modules_colors = modcol,
           within_color = modcol, 
           between_color = "lightgrey")
dev.off()
############## ANALYSIS 4.0 FIELDWORK: AGGREGATED ##############
############## 1. PREPARING THE DATA ##############

#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(bipartite)
library(igraph)
source("R/RestNullModel.R")
source("R/PosteriorProb.R")

# Load edge and node lists
nodes = read.delim("data/net1nodes_f.txt", header = T) # node traits
links = read.delim("data/list_fieldwork", header = T) # list of interactions
# Load network
data <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
data = get.adjacency(data,sparse=FALSE)
data = data[-c(0:38), ] # remove rows
data = (data[,-39:-55]) # remove columns
data <- empty(data, count = F) # remove nodes with no interactions

#Visualize the raw matrix
visweb(data)

############### 2. MODULARITY ANALYSIS ############### 

#Compute modularity
Mod <- bipartite::computeModules(data)
Mod@likelihood

#Recover the partitions
Part <- bipartite::module2constraints(Mod)
row.Part <- Part[1:nrow(data)]
col.Part <- Part[(nrow(data)+1):(nrow(data)+ncol(data))]

#Test for the significance of modularity with a Monte Carlo procedure

#Generate randomized matrices
nulls <- nullmodel(data, N=10, method="r2d")

#Calculate the modularity of the randomized matrices
mod.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(mod.nulls, function(x) x@likelihood)

#Calculate the z-score of the randomized distribution
(z <- (Mod@likelihood - mean(like.nulls))/sd(like.nulls))

#Plot the observed modularity value against the distribution of randomized values
png(filename= "figures/figS5_modularity_fieldwork_aggregated.png", res= 300, height= 3000, width= 3000)
plot(density(like.nulls), xlim=c(min((Mod@likelihood), min(like.nulls)), max((Mod@likelihood), max(like.nulls))), 
     main="Observed vs. randomized")
abline(v=(Mod@likelihood), col="red", lwd=2)
dev.off()

#Estimate the P-value
mean(like.nulls) # E (expected modularity)
z # z
sd(like.nulls)
Mod@likelihood
praw <- sum(like.nulls>(Mod@likelihood)) / length(like.nulls)
ifelse(praw > 0.5, 1-praw, praw)

############### 3a. NESTEDNESS ANALYSIS (RESTRICTED) ############### 

#Calculate the desired nestedness metric (here NODF) for the original network.
obs <- unlist(bipartite::nest.smdm(x = data, 
                                   constraints = Part, #Input the modular structured recovered from step 2
                                   weighted = F, #By considering the edge weights, you are choosing WNODA
                                   decreasing = "fill")) # needs to be fill for binary matrices
#Check the scores
obs

#Calculate constrained interaction probabilities considering the network's modular structure
Pij <- PosteriorProb(M = data, 
                     R.partitions = row.Part, C.partitions = col.Part, #Input the modular structured recovered from step 2
                     Prior.Pij = "degreeprob", #Choose the null model
                     Conditional.level = "modules") #Choose the kind of constraints

#Generate randomized networks with the null model (restricted)
set.seed(1)
nulls <- RestNullModel(M = data, 
                       Pij.Prob = Pij, #Recover the probabilities calculated in the previous command
                       Numbernulls = 100, #This step may take long, so start experimenting with low values
                       Print.null = T, allow.degeneration = F, #Choose whether you allow orphan rows and columns to be removed or not
                       return.nonrm.species = F, 
                       connectance = T, byarea = T, 
                       R.partitions = row.Part, C.partitions = col.Part)

#Calculate the same nestedness metric for all randomized networks
null <- sapply(nulls, function(x) bipartite::nest.smdm(x = x, 
                                                       constraints = Part, 
                                                       weighted = F, 
                                                       decreasing = "fill"))
NODF.null <- unlist(null[3,]) # NODF for entire matrix
NODFsm.null <- unlist(null[8,]) # NODF for nodes belonging to the same module
NODFdm.null <- unlist(null[9,]) # NODF for nodes belonging to different modules

#Plot the observed nestedness value against the distribution of randomized values
png(filename= "figures/figS5_nestedness_fieldwork_aggregated.png", res= 300, height= 1800, width= 6000)
par(mfrow = c(1,3))
plot(density(NODF.null), xlim=c(min(obs[3], min(NODF.null)), max(obs[3], max(NODF.null))), 
     main="Observed vs. randomized", xlab = "NODF matrix")
abline(v=obs[3], col="red", lwd=2)    
plot(density(NODFsm.null), xlim=c(min(obs[8], min(NODFsm.null)), max(obs[8], max(NODFsm.null))), 
     main="Observed vs. randomized", xlab = "NODFsm matrix")
abline(v=obs[8], col="red", lwd=2)    
plot(density(NODFdm.null), xlim=c(min(obs[9], min(NODFdm.null)), max(obs[9], max(NODFdm.null))), 
     main="Observed vs. randomized", xlab = "NODFdm matrix")
abline(v=obs[9], col="red", lwd=2)
dev.off()

#Estimate the E, Z-, P-values

#Nestedness in the entire network
mean(NODF.null) # Expected value
obs[3] - mean(NODF.null)/sd(NODF.null) # Z-value
praw.NODF <- sum(NODF.null>obs[3]) / length(NODF.null)
p.NODF <- ifelse(praw.NODF > 0.5, 1- praw.NODF, praw.NODF)    # P-value
p.NODF

#Nestedness within the modules
mean(NODFsm.null) # Expected value
obs[8] - mean(NODFsm.null)/sd(NODFsm.null) # Z-value
praw.NODFsm <- sum(NODFsm.null>obs[8]) / length(NODFsm.null)
p.NODFsm <- ifelse(praw.NODFsm > 0.5, 1- praw.NODFsm, praw.NODFsm)    # P-value
p.NODFsm

#Nestedness between the modules
mean(NODFdm.null) # Expected value
obs[9] - mean(NODFdm.null)/sd(NODFdm.null) # Z-value
praw.NODFdm <- sum(NODFdm.null>obs[9]) / length(NODFdm.null)
p.NODFdm <- ifelse(praw.NODFdm > 0.5, 1- praw.NODFdm, praw.NODFdm)    # P-value
p.NODFdm


############### 3b. NESTEDNESS ANALYSIS (FREE) ############### 

#Calculate the desired nestedness metric (here NODF) for the original network.
obs <- unlist(bipartite::nest.smdm(x = data, 
                                   constraints = Part, #Input the modular structured recovered from step 2
                                   weighted = F, #By considering the edge weights, you are choosing WNODA
                                   decreasing = "fill")) # needs to be fill for binary matrices
#Check the scores
obs

#Generate randomized networks with the null model (free)
set.seed(1)
nulls <- nullmodel(data, N = 100, method = "vaznull")

#Calculate the same nestedness metric for all randomized networks
null <- (sapply(nulls, function(x) bipartite::nest.smdm(x = x,
                                                        constraints = Part,
                                                        weighted = F,
                                                        decreasing = "fill")))
NODF.null <- unlist(null[3,]) # NODF for entire matrix
NODFsm.null <- unlist(null[8,]) # NODF for nodes belonging to the same module
NODFdm.null <- unlist(null[9,]) # NODF for nodes belonging to different modules

#Plot the observed nestedness value against the distribution of randomized values
png(filename= "figures/figS5_nestednessFree_fieldwork_aggregated.png", res= 300, height= 1800, width= 6000)
par(mfrow = c(1,3))
plot(density(NODF.null), xlim=c(min(obs[3], min(NODF.null)), max(obs[3], max(NODF.null))), 
     main="Observed vs. randomized", xlab = "NODF matrix")
abline(v=obs[3], col="red", lwd=2)    
plot(density(NODFsm.null), xlim=c(min(obs[8], min(NODFsm.null)), max(obs[8], max(NODFsm.null))), 
     main="Observed vs. randomized", xlab = "NODFsm matrix")
abline(v=obs[8], col="red", lwd=2)    
plot(density(NODFdm.null), xlim=c(min(obs[9], min(NODFdm.null)), max(obs[9], max(NODFdm.null))), 
     main="Observed vs. randomized", xlab = "NODFdm matrix")
abline(v=obs[9], col="red", lwd=2)    
dev.off()

#Estimate the E, Z-, P-values

#Nestedness in the entire network
mean(NODF.null) # Expected value
obs[3] - mean(NODF.null)/sd(NODF.null) # Z-value
praw.NODF <- sum(NODF.null>obs[3]) / length(NODF.null)
p.NODF <- ifelse(praw.NODF > 0.5, 1- praw.NODF, praw.NODF)    # P-value
p.NODF

#Nestedness within the modules
mean(NODFsm.null) # Expected value
obs[8] - mean(NODFsm.null)/sd(NODFsm.null) # Z-value
praw.NODFsm <- sum(NODFsm.null>obs[8]) / length(NODFsm.null)
p.NODFsm <- ifelse(praw.NODFsm > 0.5, 1- praw.NODFsm, praw.NODFsm)    # P-value
p.NODFsm

#Nestedness between the modules
mean(NODFdm.null) # Expected value
obs[9] - mean(NODFdm.null)/sd(NODFdm.null) # Z-value
praw.NODFdm <- sum(NODFdm.null>obs[9]) / length(NODFdm.null)
p.NODFdm <- ifelse(praw.NODFdm > 0.5, 1- praw.NODFdm, praw.NODFdm)    # P-value
p.NODFdm


############### 4. PLOTTING THE NETWORK ############### 
par(mfrow = c(1,1))
#Sort the matrix in a way that facilitates visualizing the compound topology
data.comp <- bipartite::sortmatrix(matrix = data, topology = "compound", sort_by = "weights", row_partitions = row.Part, col_partitions = col.Part)

#Assign colors for the modules
modcol <- rainbow((length(unique(Part))), alpha=1)

#Plot the matrix
png(filename= "figures/figSX_compound_fieldwork.png", res= 300, height= 2500, width= 3500)
plotmatrix(data.comp$matrix, 
           row_partitions = data.comp$row_partitions, 
           col_partitions = data.comp$col_partitions, 
           border = T,
           binary = F,
           modules_colors = modcol,
           within_color = modcol, 
           between_color = "lightgrey")
dev.off()
############## ANALYSIS 5.0 PALYNOLOGICAL: AGGREGATED ##############
############## 1. PREPARING THE DATA ##############

#Delete all previous objects
rm(list= ls())
#Clear the console
cat("\014") 

#Load the required packages and functions.
library(bipartite)
library(igraph)
source("R/RestNullModel.R")
source("R/PosteriorProb.R")

# Load edge and node lists
nodes = read.delim("data/net1nodes_p.txt", header = T) # node traits
links = read.delim("data/list_palynological", header = T) # list of interactions
# Load network
data <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Create an adjacent matrix for bipartite
data = get.adjacency(data,sparse=FALSE)
data = data[-c(0:14), ] # remove rows
data = data[,-15:-100] # remove columns
data <- empty(data, count = F) # remove nodes with no interactions

#Visualize the raw matrix
visweb(data)

############### 2. MODULARITY ANALYSIS ############### 

#Compute modularity
Mod <- bipartite::computeModules(data)
Mod@likelihood

#Recover the partitions
Part <- bipartite::module2constraints(Mod)
row.Part <- Part[1:nrow(data)]
col.Part <- Part[(nrow(data)+1):(nrow(data)+ncol(data))]

#Test for the significance of modularity with a Monte Carlo procedure

#Generate randomized matrices
set.seed(1)
nulls <- nullmodel(data, N=10, method="r2d")

#Calculate the modularity of the randomized matrices
mod.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(mod.nulls, function(x) x@likelihood)

#Calculate the z-score of the randomized distribution
(z <- (Mod@likelihood - mean(like.nulls))/sd(like.nulls))

#Plot the observed modularity value against the distribution of randomized values
png(filename= "figures/figS5_modularity_palynological_aggregated.png", res= 300, height= 3000, width= 3000)
plot(density(like.nulls), xlim=c(min((Mod@likelihood), min(like.nulls)), max((Mod@likelihood), max(like.nulls))), 
     main="Observed vs. randomized")
abline(v=(Mod@likelihood), col="red", lwd=2)
dev.off()

#Estimate the P-value
Mod@likelihood # Obs. value
mean(like.nulls) # E (expected modularity)
z # z
sd(like.nulls)
praw <- sum(like.nulls>(Mod@likelihood)) / length(like.nulls)
ifelse(praw > 0.5, 1-praw, praw)

############### 3a. NESTEDNESS ANALYSIS (RESTRICTED) ############### 

#Calculate the desired nestedness metric (here NODF) for the original network.
obs <- unlist(bipartite::nest.smdm(x = data, 
                                   constraints = Part, #Input the modular structured recovered from step 2
                                   weighted = F, #By considering the edge weights, you are choosing WNODA
                                   decreasing = "fill")) # needs to be fill for binary matrices
#Check the scores
obs

#Calculate constrained interaction probabilities considering the network's modular structure
Pij <- PosteriorProb(M = data, 
                     R.partitions = row.Part, C.partitions = col.Part, #Input the modular structured recovered from step 2
                     Prior.Pij = "degreeprob", #Choose the null model
                     Conditional.level = "modules") #Choose the kind of constraints

#Generate randomized networks with the null model (restricted)
set.seed(1)
nulls <- RestNullModel(M = data, 
                       Pij.Prob = Pij, #Recover the probabilities calculated in the previous command
                       Numbernulls = 100, #This step may take long, so start experimenting with low values
                       Print.null = T, allow.degeneration = F, #Choose whether you allow orphan rows and columns to be removed or not
                       return.nonrm.species = F, 
                       connectance = T, byarea = T, 
                       R.partitions = row.Part, C.partitions = col.Part)

#Calculate the same nestedness metric for all randomized networks
null <- sapply(nulls, function(x) bipartite::nest.smdm(x = x, 
                                                       constraints = Part, 
                                                       weighted = F, 
                                                       decreasing = "fill"))
NODF.null <- unlist(null[3,]) # NODF for entire matrix
NODFsm.null <- unlist(null[8,]) # NODF for nodes belonging to the same module
NODFdm.null <- unlist(null[9,]) # NODF for nodes belonging to different modules

#Plot the observed nestedness value against the distribution of randomized values
png(filename= "figures/figS5_nestedness_palynological_aggregated.png", res= 300, height= 1800, width= 6000)
par(mfrow = c(1,3))
plot(density(NODF.null), xlim=c(min(obs[3], min(NODF.null)), max(obs[3], max(NODF.null))), 
     main="Observed vs. randomized", xlab = "NODF matrix")
abline(v=obs[3], col="red", lwd=2)    
plot(density(NODFsm.null), xlim=c(min(obs[8], min(NODFsm.null)), max(obs[8], max(NODFsm.null))), 
     main="Observed vs. randomized", xlab = "NODFsm matrix")
abline(v=obs[8], col="red", lwd=2)    
plot(density(NODFdm.null), xlim=c(min(obs[9], min(NODFdm.null)), max(obs[9], max(NODFdm.null))), 
     main="Observed vs. randomized", xlab = "NODFdm matrix")
abline(v=obs[9], col="red", lwd=2)    
dev.off()

#Estimate the E, Z-, P-values

#Nestedness in the entire network
mean(NODF.null) # Expected value
obs[3] - mean(NODF.null)/sd(NODF.null) # Z-value
praw.NODF <- sum(NODF.null>obs[3]) / length(NODF.null)
p.NODF <- ifelse(praw.NODF > 0.5, 1- praw.NODF, praw.NODF)    # P-value
p.NODF

#Nestedness within the modules
mean(NODFsm.null) # Expected value
obs[8] - mean(NODFsm.null)/sd(NODFsm.null) # Z-value
praw.NODFsm <- sum(NODFsm.null>obs[8]) / length(NODFsm.null)
p.NODFsm <- ifelse(praw.NODFsm > 0.5, 1- praw.NODFsm, praw.NODFsm)    # P-value
p.NODFsm

#Nestedness between the modules
mean(NODFdm.null) # Expected value
obs[9] - mean(NODFdm.null)/sd(NODFdm.null) # Z-value
praw.NODFdm <- sum(NODFdm.null>obs[9]) / length(NODFdm.null)
p.NODFdm <- ifelse(praw.NODFdm > 0.5, 1- praw.NODFdm, praw.NODFdm)    # P-value
p.NODFdm


############### 3b. NESTEDNESS ANALYSIS (FREE) ############### 

#Calculate the desired nestedness metric (here NODF) for the original network.
obs <- unlist(bipartite::nest.smdm(x = data, 
                                   constraints = Part, #Input the modular structured recovered from step 2
                                   weighted = F, #By considering the edge weights, you are choosing WNODA
                                   decreasing = "fill")) # needs to be fill for binary matrices
#Check the scores
obs

#Generate randomized networks with the null model (free)
set.seed(1)
nulls <- nullmodel(data, N = 100, method = "vaznull")

#Calculate the same nestedness metric for all randomized networks
null <- (sapply(nulls, function(x) bipartite::nest.smdm(x = x,
                                                        constraints = Part,
                                                        weighted = F,
                                                        decreasing = "fill")))
NODF.null <- unlist(null[3,]) # NODF for entire matrix
NODFsm.null <- unlist(null[8,]) # NODF for nodes belonging to the same module
NODFdm.null <- unlist(null[9,]) # NODF for nodes belonging to different modules

#Plot the observed nestedness value against the distribution of randomized values
png(filename= "figures/figS5_nestednessFree_palynological_aggregated.png", res= 300, height= 1800, width= 6000)
par(mfrow = c(1,3))
plot(density(NODF.null), xlim=c(min(obs[3], min(NODF.null)), max(obs[3], max(NODF.null))), 
     main="Observed vs. randomized", xlab = "NODF matrix")
abline(v=obs[3], col="red", lwd=2)    
plot(density(NODFsm.null), xlim=c(min(obs[8], min(NODFsm.null)), max(obs[8], max(NODFsm.null))), 
     main="Observed vs. randomized", xlab = "NODFsm matrix")
abline(v=obs[8], col="red", lwd=2)    
plot(density(NODFdm.null), xlim=c(min(obs[9], min(NODFdm.null)), max(obs[9], max(NODFdm.null))), 
     main="Observed vs. randomized", xlab = "NODFdm matrix")
abline(v=obs[9], col="red", lwd=2)    
dev.off()

#Estimate the E, Z-, P-values

#Nestedness in the entire network
mean(NODF.null) # Expected value
obs[3] - mean(NODF.null)/sd(NODF.null) # Z-value
praw.NODF <- sum(NODF.null>obs[3]) / length(NODF.null)
p.NODF <- ifelse(praw.NODF > 0.5, 1- praw.NODF, praw.NODF)    # P-value
p.NODF

#Nestedness within the modules
mean(NODFsm.null) # Expected value
obs[8] - mean(NODFsm.null)/sd(NODFsm.null) # Z-value
praw.NODFsm <- sum(NODFsm.null>obs[8]) / length(NODFsm.null)
p.NODFsm <- ifelse(praw.NODFsm > 0.5, 1- praw.NODFsm, praw.NODFsm)    # P-value
p.NODFsm

#Nestedness between the modules
mean(NODFdm.null) # Expected value
obs[9] - mean(NODFdm.null)/sd(NODFdm.null) # Z-value
praw.NODFdm <- sum(NODFdm.null>obs[9]) / length(NODFdm.null)
p.NODFdm <- ifelse(praw.NODFdm > 0.5, 1- praw.NODFdm, praw.NODFdm)    # P-value
p.NODFdm


############### 4. PLOTTING THE NETWORK ############### 
par(mfrow = c(1,1))
#Sort the matrix in a way that facilitates visualizing the compound topology
data.comp <- bipartite::sortmatrix(matrix = data, topology = "compound", sort_by = "weights", row_partitions = row.Part, col_partitions = col.Part)

#Assign colors for the modules
modcol <- rainbow((length(unique(Part))), alpha=1)

#Plot the matrix
png(filename= "figures/figSX_compound_chemical.png", res= 300, height= 2500, width= 3500)
plotmatrix(data.comp$matrix, 
           row_partitions = data.comp$row_partitions, 
           col_partitions = data.comp$col_partitions, 
           border = T,
           binary = F,
           modules_colors = modcol,
           within_color = modcol, 
           between_color = "lightgrey")
dev.off()
