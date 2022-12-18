###########
# 1 SETUP #
###########

# 1.1 Setting directories

# Set directory where packages should be installed
.libPaths('/home/danimelsz/Downloads/R/x86_64-pc-linux-gnu-library/4.1')

# Set directory
setwd("~/Desktop/Stingless_Bees/GitHub/data")

# 1.2 Load Packages

# The packages 'car' and 'tidyverse' require the following libraries in Ubuntu: 
# $ sudo apt-get install libcurl4-openssl-dev
# $ sudo apt-get install -y libssl-dev
# $ sudo apt-get install libxml2-dev

library(AICcmodavg)
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
library(nlme)
library(phytools)
library(picante)
library(reshape2)
library(RInSp)
library(rstatix)
library(tidyverse)

################
# 2 INPUT DATA #
################

# 2.1 GLOBAL DATA
resin_family.mat = read.delim("matrix_family_full.txt", row.names=1)%>% as.matrix()
resin_family2.mat <- empty(resin_family.mat, count = F)
all = specieslevel(resin_family2.mat, index="ALL")

# 2.2 CHEMICAL DATA
chem.mat = read.delim("matrix_family_chem.txt", row.names=1)%>% as.matrix()
chem2.mat <- empty(chem.mat, count = F)
all.chem = specieslevel(chem2.mat, index="ALL")

# 2.3 FIELDWORK DATA
field.mat = read.delim("matrix_family_fieldwork.txt", row.names=1) %>% as.matrix()
field2.mat = empty(field.mat, count=F)
all.field = specieslevel(field2.mat, index="ALL")

# 2.4 PALYNOLOGICAL DATA
pal.mat = read.delim("matrix_family_pal.txt", row.names=1)%>% as.matrix()
pal2.mat = empty(pal.mat, count=F)
all.pal = specieslevel(pal2.mat, index="ALL")

# 2.5 TRAITS
traits = data.frame("bee_species" = c('Axestotrigona_ferruginea', 'Frieseomelitta_longipes', 'Frieseomelitta_silvestrii', 'Frieseomelitta_varia', 'Geniotrigona_thoracica', 'Heterotrigona_erythrogastra', 'Heterotrigona_hobbyi', 'Heterotrigona_itama', 'Homotrigona_apicalis', 'Homotrigona_binghami', 'Homotrigona_canifrons','Homotrigona_fimbriata', 'Homotrigona_haematoptera', 'Homotrigona_melanoleuca', 'Lepidotrigona_terminata', 'Lepidotrigona_ventralis', 'Lestrimelitta_limao', 'Lisotrigona_cacciae', 'Lisotrigona_furva', 'Melipona_beecheii', 'Melipona_compressipes','Melipona_fasciculata','Melipona_favosa','Melipona_flavolineata','Melipona_fuliginosa','Melipona_mandacaia','Melipona_marginata', 'Melipona_orbignyi', 'Melipona_quadrifasciata','Melipona_scutellaris','Melipona_seminigra_merrillae', 'Melipona_subnitida', 'Nannotrigona_testaceicornis', 'Paratrigona_anduzei', 'Plebeia_droryana', 'Plebeia_emerina', 'Plebeia_frontalis', 'Plebeia_lucii', 'Plebeia_remota', 'Scaptotrigona_bipunctata', 'Scaptotrigona_depilis', 'Scaptotrigona_jujuyensis', 'Scaptotrigona_postica', 'Tetragona_clavipes', 'Tetragonisca_angustula', 'Tetragonisca_fiebrigi', 'Tetragonula_atripes', 'Tetragonula_biroi','Tetragonula_carbonaria','Tetragonula_collina','Tetragonula_fuscobalteata','Tetragonula_geissleri','Tetragonula_hockingsi','Tetragonula_laeviceps','Tetragonula_melanocephala','Tetragonula_melina','Tetragonula_minor','Tetragonula_pagdeni','Tetragonula_reepeni', 'Tetragonula_rufibasalis', 'Tetragonula_sapiens', 'Trigona_spinipes', 'Trigona_corvina', 'Trigona_fulviventris', 'Trigona_recursa', 'Trigona_williana', 'Wallacetrigona_incisa'), 
                    'body_size_mean' = c(1.925, NA, 1.01, 1.74, 1.835, 1.15, NA, 1.645, 1.46, 1.495, 1.805, 1.925, NA, 1.44, 1.97, 1.145, 2.181, 0.93, 1.027, 2.6809, 3.725, 3.3, 2.633, 2.671, 3.81, 2.7083, 2.299, 3.5, 3.2105, 3.912, 3.02, 3.05, 1.624, 2.00, 1.476, 1.24075, 1.02, 1.00, 1.865, 2.32, 1.71, NA, 2.257, 1.886, 1.257, 1.145, NA, 1.21, 1.109, 1.5, 0.91, 1.915, NA, 1.2073, NA, NA, NA, 0.981, NA, NA, 3.76, 2.224, 2.125, 1.36, 1.724, 2.625, NA), 
                    'biogeographical_region' = c('African','Neotropical','Neotropical','Neotropical','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Neotropical','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Indo-Malayan-Australasian','Neotropical','Neotropical','Neotropical','Neotropical','Neotropical','Indo-Malayan-Australasian'), 
                    'phylogeny' = c('no','no','yes','yes','yes','yes','yes','yes','yes','yes','yes','yes','yes','yes','yes','yes','yes','no','yes','yes','yes','yes','yes','yes','yes','yes','yes','no','yes','yes','yes','no','yes','no','yes','no','yes','no','no','yes','yes','no','yes','yes','yes','yes','yes','no','yes','yes','yes','yes','yes','yes','yes','yes','yes','yes','no','yes','yes','no','yes','no','yes','yes','yes'),
                    'number_of_articles' = c(1,1,1,2,3,2,1,5,4,2,2,3,2,2,NA,2,1,2,2,1,3,4,2,1,1,1,1,1,8,3,1,3,1,1,1,1,1,1,2,1,2,1,4,2,7,1,1,1,5,3,2,2,2,6,2,1,1,1,1,1,3,2,1,1, 1, 1, 1), 
                    'number_of_interactions' = c(1,1,1,30,7,3,1,9,8,7,12,9,4,5,NA,4,26,5,4,2,4,114,3,74,1,24,1,1,47,4,1,84,23,1,2,1,1,1,3,1,2,9,88,2,97,9,3,1,25,18,4,13,2,19,8,1,1,1,3,3,3,3, 1, 1, 22, 1, 2), 
                    'botanical_origin_method' = c('chem','chem','chem','chem_pal','chem_field','chem_field','chem','chem_field','chem_field','chem_field','chem_field','chem_field','chem_field','chem_field', NA, 'chem_field','pal','chem','chem','chem','chem_field','chem_pal','chem_field','pal','field','pal','chem','chem','chem_field_pal','chem_field','field','chem_pal','pal','chem','field','field','field','field','chem_field','chem_field','chem','chem','chem_pal','field','chem_field_pal','chem','field','chem','chem_field','chem_field','chem_field','chem_field','field','chem_field','chem_field','field','chem','chem','field','chem','chem_field','field','field','field','pal','chem','chem'),
'module' = c('M6', NA, NA, 'M6', 'M2', 'M2', NA, 'M2', 'M2', 'M2', 'M2', 'M2', NA, 'M1', NA, 'M2', 'M7', 'M4', 'M2', 'M6', 'M1','M4','M1','M4','M4','M8','M1', NA, 'M7','M1','M4', 'M4', 'M7', 'M1', 'M2', 'M3', 'M3', 'M1', 'M1', 'M2', 'M1', NA, 'M5', 'M2', 'M7', 'M6', NA, 'M3','M1','M2','M2','M2', NA, 'M2',NA,NA,NA,'M1',NA, NA, 'M1', 'M1', 'M3', 'M3', 'M6', NA, NA),
'five.modules' = c('M3', 'M1', 'M3', 'M3', 'M1', 'M1', 'M1', 'M1', 'M1', 'M1', 'M1','M1', 'M1', 'M1', NA, 'M1', 'M3', 'M1', 'M1', 'M4', 'M4','M2','M4','M2','M1','M2','M1', 'M4', 'M5','M1','M2', 'M2', 'M3', 'M1', 'M1', 'M1', 'M1', 'M1', 'M1', 'M5', 'M1', 'M3', 'M2', 'M1', 'M3', 'M3', 'M1', 'M1','M1','M1','M1','M1','M1','M1','M1','M1','M1','M1','M1', 'M1', 'M1', 'M1', 'M1', 'M1', 'M3', 'M3', 'M1')) 

# 2.6 PHYLOGENY
tree = read.nexus(file = "tree_edited.nex")
plotTree(tree,node.numbers=T,fsize=0.35)

########################
# 3 CENTRALITY METRICS #
########################

all = specieslevel(resin_family2.mat, index="ALL")
all.f = as.data.frame(all$`higher level`)
all.f = cbind(rownames(all.f), all.f)
rownames(all.f) = NULL
colnames(all.f) = c('bee_species', 'degree', 'norm.degree', 'strength', 'push.pull', 'nestedrank', 'PDI', 'resource_range', 'specificity_index', 'PSI', 'NSI', 'betweenness', 'weighted_betweenness', 'closeness', 'weighted_closeness', 'fisher_alpha', 'partner_diversity', 'effective_partners', 'proportional_generality', 'proportional_similarity', 'specialization_d')
all.data.f = merge(all.f, traits, all=T)
all.data.f = all.data.f %>% drop_na(body_size_mean) # remove NA values in ITD
all.data.f = all.data.f %>% drop_na(specialization_d) # remove NA values in centrality-metrics

#################
# 4 PGLS MODELS #
#################

# 4.1 Prepare data
itd.data = all.data.f[!(all.data.f$phylogeny=='no'),] # remove species absent in phylogeny
itd.data = itd.data[c(1,2, 12, 14, 21, 22)] # select columns of species and ITD
log_ITD_vector = log(itd.data$body_size_mean)# calculate a vector of log ITD
itd.data$log_ITD = log_ITD_vector # add the vector as a new column in the dataframe
itd.data = data.frame(itd.data[,-1], row.names=itd.data[,1]) # convert the column "bee species" as row names
# Check if ITD data is ok
itd.data = itd.data[!(row.names(itd.data)) %in% c(""),] # remove species in the dataframe but not in the tree
ITD_setnames = setNames(itd.data$body_size_mean,rownames(itd.data)) # create a vector with log ITD associated with species names
obj.itd = name.check(tree,ITD_setnames) # check if terminals in tree and row names in dataframe are ok
tree.cut.itd = drop.tip(tree,obj.itd$tree_not_data) # prune species in tree absent in dataframe
name.check(tree.cut.itd,ITD_setnames) # check if everything is ok

# 4.1 Fitting evolutionary models (TABLE S3)

# 4.1.1 ITD
fitBM.ITD<-fitContinuous(tree,ITD_setnames)
fitOU.ITD<-fitContinuous(tree,ITD_setnames,model="OU")
fitEB.ITD<-fitContinuous(tree,ITD_setnames,model="EB")

aic.vals.ITD<-setNames(c(fitBM.ITD$opt$aicc,fitOU.ITD$opt$aicc,fitEB.ITD$opt$aicc),
    c("BM","OU","EB"))
aic.vals.ITD
aic.w(aic.vals.ITD) # OU is better for ITD

# 4.1.2 SPECIALIZATION D'
log_d_vector = log(itd.data$specialization_d) # calculate a vector of log d'
itd.data$log_d = log_d_vector # add the vector as a new column in the dataframe
itd.data.d = itd.data[-c(20,39, 42),] # remove rows with -Inf values in log d'
d_setnames = setNames(itd.data.d$log_d,rownames(itd.data.d)) # create a vector with log d' associated with species names
obj.d = name.check(tree,d_setnames) # check if terminals in tree and row names in dataframe are ok
tree.cut.d = drop.tip(tree,obj.d$tree_not_data) # prune the tree
name.check(tree.cut.d,d_setnames) # check if everything is ok

fitBM.d<-fitContinuous(tree,d_setnames)
fitOU.d<-fitContinuous(tree,d_setnames,model="OU")
fitEB.d<-fitContinuous(tree,d_setnames,model="EB")

aic.vals.d<-setNames(c(fitBM.d$opt$aicc,fitOU.d$opt$aicc,fitEB.d$opt$aicc),
    c("BM","OU","EB"))
aic.vals.d
aic.w(aic.vals.d) # OU is better for specialization d

# 4.1.2 RELATIVE DEGREE
log_r_vector = log(itd.data$degree) # calculate a vector of log relative degree
itd.data$log_r = log_r_vector # add the vector as a new column in the dataframe
itd.data.r = itd.data[-c(6,8, 13, 16, 18, 20, 26, 28, 35, 38, 40, 42),]
r_setnames = setNames(itd.data.r$log_r,rownames(itd.data.r)) # create a vector with log r associated with species names
obj.r = name.check(tree,r_setnames) # check if terminals in tree and row names in dataframe are ok
tree.cut.r = drop.tip(tree,obj.r$tree_not_data) # prune the tree
name.check(tree.cut.r,r_setnames) # check if everything is ok

fitBM.r<-fitContinuous(tree,r_setnames)
fitOU.r<-fitContinuous(tree,r_setnames,model="OU")
fitEB.r<-fitContinuous(tree,r_setnames,model="EB")

aic.vals.r<-setNames(c(fitBM.r$opt$aicc,fitOU.r$opt$aicc,fitEB.r$opt$aicc),
    c("BM","OU","EB"))
aic.vals.r
aic.w(aic.vals.r) # OU is better for degree

# 4.1.3 BETWEENNESS
log_b_vector = log(itd.data$betweenness) # calculate a vector of log relative degree
itd.data$log_b = log_b_vector # add the vector as a new column in the dataframe
itd.data.b = itd.data[-c(6, 8, 13, 18, 20, 26, 28, 35, 38, 40, 42),] # remove rows with -Inf values in log b
b_setnames = setNames(itd.data.b$log_b,rownames(itd.data.b)) # create a vector with log r associated with species names
obj.b = name.check(tree,b_setnames) # check if terminals in tree and row names in dataframe are ok
tree.cut.b = drop.tip(tree,obj.b$tree_not_data) # prune the tree
name.check(tree.cut.b,b_setnames) # check if everything is ok

fitBM.b<-fitContinuous(tree,b_setnames)
fitOU.b<-fitContinuous(tree,b_setnames,model="OU")
fitEB.b<-fitContinuous(tree,b_setnames,model="EB")

aic.vals.b<-setNames(c(fitBM.b$opt$aicc,fitOU.b$opt$aicc,fitEB.b$opt$aicc),
    c("BM","OU","EB"))
aic.vals.b
aic.w(aic.vals.b) # OU is better for betweenness

# 4.1.4 CLOSENESS
log_c_vector = log(itd.data$closeness) # calculate a vector of log relative degree
itd.data$log_c = log_c_vector # add the vector as a new column in the dataframe
itd.data.c = itd.data[-c(42),] # remove rows with -Inf values in log b
c_setnames = setNames(itd.data$log_c,rownames(itd.data.c)) # create a vector with log r associated with species names
obj.c = name.check(tree,c_setnames) # check if terminals in tree and row names in dataframe are ok
tree.cut.c = drop.tip(tree,obj.c$tree_not_data) # prune the tree
name.check(tree.cut.c,c_setnames) # check if everything is ok

fitBM.c<-fitContinuous(tree,c_setnames)
fitOU.c<-fitContinuous(tree,c_setnames,model="OU")
fitEB.c<-fitContinuous(tree,c_setnames,model="EB")

aic.vals.c<-setNames(c(fitBM.c$opt$aicc,fitOU.c$opt$aicc,fitEB.c$opt$aicc),
    c("BM","OU","EB"))
aic.vals.c
aic.w(aic.vals.c) # OU is better for closeness


# 4.2 PGLS (TABLE 1)

# 4.2.1 SPECIALIZATION D' ~

# Models
pglsOU.d1 = gls(log_d~1, correlation = corMartins(1, phy=tree.cut.d), data = itd.data.d, method="ML")
pglsOU.d2 = gls(log_d~log_ITD, correlation = corMartins(1, phy=tree.cut.d), data = itd.data.d, method="ML")

# AIC 
pgls_d = list(pglsOU.d1, pglsOU.d2)
aictab(cand.set = pgls_d)
summary(pglsOU.d2)

# Diagnostic
plot(pglsOU.d2)

# 4.2.2 Relative DEGREE ~
pglsOU.r1 = gls(log_r~1, correlation = corMartins(1, phy=tree.cut.r), data = itd.data.r, method="ML")
pglsOU.r2 = gls(log_r~log_ITD, correlation = corMartins(1, phy=tree.cut.r), data = itd.data.r, method="ML")

# AIC 
pgls_r = list(pglsOU.r1, pglsOU.r2)
aictab(cand.set = pgls_r)
summary(pglsOU.r2)

# Diagnostic
plot(pglsOU.r2)

# 4.2.3 BETWEENNESS ~
pglsOU.b1 = gls(log_b~1, correlation = corMartins(1, phy=tree.cut.b), data = itd.data.b, method="ML")
pglsOU.b2 = gls(log_b~log_ITD, correlation = corMartins(1, phy=tree.cut.b), data = itd.data.b, method="ML")

# AIC 
pgls_b = list(pglsOU.b1, pglsOU.b2)
aictab(cand.set = pgls_b)
summary(pglsOU.b2)

# Diagnostic
plot(pglsOU.b2)

# 4.2.4 CLOSENESS
pglsOU.c1 = gls(log_c~1, correlation = corMartins(1, phy=tree.cut.c), data = itd.data.c, method="ML")
pglsOU.c2 = gls(log_c~log_ITD, correlation = corMartins(1, phy=tree.cut.c), data = itd.data.c, method="ML")

# AIC 
pgls_c = list(pglsOU.c1, pglsOU.c2)
aictab(cand.set = pgls_c)
summary(pglsOU.bc)

# Diagnostic
plot(pglsOU.c1)

###################################################
# TABLE S4: LMs INCLUDING ITD AND SAMPLING EFFORT #
###################################################

# Node-level metrics
all = specieslevel(resin_family2.mat, index="ALL")
all.f = as.data.frame(all$`higher level`)
all.f = cbind(rownames(all.f), all.f)
rownames(all.f) = NULL
colnames(all.f) = c('bee_species', 'degree', 'norm.degree', 'strength', 'push.pull', 'nestedrank', 'PDI', 'resource_range', 'specificity_index', 'PSI', 'NSI', 'betweenness', 'weighted_betweenness', 'closeness', 'weighted_closeness', 'fisher_alpha', 'partner_diversity', 'effective_partners', 'proportional_generality', 'proportional_similarity', 'specialization_d')
all.data.f = merge(all.f, traits, all=T)
all.data.f = all.data.f %>% drop_na(body_size_mean) # remove NA values in ITD
all.data.f = all.data.f %>% drop_na(specialization_d) # remove NA values in centrality-metrics

# Log scale
d2 = log(all.data.f$specialization_d)
d2[is.na(d2) | d2=="-Inf"] = NA
r2 = log(all.data.f$degree)
c2 = log(all.data.f$closeness)
c2[is.na(c2) | c2=="-Inf"] = NA
b2 = log(all.data.f$betweeness)
b2[is.na(b2) | b2=="-Inf"] = NA

# LMs: Specialization d' ~ 
mod1 = lm(d2~body_size_mean*number_of_articles, data=all.data.f)
mod2 = lm(d2~body_size_mean+number_of_articles, data=all.data.f)
mod3 = lm(d2~body_size_mean, data=all.data.f)
mod4 = lm(d2~number_of_articles, data=all.data.f)
mod5 = lm(d2~1, data=all.data.f)

models_d = list(mod1, mod2, mod3, mod4, mod5)
mod.names <- c('itd*sampEff','itd+sampEff','ITD', 'sampEff', '1')
aictab(cand.set = models_d, modnames = mod.names)

summary(mod3)

# LMs: relative degree ~
modr1 = lm(r2~body_size_mean*number_of_articles, data=all.data.f)
modr2 = lm(r2~body_size_mean+number_of_articles, data=all.data.f)
modr3 = lm(r2~body_size_mean, data=all.data.f)
modr4 = lm(r2~number_of_articles, data=all.data.f)
modr5 = lm(r2~1, data=all.data.f)

models_r = list(modr1, modr2, modr3, modr4, modr5)
mod.names <- c('itd*sampEff','itd+sampEff','ITD', 'sampEff', '1')
aictab(cand.set = models_r, modnames = mod.names)

# LMs: closeness ~
modc1 = lm(r2~body_size_mean*number_of_articles, data=all.data.f)
modc2 = lm(r2~body_size_mean+number_of_articles, data=all.data.f)
modc3 = lm(r2~body_size_mean, data=all.data.f)
modc4 = lm(r2~number_of_articles, data=all.data.f)
modc5 = lm(r2~1, data=all.data.f)

models_c = list(modc1, modc2, modc3, modc4, modc5)
mod.names <- c('itd*sampEff','itd+sampEff','ITD', 'sampEff', '1')
aictab(cand.set = models_c, modnames = mod.names)

# LMs: betweenness ~
modb1 = lm(r2~body_size_mean*number_of_articles, data=all.data.f)
modb2 = lm(r2~body_size_mean+number_of_articles, data=all.data.f)
modb3 = lm(r2~body_size_mean, data=all.data.f)
modb4 = lm(r2~number_of_articles, data=all.data.f)
modb5 = lm(r2~1, data=all.data.f)

models_b = list(modb1, modb2, modb3, modb4, modb5)
mod.names <- c('itd*sampEff','itd+sampEff','ITD', 'sampEff', '1')
aictab(cand.set = models_b, modnames = mod.names)

#########################################
# TABLE S5: CORRELATION BETWEEN METRICS #
#########################################

# Select variables
cor_var = select(all.data.f, degree, betweenness, closeness, specialization_d)

# Correlation matrix
cor.mat = cor_mat(cor_var, method="pearson", conf.level=.95)
cor.mat

# Correlation matrix (p-values)
cor.pmat = cor_pmat(cor_var, method="pearson", conf.level=.95)
cor.pmat

# Correlation Plot
cor.mat %>%
  cor_reorder() %>%
  pull_lower_triangle() %>%
  cor_plot(label=T)

######################################
# FIGURE S4: ANCESTRAL STATES OF ITD #
######################################

itd.data = all.data.f[!(all.data.f$phylogeny=='no'),]
itd.data = itd.data[c(1,23)] # select columns of species and ITD
log_ITD_vector = log(itd.data$body_size_mean)# calculate a vector of log ITD
itd.data$log_ITD = log_ITD_vector # add the vector as a new column in the dataframe
itd.data = data.frame(itd.data[,-1], row.names=itd.data[,1]) # convert the column "bee species" as row names
itd.data = itd.data[!(row.names(itd.data)) %in% c(""),] # remove species in the dataframe but not in the tree
ITD_setnames = setNames(itd.data$body_size_mean,rownames(itd.data))
obj.ITD = name.check(tree,ITD_setnames)
tree.cut.ITD = drop.tip(tree,obj.ITD$tree_not_data)
name.check(tree.cut.ITD,ITD_setnames)

fit.ITD = fastAnc(tree.cut.ITD, ITD_setnames, vars=TRUE,CI=TRUE)
print(fit.ITD,printlen=10)
obj.ITD=contMap(tree.cut.ITD,ITD_setnames,plot=FALSE)
obj.ITD = setMap(obj.ITD, invert=TRUE)

pdf(file = "ASR_ITD.pdf",   # The directory you want to save the file in
    width = 30, # The width of the plot in inches
    height = 30) # The height of the plot in inches
plot(obj.ITD,fsize=c(3,3), outline=FALSE, lwd=c(15,10), leg.txt="Body size (mm)")
dev.off()

########################################################
# FIGURE S5: IS THERE BIAS ON SPECIALIATION D' #########
# FROM METHODS OF IDENTIFICATION OF BOTANICAL SOURCES? #
########################################################

# Vectors with d values
all.degree = as.numeric(all$`higher level`$d)
l.all = rep("all", length(all.degree))
chem.degree = as.numeric(all.chem$`higher level`$d)
l.chem = rep("Chem.", length(chem.degree))
field.degree = as.numeric(all.field$`higher level`$d)
l.field = rep("Field.", length(field.degree))
pal.degree = as.numeric(all.pal$`higher level`$d)
l.pal = rep("Palyn.", length(pal.degree))

# Merge data
degree = as.data.frame(c(all.degree, chem.degree, field.degree, pal.degree))
degree$method = c(l.all, l.chem, l.field, l.pal)
colnames(degree) = c("Specialization", "method")

# ANOVA
anova = (aov(data=degree, Specialization~method))
summary(anova)

# Tukey-test
TukeyHSD(anova, conf.level = .95)

# Plots
pdf(file = "figures/S5_methods.pdf",
    width = 20, # The width of the plot in inches
    height = 10) # The height of the plot in inches
par(mfrow=c(1,2))
boxplot(degree$Specialization~degree$method,
        col=c("gray", "orangered","orange", "mediumpurple"),
        ylab="Specialization d'", 
        xlab="Method to determine botanical source of resins")
plot(TukeyHSD(anova, conf.level = .95))
dev.off()
