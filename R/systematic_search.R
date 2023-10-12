###########
# 1 SETUP #
###########

# 1.1 Set directory

setwd("home/danimelsz/Desktop/Stingless_Bess/GitHub")

# 1.2 Load Packages

library(devtools)
library(dplyr)
library(ggplot2)
library(ggraph)
library(igraph)
library(litsearchr)
library(readr)
library(rworldmap)

##################
# 2 NAIVE SEARCH #
##################

# 2.1 Web Of Science

# ((propolis OR geopropolis OR resin*) AND ('resin source*' OR 'botanical source*' OR source*) AND ('stingless bee*' OR Meliponini OR Aparatrigona OR Camargoia OR Celetrigona OR Cephalotrigona OR Dolichotrigona OR Duckeola OR Friesella OR Frieseomelitta OR Geotrigona OR Lestrimelitta OR Leurotrigona OR Melipona OR Mourella OR Nannotrigona OR Nogueirapis OR Oxytrigona OR Paratrigona OR Partamona OR Plebeia OR Ptilotrigona OR Scaptotrigona OR Scaura OR Schwarziana OR Tetragona OR Tetragonisca OR Trigona OR Trigonisca))

# Load the WoS dataset

wos_naive = import_results(file = "data/naive_search/wos_naive.bib")
 
# Check the number of papers found in WoS naive search
print(nrow(wos_naive)) # 75 papers found in WoS naive search

# 2.2 Scopus

# ( TITLE-ABS-KEY ( ( propolis  OR  geopropolis  OR  resin ) ) )  AND  ( ( resin  AND source*  OR  botanical  AND source*  OR  source* ) )  AND  ( stingless  AND bee*  OR  meliponini  OR  aparatrigona  OR  camargoia  OR  celetrigona  OR  cephalotrigona  OR  dolichotrigona  OR  duckeola  OR  friesella  OR  frieseomelitta  OR  geotrigona  OR  lestrimelitta  OR  leurotrigona  OR  melipona  OR  mourella  OR  nannotrigona  OR  nogueirapis  OR  oxytrigona  OR  paratrigona  OR  partamona  OR  plebeia  OR  ptilotrigona  OR  scaptotrigona  OR  scaura  OR  schwarziana  OR  tetragona  OR  tetragonisca  OR  trigona  OR  trigonisca ) 

# Load the Scopus dataset

scopus_naive = import_results(file = "naive_search/scopus_naive.bib")

# Check the number of papers found in Scopus naive search
print(nrow(scopus_naive)) # 192 papers found in Scopus naive search

# Rename column name of the dataframe
colnames(scopus_naive)

# Check the title of the first paper in Scopus naive table
scopus_naive[1, "title"]

# 2.3 Getting potential search terms

# Extract potential terms from Scopus keywords
scopus_terms = extract_terms(keywords = scopus_naive[, "author_keywords"], method = "tagged")

# Extract new terms from WoS keywords
wos_terms = extract_terms(keywords = wos_naive[, "keywords"], method = "tagged")

# Merge scopus_terms and wos_terms
terms = unique(c(scopus_terms, wos_terms))

# The lists above indicate that some keywords might be included in a new search. We decided to add "bee-plant interactions", "botanical origin", and "plant origin". Therefore, we performed a new search in both WoS and Scopus using the following terms:
# ((propolis OR geopropolis OR resin*) AND ('resin source*' OR 'botanical source*' OR source* OR 'botanical origin*' OR 'plant origin*' OR bee-plant interaction*) AND ('stingless bee*' OR Meliponini OR Aparatrigona OR Camargoia OR Celetrigona OR Cephalotrigona OR Dolichotrigona OR Duckeola OR Friesella OR Frieseomelitta OR Geotrigona OR Lestrimelitta OR Leurotrigona OR Melipona OR Mourella OR Nannotrigona OR Nogueirapis OR Oxytrigona OR Paratrigona OR Partamona OR Plebeia OR Ptilotrigona OR Scaptotrigona OR Scaura OR Schwarziana OR Tetragona OR Tetragonisca OR Trigona OR Trigonisca))

# 2.4 Deduplicating

# First, I created the file 'merged_naive.csv' in Excel, including only relevant columns and merging both Scopus and WoS tables. Now, let's remove duplicates.

# Load merged_naive.csv
merged_naive = read.csv(file = 'data/naive_search/merged_naive.csv')

# Remove duplicates
 deduplicated_naive = litsearchr::remove_duplicates(df = merged_naive, field = 'Title', method = 'exact')

# Note that the deduplicated table still presents some duplications because exact method distinguish between uppercase and lowercase letters. So, I removed by eye the following lines:
deduplicated_naive_2 = deduplicated_naive[-c(193, 194, 195, 197, 198, 199, 200, 205, 206, 211, 214, 217, 220, 222, 225, 227, 228, 229, 234, 235, 240, 241, 244, 246, 247, 248, 249, 250, 251, 252, 253, 257, 258, 260, 261, 262, 263, 267),] 

# Export the dataframe as csv
write.csv(deduplicated_naive_2, "data/naive_search/deduplicated_merged_naive.csv", row.names = F)

#######################
# 3 SYSTEMATIC SEARCH #
#######################

# 3.1 Co-occurrence networks
 
# Merge titles and abstracts
docs = paste(deduplicated_naive_2[, "Title"], deduplicated_naive_2[, "Abstract"])

# Create a matrix of terms
dfm = create_dfm(elements = docs, features = terms)

# Remove non-important terms
dfm2 = dfm[,-6]
dfm2 = dfm2[,-22]
dfm2 = dfm2[,-36]
dfm2 = dfm2[,-42]

# Transform the matrix into a network
graph = create_network(dfm2, min_studies = 5, min_occ = 3) 

# Create a rank of node importance
strengths = strength(graph)
data.frame(term=names(strengths), strength=strengths, row.names=NULL) %>%
  mutate(rank=rank(strength, ties.method='min')) %>%
  arrange(strength) ->
  term_strengths

# The top terms are the most weakly linked to the others. Therefore, some of them should not be included in the final literature search. How to decide a cutoff to separate weak and strong terms? Based on the aforementioned rank, litsearchr uses an algorithm to find a node importance cutoff.

# Get cutoff
cutoff = find_cutoff(graph, method="cumulative", percent=0.8)

# Get a list of terms using the cutoff
g_redux = reduce_graph(graph, cutoff)
selected_terms = get_keywords(g_redux)
selected_terms

# 3.2 Final search

# After the revision of search terms through co-occurrence network, we need to group our terms in subtopics. Our naive search had three subtopics: 
# (i) bee product (i.e., resins)
# (ii) product source (i.e., plants)
# (iii) taxonomic group definition (i.e., stingless bees)
# Although there are algorithms to group terms into clusters, these are not so reliable. Therefore, litsearchr documentation recommends doing this step by eye. This is not a problem because the list of terms is not so extensive.

# Group interesting terms into clusters 
grouped_terms = list(
  product = selected_terms[12],
  source = selected_terms[8],
  group = selected_terms[c(14, 15)]
)

# 3.2.1 Final search in WoS
# ((propolis OR geopropolis OR resin* OR 'plant resin') AND ('resin source*' OR 'botanical source*' OR source* OR 'botanical origin') AND ('stingless bee*' OR Meliponini OR Aparatrigona OR Camargoia OR Celetrigona OR Cephalotrigona OR Dolichotrigona OR Duckeola OR Friesella OR Frieseomelitta OR Geotrigona OR Lestrimelitta OR Leurotrigona OR Melipona OR Mourella OR Nannotrigona OR Nogueirapis OR Oxytrigona OR Paratrigona OR Partamona OR Plebeia OR Ptilotrigona OR Scaptotrigona OR Scaura OR Schwarziana OR Tetragona OR Tetragonisca OR Trigona OR Trigonisca))

# Input data from WoS
wos_co_occurrence = read.csv(file = 'data/co_occurrence_search/wos_final_dataset_papers.csv')

# Remove duplicates from the merged table 
 deduplicated_wos_post_network = litsearchr::remove_duplicates(df = wos_co_occurrence, field = 'Title', method = 'exact')

# Write a CSV file with deduplicated papers after co-occurrence network
 write.csv(deduplicated_wos_post_network, "data/co_occurrence_search/deduplicated_wos_post_network.csv", row.names = F)

# 3.2.2 Final search in Scopus
# ( TITLE-ABS-KEY ( propolis  OR  geopropolis  OR  resin  OR  plant  AND resin ) )  AND  ( ( ( resin  AND source*  OR  botanical  AND source*  OR  source*  OR  botanical  AND origin ) )  AND  ( stingless  AND bee*  OR  meliponini ) )  AND  ( stingless  AND bee*  OR  meliponini  OR  aparatrigona  OR  camargoia  OR  celetrigona  OR  cephalotrigona  OR  dolichotrigona  OR  duckeola  OR  friesella  OR  frieseomelitta  OR  geotrigona  OR  lestrimelitta  OR  leurotrigona  OR  melipona  OR  mourella  OR  nannotrigona  OR  nogueirapis  OR  oxytrigona  OR  paratrigona  OR  partamona  OR  plebeia  OR  ptilotrigona  OR  scaptotrigona  OR  scaura  OR  schwarziana  OR  tetragona  OR  tetragonisca  OR  trigona  OR  trigonisca ) 

# Input data from Scopus
scopus_co_occurrence = read.csv(file = 'data/co_occurrence_search/scopus_post_network.csv')

# Remove duplicates from the merged table 
 deduplicated_final = litsearchr::remove_duplicates(df = scopus_co_occurrence, field = 'Title', method = 'exact')
 deduplicated_final2 = deduplicated_final[-c(170, 234, 240, 260),] 
 deduplicated_final2 = deduplicated_final[-c(259, 262),] 

# Write a CSV file with deduplicated papers after co-occurrence network
write.csv(deduplicated_wos_post_network, "data/co_occurrence_search/final_dataset_litsearchr.csv", row.names = F)

#############################
# 4 SUPPLEMENTARY FIGURE S1 #
#############################

# 4.1 Co-occurrence network

png(file="co_occurrence_search/figS2a_resins.png",
width=1000, height=1000)

ggraph(graph, layout="stress") +
  coord_fixed() +
  expand_limits(x=c(-3, 3)) +
  geom_edge_link(aes(alpha=weight)) +
  geom_node_point(shape="circle filled", fill="white") +
  geom_node_text(aes(label=name), hjust="outward", check_overlap=TRUE) +
  guides(edge_alpha=FALSE)

dev.off()

# 4.2 Rank of terms

png(file="co_occurrence_search/figS2b_resins.png",
width=800, height=800)

cutoff_fig <- ggplot(term_strengths, aes(x=rank, y=strength, label=term)) +
  geom_line() +
  geom_point() +
  geom_text(data=filter(term_strengths, rank>5), hjust="right", nudge_y=20, check_overlap=TRUE)

cutoff_fig +
  geom_hline(yintercept = cutoff, linetype="dashed")

dev.off()
