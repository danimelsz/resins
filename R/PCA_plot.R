# PCA among MBSDs

# packages
if(!require(dplyr)){install.packages("dplyr")}
if(!require(factoextra)){install.packages("factoextra")}
if(!require(ggforce)){install.packages("ggforce")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(MASS)){install.packages("MASS")}
if(!require(vegan)){install.packages("vegan")}

# data
data = read_csv("data/nmds_FAIL.csv")

# ANOSIM
data_ = data[-c(1,27,39,4,22,25),]
datab = data_[, 3:ncol(data)]
anosim(datab, data_$MBSD, distance = "bray", permutations = 1000)

# analysis
nmds_result <- metaMDS(data[, 3:ncol(data)], k = 2, trymax = 100)
nmds_data <- data.frame(nmds_result$points, Group = data$MBSD)
nmds_data <- nmds_data[-c(1,27,39,4,22,25), ]
nmds_data$NMDS1 = scale(nmds_data$MDS1)
nmds_data$NMDS2 = scale(nmds_data$MDS2)
nmds_data %>%
  group_by(Group) %>%
  summarise(variance = var(NMDS1, na.rm = TRUE)) # calculate variance per group in NMDS1
nmds_data %>%
  group_by(Group) %>%
  summarise(variance = var(NMDS2, na.rm = TRUE)) # calculate variance per group in NMDS2

# Plot NMDS solution with points colored by the grouping variable
png(filename= "figures/pca_nmds_MBSD.png", res= 300, height= 3000, width= 3600)
my_colors <- c("brown3", "darkorange", "darkslateblue")
ggplot(nmds_data, aes(NMDS1, NMDS2, color = Group)) +
  geom_point(size=5, alpha=.8) +
  scale_color_manual(values = my_colors)+
  #geom_jitter(width = 1000, height = 1000, alpha=.5, size=3) +
  theme_classic() +
  labs(x = "\nNMDS 1", y = "NMDS 2\n")+
  theme(axis.text.x = element_text(size = 16.5), # size of x and y values
        axis.text.y = element_text(size = 16.5)) +
  theme(axis.title.x = element_text(size = 25),  # Set x-axis title font size
        axis.title.y = element_text(size = 25)) +
  theme(legend.text = element_text(size = 16)) + # Set the font size (here, 12 points)
  theme(legend.title = element_blank()) # hide legend title
dev.off()

# PCA
pca_result <- prcomp(data[, 3:ncol(data)], scale. = TRUE)
pca.data = as.data.frame(pca_result$x[, 1:2])

png(filename= "figures/pca_MBSD.png", res= 300, height= 3000, width= 3600)
my_colors <- c("brown3", "darkorange", "darkslateblue")
fviz_pca_ind(pca_result, geom.ind = "point", col.ind = data$MBSD,
             #addEllipses = TRUE, ellipse.level = 0.95,
             palette = my_colors, pointsize = 7,alpha.ind = 0.8, axes.linetype=NA) +
  theme_classic()+
  xlim(-7,3) + ylim(-4.5,3)+
  xlab("\nPC1") + ylab("PC2\n")+
  theme(axis.text.x = element_text(size = 16.5), # size of x and y values
        axis.text.y = element_text(size = 16.5)) +
  theme(axis.title.x = element_text(size = 25),  # Set x-axis title font size
        axis.title.y = element_text(size = 25)) +
  theme(legend.text = element_text(size = 17)) + # Set the font size (here, 12 points)
  theme(legend.title = element_blank()) + # hide legend title
  ggtitle("") # hide plot title
dev.off()
