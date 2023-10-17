#load libraries----
library(ape)
library(phylotate)
library(phytools)
library(phangorn)
library(rworldmap)
library(ggplot2)
library(ggmap)
library(RColorBrewer)
library(PCMBase)
library(TreeTools)

#MCCT geo augmented from Koile et al. 2023 - downloaded from original paper

tree.yes.geo.annot.augmented <- read_annotated("Trees/03-linguistic_geographic_augmented/ba_ok.tree", format="nexus")
#extract lat and long----
temp <- tree.yes.geo.annot.augmented$node.comment

#tips, internal nodes, and root
tmp <- as.character(gsub(".*location=\\{(.+)\\}.*", "\\1", temp))

tmp2 <-  as.numeric(unlist(strsplit(tmp, split=",")))
tmp3 <- as.data.frame(matrix(tmp2,ncol=2,byrow=T))

colnames(tmp3) <- c("Longitude","Latitude") 
coordinates.augmented <- tmp3

tree.yes.geo.annot.augmented

temp <- tree.yes.geo.annot.augmented$node.comment

# Calculate patristic distances between languages in Table S13
library(adephylo)
tips=c("D311_Bila", "C83_Bushong", "B31_Tsogo",
       "B71aLe_Teke_Leconi", "B602_Kaningi_Nord", "Glotto_Kinyarwanda",
       "Glotto_Gyele", "B22b_Koya", "B22b_Ngom", "B52_Nzebi")
dist_lang <- as.matrix(distTips(tree.yes.geo.annot.augmented, tips=tips , method="patristic"))
dist_lang <- library(dplyr)

#bujeba = gyele, Mangbutu = Efe
dist_now <- dist_lang[tips,]

dist_now  <- as.data.frame(dist_now)
dist_now  <- dist_now  %>%
  dplyr::select(tips)

max_tree_dist <- max(dist_lang)

write.csv(dist_now, "distances_cahg_languages_patristic.csv")

library(reshape2)
data_now <- as.matrix(dist_now)
data_list <- melt(data_now)

# read in ID by group to group all the distances by group together (average across different languages)
cahg_lang <- read.csv("data_cahg_languages_koile.csv")
colnames(cahg_lang)[1] <- "Var2"
cahg_distances <- merge(data_list, cahg_lang, by="Var2")
colnames(cahg_lang)[1] <- "Var1"
cahg_distances <- merge(cahg_distances, cahg_lang, by="Var1")

# aggregate all the rows that are the same for two columns
cahg_distances$value <- as.numeric(cahg_distances$value)
cahg_dist <- cahg_distances %>%
  group_by(Group.x, Group.y) %>%
  summarise_at(vars(value),
               list(mean_val=mean))

write.csv(cahg_dist, "cahg_distance_language_matrix_koile.csv")
write.csv(max_tree_dist, "max_phylo_dist_koile.csv")
##################################################################

# read in pmi distances and patristic distances
shared_words <- read.csv("shared_words_counts_september23.csv")
shared_music <- subset(shared_words, shared_words$type=="music")
adj_pmi <- xtabs(linguistic_distance ~ ID2 + ID1, data = shared_music)
adj_patristic <- xtabs(koile_patristic ~ ID2 + ID1, data = shared_music)
adj_pmi[adj_pmi==0] <- NA
adj_patristic[adj_patristic==0] <- NA

library(gplots)
n=9
adj_pmi_rounded <- round(adj_pmi, 3)
mirrored_mat <- adj_pmi_rounded[, rev(seq_len(n))]
heatmap.2(mirrored_mat, Rowv = NA, Colv = NA, col = hcl.colors(12, "Geyser", rev=TRUE), trace = "none", dendrogram = "none",
        cellnote = mirrored_mat, notecol = "black", notecex = 0.8, cexRow = 0.8, cexCol = 0.8, main="PMI Linguistic distance")
dev.off()

adj_patristic_rounded <- round(adj_patristic, 3)
mirrored_mat <- adj_patristic_rounded[, rev(seq_len(n))]
heatmap.2(mirrored_mat, Rowv = NA, Colv = NA, col = hcl.colors(12, "Geyser", rev=TRUE), trace = "none", dendrogram = "none",
          cellnote = mirrored_mat, notecol = "black", notecex = 0.8, cexRow = 0.8, cexCol = 0.8, main="Patristic linguistic distance")
dev.off()


