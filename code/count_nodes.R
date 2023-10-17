

####################################################################################
# Calculate distance between CAHG languages with number of nodes in Glottolog tree #
####################################################################################

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
library(flextable)
library(phyloWeights)
library(glottoTrees)

# Glottolog v-4-8 tree downloaded from website
tree <- ape::read.tree("tree_glottolog_newick.txt")

# Convert multiPhylo to a list of phylogenetic trees
tree_list <- as.list(tree)

# Function to find the path between two languages in a tree
find_path <- function(tree, language_a, language_b) {
  # Find the MRCA (Most Recent Common Ancestor) of the two languages
  mrca <- getMRCA(tree, c(language_a, language_b))
  
  # Find the path from language A to MRCA
  path_a <- getDescendants(tree, mrca, "preorder")
  
  # Find the path from language B to MRCA (excluding MRCA)
  path_b <- getDescendants(tree, mrca, "preorder", right = TRUE)
  
  # Combine the paths
  path <- c(language_a, path_a, path_b[-1], language_b)
  
  return(path)
}

language_metadata_v4.8 <- get_glottolog_languages(glottolog_version = "4.8")
my_lang <- read.csv("data_cahg_languages.csv")
my_lang_data <- subset(language_metadata_v4.8, language_metadata_v4.8$glottocode %in% my_lang$Glottocode)

# I want to combine two trees where languages of cahg belong
tree_names <- c("CentralSudanic", "Atlantic-Congo")

# 
multiPhylo_cahg <- get_glottolog_trees(tree_names)
tree_cahg<- assemble_rake(multiPhylo_cahg)

# change labels to glottocodes
tree_cahg_abr <- abridge_labels(tree_cahg)

plot_glotto(tree_cahg_abr, nodelabels = FALSE)

my_lang_tips <- subset(my_lang_data, my_lang_data$position=="tip")

tree_cahg_only2 <- keep_as_tip(tree_cahg_abr, label = my_lang_data$glottocode)

tip_and_node_list <- tree_cahg_only2$tip.label

# NOW FOR ALL TIPS
calculate_nodes_between_tips <- function(tree, tip1, tip2) {
  # Find the MRCA (Most Recent Common Ancestor) of the two tips
  mrca <- getMRCA(tree, c(tip1, tip2))
  
  # Count the number of nodes (internal nodes) between the tips
  num_nodes_between_tips <- length(getDescendants(tree, node = mrca)) - 2
  
  # Subtract 2 to exclude the MRCA and the common ancestor of the root
  # (assuming the root is an internal node)
  
  return(num_nodes_between_tips)
}

# Initialize a distance matrix
num_tips <- length(tip_and_node_list)
distance_matrix <- matrix(NA, nrow = num_tips, ncol = num_tips,
                          dimnames = list(tip_and_node_list, tip_and_node_list))

# Calculate distances for all pairs of tips
for (i in 1:num_tips) {
  for (j in 1:num_tips) {
    if (i != j) {
      distance <- calculate_nodes_between_tips(tree_cahg_only2, tip_and_node_list[i], tip_and_node_list[j])
      distance_matrix[i, j] <- distance
    }
  }
}

# Print the distance matrix
print(distance_matrix)
library(reshape2)
dist_df <- as.data.frame(distance_matrix)
data_list <- melt(distance_matrix)
colnames(my_lang)[5] <- "Var2"
cahg_distances <- merge(data_list, my_lang, by="Var2")
colnames(my_lang)[5] <- "Var1"

cahg_distances <- merge(cahg_distances, my_lang, by="Var1")
cahg_distances$value <- as.numeric(cahg_distances$value)
cahg_dist <- cahg_distances %>%
  group_by(Group.x, Group.y) %>%
  summarise_at(vars(value),
               list(mean_val=mean))

write.csv(cahg_dist, "cahg_distance_language_matrix_glottolog.csv")




                                 
