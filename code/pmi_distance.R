# Process ASJP dataset and calculate pmi distances

library(dplyr)
library(reshape2)

# Read dataset as downloaded from ASJP website: https://asjp.clld.org/download
data <- read.csv("pmiWorld.csv")
data <- as.matrix(data)
rows <- data[,"X"]
data <- data[,-1]

rownames(data) <- rows
#bujeba = gyele, Mangbutu = Efe
data_now <- data[c("RWANDA_TWA_KINIGI",
                   "BUSHONG_1","BAKA_2", "TSOGO",
                    "TEKE", "KANINGI",
                   "BUJEBA", "TIKAR_AKUEN", "MANGBUTU", "NGOM_KELE",
                   "C10_GANDO",
                   "BOGONGO", "BIRA"),]

data_now <- as.data.frame(data_now)
data_now <- data_now %>%
  dplyr::select("RWANDA_TWA_KINIGI",
         "BUSHONG_1","BAKA_2", "TSOGO",
          "TEKE",  "KANINGI",
         "BUJEBA", "TIKAR_AKUEN", "MANGBUTU", "NGOM_KELE",
         "C10_GANDO",
         "BOGONGO", "BIRA")

write.csv(data_now, "distances_cahg_languages.csv")
data_now <- as.matrix(data_now)
data_list <- melt(data_now)

# read in ID by group to group all the distances by group together (average across different languages)
cahg_lang <- read.csv("data_cahg_languages.csv")
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

write.csv(cahg_dist, "cahg_distance_language_matrix.csv")
