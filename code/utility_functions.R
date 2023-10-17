# Functions required to run the code from Padilla-Iglesias et al. 2023

# THESE FUNCTIONS ARE ADAPTED FROM SHENNAN ET AL. 2015 EVOL.HUM.BEH.

#Compute Great Arc Distances
distMat<-function(data)
  {
require(argosfilter)
res<-matrix(0,nrow=dim(data)[1],ncol=dim(data)[1])
for (i in 1:dim(data)[1])
  {
    for (j in 1:dim(data)[1])
      {
        if(i!=j){
        res[i,j]<-distance(lat1=data$y[i],lat2=data$y[j],lon1=data$x[i],lon2=data$x[j])
      }
      }
  }
res[is.nan(res)]=0
return(as.dist(res))
  }

#
    
BinaryCulture<-function(cultureList)
    {
res<-matrix(NA,ncol=length(cultureList),nrow=length(cultureList))

for (i in 1:length(cultureList))
    {
        for (j in 1:length(cultureList))
            {
                if (i>j)
                    {
                        res[i,j]=as.numeric(cultureList[i]!=cultureList[j])
                    }
            }
    }
return(res)
}

GetFst <-function(cultureList)
    {
res<-matrix(NA,ncol=length(cultureList),nrow=length(cultureList))

for (i in 1:length(cultureList))
    {
        for (j in 1:length(cultureList))
            {
                if (i!=j & i>j)
                    {
                        ind_i <- which(rownames(Fst_mat) == cultureList[i])
                        ind_j <- which(colnames(Fst_mat) == cultureList[j])
                                #print(cultureList)
                                #print(ind_i)
                                #print(ind_j)
                                fst_val <- Fst_mat[ind_i, ind_j]
                                res[i,j] <- fst_val
                   
                    }
            }
    }
res[is.na(res)] <- 0
return(res)
}

## FROM HERE ONWARDS FUNCTIONS ADAPTED FROM MATSUMAE ET AL. 2021 SCIENCE ADVANCES

# Scree plot
# Theme
theme_scree <- function(...) {
  theme_minimal() +
  theme(
    text = element_text(color = "#22211d", family=my_font),
    axis.line = element_blank(),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.border = element_blank(),
    ...)}

# Plot function
plot_pc_scree <- function(eigenval, title, type="PC") {

  if (type == "PC") {
    pc_index = seq(1:nrow(eigenval))
    eigenval = data.frame(val = eigenval[, 1],
                          idx = pc_index,
                          rel = eigenval[, 2],
                          cum = eigenval[, 3])
    pc_legend = "Principal Components"
    subtitle= paste("Explained variance by", pc_legend)
    pc_labels = paste(rep(c("PC"), nrow(eigenval)), seq(1:nrow(eigenval)))
    ann_offset_y = 6}
  
  else if (type == "PCo"){
    pc_index = seq(1:length(eigenval))
    eigenval = data.frame(val = eigenval,
                          idx = pc_index,
                          rel = eigenval/sum(eigenval)*100,
                          cum = cumsum(eigenval/sum(eigenval)*100))
    pc_legend = "Principal Coordinates"
    subtitle= paste("Explained variance by", pc_legend)
    pc_labels = paste(rep(c("PCo"), nrow(eigenval)), seq(1:nrow(eigenval)))
    ann_offset_y = 10}
  
  else stop("Type must be either PC or PCo")
  p <- ggplot(data=eigenval,aes(x=idx, y=rel))
  p <- p + geom_bar(stat="identity")
  p <- p + xlab(pc_legend) + ylab("Eigenvalues (%)")
  p <- p + geom_line(data=eigenval, aes(x=idx, y=cum), colour="#4fb6ca")
  p <- p + geom_point(data=eigenval, aes(x=idx, y=cum), colour="#4fb6ca")
  p <- p + scale_x_discrete(limits=pc_labels)
  p <- p + labs(title=title, subtitle=subtitle)
  p <- p + annotate ("label", x = tail(pc_index, n=1)-2.3,
                     y = tail(eigenval$cum, n=1) - ann_offset_y,
                     label="Cumulative eigenvalues",
                     colour ="#4fb6ca", label.size=NA, size=3.5)
  p <- p + theme_scree()
  return(p)}

# Heat maps
# Theme
theme_heat_maps <- function(...) {
  theme_minimal() +
  theme(
    text = element_text(color = "#22211d", family=my_font),
    legend.background = element_rect(fill = NA, color = NA),
    legend.direction = "horizontal",
    legend.position = "bottom",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size=15),
    axis.text.x = element_text(vjust = 0, angle = 90, hjust=0),
    plot.margin=unit(c(3,1,1.5,1.2),"cm"),
   
    ...)}

# Plot function
heat_map_pc <- function(pcs, type="PCo"){
  #' Plots a heat map of the PCs/PCos
  #' @param pcs: the principal components/coordinates
  #' @param type: either PC or PCo
  #' @return the heat map (ggplot)
  

  cnames <- colnames(pcs)
  
  # Normalize the PCs/PCos for each factor
  pcs_norm <- sapply(cnames, function (y) {
    x <- pcs[, y]
    x_norm <- (x-min(x))/(max(x)-min(x))
    return (x_norm)}, USE.NAMES = T)
  
  if (type == "PC"){
    primary_pc_name <- "PC 1"
    colnames(pcs_norm) <- paste("PC", seq(1, ncol(pcs_norm)))}
  if (type == "PCo"){
    primary_pc_name <- "PCo 1"
    colnames(pcs_norm) <- paste("PCo", seq(1, ncol(pcs_norm)))}
  
  rownames(pcs_norm) <- rownames(pcs)
  melted <- melt(pcs_norm)
  
  order_x <- melted %>%
  filter(Var2 == primary_pc_name) %>%
  arrange(desc(value)) %>%
  pull(Var1)
  
  # Create heat map
  h <- ggplot(data = melted, aes(ordered(Var1, levels=order_x),
                                 ordered(Var2, levels=rev(levels(Var2))),
                                 fill=value))
  
  h <- h + geom_tile(colour="grey")
  h <- h + scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                limits = c(0,1), midpoint=0.5, name = NULL)
  h <- h + theme_heat_maps()
  h <- h + coord_fixed()
  h <- h + scale_x_discrete(position = "top")
  return (h)}

print_dist <- function(d, caption) {
  #' Prints the distance matrices in a table
  #' @param d: the input distance matrix
  #' @param caption: the caption of the table
  d.m <- as.matrix(sort_dist_mat(d))
  d.m[upper.tri(d.m, diag=T)] <- NA
  colnames(d.m) <- abbreviate(colnames(d.m), minlength=8)
  rownames(d.m) <- abbreviate(rownames(d.m), minlength=8)
  options(knitr.kable.NA = '')
  d.m <- d.m[2:nrow(d.m),1:ncol(d.m)-1]
  
  kable(d.m, digits=3, format = 'latex', caption=caption) %>%
    kable_styling(latex_options = c("scale_down", "hold_position")) %>%
                  column_spec(1, border_left=T) %>%
                  column_spec(ncol(d.m)+1, border_right=T)}

# Plot function
plot_rda_as_matrix <- function(rda_matrix){
  #' Plots the results from the redundancy analysis
  #' @param rda_matrix: a matrix comprising the rda results
  #' @return the rda_matrix plot
  c <- ggplot(data = rda_matrix, aes(response, explanatory, fill = r2_adj))
  c <- c + geom_tile(color = "white")
  c <- c + scale_fill_gradient(low = "white", high = "red", na.value = "white",
                               limits = c(0,1),
  name=expression(paste("Explained variance \nin response (adjusted ",R^{2},")    ")))
  
  c <- c + geom_text(aes(response, explanatory,
                         label = paste(round(r2_adj,2), sig_level)),
                     color = "black", size = 3)

  c <- c + xlab("\nresponse") + ylab("explanatory\n")
  c <- c + theme_corr_rda()
  c <- c + coord_fixed()
  return(c)}

