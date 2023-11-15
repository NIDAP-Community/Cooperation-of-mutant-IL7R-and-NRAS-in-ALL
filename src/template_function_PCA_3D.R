# PCA 3D [CCBR] (0309860a-90f0-490d-8e61-69b483923da1): v37
PCA_3D <- function(Batch_Corrected_Counts, ccbr991_metadata_nidap) {

## This function generates an interactive 3D PCA plot

## --------- ##
## Libraries ##
## --------- ##

library(edgeR)
library(tidyverse)
library(colorspace)
library(plotly)
library(htmlwidgets)
    
## -------------------------------- ##
## User-Defined Template Parameters ##
## -------------------------------- ##

# Basic parmaters
Counts_Matrix <- Batch_Corrected_Counts
Sample_Metadata_Table <- ccbr991_metadata_nidap
FeatureID_Name_Column <- "Gene"
Sample_Names_Column <- "Sample"
Samples_to_Include <- c("control1","control2","control3","IL7R1","IL7R2","IL7R3","NRAS1","NRAS2","NRAS3","Both1","Both2","Both3")
Group_Column <- "Group"
Plot_Labels_Column <-"Label"
Samples_to_Rename_Manually <- c()
# Advanced parameters
Outlier_Samples_to_Remove <- c()
Use_CPM <- FALSE  
# Visualization parameters
Points_Size <-8
Point_Color <- c()
Plot_Title <- "PCA 3D"    

##--------------- ##
## Error Messages ##
## -------------- ##

## --------- ##
## Functions ##
## --------- ##

# generate random colors
getourrandomcolors <- function(k) {
  seed = 10
  n <- 2e3
  ourColorSpace <- colorspace::RGB(runif(n), runif(n), runif(n))
  ourColorSpace <- as(ourColorSpace, "LAB")
  currentColorSpace <- ourColorSpace@coords
  # Set iter.max to 20 to avoid convergence warnings.
  set.seed(seed)  
  km <- kmeans(currentColorSpace, k, iter.max = 20)
  return(unname(hex(LAB(km$centers))))
}

## --------------- ##
## Main Code Block ##
## --------------- ##

###########################################
#This code block does input data validation

# get the sample data

# remove outliers
Samples_to_Include <-
  Samples_to_Include[!Samples_to_Include %in% Outlier_Samples_to_Remove]
# include samples
Samples_to_Include <-
  Samples_to_Include[Samples_to_Include != FeatureID_Name_Column]
  cat(sprintf("Total number of samples included in the PCA plot: %g", length(Samples_to_Include)))
Counts_Matrix[, append(FeatureID_Name_Column, Samples_to_Include)] -> df
Sample_Metadata_Table <-
  Sample_Metadata_Table[Sample_Metadata_Table[, Sample_Names_Column]  %in% Samples_to_Include, ]
Sampinfo <- Sample_Metadata_Table
colnames(df)[colnames(df) == FeatureID_Name_Column] <- "Gene"
df -> edf.orig

###### PCA plot ##############

# evaluate and display PCA plot
rownames(Sampinfo) <- Sampinfo[, Sample_Names_Column]
Sampinfo <-
  Sampinfo[match(colnames(edf.orig[, -1]), Sampinfo[, Sample_Names_Column]),]
Sampinfo = Sampinfo[complete.cases(Sampinfo[, Sample_Names_Column]), ]
cat(paste0(
  "\nTotal number of genes included in the PCA plot: ",
  nrow(edf.orig)
))
edf <- edf.orig[, match(Sampinfo[, Sample_Names_Column], colnames(edf.orig))]
idx = !rowMeans(edf) == 0
if (Use_CPM) {
  edf = edgeR::cpm(edf[idx, ])
}
edf.orig = edf.orig[idx, ]
tedf <- t(edf)
colnames(tedf) <- edf.orig[, 1]
tedf <- tedf[, colSums(is.na(tedf)) != nrow(tedf)]
tedf <- tedf[, apply(tedf, 2, var) != 0]
pca <- prcomp(tedf, scale. = T)
pca.df <- as.data.frame(pca$x)
pca.df$group <- Sampinfo[, Group_Column]
pca.df$sample <- Sampinfo[, Plot_Labels_Column]

# pca stats
stats <-
  data.frame(id = paste0("PC", 1:length(pca$sdev)),
             eigenvalue = (pca$sdev) ^ 2) %>% mutate(variance = eigenvalue * 100 / sum(eigenvalue)) %>% mutate(cumvariance =
                                                                                                                 cumsum(variance)) %>% mutate(
                                                                                                                   variance_label = sprintf("%s (%.1f%% variance)", id, variance),
                                                                                                                   cumvariance_label = sprintf("%s (%.1f%% cumulative variance)", id, cumvariance)
                                                                                                                 )
axis_title = sub(" variance", "", stats$variance_label)
# if rename samplea
if (!is.null(Samples_to_Rename_Manually)) {
  if (Samples_to_Rename_Manually != c("")) {
    for (x in Samples_to_Rename_Manually) {
      old <- strsplit(x, ": ?")[[1]][1]
      new <- strsplit(x, ": ?")[[1]][2]
      pca.df$sample <-
        ifelse(pca.df$sample == old, new, pca.df$sample)
    }
  }
}
# set the colors to be used in the plot
if (length(unique(Sampinfo[, Group_Column])) > length(Point_Color)) {
  ## Original color-picking code.
  k = length(unique(Sampinfo[, Group_Column])) - length(Point_Color)
  more_cols <- getourrandomcolors(k)
  Point_Color <- c(Point_Color , more_cols)
} else {
  Point_Color <- Point_Color[1:length(unique(Sampinfo[, Group_Column]))]
}
names(Point_Color) <- unique(Sampinfo[, Group_Column])
# set label size (may add to user-defined parameters in the future)
Label_Font_Size <- 24

# plot PCA
cat("\nRunning PCA...")
fig <- plot_ly(
    pca.df,
    x = ~ PC1,
    y = ~ PC2,
    z = ~ PC3,
    color = ~ group,
    colors = Point_Color,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = Points_Size),
    hoverinfo = 'text',
    text = ~ sample,
    size = Label_Font_Size
  )
legend = TRUE
if (legend == FALSE) {
    fig <- hide_legend(fig)
}
fig <-
    fig %>% layout(
      title = list(text = Plot_Title, size = 5),
      scene = list(
        xaxis = list(title = axis_title[1]),
        yaxis = list(title = axis_title[2]),
        zaxis = list(title = axis_title[3])
    ),
      legend = list(
        itemsizing = 'constant',
        size = 12,
        y = 0.5
    )
)
  
###### OUTPUT ##############
# save in dataset container
# auto removed: output <- new.output()
# auto removed: output_fs <- output$fileSystem()
# html widget
html_file <- sprintf("Plot_%s.html", gsub(" ", "_", Plot_Title))
htmlwidgets::saveWidget(fig, html_file)
# output_fs$upload(html_file, html_file)
cat(sprintf("\nFiles saved in the output dataset:\n• %s", html_file))
# components table
component_file <-
    sprintf("Components_%s.csv", gsub(" ", "_", Plot_Title))
    write.csv(stats, row.names = FALSE, file  (component_file, 'w'))
cat(sprintf("\n• %s", component_file))
# coordinates table
coordinate_file <-
    sprintf("Coordinates_%s.csv", gsub(" ", "_", Plot_Title))
write.csv(pca.df,
        row.names = FALSE,
        file(coordinate_file, 'w'))
cat(sprintf("\n• %s", coordinate_file))
  
  ## keep for later solution to quality of the downloaded static image
  #fig = config(fig, toImageButtonOptions = list(format= 'svg', # one of png, svg, jpeg, webp
  #filename= sub(".html","", html_file), height= 300, width= 500,           scale= 1))  
  
# display visualization
print(fig)

}

## ---------------------------- ##
## Global Imports and Functions ##
## ---------------------------- ##

## Functions defined here will be available to call in
## the code for any table.

## --------------- ##
## End of Template ##
## --------------- ##

print("template_function_PCA_3D.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Batch_Corrected_Counts<-readRDS(paste0(rds_output,"/var_Batch_Corrected_Counts.rds"))
Input_is_Seurat_count <- 0
for(item in var_Batch_Corrected_Counts){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Batch_Corrected_Counts<-as.data.frame(var_Batch_Corrected_Counts)}else{var_Batch_Corrected_Counts <- var_Batch_Corrected_Counts}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_ccbr991_metadata_nidap<-readRDS(paste0(rds_output,"/var_ccbr991_metadata_nidap.rds"))
Input_is_Seurat_count <- 0
for(item in var_ccbr991_metadata_nidap){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_ccbr991_metadata_nidap<-as.data.frame(var_ccbr991_metadata_nidap)}else{var_ccbr991_metadata_nidap <- var_ccbr991_metadata_nidap}
invisible(graphics.off())
var_PCA_3D<-PCA_3D(var_Batch_Corrected_Counts,var_ccbr991_metadata_nidap)
invisible(graphics.off())
saveRDS(var_PCA_3D, paste0(rds_output,"/var_PCA_3D.rds"))
