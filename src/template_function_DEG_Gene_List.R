# DEG Gene List [CCBR] (f7fc094e-d625-4053-8f61-b44df3178260): v42
DEG_Gene_List <- function(DEG_Analysis) {

## This function filters DEG table

## --------- ##
## Libraries ##
## --------- ##

library(tidyverse)
library(dplyr)
library(tidyselect)
library(tibble)
library(ggplot2)
library(plotrix)
    
## -------------------------------- ##
## User-Defined Template Parameters ##
## -------------------------------- ##

# Basic parameters
deg_table <- DEG_Analysis 
gene_names_column <- "Gene"
significance_column <- "adjpval"
significance_cutoff <- 0.001
change_column <- "logFC"
change_cutoff <- 1
filtering_mode <- "in any contrast"

# Advanced parameters
include_estimates <- c("FC","logFC","tstat","pval","adjpval")
round_estimates <- TRUE
contrast_filter = "none"
contrasts = c()
groups = c()
groups_filter = "none"

# Visualization
label_font_size = 6
label_distance = 1
y_axis_expansion = 0.08
fill_colors =c("steelblue1","whitesmoke")
pie_chart_in_3d = TRUE
bar_width = 0.4
draw_bar_border = TRUE
force_barchart = TRUE
rounding_decimal_for_percent_cells = 0

## -------------------------------- ##
## Parameter Misspecifation Errors  ##
## -------------------------------- ##

## -------------------------------- ##
## Functions                        ##
## -------------------------------- ##

## --------------- ##
## Main Code Block ##
## --------------- ##

## If include_estimates param is empty due to template upgrade,
## then fill it with default values.
if (length(include_estimates) == 0){
    include_estimates <- c("FC","logFC","tstat","pval","adjpval")
}
## select DEG stat columns
estimates = paste0("_",include_estimates)
signif = paste0("_", significance_column)
change = paste0("_", change_column)
deg_table <- deg_table %>% dplyr::select(gene_names_column, ends_with(c(estimates, signif, change)))

contrasts_name = deg_table %>% dplyr::select(ends_with(signif)) %>% colnames() 
contrasts_name = unlist(strsplit(contrasts_name, signif))
if ( contrast_filter == "keep") {
    contrasts_name = intersect(contrasts_name, contrasts)    
} else if ( contrast_filter == "remove") {
    contrasts_name = setdiff(contrasts_name, contrasts)    
}
contrasts_name = paste0(contrasts_name, "_")

groups_name = deg_table %>% dplyr::select(ends_with(c("_mean","_sd"))) %>% colnames() 
groups_name = unique(gsub("_mean|_sd", "", groups_name))
if ( groups_filter == "keep") {
    groups_name = intersect(groups_name, groups)    
} else if ( contrast_filter == "remove") {
    groups_name = setdiff(groups_name, groups)    
}
groups_name = paste0(groups_name,"_")

deg_table <- deg_table %>% dplyr::select(gene_names_column, starts_with(c(groups_name,contrasts_name)))

## select filter variables
datsignif <- deg_table %>% 
    dplyr::select(gene_names_column, ends_with(signif)) %>%
    tibble::column_to_rownames(gene_names_column)
datchange <- deg_table %>%
    dplyr::select(gene_names_column, ends_with(change)) %>% tibble::column_to_rownames(gene_names_column)
genes <- deg_table[,gene_names_column]

## filter genes
if (filtering_mode == "in any contrast") {       
    significant <- apply(datsignif, 1, function(x) any(x <= significance_cutoff))
    changed <- apply(datchange, 1, function(x) any(abs(x) >= change_cutoff))
    select_genes <- genes[significant & changed]
    } else {
    significant <- apply(datsignif, 1, function(x) all(x <= significance_cutoff))
    changed <- apply(datchange, 1, function(x) all(abs(x) >= change_cutoff))
    select_genes <- genes[significant & changed]
}
# stop if 0 genes selected with the selection criteria
print(sprintf("Total number of genes selected with %s ≤ %g is %g", significance_column, significance_cutoff, sum(significant)))
print(sprintf("Total number of genes selected with %s ≤ %g and |%s| >= %g is %g", significance_column, significance_cutoff, change_column, change_cutoff, sum(significant&changed)))

if (length(select_genes) == 0) {
    stop("ERROR: Selection criteria select no genes - change stringency of the Significance cutoff and/or Change cutoff parameters")
}

##.output dataset
out <- deg_table %>% dplyr::filter(get(gene_names_column) %in% select_genes)
if (round_estimates) {
    out <- out %>% mutate_if(is.numeric, ~signif(., 3))
}

## do plot
significant <- apply( datsignif, 2, function(x) x <= significance_cutoff )  
changed <- apply( datchange, 2, function(x) abs(x) >= change_cutoff )
dd <- significant & changed
if (draw_bar_border){
    bar_border = 'black'
} else {
    bar_border = NA
}

## If fill_colors is blank due to template upgrade, then
## give it default values.
if(length(fill_colors) == 0){
    fill_colors <- c("steelblue1","whitesmoke")
}

if (filtering_mode == "in any contrast") {
    say_contrast = paste(colnames(dd), collapse=" | ")
    say_contrast = gsub("_pval|_adjpval","", say_contrast)
    tab <- reshape2::melt(apply(dd, 2, table))  %>% 
        dplyr::mutate(Significant=ifelse(Var1, "TRUE", "FALSE")) %>% 
        dplyr::mutate(Significant=factor(Significant, levels=c("TRUE","FALSE")), Count=value, Count_format=format(round(value, 1), nsmall=0, big.mark=",")) %>% 
        dplyr::mutate(Var2=gsub("_pval|_adjpval", "", Var2)) %>%
        group_by(Var2) %>% 
        dplyr::mutate(Percent=round(Count / sum(Count)*100,rounding_decimal_for_percent_cells)) %>% 
        dplyr::mutate(Label=sprintf("%s (%g%%)", Count_format, Percent))

    pp <- ggplot(tab, aes(x="", y=Count, labels=Significant, fill=Significant)) +
        geom_col(width=bar_width, position="dodge", col=bar_border) + facet_wrap(~Var2) +
        scale_fill_manual(values=fill_colors) + theme_bw(base_size=20) +
        xlab("Contrast") + ylab("Number of Genes") +
        geom_text(aes(label=Label), color=c("black"), size=label_font_size, position=position_dodge(width=bar_width), vjust=-label_distance) +
        theme(axis.ticks.x=element_blank(), axis.text.x=element_blank()) +
        ggtitle(sprintf("%s<%g & |%s|>%g %s: %s", significance_column, significance_cutoff, change_column, change_cutoff, filtering_mode, say_contrast)) + 
        theme(legend.key.size = unit(3,"line"), legend.position='top') +
        theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
        theme(strip.background = element_blank(), strip.text = element_text(size=26)) +
        xlab("") + 
        scale_y_continuous(name="", expand=c(y_axis_expansion, 0))
print(pp) 

} else {
    
    say_contrast = paste(colnames(dd), collapse=" & ")
    say_contrast = gsub("_pval|_adjpval","", say_contrast)
    dd <- apply(dd, 1, function(x) all(x==TRUE))
            
    if (force_barchart) {
        dd = data.frame(dd)   
        colnames(dd) = say_contrast                     
        tab <- reshape2::melt(apply(dd, 2, table))  %>% 
        dplyr::mutate(Significant=ifelse(Var1, "TRUE", "FALSE")) %>% 
        dplyr::mutate(Significant=factor(Significant, levels=c("TRUE","FALSE")), Count=value, Count_format=format(round(value, 1), nsmall=0, big.mark=",")) %>% 
        group_by(Var2) %>% 
        dplyr::mutate(Percent=round(Count / sum(Count)*100,rounding_decimal_for_percent_cells)) %>% 
        dplyr::mutate(Label=sprintf("%s (%g%%)", Count_format, Percent))

    pp <- ggplot(tab, aes(x="", y=Count, labels=Significant, fill=Significant)) +
        geom_col(width=bar_width, position="dodge", col=bar_border) + facet_wrap(~Var2) +
        scale_fill_manual(values=fill_colors) + theme_bw(base_size=20) +
        xlab("Contrast") + ylab("Number of Genes") +
        geom_text(aes(label=Label), color=c("black"), size=label_font_size, position=position_dodge(width=bar_width), vjust=-label_distance) +
        theme(axis.ticks.x=element_blank(), axis.text.x=element_blank()) +
        ggtitle(sprintf("%s<%g & |%s|>%g %s: %s", significance_column, significance_cutoff, change_column, change_cutoff, filtering_mode, say_contrast)) + 
        theme(legend.key.size = unit(3,"line"), legend.position='top') +
        theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
        theme(strip.background = element_blank(), strip.text = element_text(size=26)) +
        xlab("") + 
        scale_y_continuous(name="", expand=c(y_axis_expansion, 0))
print(pp) 
    } else {
    
    N = c( sum(dd), length(dd)-sum(dd))
    Nk = format(round(as.numeric(N), 1), nsmall=0, big.mark=",")
    P = round(N/sum(N)*100,rounding_decimal_for_percent_cells)
    if (label_font_size > 0) {
        labs = c(sprintf("Significant\n%s (%g%%)", Nk[1], P[1]) , sprintf("Non-Significant\n%s (%g%%)", Nk[2], P[2]))
    } else { 
        labs = NULL
    }
    if (pie_chart_in_3d) {
        pie3D(N, radius=0.8, height=0.06, col = fill_colors, theta=0.9, start=0, explode=0, labels=labs, labelcex=label_font_size, shade=0.7, sector.order=1:2, border=FALSE)
    title(main=sprintf("%s<%g & |%s|>%g %s: %s", significance_column, significance_cutoff, change_column, change_cutoff, filtering_mode, say_contrast), cex.main=4, line=-2)
    } else {
        labs = gsub("\n", ": ", labs)
        pie3D(N, radius=0.8, height=0.06, col = fill_colors, theta=0.9, start=45, explode=0, labels=labs, labelcex=label_font_size, shade=0.7, sector.order=1:2, border=NULL)
    title(main=sprintf("%s<%g & |%s|>%g %s: %s", significance_column, significance_cutoff, change_column, change_cutoff, filtering_mode, say_contrast), cex.main=4, line=-2)
    }
}

} 

return(out)
    
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in the code for any table.

#######################
## End of Template   ##
#######################

print("template_function_DEG_Gene_List.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_DEG_Analysis<-readRDS(paste0(rds_output,"/var_DEG_Analysis.rds"))
Input_is_Seurat_count <- 0
for(item in var_DEG_Analysis){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_DEG_Analysis<-as.data.frame(var_DEG_Analysis)}else{var_DEG_Analysis <- var_DEG_Analysis}
invisible(graphics.off())
var_DEG_Gene_List<-DEG_Gene_List(var_DEG_Analysis)
invisible(graphics.off())
saveRDS(var_DEG_Gene_List, paste0(rds_output,"/var_DEG_Gene_List.rds"))
