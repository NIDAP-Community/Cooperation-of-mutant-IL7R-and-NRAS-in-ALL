# Volcano Plot - Summary [CCBR] (eeced39d-ed52-4b16-9847-4282971775e6): v387
Volcano_Summary_MS_Contrasts_Only <- function(DEG_Analysis_MS_Contrasts_Only) {
    # image: png
    
   stattype<-"pval"
   add_deg_columns<-c("FC","logFC","tstat","pval","adjpval")
   image_width = 15
   image_height = 15
   image_resolution = 300
   aspect_ratio = 0

    if (FALSE & (("png" =='svg')))  {
        library(svglite)
        svglite::svglite(
        file="Volcano_Summary_MS_Contrasts_Only.png",
        width=image_width,
        height=image_height,
        pointsize=1,
        bg="white",
    )} else if ("png" == 'png') {
      png(
      filename="Volcano_Summary_MS_Contrasts_Only.png",
      width=image_width,
      height=image_height,
      units="in",
      pointsize=4,
      bg="white",
      res=image_resolution,
      type="cairo")
    }

    suppressMessages(library(ggplot2))
    suppressMessages(library(dplyr))
    suppressMessages(library(ggrepel))

    genesmat <- DEG_Analysis_MS_Contrasts_Only
    value_to_sort_the_output_dataset = "t-statistic"
    
    volcols<-colnames(genesmat)
     print(volcols)
    statcols<-volcols[grepl("logFC",volcols)]
    contrasts<-unique(gsub("_logFC","",statcols))   
   
    Plots <- list()
    df_outs <- list()
    for(contrast in contrasts){
        print(paste0("Doing contrast: ",contrast))
        lfccol=paste0(contrast,"_logFC")
        pvalcol=paste0(contrast,"_",stattype)
        tstatcol=paste0(contrast,"_","tstat")

        print(paste0("Fold change column: ",lfccol))
        print(paste0(stattype," column: ",pvalcol))

        no_genes_to_label <- 30
    if (value_to_sort_the_output_dataset=="fold-change") {
        genesmat %>% dplyr::arrange(desc(abs(genesmat[,lfccol]))) -> genesmat
    } else if (value_to_sort_the_output_dataset=="p-value") {
        genesmat %>% dplyr::arrange(genesmat[,pvalcol]) -> genesmat
    } else if (value_to_sort_the_output_dataset == "t-statistic") {
        genesmat %>% dplyr::arrange(desc(abs(genesmat[,tstatcol]))) -> genesmat
    }
    print(paste0("Total number of genes included in volcano plot: ", nrow(genesmat)))
    if (TRUE){
        negative_log10_p_values <- -log10(genesmat[,pvalcol])
        ymax <- ceiling(max(negative_log10_p_values[is.finite(negative_log10_p_values)]))
    } else {
        ymax = 10
    }
    if (TRUE){
        xmax1 = ceiling(max(genesmat[,lfccol]))
        xmax2 = ceiling(max(-genesmat[,lfccol]))
        xmax=max(xmax1,xmax2)
    } else {
        xmax = 5
    }
   

    ## work with a list of genes
if (FALSE){
    gl <- trimws(unlist(strsplit(c("Provide list of genes-comma separated"), ",")), which=c("both"))
        ind <- match(gl, genesmat$Gene) # get the indices of the listed genes
        gene_list_ind <- c(1:no_genes_to_label,ind) # when list provided
        color_gene_label <- c(rep(c("black"), no_genes_to_label), rep(c("green3"),length(ind)))
   }else if (FALSE){
        gl <- trimws(unlist(strsplit(c("Provide list of genes-comma separated"), ",")), which=c("both")) # unpack the gene list provided by the user and remove white spaces
        ind <- match(gl, genesmat$Gene) # get the indices of the listed genes
        gene_list_ind <- ind # when list provided
        color_gene_label <- rep(c("green3"), length(ind))
    } else {
        if (no_genes_to_label>0) {
        gene_list_ind <- 1:no_genes_to_label # if no list provided label the number of genes given by the user
        color_gene_label <- rep(c("black"), no_genes_to_label)
        } else if (no_genes_to_label ==0) {
        gene_list_ind <-0
        }
        }   

## special nudge/repel of specific genes
if (FALSE){
    gn <- trimws(unlist(strsplit(c("Provide list of genes-comma separated"), ",")), which=c("both"))
    ind_gn <- match(gn, genesmat$Gene[gene_list_ind]) # get the indices of the listed genes
    nudge_x_all <- rep(c(0.2), length(genesmat$Gene[gene_list_ind]))
    nudge_y_all <- rep(c(0.2), length(genesmat$Gene[gene_list_ind]))
    nudge_x_all[ind_gn] <- c(2)
    nudge_y_all[ind_gn] <- c(2)
} else {
    nudge_x_all <- 0.2
    nudge_y_all <- 0.2
}
   
    ## flip contrast section
    flipVplot <- FALSE
        indc <- which(colnames(genesmat) == lfccol) # get the indice of the column that contains the contrast_logFC data

        if (length(indc)==0){
            print("Please rename the logFC column to include the contrast evaluated.")
        } else{
        old_contrast <- colnames(genesmat)[indc]
        }  
    # actually flip contrast
    if (flipVplot){
        # get the indice of the contrast to flip
        indcc <- match(old_contrast,colnames(genesmat)) 
        # create flipped contrast label
        splt1 <- strsplit(old_contrast, "_") # split by underline symbol to isolate the contrast name
        splt2 <- strsplit(splt1[[1]][1],"-") # split the contrast name in the respective components
        flipped_contrast <- paste(splt2[[1]][2], splt2[[1]][1],sep="-") #flip contrast name
        new_contrast_label <- paste(flipped_contrast, c("logFC"), sep = "_") 
        # rename contrast column to the flipped contrast
        colnames(genesmat)[indcc] <- new_contrast_label
        # flip the contrast data around y-axis
        genesmat[,indcc] <- -genesmat[indcc]
    } else{ new_contrast_label<- old_contrast}

    grm<-genesmat[,c(new_contrast_label,pvalcol)]
    grm[,"neglogpval"]<- -log10(genesmat[,pvalcol])
    colnames(grm)=c("FC","pval","neglogpval")
    print(grm[gene_list_ind,])
    p <- ggplot(grm,
        aes_string(x = "FC", y ="neglogpval" ))+ # modified by RAS
        theme_classic() +
        geom_point(
            color='black',
            size = 2) +
        geom_vline(xintercept=c(-1,1), color='red', alpha=1.0) + 
        geom_hline(yintercept=-log10(0.001), color='blue', alpha=1.0) +  
        geom_point(
            data = grm[genesmat[,pvalcol] < 0.001,],
            color = 'lightgoldenrod2',
            size = 2) +
        geom_point(
            data = grm[genesmat[,pvalcol] < 0.001 & abs(grm[,"FC"])>1,], 
            color = 'red',
            size = 2) +
        geom_text_repel(
            data = grm[gene_list_ind,], 
           label = genesmat$Gene[gene_list_ind], 
            color = color_gene_label,
            fontface = 1,
            nudge_x = nudge_x_all,
            nudge_y = nudge_y_all,
            size = 4,
            segment.size = 0.5) +
        xlim(-xmax,xmax) +
        ylim(0,ymax) + xlab(new_contrast_label) + ylab(pvalcol)

        if (aspect_ratio > 0){
            p <- p + coord_fixed(ratio=aspect_ratio)
        }
        
        #print(p)
        Plots[[contrast]]=p
        print(head(genesmat[,pvalcol] ))
         filtered_genes =   genesmat$Gene[genesmat[,pvalcol] < 0.001 & abs(grm[,"FC"])>1]
        #print(filtered_genes)
        repeated_column = rep(contrast, length(filtered_genes))
        ## If param empty upon template upgrade, fill it with default value.
        if (length(add_deg_columns) == 0) {
            add_deg_columns <- c("FC","logFC","tstat","pval","adjpval")
        }
        ## Get columns for output table.
        if (add_deg_columns == "none") {
            new_df <- data.frame(filtered_genes, repeated_column)
            names(new_df) <- c("Gene", "Contrast")
        } else {
            add_deg_columns = setdiff(add_deg_columns, "none")
            out_columns = paste(contrast, add_deg_columns, sep="_")
            deg = genesmat[,c("Gene", out_columns)]
            names(deg)[1] = "Gene"
            new_df <- data.frame(filtered_genes, repeated_column) %>% dplyr::left_join(deg, by=c("filtered_genes"="Gene"))
            names(new_df) <- c("Gene", "Contrast", add_deg_columns)
        }

        df_out1 <- new_df
        df_outs[[contrast]]=df_out1
   }

   Use_default_grid_layout = TRUE
   require(gridExtra)
   nplots=length(Plots)
   if (Use_default_grid_layout) {
    nrows=ceiling(nplots/ceiling(sqrt(nplots)))
   } else {
    nrows = 1
   }

 do.call("grid.arrange", c(Plots, nrow=nrows))
   print("done plotting")

    df_out <- unique(do.call("rbind", df_outs))
     print(head(df_out))
    print(colnames(df_out))
    return(df_out)
}

 

print("template_function_Volcano_Summary_MS_Contrasts_Only.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_DEG_Analysis_MS_Contrasts_Only<-readRDS(paste0(rds_output,"/var_DEG_Analysis_MS_Contrasts_Only.rds"))
Input_is_Seurat_count <- 0
for(item in var_DEG_Analysis_MS_Contrasts_Only){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_DEG_Analysis_MS_Contrasts_Only<-as.data.frame(var_DEG_Analysis_MS_Contrasts_Only)}else{var_DEG_Analysis_MS_Contrasts_Only <- var_DEG_Analysis_MS_Contrasts_Only}
invisible(graphics.off())
var_Volcano_Summary_MS_Contrasts_Only<-Volcano_Summary_MS_Contrasts_Only(var_DEG_Analysis_MS_Contrasts_Only)
invisible(graphics.off())
saveRDS(var_Volcano_Summary_MS_Contrasts_Only, paste0(rds_output,"/var_Volcano_Summary_MS_Contrasts_Only.rds"))
