# Normalization [CCBR] (afc2524c-9bae-4873-98c1-e06ef8f4632b): v260
Normalized_Counts <- function(Filtered_Counts, ccbr991_metadata_nidap) {
    
    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(limma)
    library(tidyverse)
    library(edgeR)
    library(ggplot2)
    library(plotly)
    library(dplyr)
    library(RColorBrewer)
    library(colorspace)
    library(stringr)
    library(RCurl)
    library(reshape2)
    library(gridExtra)
    library(amap)
    library(lattice)
    library(gplots)
    library(gridGraphics)
    library(dendsort)

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##
    
    #Basic Parameters:
    counts_matrix <- Filtered_Counts
    sample_metadata <- ccbr991_metadata_nidap
    gene_names_column <- "Gene"
    columns_to_include = c("Gene","control1","control2","control3","IL7R1","IL7R2","IL7R3","NRAS1","NRAS2","NRAS3","Both1","Both2","Both3")
    sample_names_column <- "Sample"
    groups_column <- "Group"
    labels_column <- "Label"

    #Normalization Parameters:
    input_in_log_counts <- FALSE
    normalization_method <- "quantile"
    
    #PCA parameters:
    samples_to_rename_manually_on_pca <- c("")
    add_labels_to_pca <- TRUE
    principal_component_on_x_axis_for_pca <- 1
    principal_component_on_y_axis_for_pca <- 2
    legend_position_for_pca <- "top"
    label_offset_x_for_pca <- 2
    label_offset_y_for_pca <- 2
    label_font_size_for_pca <- 3
    point_size_for_pca <- 2

    #Histogram parameters:
    color_histogram_by_group <- FALSE
    set_min_max_for_x_axis_for_histogram <- FALSE 
    minimum_for_x_axis_for_histogram <- -1
    maximum_for_x_axis_for_histogram <- 1
    legend_font_size_for_histogram <- 10
    legend_position_for_histogram <- "top"
    number_of_histogram_legend_columns <- 6
    
    #Visualization Parameters:
    number_of_image_rows <- 2
    colors_for_plots <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
    make_plots_interactive <- FALSE
    plot_correlation_matrix_heatmap <- TRUE
    

    ##--------------- ##
    ## Error Messages ##
    ## -------------- ##

    
    ## --------- ##
    ## Functions ##
    ## --------- ##

    getourrandomcolors<-function(k){
        seed=10
        n <- 2e3
        ourColorSpace <- colorspace::RGB(runif(n), runif(n), runif(n))
        ourColorSpace <- as(ourColorSpace, "LAB")
        currentColorSpace <- ourColorSpace@coords
        # Set iter.max to 20 to avoid convergence warnings.
        set.seed(seed)
        km <- kmeans(currentColorSpace, k, iter.max=20)
        return( unname(hex(LAB(km$centers))))
    }

    make_heatmap <- function(counts_matrix, metadata,colorval) {
        mat <- as.matrix(counts_matrix) 
        tcounts=t(mat)
        tcounts=merge(metadata,tcounts,by.x=sample_names_column,by.y='row.names')
        rownames(tcounts)=tcounts[,labels_column]
        tcounts=tcounts[,!colnames(tcounts)%in%colnames(metadata)]
        d=Dist(tcounts,method="correlation",diag=TRUE)
        dend = rev(dendsort(as.dendrogram(hclust( d,method="average"))))
        m=as.matrix(d)
        sample_metadata <- metadata
        rownames(sample_metadata) = sample_metadata[[labels_column]]
        idx = as.factor(sample_metadata[rownames(m),groups_column])
        col = colorval
        cols <- col[idx]
        new.palette=colorRampPalette(c("blue","green","yellow"),space="rgb")
  
    mk<-function(){
        if(length(colnames(m))>20){
            par(mar=c(0,0,0,0))
            heatmap.2(m,
                        labRow = NA, 
                        labCol = NA,
                        col=new.palette(20),
                        trace="none",
                        colRow = col[idx], 
                        colCol = col[idx],
                        rowDendrogram=dend,
                        colDendrogram=dend,
                        RowSideColors = col[idx],
                        ColSideColors = col[idx],
                        dendrogram = "row",
                        cexRow=3,
                        cexCol=3,
                        margins=c(0,0),   
                        lmat=rbind( c(0,0,2),c(4,1,3) ,c(0,5,6) ), 
                        lhei=c(.2,4,2), 
                        lwid=c(1, .2,4 ), 
                        key.par=list(mgp=c(1.75, .5, 0), 
                        mar=c(7, 2, 3.5, 0), 
                        cex.axis=.1, 
                        cex.lab=3, 
                        cex.main=1, 
                        cex.sub=1),
                        key.xlab = "Correlation",
                        key.ylab="Count",
                        key.title=" ")       
        } else {
            heatmap.2(m,col=new.palette(20),
                        trace="none",
                        colRow = col[idx], 
                        colCol = col[idx],
                        rowDendrogram=dend,
                        colDendrogram=dend,
                        RowSideColors = col[idx],
                        ColSideColors = col[idx],
                        dendrogram = "row",
                        cexRow=3,cexCol=3,margins=c(4,1),  
                        lmat=rbind( c(0,0,2),c(4,1,3) ,c(0,5,6) ), 
                        lhei=c( .2,4,2), 
                        lwid=c(1, .2,4),
                        key.par=list(mgp=c(1.75, .5, 0), mar=c(7, 2, 3.5, 0), cex.axis=.1, cex.lab=3, cex.main=1, cex.sub=1),
                        key.xlab = "Correlation",
                        key.ylab="Count",
                        key.title=" ")
            }
        }
  
            tg<-mk()
            grid.echo(mk)
            gh1<-grid.grab()
            mklegend<-function(){
            plot.new()
            legend(x="top", legend=levels(idx), col=col[as.factor(levels(idx))],pch=15,x.intersp=3,bty ="n",cex=2)
            }
        grid.echo(mklegend )
        gh2<-grid.grab()
        lay <- c(1,3)
        grid.newpage()
        grid.arrange(gh1,gh2,nrow=1,widths=c(unit(1000, "bigpts"),unit(300, "bigpts")))
        gh<-grid.grab()
        return(gh)
    }

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##
samples_to_include=columns_to_include[columns_to_include%in%sample_metadata[,sample_names_column,drop=T]]
    anno_col=columns_to_include[columns_to_include%in%sample_metadata[,sample_names_column,drop=T]==F]

    samples_to_include <- samples_to_include[samples_to_include != gene_names_column]
    samples_to_include <- samples_to_include[samples_to_include != "Gene"]
    samples_to_include <- samples_to_include[samples_to_include != "GeneName"]
    
##create unique rownames to correctly add back Annocolumns at end of template
counts_matrix[,gene_names_column]=paste0(counts_matrix[,gene_names_column],'_',1:nrow(counts_matrix))

    anno_col=c(anno_col,gene_names_column)%>%unique
    anno_tbl=counts_matrix[,anno_col,drop=F]%>%as.data.frame
 
    
    df.filt <- counts_matrix[,samples_to_include]
    gene_names <- NULL
    gene_names$GeneID <- counts_matrix[,1]

    
    
    sample_metadata <- sample_metadata[match(colnames(df.filt),sample_metadata[[sample_names_column]]),] #First match sample metadata to counts matrix
    sample_metadata <- sample_metadata[rowSums(is.na(sample_metadata)) != ncol(sample_metadata), ] # Remove empty rows
    sample_metadata <- sample_metadata[, colSums(is.na(sample_metadata)) == 0] #Remove empty columns
    rownames(sample_metadata) <- sample_metadata[[sample_names_column]]

    df.filt <- df.filt[,match(sample_metadata[[sample_names_column]],colnames(df.filt))] #Match counts matrix columns to sample metadata
    
    #If input is in log space, linearize
    if(input_in_log_counts == TRUE){
        x <- DGEList(counts=2^df.filt, genes=gene_names)
    } else {
        x <- DGEList(counts=df.filt, genes=gene_names)     
    }

    v <- voom(x,normalize=normalization_method)
    rownames(v$E) <- v$genes$GeneID
    as.data.frame(v$E) %>% rownames_to_column(gene_names_column) -> df.voom
    print(paste0("Total number of features included: ", nrow(df.voom)))

    #Start PCA Plot:
    
    edf <- v$E
    tedf <- t(edf)
    tedf <- tedf[, colSums(is.na(tedf)) != nrow(tedf)]
    tedf <- tedf[, apply(tedf, 2, var) != 0]
    pca <- prcomp(tedf, scale.=T)
    
    pcx <- paste0("PC",principal_component_on_x_axis_for_pca)
    pcy <- paste0("PC",principal_component_on_y_axis_for_pca)
    pca.df <- as.data.frame(pca$x) %>% dplyr::select(.data[[pcx]], .data[[pcy]])
    pca.df$group <- sample_metadata[[groups_column]]
    pca.df$sample <- sample_metadata[[labels_column]]
    perc.var <- (pca$sdev^2/sum(pca$sdev^2))*100
    perc.var <- formatC(perc.var,format = "g",digits=4)
    pc.x.lab <- paste0(pcx," ", perc.var[principal_component_on_x_axis_for_pca],"%")
    pc.y.lab <- paste0(pcy," ", perc.var[principal_component_on_y_axis_for_pca],"%")
    labelpos <- pca.df
    labelpos$mean_y <- pca.df[[pcy]]+label_offset_y_for_pca
    labelpos$mean_x <- pca.df[[pcx]]+label_offset_x_for_pca
    pca.df$xdata <- pca.df[[pcx]]
    pca.df$ydata <- pca.df[[pcy]]

    # Manual changes to sample names
    replacements = samples_to_rename_manually_on_pca

    if (!is.null(replacements)) {
        if (replacements != c("")) {
            for (x in replacements) {
                old <- strsplit(x, ": ?")[[1]][1]
                new <- strsplit(x, ": ?")[[1]][2]
                pca.df$sample <- ifelse(pca.df$sample==old, new, pca.df$sample)
            }
        }
    }

    colorlist <- c("#5954d6","#e1562c","#b80058","#00c6f8","#d163e6","#00a76c","#ff9287","#008cf9","#006e00","#796880","#FFA500","#878500")
    names(colorlist) <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
    if(length(colors_for_plots) == 0){
        colors_for_plots <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
    }
    colorval <- colorlist[colors_for_plots]
    colorval <- unname(colorval) #remove names which affect ggplot

    if (length(unique(sample_metadata[[groups_column]])) > length(colorval)) {
        ## Original color-picking code.
        k=length(unique(sample_metadata[[groups_column]]))-length(colorval)
        more_cols<- getourrandomcolors(k) 
        colorval <- c(colorval , more_cols)
    }

    if (add_labels_to_pca == TRUE){
    g <- ggplot(pca.df, aes(x=xdata, y=ydata)) +
        theme_bw() +
        theme(legend.title=element_blank()) +
        theme(legend.position=legend_position_for_pca) +
        geom_point(aes(color=group), size=point_size_for_pca) +
        geom_text(data=labelpos, aes(x=labelpos$mean_x, y=labelpos$mean_y, 
            label=sample, color=group, vjust="inward", hjust="inward"), size=label_font_size_for_pca, show.legend=FALSE) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
        scale_colour_manual(values = colorval) +
        xlab(pc.x.lab) + ylab(pc.y.lab)
    } else {
    g <- ggplot(pca.df, aes(x=xdata, y=ydata)) +
        theme_bw() +
        theme(legend.title=element_blank()) +
        theme(legend.position=legend_position_for_pca) +
        geom_point(aes(color=group,text=sample), size=point_size_for_pca) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
        scale_colour_manual(values = colorval) +
        xlab(pc.x.lab) + ylab(pc.y.lab)    
    }

    par(mfrow = c(2,1))

    #Histogram Plot:
    
    df.m <- melt(edf,id.vars=c(gene_names_column))
    df.m = dplyr::rename(df.m,sample=Var2)

    if(set_min_max_for_x_axis_for_histogram == TRUE){
        xmin = minimum_for_x_axis_for_histogram
        xmax = maximum_for_x_axis_for_histogram
    } else {
        xmin = min(df.m$value)
        xmax = max(df.m$value)
    }

    if(color_histogram_by_group == TRUE){
        df.m %>% mutate(colgroup = sample_metadata[sample,groups_column]) -> df.m
        df.m = df.m[complete.cases(df.m[, "colgroup"]),]
        df.m$colgroup = gsub("\\s","_",df.m$colgroup)
        df.m$colgroup = factor(df.m$colgroup, levels=unique(df.m$colgroup))
        #print(unique(df.m$sample))

        # plot Density 
        g2 = ggplot(df.m, aes(x=value, group=sample)) + 
            geom_density(aes(colour = colgroup)) +
            xlab("Filtered Counts") + ylab("Density") +
            theme_bw() +
            theme(legend.position=legend_position_for_histogram,legend.text = element_text(size = legend_font_size_for_histogram)) + 
            ggtitle("Frequency Histogram") +
            xlim(xmin,xmax) +
            #scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),40)) +
            scale_colour_manual(values=colorval)
    } else {
        
        df.m$sample = sample_metadata[df.m$sample,labels_column]
        n=length(unique(df.m$sample))
        cols<- getourrandomcolors(n) 
        
        g2 = ggplot(df.m, aes(x=value, group=sample)) + 
            geom_density(aes(colour = sample)) +
            xlab("Filtered Counts") + ylab("Density") +
            theme_bw() +
            theme(legend.position=legend_position_for_histogram,legend.text = element_text(size = legend_font_size_for_histogram)) +  
            ggtitle("Frequency Histogram") +
            xlim(xmin,xmax) +
            #scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),n)) +
            scale_colour_manual(values=cols)#+
            guides(linetype = guide_legend(ncol = number_of_histogram_legend_columns))
    }

    # dev.off()

    imageWidth = 3000
    imageHeight = 1500*2
    dpi = 300

    png(
      filename="Normalized_Counts.png",
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")

    if(plot_correlation_matrix_heatmap == TRUE){
        if(make_plots_interactive == TRUE){
            p1=(g)%>%ggplotly(tooltip = c("sample","group"))
            p2=(g2+theme(legend.position = "none")) %>%ggplotly(tooltip = c("sample"))
            fig=subplot(p1,p2,which_layout = 'merge',margin=.05,shareX = F,shareY = F,titleY = T,titleX = T,widths=c(.5,.5),nrows = 1)
            fig=fig %>% layout(title = 'Interactive PCA and Histogram')
            print(fig)
        } else {
            require(gridExtra)
            gh<-make_heatmap(df.filt,sample_metadata,colorval)
            grid.arrange(g,g2,gh, nrow=number_of_image_rows)
            # dev.off()
        }  
    } else {
        if(make_plots_interactive == TRUE){
            p1=(g)%>%ggplotly(tooltip = c("sample","group"))
            p2=(g2+theme(legend.position = "none")) %>%ggplotly(tooltip = "sample" )
            fig=subplot(p1,p2,which_layout = 'merge',margin=.05,shareX = F,shareY = F,titleY = T,titleX = T,widths=c(.5,.5),nrows = 1)
            fig=fig %>% layout(title = 'Interactive PCA and Histogram')
            print(fig)
        } else {
            grid.arrange(g,g2, nrow=number_of_image_rows)
            # dev.off()
        }
        }    

print("Sample columns")
print(colnames(df.voom)[!colnames(df.voom)%in%gene_names_column])
print("Feature Columns")
print(colnames(anno_tbl))

    df.voom=merge(anno_tbl,df.voom,by=gene_names_column,all.y=T)
   df.voom[,gene_names_column]=gsub('_[0-9]+$',"",df.voom[,gene_names_column])

    return(df.voom)
}

#################################################
## Global imports and functions included below ##
#################################################

print("template_function_Normalized_Counts.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Filtered_Counts<-readRDS(paste0(rds_output,"/var_Filtered_Counts.rds"))
Input_is_Seurat_count <- 0
for(item in var_Filtered_Counts){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Filtered_Counts<-as.data.frame(var_Filtered_Counts)}else{var_Filtered_Counts <- var_Filtered_Counts}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_ccbr991_metadata_nidap<-readRDS(paste0(rds_output,"/var_ccbr991_metadata_nidap.rds"))
Input_is_Seurat_count <- 0
for(item in var_ccbr991_metadata_nidap){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_ccbr991_metadata_nidap<-as.data.frame(var_ccbr991_metadata_nidap)}else{var_ccbr991_metadata_nidap <- var_ccbr991_metadata_nidap}
invisible(graphics.off())
var_Normalized_Counts<-Normalized_Counts(var_Filtered_Counts,var_ccbr991_metadata_nidap)
invisible(graphics.off())
saveRDS(var_Normalized_Counts, paste0(rds_output,"/var_Normalized_Counts.rds"))
