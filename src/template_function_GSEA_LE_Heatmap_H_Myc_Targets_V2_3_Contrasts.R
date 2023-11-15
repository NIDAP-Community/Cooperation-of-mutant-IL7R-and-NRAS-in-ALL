# Expression Heatmap [CCBR] (89a32987-10d9-4233-91c3-e9adf3dcc517): v548
GSEA_LE_Heatmap_H_Myc_Targets_V2_3_Contrasts <- function(Batch_Corrected_Counts,ccbr991_metadata_nidap) {
    ## This function uses pheatmap to draw a heatmap, scaling first by rows
    ## (with samples in columns and genes in rows)

    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(colorspace)
    library(dendsort)
    library(ComplexHeatmap)
    library(dendextend)
    library(tibble)
    library(stringr)
    library(RColorBrewer)
    library(dplyr)
    library(grid)
    library(gtable)
    library(gridExtra)
    library(gridGraphics)
    library(circlize)

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##

    #Basic Parameters:
    counts_matrix <- Batch_Corrected_Counts
    sample_metadata <- ccbr991_metadata_nidap
    gene_column_name <- "Gene"
    group_columns <- c("Group")
    sample_name_column <- "Sample"
    samples_to_include = c("control1","control2","control3","IL7R1","IL7R2","IL7R3","NRAS1","NRAS2","NRAS3","Both1","Both2","Both3")
    include_all_genes <- FALSE
    filter_top_genes_by_variance = FALSE
    top_genes_by_variance_to_include <- 500
    specific_genes_to_include_in_heatmap = "Mcm4 Noc4l Ung Nolc1 Rrp9 Tcof1 Pprc1 Rrp12 Utp20 Plk1 Mcm5 Srm Myc Prmt3 Gnl3 Pa2g4 Nip7 Ddx18 Bysl Nop56 Mrto4 Mybbp1a Las1l Wdr43 Hk2 Grwd1 Ipo4 Ndufaf4 Tmem97 Tbrg4 Farsa Wdr74 Slc29a2 Mphosph10 Rcl1 Phb Aimp2 Cdk4 Dctpp1 Dusp2 Exosc5 Hspe1 Imp4 Map3k6 Pes1 Plk4 Ppan Pus1 Rabepk Slc19a1 Sord Supv3l1 Tfb2m" 
    
    #Visual Parameters:
    heatmap_color_scheme <- "Default"
    autoscale_heatmap_color <- TRUE
    set_min_heatmap_color <- -2
    set_max_heatmap_color <- 2
    group_colors <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
    assign_group_colors <- FALSE
    assign_color_to_sample_groups <- c()
    legend_font_size <- 10 
    display_gene_names <- TRUE
    gene_name_font_size <- 8
    display_sample_names <- TRUE
    sample_name_font_size <- 8
    display_dendrograms <- TRUE
    reorder_dendrogram <- TRUE
    reorder_dendrogram_order <- c("control3","control2","control1","IL7R1","IL7R3","IL7R2","NRAS1","NRAS3","NRAS2","Both1","Both3","Both2")
    manually_rename_samples <- FALSE
    samples_to_rename <- c("")
    display_numbers <- FALSE
    aspect_ratio <- "Auto"

    #Advanced Parameters
    distance_metric <- "euclidean"
    clustering_method <- "average"
    center_and_rescale_expression <- TRUE
    cluster_genes <- TRUE
    cluster_samples <- TRUE
    arrange_sample_columns <- FALSE
    order_by_gene_expression <- FALSE
    gene_to_order_columns <- " "
    gene_expression_order <- "low_to_high"
    return_z_scores <- FALSE

    ##--------------- ##
    ## Error Messages ##
    ## -------------- ##

    if(include_all_genes == TRUE && filter_top_genes_by_variance == TRUE){
        stop("ERROR: Choose only one of 'Include all genes' or 'Filter top genes by variance' as TRUE")
    }

    if((cluster_samples == TRUE && arrange_sample_columns == TRUE) | (arrange_sample_columns == TRUE && order_by_gene_expression == TRUE) | 
    (arrange_sample_columns == TRUE && cluster_samples == TRUE) | (cluster_samples == FALSE && arrange_sample_columns == FALSE && order_by_gene_expression == FALSE)) {
     stop("ERROR: Choose only one of 'Cluster Samples', 'Arrange sample columns', or 'order by gene expression' as TRUE")   
    }

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

    ## Begin pal() color palette function∂:
    pal = function (n, h=c(237, 43), c=100, l=c(70, 90), power=1, fixup=TRUE, gamma=NULL, alpha=1, ...) {
        if (n < 1L) {
            return(character(0L))
        }
        h <- rep(h, length.out = 2L)
        c <- c[1L]
        l <- rep(l, length.out = 2L)
        power <- rep(power, length.out = 2L)
        rval <- seq(1, -1, length = n)
        rval <- hex(
            polarLUV(
                L = l[2L] - diff(l) * abs(rval)^power[2L], 
                C = c * abs(rval)^power[1L],
                H = ifelse(rval > 0, h[1L], h[2L])
            ),
            fixup=fixup, ...
        )
        if (!missing(alpha)) {
            alpha <- pmax(pmin(alpha, 1), 0)
            alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
                width = 2L, upper.case = TRUE)
            rval <- paste(rval, alpha, sep = "")
        }
        return(rval)
    } 
    # End pal() color palette function:

    ## Begin doheatmap() function:
    doheatmap <- function(dat, clus, clus2, ht, rn, cn, col, dispnum) {
        #require(pheatmap)
        #require(dendsort)
        col.pal <- np[[col]]
        if (FALSE) {
            col.pal = rev(col.pal)
        }
        # Define metrics for clustering
        drows1 <- distance_metric
        dcols1 <- distance_metric
        minx = min(dat)
        maxx = max(dat)
        if (autoscale_heatmap_color) {
            breaks = seq(minx, maxx, length=100)
            legbreaks = seq(minx, maxx, length=5)
        } else {
            breaks = seq(set_min_heatmap_color, set_max_heatmap_color, length=100)
            legbreaks = seq(set_min_heatmap_color, set_max_heatmap_color, length=5)
        }
        breaks = sapply(breaks, signif, 4)
        legbreaks = sapply(legbreaks, signif, 4)
        # Run cluster method using 
        if(distance_metric != "correlation"){
            hc = hclust(dist(t(dat),method=distance_metric), method=clustering_method)
            hcrow = hclust(dist(dat,method=distance_metric), method=clustering_method)
            if (FALSE) {
                sort_hclust <- function(...) as.hclust(rev(dendsort(as.dendrogram(...))))
            } else {
                sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
            }
            if (clus) {
                colclus <- sort_hclust(hc)
            } else {
                colclus = FALSE
            }
            if (clus2) {
                rowclus <- sort_hclust(hcrow)
            } else {
                rowclus = FALSE
            }
        } else {
            clus <- clus
            rowclus <- clus2
        }
        if (display_dendrograms) {
            treeheight <- 25
        } else {
            treeheight <- 0
        }

        hm.parameters <- list(
            dat, 
            color=col.pal,
            legend_breaks=legbreaks,
            legend=TRUE,
            scale="none",
            treeheight_col=treeheight,
            treeheight_row=treeheight,
            kmeans_k=NA,
            breaks=breaks,
            display_numbers=dispnum,
            number_color = "black",
            fontsize_number = 8,
            height=80,
            cellwidth = NA, 
            cellheight = NA, 
            fontsize= legend_font_size,   
            fontsize_row=gene_name_font_size,
            fontsize_col=sample_name_font_size,
            show_rownames=rn, 
            show_colnames=cn,
            cluster_rows=rowclus, 
            cluster_cols=clus,
            #clustering_distance_rows=drows1, 
            #clustering_distance_cols=dcols1,
            annotation_col = annotation_col,
            annotation_colors = annot_col,
            labels_col = labels_col
        )
        mat = t(dat)
        callback = function(hc, mat) {
            dend = rev(dendsort(as.dendrogram(hc)))
            if(reorder_dendrogram == TRUE) {
                dend %>% dendextend::rotate(reorder_dendrogram_order) -> dend
            } else {
                dend %>% dendextend::rotate(c(1:nobs(dend))) 
            }
            as.hclust(dend)
        }

        #plot(dend)
        #as.hclust(dend)

        ## Make Heatmap
        if(distance_metric == "correlation"){
            phm <- do.call("pheatmap", c(hm.parameters, list(clustering_distance_rows=drows1,clustering_distance_cols=dcols1)))    
        } else {
            phm <- do.call("pheatmap", c(hm.parameters, list(clustering_callback=callback)))
        }
        
    }
    # End doheatmap() function.

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##

    ## Build different color spectra options for heatmap:
    np0 = pal(100) 
    np1 = diverge_hcl(100, c=100, l=c(30, 80), power=1) # Blue to Red
    np2 = heat_hcl(100, c=c(80, 30), l=c(30, 90), power=c(1/5, 2)) # Red to Vanilla
    np3 = rev(heat_hcl(100, h=c(0, -100), c=c(40, 80), l=c(75, 40), power=1)) # Violet to Pink
    np4 = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100)) #Red to yellow to blue
    np5 = colorRampPalette(c("steelblue","white", "red"))(100)  # Steelblue to White to Red

    ## Gather list of color spectra and give them names for the GUI to show.
    np = list(np0, np1, np2, np3, np4, np5)
    names(np) = c("Default","Blue to Red","Red to Vanilla","Violet to Pink","Bu Yl Rd","Bu Wt Rd")

    ## Parse input counts matrix. Subset by samples.
    df1 <- counts_matrix
    # Swap out Gene Name column name, if it's not 'Gene'.
    if(gene_column_name != "Gene"){
        # Drop original Gene column
        df1 = df1[,!(colnames(df1)%in% c("Gene")) ]
        # Rename column to Gene
        colnames(df1)[which(colnames(df1) == gene_column_name)] <- 'Gene'
    }
    # Get sample columns
    samples_to_include <- samples_to_include[samples_to_include != gene_column_name]
    samples_to_include <- samples_to_include[samples_to_include != "Gene"]
    samples_to_include <- samples_to_include[samples_to_include != "GeneName"]

    # Build new counts matrix containing only sample subset chosen by user.
    df1 <- df1[,append("Gene", samples_to_include)]
    df.orig = df1
    df.orig %>% dplyr::group_by(Gene) %>% summarise_all(funs(mean)) -> df
    df.mat = df[ , (colnames(df) != "Gene" )] %>% as.data.frame
    df %>% dplyr::mutate(Gene = stringr::str_replace_all(Gene, "_", " ")) -> df
    df %>% dplyr::mutate(Gene = stringr::str_wrap(Gene,50)) -> df            ##### MC: added wrapper for looooonnggg names
    row.names(df.mat) <- df$Gene
    rownames(df.mat) <- str_wrap(rownames(df.mat),10)
    df.mat <- as.data.frame(df.mat)

    ## Subset counts matrix by genes.
    # Toggle to include all genes in counts matrix (in addition to any user-submitted gene list).
    if (include_all_genes == FALSE) {
        # Add user-submitted gene list (optional).
        genes_to_include_parsed = c()
        genes_to_include_parsed = strsplit(specific_genes_to_include_in_heatmap, " ")[[1]]
        df.mat[genes_to_include_parsed,] -> df.final.extra.genes
        if(filter_top_genes_by_variance == TRUE) {
            # Want to filter all genes by variance.
            df.final = as.matrix(df.mat)
            var <- matrixStats::rowVars(df.final)
            df <- as.data.frame(df.final)
            rownames(df) <- rownames(df.final)
            df.final <- df
            df.final$var <- var
            df.final %>% rownames_to_column("Gene") -> df.final 
            df.final %>% dplyr::arrange(desc(var)) -> df.final
            df.final.extra.genes = dplyr::filter(df.final, Gene %in% genes_to_include_parsed)
            df.final = df.final[1:top_genes_by_variance_to_include,]
            df.final = df.final[complete.cases(df.final),]
            # Rbind user gene list to variance-filtered gene list and deduplicate.
            df.final <- rbind(df.final, df.final.extra.genes)
            df.final <- df.final[!duplicated(df.final),] 
            rownames(df.final) <- df.final$Gene
            df.final$Gene <- NULL
            df.final$var <- NULL
        } else {
            # Want to use ONLY user-provided gene list.
            df.final <- df.final.extra.genes
            df.final <- df.final[!duplicated(df.final),]
            # Order genes in heatmap by user-submitted order of gene names.
            df.final <- df.final[genes_to_include_parsed,]
            #df.final$Gene <- NULL
        }
    } else {
        df.final <- df.mat
        df.final$Gene <- NULL
    }
    
        ## Optionally apply centering and rescaling (default TRUE).
    if (center_and_rescale_expression == TRUE) {
            tmean.scale = t(scale(t(df.final)))
            tmean.scale = tmean.scale[!is.infinite(rowSums(tmean.scale)),]
            tmean.scale = na.omit(tmean.scale)
    } else {
            tmean.scale = df.final
    }

    if(order_by_gene_expression == TRUE){
        gene_to_order_columns <- gsub(" ","",gene_to_order_columns)
        if(gene_expression_order == "low_to_high"){
        tmean.scale <- tmean.scale[,order(tmean.scale[gene_to_order_columns,])] #order from low to high 
        } else{
        tmean.scale <- tmean.scale[,order(-tmean.scale[gene_to_order_columns,])] #order from high to low  
        }
    }

    df.final <- as.data.frame(tmean.scale)

    ## Parse input sample metadata and add annotation tracks to top of heatmap.
    annot <- sample_metadata
    # Filter to only samples user requests.
    annot %>% dplyr::filter(.data[[sample_name_column]] %in% samples_to_include) -> annot
    annot %>% arrange(match(.data[[sample_name_column]], colnames(df.final))) -> annot
    # Arrange sample options.
    if(arrange_sample_columns) {
      annot %>% dplyr::arrange_(.dots=group_columns) -> annot
      df.final <- df.final[,match(annot[[sample_name_column]],colnames(df.final))] 
    }
    # Build subsetted sample metadata table to use for figure.

    colorlist <- c("#5954d6","#e1562c","#b80058","#00c6f8","#d163e6","#00a76c","#ff9287","#008cf9","#006e00","#796880","#FFA500","#878500")
    names(colorlist) <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
    group_colors <- colorlist[group_colors]

    annot %>% dplyr::select(group_columns) -> annotation_col    
    annotation_col = as.data.frame(unclass(annotation_col))
    annotation_col[] <- lapply(annotation_col,factor)
    x <- length(unlist(lapply(annotation_col,levels)))
    if(x>length(group_colors)){
        k=x-length(group_colors)
        more_cols<- getourrandomcolors(k) 
        group_colors <- c(group_colors, more_cols)
    }
    rownames(annotation_col) <- annot[[sample_name_column]]
    annot_col = list()
    b=1
    i=1
    while (i <= length(group_columns)){
        nam <- group_columns[i]
        grp <- as.factor(annotation_col[,i])
        c <- b+length(levels(grp))-1
        col = group_colors[b:c]
        names(col) <- levels(grp)
        assign(nam,col)
        annot_col = append(annot_col,mget(nam))
        b = c+1
        i=i+1
    }

    if(assign_group_colors == TRUE){
            colassign <- assign_color_to_sample_groups
            groupname <- c()
            groupcol <- c() 
            for (i in 1:length(colassign)) {
                groupname[i] <- strsplit(colassign[i], ": ?")[[1]][1]
                groupcol[i] <- strsplit(colassign[i], ": ?")[[1]][2]
            }
            annot_col[[1]][groupname] <- groupcol
    }

    ## Setting labels_col for pheatmap column labels.
    if (manually_rename_samples == TRUE) {
        # Use user-provided names to rename samples.
        replacements = samples_to_rename
        old <- c()
        new <- c()
        labels_col <- colnames(df.final)
        for (i in 1:length(replacements)) {
            old <- strsplit(replacements[i], ": ?")[[1]][1]
            new <- strsplit(replacements[i], ": ?")[[1]][2]
            old=gsub("^[[:space:]]+|[[:space:]]+$","",old)
            new=gsub("^[[:space:]]+|[[:space:]]+$","",new)
            labels_col[labels_col==old]=new           
        }
    } else {
        ## Use original column names for samples.
        labels_col <- colnames(df.final)
    }

    ## Print number of genes to log.
    print(paste0("The total number of genes in heatmap: ", nrow(df.final)))

    ## Make the final heatmap.
    p <- doheatmap(dat=df.final, clus=cluster_samples, clus2=cluster_genes, ht=50, rn=display_gene_names, cn=display_sample_names, col=heatmap_color_scheme, dispnum=display_numbers)
    p@matrix_color_mapping@name <- " "
    p@matrix_legend_param$at <- as.numeric(formatC(p@matrix_legend_param$at, 2))
    p@column_title_param$gp$fontsize <- 10
    print(p)

    ## If user sets toggle to TRUE, return Z-scores.
    ## Else return input counts matrix by default (toggle FALSE).
    ## Returned matrix includes only genes & samples used in heatmap.
    if(return_z_scores){
        df.new <- data.frame(tmean.scale) # Convert to Z-scores.
        df.new %>% rownames_to_column("Gene") -> df.new
        return(df.new)
    } else {
        df.final %>% rownames_to_column("Gene") -> df.new
        return(df.new)
    }
}

## --------------- ##
## End of Template ##
## --------------- ##

###Global code ####

pheatmap = function(mat, 
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), 
    kmeans_k = NA, 
    breaks = NA, 
    border_color = ifelse(nrow(mat) < 100 & ncol(mat) < 100, "grey60", NA),
    cellwidth = NA, 
    cellheight = NA, 
    scale = "none", 
    cluster_rows = TRUE,
    cluster_cols = TRUE, 
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean", 
    clustering_method = "complete",
    clustering_callback = NA, 
    cutree_rows = NA, 
    cutree_cols = NA,
    treeheight_row = ifelse(class(cluster_rows) == "hclust" || cluster_rows, 50, 0), 
    treeheight_col = ifelse(class(cluster_cols) == "hclust" || cluster_cols, 50, 0), 
    legend = TRUE, 
    legend_breaks = NA,
    legend_labels = NA, 
    annotation_row = NA, 
    annotation_col = NA,
    annotation = NA, 
    annotation_colors = NA, 
    annotation_legend = TRUE,
    annotation_names_row = TRUE, 
    annotation_names_col = TRUE,
    drop_levels = TRUE, 
    show_rownames = TRUE, 
    show_colnames = TRUE, 
    main = NA,
    fontsize = 10, 
    fontsize_row = fontsize, 
    fontsize_col = fontsize,
    angle_col = c("270", "0", "45", "90", "315"), 
    display_numbers = FALSE,
    number_format = "%.2f", 
    number_color = "grey30", 
    fontsize_number = 0.8 * fontsize, 
    gaps_row = NULL, 
    gaps_col = NULL, 
    labels_row = NULL,
    labels_col = NULL, 
    filename = NA, 
    width = NA, 
    height = NA,
    silent = FALSE, 
    na_col = "#DDDDDD", 
    name = NULL,

    # other graphic parameters for fonts
    fontfamily = "",
    fontfamily_row = fontfamily,
    fontfamily_col = fontfamily,
    fontface = 1,
    fontface_row = fontface,
    fontface_col = fontface,

    # argument specific for Heatmap()
    heatmap_legend_param = list(),
    ...,
    run_draw = FALSE
) {
  
   # if(is.data.frame(mat)) {
   #     warning("The input is a data frame, convert it to the matrix.")
        mat = as.matrix(mat)
   # }

    if(!identical(kmeans_k, NA)) {
        warning("argument `kmeans_k` is not suggested to use in pheatmap -> Heatmap translation because it changes the input matrix. You might check `row_km` and `column_km` arguments in Heatmap().")
        km = kmeans(mat, centers = kmeans_k)
        mat = km$centers
        rownames(mat) = paste0("Cluster: ", seq_along(km$size), ", Size: ", km$size)
    }

    if("row" %in% scale) {

        if(any(is.na(mat))) {
            mat = (mat - rowMeans(mat, na.rm = TRUE))/rowSds(mat, na.rm = TRUE)
        } else {
            mat = t(scale(t(mat)))
        }
    } else if("column" %in% scale) {
        if(any(is.na(mat))) {
            mat = t((t(mat) - colMeans(mat, na.rm = TRUE))/colSds(mat, na.rm = TRUE))
        } else {
            mat = scale(mat)
        }
    }

    ht_param = list(matrix = mat)

    if(!identical(scale, "none") && !identical(breaks, NA)) {
        warning("It not suggested to both set `scale` and `breaks`. It makes the function confused.")
    }

    # if color is a color mapping function
    if(is.function(color)) {
        ht_param$col = color
        if(!identical(breaks, NA)) {
            warning("`breaks` is ignored when `color` is set as a color mapping function.")
        }
    } else {
        if(identical(breaks, NA)) {
            n_col = length(color)
            if(identical(scale, "row") || identical(scale, "column")) {
                lim = max(abs(mat), na.rm = TRUE)
                ht_param$col = colorRamp2(seq(-lim, lim, length.out = n_col), color)
            } else {
                ht_param$col = colorRamp2(seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out = n_col), color)
            }
        } else  {
            if(length(breaks) == length(color) + 1) {
                ht_param$col = local({
                    breaks = breaks
                    color = color
                    fun = function(x) {
                        n = length(color)
                        df = data.frame(start = c(-Inf, breaks[seq_len(n)], breaks[n+1]), 
                                        end = c(breaks[1], breaks[1+seq_len(n)], Inf))
                        # tell which interval x is in
                        ind = numeric(length(x))
                        for(i in seq_along(x)) {
                            ind[i] = which(df$start <= x[i] & df$end > x[i])
                        }
                        ind = ind - 1
                        ind[ind < 1] = 1
                        ind[ind > n] = n
                        color[ind]
                    }
                    attr(fun, "breaks") = breaks
                    fun
                })
            } else if(length(breaks) == length(color)) {
                ht_param$col = colorRamp2(breaks, color)
            } else {
                n_col = length(color)
                ht_param$col  = colorRamp2(seq(min(breaks), max(breaks), length.out = n_col), color)
                warning("`breaks` does not have the same length as `color`. The colors are interpolated from the minimal to the maximal of `breaks`.")
            }
        }
    }
    
    if(!identical(filename, NA)) {
        warning("argument `filename` is not supported in pheatmap -> Heatmap translation, skip it.")
    }

    if(!identical(width, NA)) {
        warning("argument `width` is not supported in pheatmap -> Heatmap translation, skip it.")
    }
    
    if(!identical(height, NA)) {
        warning("argument `height` is not supported in pheatmap -> Heatmap translation, skip it.")
    }

    if(!identical(silent, FALSE)) {
        warning("argument `silent` is not supported in pheatmap -> Heatmap translation, skip it.")
    }

    ht_param$rect_gp = gpar(col = border_color)

    if(nrow(mat) > 1000 || ncol(mat) > 1000) {
        if(!is.na(border_color)) {
            warning("border color is set for the matrix with large numbers of rows or columns. You might only be able to see the border colors in the plot. Set `border_color = NA` to get rid of it.")
        }
    }

    if(!identical(cellwidth, NA)) {
        ht_param$width = ncol(mat)*unit(cellwidth, "pt")
    }

    if(!identical(cellheight, NA)) {
        ht_param$height = nrow(mat)*unit(cellheight, "pt")
    }

    if(identical(clustering_distance_rows, "correlation")) clustering_distance_rows = "pearson"
    if(identical(clustering_distance_cols, "correlation")) clustering_distance_cols = "pearson"

    ht_param$cluster_rows = cluster_rows
    ht_param$cluster_columns = cluster_cols 
    ht_param$clustering_distance_rows = clustering_distance_rows
    ht_param$clustering_distance_columns = clustering_distance_cols 
    ht_param$clustering_method_rows = clustering_method
    ht_param$clustering_method_columns = clustering_method

    if(!is.na(cutree_rows)) {
        if(inherits(cluster_rows, c("logical", "hclust", "dendrogram"))) {
            ht_param$row_split = cutree_rows
            ht_param$row_gap = unit(4, "bigpts")
            ht_param["row_title"] = list(NULL)
        }
    }
    if(!is.na(cutree_cols)) {
        if(inherits(cluster_cols, c("logical", "hclust", "dendrogram"))) {
            ht_param$column_split = cutree_cols
            ht_param$column_gap = unit(4, "bigpts")
            ht_param["column_title"] = list(NULL)
        }
    }
    
    ht_param$row_dend_width = unit(treeheight_row, "pt")
    ht_param$column_dend_height = unit(treeheight_col, "pt")

    ht_param$show_heatmap_legend = legend

    if(identical(scale, "row") || identical(scale, "column")) {
        if(identical(legend_breaks, NA)) {
            lim = quantile(abs(mat), 0.975, na.rm = TRUE)

            le = pretty(c(-lim, lim), n = 3)
            if(length(le) == 7 && le[1] == -3) {
                le = c(-3, -1.5, 0, 1.5, 3)
            } else if(! 0 %in% le) {
                le = c(le[1], le[1]/2, 0, le[length(le)]/2, le[length(le)])
            }
            legend_breaks = le
        }
    }
    if(!identical(legend_breaks, NA)) {
        heatmap_legend_param$at = legend_breaks
    }
    if(!identical(legend_labels, NA)) {
        heatmap_legend_param$labels = legend_labels
    }
    ht_param$heatmap_legend_param = heatmap_legend_param

    if(identical(annotation_colors, NA)) {
        annotation_colors = list()
    }
    if(!identical(annotation_col, NA)) {
        acn = rownames(annotation_col)
        mcn = colnames(mat)
        if(!is.null(acn)) {
            if(acn[1] %in% mcn) {
                if(length(union(acn, mcn)) == length(mcn)) {
                    if(!identical(acn, mcn)) {
                        warning("Column annotation has different order from matrix columns. Adjust the column annotation based on column names of the matrix.")
                    }
                    annotation_col = annotation_col[mcn, , drop = FALSE]
                }
            }
        }
        for(nm in colnames(annotation_col)) {
            if(nm %in% names(annotation_colors)) {
                if(is.null(names(annotation_colors[[nm]])) && is.numeric(annotation_col[, nm])) {
                    foo_x = annotation_col[, nm]
                    foo_n_col = length(annotation_colors[[nm]])
                    annotation_colors[[nm]] = colorRamp2(seq(min(foo_x), max(foo_x), length.out = foo_n_col), annotation_colors[[nm]])
                }
            }
        }
        ht_param$top_annotation = HeatmapAnnotation(df = annotation_col[, ncol(annotation_col):1, drop = FALSE], 
            col = annotation_colors, show_legend = annotation_legend,
            show_annotation_name = annotation_names_col, gp = gpar(col = border_color),
            annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
            simple_anno_size = unit(10, "bigpts"), gap = unit(2, "bigpts"))
    }
    if(!identical(annotation_row, NA)) {
        arn = rownames(annotation_row)
        mrn = rownames(mat)
        if(!is.null(arn)) {
            if(arn[1] %in% mrn) {
                if(length(union(arn, mrn)) == length(mrn)) {
                    if(!identical(arn, mrn)) {
                        warning("Row annotation has different order from matrix rows. Adjust the row annotation based on row names of the matrix.")
                    }
                    annotation_row = annotation_row[mrn, , drop = FALSE]
                }
            }
        }
        for(nm in colnames(annotation_row)) {
            if(nm %in% names(annotation_colors)) {
                if(is.null(names(annotation_colors[[nm]])) && is.numeric(annotation_row[, nm])) {
                    foo_x = annotation_row[, nm]
                    foo_n_col = length(annotation_colors[[nm]])
                    annotation_colors[[nm]] = colorRamp2(seq(min(foo_x), max(foo_x), length.out = foo_n_col), annotation_colors[[nm]])
                }
            }
        }
        ht_param$left_annotation = rowAnnotation(df = annotation_row[, ncol(annotation_row):1, drop = FALSE], 
            col = annotation_colors, show_legend = annotation_legend,
            show_annotation_name = annotation_names_row, gp = gpar(col = border_color),
            annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
            simple_anno_size = unit(10, "bigpts"), gap = unit(2, "bigpts"))
    }

    if(!identical(annotation, NA)) {
        warning("argument `annotation` is not supported in pheatmap -> Heatmap translation, skip it.")
    }

    if(identical(drop_levels, FALSE)) {
        warning("argument `drop_levels` is enfored to be TRUE, skip it.")
    }

    ht_param$show_row_names = show_rownames
    ht_param$show_column_names = show_colnames
    
    ht_param$row_names_gp = gpar(fontsize = fontsize_row, fontfamily = fontfamily_row, fontface = fontface_row)
    ht_param$column_names_gp = gpar(fontsize = fontsize_col, fontfamily = fontfamily_col, fontface = fontface_col)

    angle_col = match.arg(angle_col)[1]
    angle_col = switch(angle_col, 
                        "0" = 0,
                        "45" = 45,
                        "90" = 90,
                        "270" = 90,
                        "315" = -45)
    ht_param$column_names_rot = angle_col
    if(angle_col == 0) {
        ht_param$column_names_centered = TRUE
    }

    if(is.logical(display_numbers)) {
        if(display_numbers) {
            ht_param$layer_fun = local({
                number_format = number_format
                number_color = number_color
                fontsize_number = fontsize_number
                mat = mat
                function(j, i, x, y, w, h, fill) {
                    grid.text(sprintf(number_format, pindex(mat, i, j)), x = x, y = y, gp = gpar(col = number_color, fontsize = fontsize_number))
                }
            })
        }
    } else if(is.matrix(display_numbers)) {
        if(!identical(dim(display_numbers), dim(mat))) {
            stop_wrap("dimension of `display_numbers` should be the same as the input matrix.")
        }
        ht_param$layer_fun = local({
            number_color = number_color
            fontsize_number = fontsize_number
            mat = display_numbers
            function(j, i, x, y, w, h, fill) {
                grid.text(pindex(mat, i, j), x = x, y = y, gp = gpar(col = number_color, fontsize = fontsize_number))
            }
        })
    }
    
    if(!is.null(labels_row)) {
        ht_param$row_labels = labels_row
    }

    if(!is.null(labels_col)) {
        ht_param$column_labels = labels_col
    }

    if(!is.null(gaps_row)) {
        if(inherits(cluster_rows, c("hclust", "dendrogram"))) {
            stop_wrap("`gaps_row` should not be set when `cluster_rows` is set as a clustering object.")
        }
        if(identical(cluster_rows, TRUE)) {
            stop_wrap("`gaps_row` should not be set when `cluster_rows` is set to TRUE.")
        }
        slices = diff(c(0, gaps_row, nrow(mat)))
        ht_param$row_split = rep(seq_along(slices), times = slices)
        ht_param$row_gap = unit(4, "bigpts")
        ht_param["row_title"] = list(NULL)
    }
    if(!is.null(gaps_col)) {
        if(inherits(cluster_cols, c("hclust", "dendrogram"))) {
            stop_wrap("`gaps_col` should not be set when `cluster_cols` is set as a clustering object.")
        }
        if(identical(cluster_cols, TRUE)) {
            stop_wrap("`gaps_col` should not be set when `cluster_cols` is set to TRUE.")
        }
        slices = diff(c(0, gaps_col, ncol(mat)))
        ht_param$column_split = rep(seq_along(slices), times = slices)
        ht_param$column_gap = unit(4, "bigpts")
        ht_param["column_title"] = list(NULL)
    }

    #if(!identical(clustering_callback, NA)) {
    #    if(!identical(ht_param$cluster_rows, FALSE)) {
    #        row_hclust = hclust(get_dist(mat, ht_param$clustering_distance_rows), ht_param$clustering_method_rows)
    #        row_hclust = clustering_callback(row_hclust, ...)
    #        ht_param$cluster_rows = row_hclust
    #    }
    #    if(!identical(ht_param$cluster_columns, FALSE)) {
            column_hclust = #hclust(get_dist(t(mat), 
            hclust(dist(t(mat), ht_param$clustering_distance_columns), ht_param$clustering_method_columns) 
            column_hclust = clustering_callback(column_hclust, ...)
            ht_param$cluster_columns = column_hclust
    #    }
    #}

    ### default from pheatmap
    ht_param$name = name
    ht_param$row_dend_reorder = FALSE
    ht_param$column_dend_reorder = FALSE

    if(!identical(main, NA)) {
        ht_param$column_title = main
        ht_param$column_title_gp = gpar(fontface = "bold", fontsize = 1.3*fontsize)
    }
    ht_param = c(ht_param, list(...))
    ht = do.call(Heatmap, ht_param)
    attr(ht, "translate_from") = "pheatmap"

    if(run_draw) {
        draw(ht)
    } else {
        ht
    }
}

print("template_function_GSEA_LE_Heatmap_H_Myc_Targets_V2_3_Contrasts.R #########################################################################")
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
var_GSEA_LE_Heatmap_H_Myc_Targets_V2_3_Contrasts<-GSEA_LE_Heatmap_H_Myc_Targets_V2_3_Contrasts(var_Batch_Corrected_Counts,var_ccbr991_metadata_nidap)
invisible(graphics.off())
saveRDS(var_GSEA_LE_Heatmap_H_Myc_Targets_V2_3_Contrasts, paste0(rds_output,"/var_GSEA_LE_Heatmap_H_Myc_Targets_V2_3_Contrasts.rds"))
