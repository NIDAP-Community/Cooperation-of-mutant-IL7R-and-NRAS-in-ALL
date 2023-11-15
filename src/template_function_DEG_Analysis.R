# DEG Analysis [CCBR] (de953b9c-b7f3-49ac-b883-55070d759e5d): v189
DEG_Analysis <- function(Filtered_Counts, ccbr991_metadata_nidap) {
    
    ## --------- ##
    ## Libraries ##
    ## --------- ##
    
    library(limma)
    library(tidyverse)
    library(edgeR)
    library(stringr)
    library(grid)
    library(gridExtra)
    

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##

    #Basic Parameters:
    counts_matrix <- Filtered_Counts 
    sample_metadata <- ccbr991_metadata_nidap
    gene_names_column="Gene"
    sample_name_column<-"Sample"
    columns_to_include = c("Gene","control1","control2","control3","IL7R1","IL7R2","IL7R3","NRAS1","NRAS2","NRAS3","Both1","Both2","Both3")
    contrast_variable_column<-c("Group")
    contrasts<-c("IL7R-Control","NRAS-Control","Both-Control","Both-IL7R","Both-NRAS")
    covariates_columns=c("Group")

    #Advanced Parameters:
    input_in_log_counts <- FALSE
    return_mean_and_sd<-FALSE
    return_normalized_counts<-TRUE
    normalization_method<-"quantile"
    
    ##--------------- ##
    ## Error Messages ##
    ## -------------- ##

    if(make.names(colnames(counts_matrix))!=colnames(counts_matrix)){
        print("Error: The following counts matrix column names are not valid:\n")
        print(colnames(counts_matrix)[make.names(colnames(counts_matrix))!=colnames(counts_matrix)])

        print("Likely causes are columns starting with numbers or other special characters eg spaces.")
        stop("Bad column names.")
    }
    
    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##
    
    samples_to_include=columns_to_include[columns_to_include%in%sample_metadata[,sample_name_column,drop=T]]
    anno_col=columns_to_include[columns_to_include%in%sample_metadata[,sample_name_column,drop=T]==F]

    samples_to_include <- samples_to_include[samples_to_include != gene_names_column]
    samples_to_include <- samples_to_include[samples_to_include != "Gene"]
    samples_to_include <- samples_to_include[samples_to_include != "GeneName"]
    

    ##create unique rownames to correctly add back Annocolumns at end of template
    counts_matrix[,gene_names_column]=paste0(counts_matrix[,gene_names_column],'_',1:nrow(counts_matrix))

    anno_col=c(anno_col,gene_names_column)%>%unique
    anno_tbl=counts_matrix[,anno_col,drop=F]%>%as.data.frame

    df.m <- counts_matrix[,c(gene_names_column,samples_to_include)]
    gene_names <- NULL
    gene_names$GeneID <- counts_matrix[,gene_names_column]
    
    ### This code block does input data validation
    sample_metadata <- sample_metadata[match(colnames(df.m),sample_metadata[,sample_name_column]),]
    sample_metadata <- sample_metadata[rowSums(is.na(sample_metadata)) != ncol(sample_metadata), ]
    df.m <- df.m[,match(sample_metadata[,sample_name_column],colnames(df.m))]
    
    #Create DGEList object from counts
    if(input_in_log_counts == TRUE){
        x <- DGEList(counts=2^df.m, genes=gene_names)
    } else {
        x <- DGEList(counts=df.m, genes=gene_names) 
    }
    
    #Put covariates in order 
    covariates_columns=covariates_columns[order(covariates_columns!=contrast_variable_column)]
    
    for(ocv in covariates_columns){
        sample_metadata[,ocv]=gsub(" ","_",sample_metadata[,ocv])
    }

    contrasts=gsub(" ","_",contrasts)
    cov <- covariates_columns[!covariates_columns %in% contrast_variable_column]

    #Combine columns if 2-factor analysis
    if(length(contrast_variable_column)>1){
        sample_metadata %>% dplyr::mutate(contmerge = paste0(.data[[contrast_variable_column[1]]],".",.data[[contrast_variable_column[2]]])) -> sample_metadata
    } else {
        sample_metadata %>% dplyr::mutate(contmerge = .data[[contrast_variable_column]]) -> sample_metadata
    }

    contrast_var <- factor(sample_metadata$contmerge)

    if(length(cov) >0){
        dm.formula <- as.formula(paste("~0 +", paste("contmerge", paste(cov, sep="+", collapse="+"),sep="+")))
        design=model.matrix(dm.formula, sample_metadata)
        colnames(design) <- gsub("contmerge","",colnames(design))
    } else {
        dm.formula <- as.formula(~0 + contmerge)
        design=model.matrix(dm.formula, sample_metadata)
        colnames(design) <- levels(contrast_var)
    }

    #colnames(design) <- str_replace_all(colnames(design), contrast_variable_column, "")
    
    if (normalization_method %in% c("TMM","TMMwzp","RLE","upperquartile")){
        x <- calcNormFactors(x, method = normalization_method) 
        rownames(x) <- x$genes$GeneID
        v <- voom(x,design=design,normalize="none")
    } else {
        v <- voom(x,design=design,normalize=normalization_method,save.plot = TRUE)
    }
    
    rownames(v$E) <- v$genes$GeneID
    as.data.frame(v$E) %>% rownames_to_column("Gene") -> df.voom
    fit <- lmFit(v, design)
    cm <- makeContrasts(contrasts = contrasts, levels=design)

    #Print Mean-variance Plot
    sx <- v$voom.xy$x
    sy <- v$voom.xy$y
    xyplot <- as.data.frame(cbind(sx,sy))
    voomline <- as.data.frame(cbind(x=v$voom.line$x,y=v$voom.line$y))
    
    g <- ggplot() +
        geom_point(data=xyplot, aes(x=sx,y=sy),size=1) +
        theme_bw() +
        geom_smooth(data=voomline, aes(x=x,y=y),color = "red") +
        ggtitle("voom: Mean-variance trend") +
        xlab(v$voom.xy$xlab) + ylab(v$voom.xy$ylab) + 
        theme(axis.title=element_text(size=12),
        plot.title = element_text(size = 14, face = "bold",hjust = 0.5))

    #Print out sample numbers:
    #
    sampsize <- colSums(design)
    titleval <- "Please note Sample size:"
    titletext <- paste(names(sampsize), sampsize, sep = "=", collapse = " \n ") 
    titleall <- paste(titleval,"\n",titletext,"\n\n\n")

    contrast <- colnames(cm)
    connames <- strsplit(contrast,"-")
    connames <- lapply(connames,function(x) {gsub("\\(","",gsub("\\)","",x))})
    contrastsize <- lapply(connames,function(x) sampsize[unlist(x)])
    footnotetext <- paste(contrast, contrastsize, sep = " : ", collapse = "\n") 
    footnotetext <- paste("\n\n\nContrasts:\n",footnotetext)

    textall <- textGrob(paste0(titleall, footnotetext),gp=gpar(fontsize=10))

    #Run Contrasts
    fit2 <- contrasts.fit(fit, cm)
    fit2 <- eBayes(fit2)
    logFC = fit2$coefficients
    colnames(logFC)=paste(colnames(logFC),"logFC",sep="_")
    tstat = fit2$t
    colnames(tstat)=paste(colnames(tstat),"tstat",sep="_")
    FC = 2^fit2$coefficients
    FC = ifelse(FC<1,-1/FC,FC)
    colnames(FC)=paste(colnames(FC),"FC",sep="_")
    pvalall=fit2$p.value
    colnames(pvalall)=paste(colnames(pvalall),"pval",sep="_")
    pvaladjall=apply(pvalall,2,function(x) p.adjust(x,"BH"))
    colnames(pvaladjall)=paste(colnames(fit2$coefficients),"adjpval",sep="_")

    
    if(return_mean_and_sd == TRUE){
        tve <- t(v$E)        
        mean.df <- as.data.frame(tve) %>% rownames_to_column("Sample") %>% dplyr::mutate(group=sample_metadata[sample_metadata[,sample_name_column]==Sample,contrast_variable_column]) %>% group_by(group) %>% summarise_all(funs(mean)) %>% as.data.frame()
        mean.df[,-c(1,2)] %>% as.matrix() %>% t() -> mean
        colnames(mean) <- mean.df[,1]
        colnames(mean)=paste(colnames(mean),"mean", sep="_")
        colnames(mean) = gsub("\\.", "_", colnames(mean))
        
        sd.df <- as.data.frame(tve) %>% rownames_to_column("Sample") %>% dplyr::mutate(group=sample_metadata[sample_metadata[,sample_name_column]==Sample,contrast_variable_column]) %>% group_by(group) %>% summarise_all(funs(sd)) %>% as.data.frame()
        sd.df[,-c(1,2)] %>% as.matrix() %>% t() -> sd
        colnames(sd) <- sd.df[,1]
        colnames(sd)=paste(colnames(sd), "sd",sep="_")
        colnames(sd) = gsub("\\.", "_", colnames(sd))
    finalres=as.data.frame(cbind(mean, sd,  FC, logFC, tstat, pvalall, pvaladjall)) 
    } else {
        finalres=as.data.frame(cbind(FC, logFC, tstat, pvalall, pvaladjall))
    }

    if(return_normalized_counts == TRUE){
        finalres = as.data.frame(cbind(finalres, v$E))
    }

    finalres %>% rownames_to_column("Gene") -> finalres
    print(paste0("Total number of genes included: ", nrow(finalres)))

    getgenelists <- function(FClimit,pvallimit,pval){
        upreggenes <- list()
        downreggenes <- list()
        for(i in 1:length(contrasts)){
            if(pval == "pval"){
            finalres %>% dplyr::filter(.data[[colnames(FC)[i]]] > FClimit & .data[[colnames(pvalall)[i]]] < pvallimit) %>% pull(Gene) %>% length() -> upreggenes[[i]] 
            finalres %>% dplyr::filter(.data[[colnames(FC)[i]]] < -FClimit & .data[[colnames(pvalall)[i]]] < pvallimit) %>% pull(Gene) %>% length() -> downreggenes[[i]]        
        } else {
            finalres %>% dplyr::filter(.data[[colnames(FC)[i]]] > FClimit & .data[[colnames(pvaladjall)[i]]] < pvallimit) %>% pull(Gene) %>% length() -> upreggenes[[i]] 
            finalres %>% dplyr::filter(.data[[colnames(FC)[i]]] < -FClimit & .data[[colnames(pvaladjall)[i]]] < pvallimit) %>% pull(Gene) %>% length() -> downreggenes[[i]] 
        }
        }
        names(upreggenes) <- contrasts
        names(downreggenes) <- contrasts
        allreggenes <- rbind(unlist(upreggenes),unlist(downreggenes))
        rownames(allreggenes) <- c(paste0("upreg>",FClimit, ", ",pval,"<",pvallimit),paste0("downreg<-",FClimit, ", ",pval,"<",pvallimit))
        return(allreggenes)
    }

    FCpval1 <- getgenelists(FClimit = 1.2, pvallimit = 0.05,"pval")
    FCpval2 <- getgenelists(FClimit = 1.2, pvallimit = 0.01,"pval")
    FCadjpval1 <- getgenelists(FClimit = 1.2, pvallimit = 0.05,"adjpval")
    FCadjpval2 <- getgenelists(FClimit = 1.2, pvallimit = 0.01,"adjpval")

    wraplines <- function(y){
        j = unlist(strsplit(y,"-"))
        k = strwrap(j, width = 10)
        l = paste(k,collapse="\n-")
        return(l)
    }
    
    pvaltab <- rbind(FCpval1,FCpval2,FCadjpval1,FCadjpval2)
    colnames(pvaltab) <- sapply(colnames(pvaltab), function(x) wraplines(x))
    table2 <- tableGrob(pvaltab, theme=ttheme_default(base_size = 10))
    table2$layout$clip <- "off"

    layout <- rbind(c(1,2),
                    c(1,2),
                   c(3,3))
    

    #Printing all together (tables and plot)
    grid.newpage()
    grid.arrange(textall, g, table2, layout_matrix=layout)

    #Printing in brand new multiviz
    grid.newpage()
    print(g)
    grid.newpage()
    grid.draw(textall)
    grid.newpage()
    grid.draw(table2)

### add back Anno columns and Remove row number from Feature Column
    colnames(finalres)[colnames(finalres)%in%"Gene"]=gene_names_column

    finalres=merge(anno_tbl,finalres,by=gene_names_column,all.y=T)
   finalres[,gene_names_column]=gsub('_[0-9]+$',"",finalres[,gene_names_column])

    call_me_alias<-colnames(finalres)
    colnames(finalres)<-gsub("\\(|\\)","", call_me_alias)
    df.final<-(finalres)
   
    return(df.final) 
}

print("template_function_DEG_Analysis.R #########################################################################")
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
var_DEG_Analysis<-DEG_Analysis(var_Filtered_Counts,var_ccbr991_metadata_nidap)
invisible(graphics.off())
saveRDS(var_DEG_Analysis, paste0(rds_output,"/var_DEG_Analysis.rds"))
