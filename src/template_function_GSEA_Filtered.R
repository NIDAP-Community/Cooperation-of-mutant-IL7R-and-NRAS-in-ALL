# GSEA Filters [CCBR] (d3c1d012-4ecd-4ac3-9e85-fd68cda28ba0): v108
GSEA_Filtered <- function(GSEA_Preranked, msigdb_v6_2_with_orthologs ) {

## This function filters GSEA Table

## --------- ##
## Libraries ##
## --------- ##

library(dplyr)
library(ggplot2)
library(plotly)
library(RColorBrewer)
library(tidyverse)    

## -------------------------------- ##
## User-Defined Template Parameters ##
## -------------------------------- ##

# Primary inputs
gsea_table <- GSEA_Preranked
pathways_database <- msigdb_v6_2_with_orthologs

# Basic parameters
species <- "Mouse"
p_value_filter = 'adjusted p-value'
p_value_threshold = 0.05
enrichment_score_filter = 'NES (Normalized Enrichment Score)'
enrichment_score_threshold = 0
enrichment_score_sign = "+/-"
size_filter = "Pathway size"
size_cutoff = 0
top_rank_filter = "all"

# Advanced parameters
columns_to_sort_output_by = c()
sort_output_in_decreasing_order = FALSE
collections_to_include =  c()
pathways_to_include = c()
gene_filter_universe = "Leading Edge (LE)"
genes_to_include = c()
contrast_filter = "none"
contrasts = c()

#Visualization parameters
bubble_color_variable = "collection" # collection
bubble_color_opacity = 0.95
bubble_maximal_size = 2
x_axis_minimum = c()
x_axis_maximum = c()
y_axis_minimum = c()
y_axis_maximum = c()

# Legacy parameters
tested_contrast = c()
gene_score = c()

##--------------- ##
## Error Messages ##
## -------------- ##

## --------- ##
## Functions ##
## --------- #

filter.message <-
  function(filter, warn = FALSE, condition, output) {
    n = length(unique(output$pathway))
    if (!warn) {
      if (n == 0) {
        cat(
          sprintf(
            "ERROR: Filter by %s (%s) returned %g unique pathway(s)\n",
            filter,
            condition,
            n
          )
        )
        stop("Filter condition error\n")
      } else {
        cat(
          sprintf(
            "OK: Filter by %s (%s) returned %g unique pathway(s)\n",
            filter,
            condition,
            n
          )
        )
      }
    } else {
      cat(
        sprintf(
          "WARNING: Filter by %s (%s) not specified correctly; this filter is not applied\n",
          filter,
          condition
        )
      )
    }
  }

# adjust old GSEA table (v67 or lower)
adjust.v67 <-
  function(input,
           gS,
           tC,
           sp,
           db,
           required_columns = c(
             "contrast",
             "geneScore",
             "fdr_correction_mode",
             "collection",
             "pathway",
             "pval",
             "padj",
             "ES",
             "NES",
             "nMoreExtreme",
             "size",
             "leadingEdge",
             "size_leadingEdge",
             "inPathway"
           )) {
    missing_columns = setdiff(required_columns, colnames(input))
    
    if (length(missing_columns) > 0) {
      cat("WARNING: Output from outdated 'Preranked GSEA [CCBR]' detected.\n")
      cat(
        sprintf(
          "\nWARNING: Missing columns are added (%s)\n\t 'inPathway' column includes all genes annotated to a gene set - run the latest released version of 'Preranked GSEA [CCBR]'\n\t to return only genes mapped in the dataset.\n",
          paste(missing_columns, collapse = ", ")
        )
      )
      
      if (is.null(gS)) {
        gS = "NULL"
        cat(
          "\nWARNING: 'Gene score' parameter not provided therefore 'geneScore' column is assigned the default NULL value\n\t if downstream 'GSEA Running Score Diagram & Leading Edge Heatmap [Bulk] [CCBR]' will be linked to this output,\n\t it will fail due to this specification of the 'Gene score' parameter.\n"
        )
      } else {
        if (length(gS) > 1) {
          gS = gS[1]
          cat(
            sprintf(
              "\nWARNING: too many values for 'Gene score' provided; only the first one used, '%s'\n",
              gS
            )
          )
        }
        has.underscore = substring(gS, 1, 1) == "_"
        if (!has.underscore) {
          gS = paste0("_", gS[1])
          cat(
            sprintf(
              "\nWARNING: 'Gene score' should start with underscore to match column naming convention in a DEG table used for GSEA run (e.g. '_tstat');\n\t column 'geneScore' is assigned the provided value with '_' added ('%s');\n\t if downstream 'GSEA Running Score Diagram & Leading Edge Heatmap [Bulk] [CCBR]' will be linked to this output,\n\t it may fail if this is not an adequate specification of the 'Gene score' parameter.\n",
              gS
            )
          )
        }
      }
      
      if (is.null(tC)) {
        tC = "NULL"
        cat(
          "\nWARNING: 'Tested contrast' parameter not provided therefore 'contrast' column is assigned the default NULL value;\n\t if downstream 'GSEA Running Score Diagram & Leading Edge Heatmap [Bulk] [CCBR]' will be linked to this output,\n\t it will fail due to this specification of the 'Tested contrast' parameter.\n"
        )
      } else {
        if (length(tC) > 1) {
          tC = tC[1]
          cat(
            sprintf(
              "\nWARNING: too many values for 'Tested contrast' provided; only the first one used, '%s'\n",
              tC
            )
          )
        }
        is.contrast = any(grepl("-", tC))
        if (!is.contrast) {
          cat(
            sprintf(
              "\nWARNING: 'Tested contrast' should be specified exactly the same way as when Preranked GSEA was run (e.g. treated-control);\n\t column 'contrast' is assigned the provided '%s' value;\n\t if downstream 'GSEA Running Score Diagram & Leading Edge Heatmap [Bulk] [CCBR]' will be linked to this output,\n\t it may fail if this is not an adequate specification of the 'Tested contrast' parameter.\n",
              tC[1]
            )
          )
        }
      }
      
      input <-
        input %>% dplyr::mutate(
          contrast = tC,
          geneScore = gS,
          fdr_correction_mode = "over all collections",
          size_leadingEdge = sapply(strsplit(leadingEdge, ","), length)
        )
      input <-
        input[, match(required_columns[required_columns != 'inPathway'], colnames(input))]
      db <-
        db %>% filter(db[["species"]] == sp) %>% filter(`%in%`(db[["collection"]] , unique(input$collection))) %>% filter(`%in%`(db[["gene_set_name"]], unique(input$pathway))) %>% collect()
      db <-
        db %>% dplyr::group_by(collection, gene_set_name, species) %>% dplyr::summarize(inPathway = paste(unique(gene_symbol), collapse = ",")) %>% dplyr::ungroup() %>% data.frame()
      input <-
        input %>% dplyr::left_join(db,
                                   by = c("collection" = "collection", "pathway" = "gene_set_name"))
      input <-
        input[, na.omit(match(c(required_columns, "species"), colnames(input)))]
      return(input)
      
    } else {
      cat("Filtering Preranked GSEA table.\n")
      return(input)
    }
  }

    
## --------------- ##
## Main Code Block ##
## --------------- ##

# translate filters to column names
filterInput_byScore = switch(
  enrichment_score_filter,
  'NES (Normalized Enrichment Score)' = 'NES',
  'ES (Enrichment Score)' = 'ES'
)
filterInput_byPvalue = switch(p_value_filter,
                              'p-value' = 'pval',
                              'adjusted p-value' = 'padj')

## adjust old Preranked GSEA output (release v67 or lower)
gsea_table <-
  adjust.v67(
    input = gsea_table,
    gS = gene_score,
    tC = tested_contrast,
    sp = species,
    db = pathways_database
  )

## apply filters
cat("\n\nFiltering steps\n\n")

gsea_filtered <- gsea_table %>%
  dplyr::filter(get(filterInput_byPvalue) <= p_value_threshold)
filter.message(filter = p_value_filter,
               condition = p_value_threshold,
               output = gsea_filtered)

if (enrichment_score_sign == "+") {
  gsea_filtered <- gsea_filtered %>%
    dplyr::filter(get(filterInput_byScore) >= enrichment_score_threshold)
  filter.message(
    filter = "value of GSEA score",
    condition = sprintf("%s > %g", filterInput_byScore, enrichment_score_threshold),
    output = gsea_filtered
  )
  
} else if (enrichment_score_sign == "-") {
  gsea_filtered <- gsea_filtered %>%
    dplyr::filter(get(filterInput_byScore) <= -enrichment_score_threshold)
  filter.message(
    filter = "value of GSEA score",
    condition = sprintf(
      "%s < %s%g",
      filterInput_byScore,
      ifelse(enrichment_score_threshold == 0, "", "-"),
      enrichment_score_threshold
    ),
    output = gsea_filtered
  )
  
} else {
  gsea_filtered <-
    gsea_filtered %>% dplyr::filter(abs(get(filterInput_byScore)) >= enrichment_score_threshold)
  filter.message(
    filter = "value of GSEA score",
    condition = sprintf("|%s| > %g", filterInput_byScore, enrichment_score_threshold),
    output = gsea_filtered
  )
}

if (size_cutoff > 0) {
  if (size_filter == "Pathway size") {
    gsea_filtered <-
      gsea_filtered %>% dplyr::filter(size >= size_cutoff)
    filter.message(filter = size_filter,
                   condition = size_cutoff,
                   output = gsea_filtered)
    
  } else {
    gsea_filtered <-
      gsea_filtered %>% dplyr::filter(size_leadingEdge >= size_cutoff)
    filter.message(filter = size_filter,
                   condition = size_cutoff,
                   output = gsea_filtered)
  }
  
} else {
  filter.message(
    filter = size_filter,
    condition = paste0(">", size_cutoff),
    output = gsea_filtered
  )
}

if (!is.null(collections_to_include)) {
  gsea_filtered <-
    gsea_filtered %>% dplyr::filter(collection %in% collections_to_include)
  filter.message(filter = "Collection filter ON",
                 condition = "Collections to include",
                 output = gsea_filtered)
  cat(sprintf("    Found collections: %s\n", paste(
    unique(gsea_filtered$collection), collapse = ", "
  )))
}

if (!is.null(pathways_to_include)) {
  gsea_filtered <-
    gsea_filtered %>% dplyr::filter(pathway %in% pathways_to_include)
  filter.message(filter = "Pathway filter ON",
                 condition = "Pathways to include",
                 output = gsea_filtered)
  cat(sprintf("    Found pathways: %s\n", paste(
    unique(gsea_filtered$pathway), collapse = ", "
  )))
}

if (!is.null(genes_to_include)) {
  if (gene_filter_universe == 'Leading Edge (LE)') {
    index_pathway = lapply(genes_to_include, function(x)
      grep(paste0("\\b\\Q", x, "\\E\\b"), gsea_filtered$leadingEdge))
    found = sapply(index_pathway, function(x) {
      length(x) > 0
    })
    index_pathway = Reduce(union, index_pathway)
    gsea_filtered <- gsea_filtered %>% dplyr::slice(index_pathway)
    
    filter.message(filter = "Gene filter ON",
                   condition = "Leading Edge",
                   output = gsea_filtered)
    cat(sprintf(
      "    Found genes: %s\n",
      paste(genes_to_include[which(found)], collapse = ", ")
    ))
    if (sum(!found) > 0) {
      cat(sprintf(
        "    Missing genes: %s\n",
        paste(genes_to_include[which(!found)], collapse = ", ")
      ))
    }
    
  } else if (gene_filter_universe == "Pathway") {
    index_pathway = lapply(genes_to_include, function(x)
      grep(paste0("\\b\\Q", x, "\\E\\b"), gsea_filtered$inPathway))
    found = sapply(index_pathway, function(x) {
      length(x) > 0
    })
    index_pathway = Reduce(union, index_pathway)
    gsea_filtered <- gsea_filtered %>% dplyr::slice(index_pathway)
    
    filter.message(filter = "Gene filter ON",
                   condition = "in Pathway",
                   output = gsea_filtered)
    cat(sprintf(
      "    Found genes: %s\n",
      paste(genes_to_include[which(found)], collapse = ", ")
    ))
    if (sum(!found) > 0) {
      cat(sprintf(
        "    Missing genes: %s\n",
        paste(genes_to_include[which(!found)], collapse = ", ")
      ))
    }
    
  }
}

if (contrast_filter == "remove") {
  if (!is.null(contrasts)) {
    all_contrasts = unique(gsea_filtered$contrast)
    gsea_filtered <-
      gsea_filtered %>% dplyr::filter(!contrast %in% contrasts)
    removed = setdiff(all_contrasts, unique(gsea_filtered$contrast))
    if (length(removed) < 1) {
      filter.message(
        warn = TRUE,
        filter = "Contrast filter",
        condition = sprintf("remove %s missing", paste(contrasts, collapse = ", ")),
        output = gsea_filtered
      )
    } else {
      filter.message(filter = "Contrast filter",
                     condition = contrast_filter,
                     output = gsea_filtered)
      cat(sprintf(
        "    Removed contrast(s): %s\n",
        paste(removed, collapse = ", ")
      ))
      cat(sprintf("    Keep contrast(s): %s\n", paste(
        unique(gsea_filtered$contrast), collapse = ", "
      )))
    }
    
  } else if (is.null(contrasts)) {
    filter.message(
      warn = TRUE,
      filter = "Contrast filter",
      condition = class(contrasts),
      output = gsea_filtered
    )
  }
  
} else if (contrast_filter == "keep") {
  if (!is.null(contrasts)) {
    all_contrasts = unique(gsea_filtered$contrast)
    gsea_filtered <-
      gsea_filtered %>% dplyr::filter(contrast %in% contrasts)
    kept = intersect(all_contrasts, unique(gsea_filtered$contrast))
    removed = setdiff(all_contrasts, unique(gsea_filtered$contrast))
    if (length(kept) < 1) {
      filter.message(
        warn = TRUE,
        filter = "Contrast filter",
        condition = sprintf("keep %s missing", paste(contrasts, collapse = ", ")),
        output = gsea_filtered
      )
    } else {
      filter.message(filter = "Contrast filter",
                     condition = contrast_filter,
                     output = gsea_filtered)
      cat(sprintf(
        "    Removed contrast(s): %s\n",
        paste(removed, collapse = ", ")
      ))
      cat(sprintf("    Kept contrast(s): %s\n", paste(
        unique(gsea_filtered$contrast), collapse = ", "
      )))
    }
    
  } else if (is.null(contrasts)) {
    filter.message(
      warn = TRUE,
      filter = "Contrast filter",
      condition = class(contrasts),
      output = gsea_filtered
    )
  }
  
} else if (contrast_filter == "none") {
  if (!is.null(contrasts)) {
    filter.message(
      warn = TRUE,
      filter = "Contrast filter",
      condition = sprintf("none; %s", paste(contrasts, collapse = ", ")),
      output = gsea_filtered
    )
  }
}

top_rank_filter = tolower(top_rank_filter)

if (top_rank_filter == "all") {
  filterInput_keep = paste("Inf (ALL)")
  top_rank_filter = Inf
  
} else {
  top_rank_filter = as.numeric(top_rank_filter)
  
  if (is.na(top_rank_filter)) {
    stop("ERROR in Top rank filter; enter ALL (case insensitive) or a numeric rank.\n")
  }
  else if (top_rank_filter <= 0) {
    top_rank_filter = 1
    cat("WARNING: 'Top rank filter' cannot be 0 or less; its value was changed to 1.\n")
  }
  filterInput_keep = top_rank_filter
}

gsea_filtered <-
  gsea_filtered %>% dplyr::group_by(contrast, collection) %>% dplyr::mutate(p_rank = rank(pval, ties.method =
                                                                                            "min")) %>% dplyr::filter(p_rank <= top_rank_filter) %>% dplyr::select(-p_rank)
filter.message(
  filter = "Top significant filter",
  condition = sprintf(
    "up to p-value rank of %s per contrast and collection",
    filterInput_keep
  ),
  output = gsea_filtered
)

# OUTPUT ====

## sort output
if (sort_output_in_decreasing_order) {
  columns_to_sort_output_by = sapply(columns_to_sort_output_by, function(x)
    sprintf("desc(%s)", x))
}
gsea_filtered <-
  gsea_filtered %>% dplyr::arrange_(.dots = columns_to_sort_output_by)

## do plot and return dataset

if (nrow(gsea_filtered) == 0) {
  stop("ERROR: filtering returned 0 pathways")
  
} else {
  cat("\n\nFiltered pathways\n")
  tab = table(gsea_filtered$collection, gsea_filtered$contrast) %>% addmargins(margin =
                                                                                 c(1, 2))
  print(tab)
  
  cat("\n\nGSEA statistics (filtered pathways)\n\n")
  print(tibble(gsea_filtered))
  
  df <-
    gsea_filtered %>% dplyr::select(-leadingEdge, -inPathway) %>% dplyr::mutate(textContrast =
                                                                                  sprintf(
                                                                                    "%s: NES = %g, P-value = %g",
                                                                                    contrast,
                                                                                    signif(NES, 2),
                                                                                    signif(pval, 1)
                                                                                  )) %>% dplyr::group_by(pathway, collection) %>%
    dplyr::summarize(
      mean_pval = mean(pval),
      n_contrast = length(contrast),
      Pathway_size = mean(size),
      mean_NES = mean(NES),
      individual_values = paste(textContrast, collapse = "\n")
    ) %>%
    dplyr::mutate(
      textPathway = sprintf(
        "%s<br>%s<br>mean NES = %g, mean P-value = %g<br>Pathway size = %g<br><br>N contrasts = %g\n%s",
        collection,
        pathway,
        signif(mean_NES, 2),
        signif(mean_pval, 1),
        Pathway_size,
        n_contrast,
        individual_values
      )
    ) %>%
    dplyr::select(
      collection,
      pathway,
      Pathway_size,
      n_contrast,
      mean_NES,
      mean_pval,
      individual_values,
      textPathway
    )
  
  find_sort = grepl(paste(columns_to_sort_output_by, collapse = "|"),
                    colnames(df))
  if (sum(find_sort) == 0) {
    sort_by = c("n_contrast", "collection", "mean_pval")
  } else {
    sort_by = colnames(df)[which(find_sort)]
  }
  
  if (sort_output_in_decreasing_order) {
    sort_by = sapply(sort_by, function(x)
      sprintf("desc(%s)", x))
  }
  df <- df %>% dplyr::arrange_(.dots = sort_by)
  
  #cat("\n\nGSEA contrast summary (filtered pathways)\n\n")
  #print(tibble(df %>% dplyr::select(-textPathway) %>% dplyr::rename("mean_size"="Pathway_size")))
  
  if (bubble_color_variable == "pathway size") {
    ggp <-
      plot_ly(
        data = df,
        x = ~ mean_NES,
        y = ~ -log10(mean_pval),
        color = ~ Pathway_size,
        size = ~ n_contrast,
        text = ~ textPathway,
        hoverinfo = "text",
        opacity = bubble_color_opacity,
        marker = list(sizeref = ~ 2 * max(n_contrast) / bubble_maximal_size ** 2)
      ) %>%
      
      plotly::layout(
        xaxis = list(
          title = "Normalized Enrichment Score (mean)",
          tickfont = list(size = 15),
          titlefont = list(size = 15),
          showgrid = TRUE
        ),
        yaxis = list(
          title = "-log10 P-value (mean)",
          tickfont = list(size = 15),
          titlefont = list(size = 15),
          showgrid = TRUE
        ),
        legend = list(itemsizing = 'constant'),
        title = list(
          text = "Filtered pathways; bubble area proportional to the number of filtered contrasts per pathway",
          x = 0,
          font = list(size = 15)
        ),
        showlegend = TRUE
      )
    
  } else {
    ggp <-
      plot_ly(
        data = df,
        x = ~ mean_NES,
        y = ~ -log10(mean_pval),
        color = ~ collection,
        size = ~ n_contrast,
        text = ~ textPathway,
        hoverinfo = "text",
        opacity = bubble_color_opacity,
        marker = list(sizeref = ~ 2 * max(n_contrast) / bubble_maximal_size ** 2)
      ) %>%
      
      plotly::layout(
        xaxis = list(
          title = "Normalized Enrichment Score (mean)",
          tickfont = list(size = 15),
          titlefont = list(size = 15),
          showgrid = TRUE
        ),
        yaxis = list(
          title = "-log10 P-value (mean)",
          tickfont = list(size = 15),
          titlefont = list(size = 15),
          showgrid = TRUE
        ),
        legend = list(itemsizing = 'constant'),
        title = list(
          text = "Filtered pathways; bubble area proportional to the number of filtered contrasts per pathway",
          x = 0,
          font = list(size = 15)
        ),
        showlegend = TRUE
      )
  }
  # custom axis range?
  
  if (!(is.null(x_axis_minimum) & is.null(x_axis_maximum))) {
    if (is.null(x_axis_minimum))
      x_axis_minimum = floor(min(df$mean_NES))
    if (is.null(x_axis_maximum))
      x_axis_maximum = ceiling(max(df$mean_NES))
    ggp <-
      ggp %>% plotly::layout(xaxis = list(range = list(x_axis_minimum, x_axis_maximum)))
  }
  
  if (!(is.null(y_axis_minimum) & is.null(y_axis_maximum))) {
    if (is.null(y_axis_minimum))
      y_axis_minimum = floor(min(-log10(df$mean_pval)))
    if (is.null(y_axis_maximum))
      y_axis_maximum = ceiling(max(-log10(df$mean_pval)))
    ggp <-
      ggp %>% plotly::layout(yaxis = list(range = list(y_axis_minimum, y_axis_maximum)))
  }
  
  print(ggp)
  return(gsea_filtered)
  
}

}

## ---------------------------- ##
## Global Imports and Functions ##
## ---------------------------- ##

## Functions defined here will be available to call in
## the code for any table.

## --------------- ##
## End of Template ##
## --------------- ##

print("template_function_GSEA_Filtered.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_GSEA_Preranked<-readRDS(paste0(rds_output,"/var_GSEA_Preranked.rds"))
Input_is_Seurat_count <- 0
for(item in var_GSEA_Preranked){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_GSEA_Preranked<-as.data.frame(var_GSEA_Preranked)}else{var_GSEA_Preranked <- var_GSEA_Preranked}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_msigdb_v6_2_with_orthologs<-readRDS(paste0(rds_output,"/var_msigdb_v6_2_with_orthologs.rds"))
Input_is_Seurat_count <- 0
for(item in var_msigdb_v6_2_with_orthologs){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_msigdb_v6_2_with_orthologs<-as.data.frame(var_msigdb_v6_2_with_orthologs)}else{var_msigdb_v6_2_with_orthologs <- var_msigdb_v6_2_with_orthologs}
invisible(graphics.off())
var_GSEA_Filtered<-GSEA_Filtered(var_GSEA_Preranked,var_msigdb_v6_2_with_orthologs)
invisible(graphics.off())
saveRDS(var_GSEA_Filtered, paste0(rds_output,"/var_GSEA_Filtered.rds"))
