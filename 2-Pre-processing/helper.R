# helper functions
library(Seurat)
library(dplyr)

# function to shuffle data points for plotting
shuf <- function(df){
  return(df[sample(1:dim(df)[1], dim(df)[1]),])
}

# function to perform dimensionality reductions on default assay of seurat_obj and normalization on ADT assay, if ADT_norm=TRUE
dim_reductions_RNA_ADT <- function(s_obj, ADT_norm=TRUE){
  set.seed(1989)
  s_obj@reductions <- list()
  all.genes <- rownames(s_obj)
  s_obj <- s_obj %>% 
    NormalizeData() %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() 
  
  s_obj <- s_obj %>%
    RunPCA(features = VariableFeatures(object = s_obj), npcs=30) %>%
    RunUMAP(dims = 1:20, n.neighbors=20L, reduction.name = "umapndim20nn20lognorm", reduction.key = "UMAPndim20nn20lognorm_", verbose = FALSE) %>%
    RunUMAP(dims = 1:25, n.neighbors=20L, reduction.name = "umapndim25nn20lognorm", reduction.key = "UMAPndim25nn20lognorm_", verbose = FALSE) %>%
    RunUMAP(dims = 1:30, n.neighbors=20L, reduction.name = "umapndim30nn20lognorm", reduction.key = "UMAPndim30nn20lognorm_", verbose = FALSE) %>%
    RunUMAP(dims = 1:20, n.neighbors=25L, reduction.name = "umapndim20nn25lognorm", reduction.key = "UMAPndim20nn25lognorm_", verbose = FALSE) %>%
    RunUMAP(dims = 1:25, n.neighbors=25L, reduction.name = "umapndim25nn25lognorm", reduction.key = "UMAPndim25nn25lognorm_", verbose = FALSE) %>%
    RunUMAP(dims = 1:30, n.neighbors=25L, reduction.name = "umapndim30nn25lognorm", reduction.key = "UMAPndim30nn25lognorm_", verbose = FALSE) %>%
    RunUMAP(dims = 1:20, n.neighbors=30L, reduction.name = "umapndim20nn30lognorm", reduction.key = "UMAPndim20nn30lognorm_", verbose = FALSE) %>%
    RunUMAP(dims = 1:25, n.neighbors=30L, reduction.name = "umapndim25nn30lognorm", reduction.key = "UMAPndim25nn30lognorm_", verbose = FALSE) %>%
    RunUMAP(dims = 1:30, n.neighbors=30L, reduction.name = "umapndim30nn30lognorm", reduction.key = "UMAPndim30nn30lognorm_", verbose = FALSE) 
  
  if(ADT_norm){
    s_obj <- s_obj %>%
      NormalizeData(assay = "ADT", normalization.method = "CLR",margin = 2) %>%
      ScaleData(assay = "ADT")
    
    return(s_obj)
  }
}

## geneset score function adjustment:
geneset_score = function(module_tbl, counts_raw, min_cpm = 0, limit_pct = 1) {
  # perform the cell type enrichment calculation based on rescaled values

  module_list <- module_tbl %>%
    filter(.$gene %in% rownames(counts_raw)) %>%
    with(split(.$gene, celltype))

  if (class(counts_raw) != "matrix") { stop("expression matrix is not a matrix") }
  if (max(counts_raw) < 100) { stop("expression values appear to be log-scaled") }

  # filter matrix for expressed genes only
  counts_raw = filter_mat_by_cpm(counts_raw = counts_raw, min_cpm = min_cpm)

  # rescale matrix for expressed genes only
  counts_raw_subs = normalize_mat_by_gene(counts_raw = subset(counts_raw, rownames(counts_raw) %in% unlist(module_list)), limit_pct = limit_pct)


  # check if enough genes pass filter
  if (min(lengths(module_list)) < 3) { stop("too few genes per celltype") }

  # calculate average z-score per celltype
  celltype_scores_tbl = tibble()
  for (ct in names(module_list)) {
    celltype_scores_tbl =
      bind_rows(
        celltype_scores_tbl,
        tibble(
          cell = colnames(counts_raw_subs),
          celltype = ct,
          score = colMeans(subset(counts_raw_subs, rownames(counts_raw_subs) %in% module_list[[ct]]))
        )
      )
    ct_scores = colnames(counts_raw_subs)
  }

  celltype_scores_tbl <- celltype_scores_tbl %>%
    spread(celltype, score) %>%
    rename_at(vars(-contains("cell")), list(~paste0(., ".score")))

  return(celltype_scores_tbl)
}


## calculate cluster averages for large objects:

calc_clust_averages_large <- function(metadata, data, group){
  
  # get relevant metadata
  metadata <- metadata %>%
    select("cell", group)
  
  # manipulate data to merge with metadata
  # genes are already columns, transpose, so genes are rows
  # remove all genes with less than 50 cells
  xx <- tabulate(data@i + 1)
  rows_to_keep <- rownames(data)[xx > 49]
  data <- data[rows_to_keep,]
  
  # split up sparse matrix by clusters and then calculate average of each df individually
  current_mean <- data.frame(rownames(data))
  colnames(current_mean) <- "gene"
  
  for (l in 1:length(unique(metadata[,eval(group)]))) {
    l <- l-1
    #determine cells belonging to each cluster
    cells <- subset(metadata, metadata[,eval(group)] == l)[,1]
    # subset for each cluster
    data1 <- as.data.frame(data[,cells])
    current_mean$cluster <- Matrix::rowMeans(data1)
    colnames(current_mean)[2+l] <- paste0("cluster_", l)
  }
  
  return(current_mean)
}

