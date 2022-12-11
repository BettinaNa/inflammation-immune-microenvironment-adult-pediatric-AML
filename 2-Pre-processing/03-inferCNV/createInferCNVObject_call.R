createInferCNVObject_call <- function(seurat_obj=seurat_obj, s=sample, ...){

    #define cells of sample and control cells
    sample_cells <- seurat_obj@meta.data %>% rownames_to_column("cell") %>% filter(samples == eval(s)) %>% .$cell
    control_cells <- seurat_obj@meta.data %>% rownames_to_column("cell") %>% filter(samples %in% str_subset(unique(seurat_obj@meta.data$samples), "Control")) %>% .$cell
     
    # subsample control cells to 400 cells/cell type - 40 from each patient
    set.seed(4561)
    subsample_control <- seurat_obj@meta.data[control_cells,] %>% rownames_to_column("cell") %>% select(cell,broad_cell_type, samples) %>% group_by(broad_cell_type, samples) %>% sample_n(40, replace = TRUE) %>% .$cell %>% unique()
    
    if(str_detect(eval(s),"Control")){
      set.seed(4561)
      subsample_control <- seurat_obj@meta.data[control_cells,] %>% rownames_to_column("cell") %>% select(cell,broad_cell_type, samples) %>% group_by(broad_cell_type, samples) %>% sample_n(40, replace = TRUE) 
      subsample_control <- subsample_control %>% filter(samples != eval(s)) %>% .$cell %>% unique()
    }


    cnv_analysis_set <- c(subsample_control,sample_cells)

    # extract counts from seurat
    exp_mat1 <- GetAssayData(seurat_obj, assay = "RNA_soupx", slot = "counts")
    exp_mat1 <- exp_mat1[,cnv_analysis_set]

    # Ensembl 93 - July 2018 - used by Cell Ranger GRCh38 3.0.0
    ensembl <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                  dataset = "hsapiens_gene_ensembl", host = "jul2018.archive.ensembl.org")
    genes_tbl <- biomaRt::getBM(attributes =
                  c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position"),
                  filters = "hgnc_symbol", values = rownames(exp_mat1), mart = ensembl)

    genes_tbl <- genes_tbl %>%
                    filter(!chromosome_name %in% c("MT", "X", "Y")) %>%
                    filter(!str_detect(chromosome_name, "CHR_")) %>%
                    filter(!str_detect(chromosome_name, "GL")) %>%
                    filter(hgnc_symbol %in% rownames(exp_mat1)) %>%
                    arrange(hgnc_symbol, ensembl_gene_id) %>%
                    group_by(hgnc_symbol) %>% slice(1) %>% ungroup() %>%
                    mutate(chromosome_name = as.integer(chromosome_name)) %>%
                    arrange(chromosome_name, start_position) %>%
                    select(-ensembl_gene_id)

    write_tsv(genes_tbl, paste0(out_path, "infercnv-genes-",eval(s),".txt"), col_names = FALSE)


    if(length(sample_cells)>6500){
      subsets <- list()
      # divide data set in 2 groups:
      if(length(sample_cells)>13000){
        sets=3
        subsets[["set1"]] <- sample_cells[seq(1, length(sample_cells), 3)]
        subsets[["set2"]] <- sample_cells[seq(2, length(sample_cells), 3)]
        subsets[["set3"]] <- sample_cells[seq(3, length(sample_cells), 3)]
      } else {
        sets=2
        subsets[["set1"]] <- sample_cells[seq(1, length(sample_cells), 2)]
        subsets[["set2"]] <- sample_cells[seq(2, length(sample_cells), 2)]
      }

      for (i in names(subsets)){
          cnv_analysis_set <- c(subsets[[eval(i)]], subsample_control)

          exp_mat <- exp_mat1[,cnv_analysis_set]

          # import metadata
          metadata_tbl <- seurat_obj@meta.data %>% as_tibble(rownames="cell")
          # using sample names as classification
          metadata_tbl <- metadata_tbl %>% mutate(samples1=ifelse(str_detect(samples, "Control"), "Control", samples)) %>% select(cell, samples, samples1, broad_cell_type) %>% mutate(cluster = paste0(broad_cell_type, "_", samples1))

          if(str_detect(eval(s),"Control")){
              metadata_tbl$cluster <- ifelse(metadata_tbl$samples == eval(s), paste0(metadata_tbl$broad_cell_type, "_", metadata_tbl$samples), metadata_tbl$cluster)
          }
          
          metadata_tbl = metadata_tbl %>% filter(cell %in% colnames(exp_mat)) %>% select(cell, cluster) %>% arrange(cluster)

          # if there is only 1 cell of a given label, inferCNV fails
          if(length(names(which(table(metadata_tbl$cluster)==1)))>0 ){

                  # duplicate entry for given cell
                  rows_add <- metadata_tbl %>% filter(cluster %in% names(which(table(metadata_tbl$cluster)==1)))

                  # get expr data to add as well
                  exp_mat_add <- as.matrix(exp_mat[,rows_add$cell])

                  # rename 
                  rows_add$cell <- paste0(rows_add$cell, "_dup")
                  colnames(exp_mat_add) <- rows_add$cell

                  # add to original object
                  metadata_tbl <- rbind(metadata_tbl, rows_add)
                  exp_mat <- cbind(as.matrix(exp_mat), exp_mat_add)

          }
          
          write_tsv(metadata_tbl, paste0(out_path, "infercnv-annot-",eval(s),"_",i,".txt"), col_names = FALSE)

          # control group names:
          controls <- str_subset(unique(metadata_tbl$cluster), "_Control$")

          # create infercnv object
          infercnv_obj  <- infercnv::CreateInfercnvObject(
                              raw_counts_matrix = exp_mat,
                              gene_order_file = paste0(out_path,"infercnv-genes-",eval(s),".txt"), 
                              annotations_file = paste0(out_path, "infercnv-annot-",eval(s),"_",i,".txt"),
                              ref_group_names = controls
                              )

          saveRDS(infercnv_obj, paste0(out_path, "infercnv_obj_", eval(s),"-",i,".rds"))
    
      }

    } else {

      # import metadata
      metadata_tbl <- seurat_obj@meta.data %>% as_tibble(rownames="cell")
      metadata_tbl <- metadata_tbl %>% mutate(samples1=ifelse(str_detect(samples, "Control"), "Control", samples)) %>% select(cell, samples, samples1, broad_cell_type) %>% mutate(cluster = paste0(broad_cell_type, "_", samples1))

      if(str_detect(eval(s),"Control")){
          metadata_tbl$cluster <- ifelse(metadata_tbl$samples == eval(s), paste0(metadata_tbl$broad_cell_type, "_", metadata_tbl$samples), metadata_tbl$cluster)
      }
      
      metadata_tbl <- metadata_tbl %>% filter(cell %in% colnames(exp_mat1)) %>% select(cell, cluster) %>% arrange(cluster)

      # if there is only 1 cell of a given label, inferCNV fails
      if(length(names(which(table(metadata_tbl$cluster)==1)))>0 ){

              # duplicate entry for given cell
              rows_add <- metadata_tbl %>% filter(cluster %in% names(which(table(metadata_tbl$cluster)==1)))

              # get expr data to add as well
              exp_mat_add <- as.matrix(exp_mat1[,rows_add$cell])

              # rename 
              rows_add$cell <- paste0(rows_add$cell, "_dup")
              colnames(exp_mat_add) <- rows_add$cell

              # add to original object
              metadata_tbl <- rbind(metadata_tbl, rows_add)
              exp_mat1 <- cbind(as.matrix(exp_mat1), exp_mat_add)
      }
      
      write_tsv(metadata_tbl, paste0(out_path, "infercnv-annot-",eval(s),".txt"), col_names = FALSE)

      # control group names:
      controls <- str_subset(unique(metadata_tbl$cluster), "_Control$")

      # create infercnv object
      infercnv_obj <- infercnv::CreateInfercnvObject(
                          raw_counts_matrix = exp_mat1,
                          gene_order_file = paste0(out_path,"infercnv-genes-",eval(s),".txt"), 
                          annotations_file = paste0(out_path,"infercnv-annot-",eval(s),".txt"),
                          ref_group_names = controls
                         )
      
      saveRDS(infercnv_obj, paste0(out_path, "infercnv_obj_", eval(s),".rds"))
    }
}


# adjust function to run T/B/NK cells:
createInferCNVObject_call_BTNK <- function(seurat_obj=seurat_obj, s=broad_DE, ...){

    sample_cells = seurat_obj@meta.data %>% rownames_to_column("cell") %>% filter(broad_DE == eval(s)) %>% filter(!(samples %in% str_subset(unique(seurat_obj@meta.data$samples), "Control"))) %>% .$cell

    subsets <- list()
    # divide data set in 2 groups:
    if (length(sample_cells)>12000) {
        subsets[["set1"]] <- sample_cells[seq(1, length(sample_cells), 3)]
        subsets[["set2"]] <- sample_cells[seq(2, length(sample_cells), 3)]
        subsets[["set3"]] <- sample_cells[seq(3, length(sample_cells), 3)]

    } else if (length(sample_cells)>10000){
        subsets[["set1"]] <- sample_cells[seq(1, length(sample_cells), 2)]
        subsets[["set2"]] <- sample_cells[seq(2, length(sample_cells), 2)]
    }  else {
        subsets[["set1"]] <- sample_cells
    }


    control_cells=seurat_obj@meta.data %>% rownames_to_column("cell") %>% filter(broad_DE == eval(s)) %>% filter(samples %in% str_subset(unique(seurat_obj@meta.data$samples), "Control")) %>% .$cell
    
    
    # subsample control cells to 400 cells/cell type - 40 from each patient
    set.seed(4561)
    subsample_control <- seurat_obj@meta.data[control_cells,] %>% rownames_to_column("cell") %>% select(cell,Cell_type_identity, samples) %>% group_by(Cell_type_identity, samples) %>% sample_n(100, replace = TRUE) %>% .$cell %>% unique()
    
    cnv_analysis_set <- c(subsample_control,sample_cells)

    # extract counts from seurat
    exp_mat1 = GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
    exp_mat1 = exp_mat1[,cnv_analysis_set]

    # Ensembl 93 - July 2018 - used by Cell Ranger GRCh38 3.0.0
    ensembl = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
      dataset = "hsapiens_gene_ensembl", host = "jul2018.archive.ensembl.org")
    genes_tbl = biomaRt::getBM(attributes =
        c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position"),
        filters = "hgnc_symbol", values = rownames(exp_mat1), mart = ensembl)

    genes_tbl = genes_tbl %>%
      filter(!chromosome_name %in% c("MT", "X", "Y")) %>%
      filter(!str_detect(chromosome_name, "CHR_")) %>%
      filter(!str_detect(chromosome_name, "GL")) %>%
      filter(hgnc_symbol %in% rownames(exp_mat1)) %>%
      arrange(hgnc_symbol, ensembl_gene_id) %>%
      group_by(hgnc_symbol) %>% slice(1) %>% ungroup() %>%
      mutate(chromosome_name = as.integer(chromosome_name)) %>%
      arrange(chromosome_name, start_position) %>%
      select(-ensembl_gene_id)

    write_tsv(genes_tbl, paste0(out_path, "infercnv-genes-",eval(s),".txt"), col_names = FALSE)


    for (i in names(subsets)){
        cnv_analysis_set <- c(subsets[[eval(i)]], subsample_control)

        exp_mat = exp_mat1[,cnv_analysis_set]
        dim(exp_mat)

        # import metadata
        metadata_tbl = seurat_obj@meta.data %>% as_tibble(rownames="cell")
        # using sample names as classification
        #metadata_tbl = metadata_tbl %>% mutate(broad_cell_type1=gsub(" ","_", gsub("\\+","pos",broad_cell_type)) )
        metadata_tbl = metadata_tbl %>% mutate(samples1=ifelse(str_detect(samples, "Control"), "Control", samples)) %>% select(cell, samples, samples1, Cell_type_identity) %>% mutate(cluster = paste0(samples1, "_", Cell_type_identity)) %>% mutate(cluster=gsub(" ", "_", gsub("\\+", "_pos", gsub("-", "_", cluster))))

        # if(str_detect(eval(s),"Control")){
        #     metadata_tbl$cluster <- ifelse(metadata_tbl$samples == eval(s), paste0(metadata_tbl$Broad_cell_identity, "_", metadata_tbl$samples), metadata_tbl$cluster)
        # }
        
        metadata_tbl = metadata_tbl %>% filter(cell %in% colnames(exp_mat)) %>% select(cell, cluster) %>% arrange(cluster)

        if(length(names(which(table(metadata_tbl$cluster)==1)))>0 ){

                # duplicate entry for given cell
                rows_add <- metadata_tbl %>% filter(cluster %in% names(which(table(metadata_tbl$cluster)==1)))

                # get expr data to add as well
                exp_mat_add <- as.matrix(exp_mat[,rows_add$cell])

                # rename 
                rows_add$cell <- paste0(rows_add$cell, "_dup")
                colnames(exp_mat_add) <- rows_add$cell

                # add to original object
                metadata_tbl <- rbind(metadata_tbl, rows_add)
                exp_mat <- cbind(as.matrix(exp_mat), exp_mat_add)


        }
        
        metadata_tbl
        #table(metadata_tbl$cluster)
        write_tsv(metadata_tbl, paste0(out_path, "infercnv-annot-",eval(s),"_",i,".txt"), col_names = FALSE)

        # control group names:
        controls <- str_subset(unique(metadata_tbl$cluster), "Control_")

        # create infercnv object
        infercnv_obj =
        infercnv::CreateInfercnvObject(
            raw_counts_matrix = exp_mat,
            gene_order_file = paste0(out_path,"infercnv-genes-",eval(s),".txt"), 
            annotations_file = paste0(out_path, "infercnv-annot-",eval(s),"_",i,".txt"),
            ref_group_names = controls
        )

        saveRDS(infercnv_obj, paste0(out_path, "infercnv_obj_", eval(s),"-",i,".rds"))
   

    }

}
