## SHUFFLE ###
shuf <- function(df){
  return(df[sample(1:dim(df)[1], dim(df)[1]),])
}

## GO ####
# make function to run GO on DE gene list seurat output:
run_GO <- function(s_obj, DE_genes, s_name,... ){
    out <- paste0(temp_dir,"GO_", s_name,"/")
    dir.create(out)

    # define background
    universe <- rownames(s_obj@assays$RNA@counts)[which(rowSums(as.matrix(s_obj@assays$RNA@counts))>0)]
    # convert gene names to entrez gene id
    universe_entrez <- mapIds(org.Hs.eg.db, universe, 'ENTREZID', 'SYMBOL')
    universe_entrez <- universe_entrez[!is.na(universe_entrez)]

    gene_list <- list()
    up_genes <- DE_genes %>% dplyr::filter(avg_logFC>0)
    up_genes_entrez <- mapIds(org.Hs.eg.db, rownames(up_genes), 'ENTREZID', 'SYMBOL')
    gene_list[["up_genes"]] <- up_genes_entrez[!is.na(up_genes_entrez)]

    down_genes <- DE_genes %>% dplyr::filter(avg_logFC<0)
    down_genes_entrez <- mapIds(org.Hs.eg.db, rownames(down_genes), 'ENTREZID', 'SYMBOL')
    gene_list[["down_genes"]] <- down_genes_entrez[!is.na(down_genes_entrez)]

    GO_list <- list()
    for (i in names(gene_list)) {
        for (o in c("BP", "CC", "MF")){

        GO_list[[paste0(eval(i),"_",eval(o))]] <- enrichGO(gene=gene_list[[eval(i)]],OrgDb = org.Hs.eg.db, ont=eval(o), universe=universe_entrez, pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
        write.csv(GO_list[[paste0(eval(i),"_",eval(o))]], paste0(out, "GO_",eval(i),"_",eval(o),"_",s_name,"_before_filtering.csv"))

        GO_list[[paste0(eval(i),"_",eval(o))]] <- simplify(GO_list[[paste0(eval(i),"_",eval(o))]], cutoff=0.7, by="p.adjust", select_fun=min)

        write.csv(GO_list[[paste0(eval(i),"_",eval(o))]], paste0(out, "GO_",eval(i),"_",eval(o),"_",s_name,".csv"))
        
        pdf(paste0(out, "Dotplot_top_30_",eval(i),"_",eval(o),"_",s_name,".pdf"), onefile=FALSE, useDingbats=FALSE, height=6, width=10)
        print(dotplot(GO_list[[paste0(eval(i),"_",eval(o))]], showCategory=30))
        dev.off()
        }
    }

  return(GO_list)
}

## GSEA ####
# run GSEA on seurat gene list output
run_GSEA <- function(DE_gene_table, s_name){
    m_df = msigdbr::msigdbr(species = "Homo sapiens")
    m_t2g = m_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
    m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

    # make sub-folders for each guide
    out <- paste0(temp_dir,"GSEA_", s_name,"/")
    dir.create(out)
    # make sorted data.frame
    df1 <- DE_gene_table %>% rownames_to_column("gene")
    df_rank <- rank(df1[,"avg_logFC"], na.last=T, ties.method="average")
    df_rank_lh <- rank(-(df1[,"avg_logFC"]), na.last=T, ties.method="average")
    names(df_rank) <- df1[,"gene"]
    names(df_rank_lh) <- df1[,"gene"]

    print(paste("Running fgsea for", eval(s_name)))
    fgsea_lst <- list()
    fgsea_lst[["enr_data_fgsea"]] <- fgsea::fgsea(pathways = m_list, stats = df_rank,minSize=10,maxSize=500,nperm=10000)
    fgsea_lst[["depl_data_fgsea"]] <- fgsea::fgsea(pathways = m_list, stats = df_rank_lh,minSize=10,maxSize=500,nperm=10000)
    topPathwaysUp <- fgsea_lst[["enr_data_fgsea"]][ES > 0, ][head(order(pval), n=15), pathway]
    topPathwaysDown <- fgsea_lst[["depl_data_fgsea"]][ES > 0, ][head(order(pval), n=15), pathway]
    topPathways <- c(topPathwaysUp, topPathwaysDown)
    
    #plot top pathways
    pdf(paste0(out, "top_pathways_up_", s_name,".pdf"), onefile=F, useDingbats=F, width=14, height=8)
    print(fgsea::plotGseaTable(m_list[topPathwaysUp], df_rank, fgsea_lst[["enr_data_fgsea"]] , 
              gseaParam = 0.5))
    dev.off()

    pdf(paste0(out, "top_pathways_down_", s_name,".pdf"), onefile=F, useDingbats=F, width=14, height=8)
    print(fgsea::plotGseaTable(m_list[topPathwaysDown], df_rank_lh, fgsea_lst[["depl_data_fgsea"]], 
              gseaParam = 0.5))
    dev.off()

    data.table::fwrite(fgsea_lst[["enr_data_fgsea"]], file=paste0(out,"fGSEA_enrichment_", s_name, ".txt"), sep="\t", sep2=c("", " ", ""))
    data.table::fwrite(fgsea_lst[["depl_data_fgsea"]], file=paste0(out,"fGSEA_depletion_", s_name, ".txt"), sep="\t", sep2=c("", " ", ""))

    return(fgsea_lst)
}
