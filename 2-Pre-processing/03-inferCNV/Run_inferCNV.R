args <- commandArgs(trailingOnly = TRUE)

if (length(args)==2){
    sample <- as.character(args[1])
    num_threads <- as.numeric(args[2])
    out_path <- "/gpfs/data/aifantislab/home/nadorb01/projects/ad_ped_AML_v5/results/infercnv-broad-cell-type-per-patient-res.2/"
} else if(length(args)==3){
    sample <- as.character(args[1])
    num_threads <- as.numeric(args[2])
    out_path <- as.character(args[3])
    } else {
        stop ("Please provide the following arguments: sample_name, num_threads, (out_path, default set to 'ad_ped_AML_v5/results/infercnv-broad-cell-type-per-patient/')")
    }

library(infercnv)

infercnv_obj <- readRDS(paste0(out_path, "infercnv_obj_", sample,".rds"))# cutoff=0.1 works well for 10x Genomics
# analysis_mode 'subclusters' partition cells into groups having consistent patterns of CNV
out <- paste0(out_path,"infercnv-out-",sample, "/")
infercnv_obj <- infercnv::run(
  infercnv_obj = infercnv_obj,
  cutoff = 0.1, min_cells_per_gene = 10,
  out_dir = out,
  num_threads = num_threads,
  analysis_mode = "subclusters",
  denoise = TRUE,
  HMM = TRUE,
  plot_steps = FALSE, no_prelim_plot = TRUE, png_res = 300
)

# writes metadata in the output dir
infercnv::add_to_seurat(infercnv_output_path = paste0(out_path, "/infercnv-out-", sample))
saveRDS(infercnv_obj, paste0(out_path, "infercnv_obj_", sample,".rds"))
