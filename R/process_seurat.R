#' Basic processing of scRNA-Seq data in Seurat objects
#'
#' This function normalizes, scales data, and creates QC plots. Run this downstream of UMI and gene filtering + doublet detection
#'
#' @importFrom utils head write.csv write.table
#' @importFrom grDevices pdf dev.off
#' @param obj A seurat object
#' @param num_most_variable_features (Optional) input to FindVariableFeatures
#' @param logfile File to write console output to
#' @param outdir Output directory
#' @param vars_to_color Columns of the Seurat metadata to color the UMAP by
#' @param image_save_mode (Optional) Boolean indicating whether to save images
#' @param table_save_mode (Optional) Boolean indicating whether to save tables
#'
#' @returns A processed seurat object
#' @export
#'
process_seurat <- function(obj, num_most_variable_features = 2000, logfile, outdir, vars_to_color = c("seurat_clusters"),
                           image_save_mode = T, table_save_mode = T) {

  logger::log_info("Starting the pipeline. Check logfile at ", logfile, " for more info")

  # Writing log to logfile
  sink(logfile)

  # Normalize
  obj <- Seurat::NormalizeData(obj)

  # Finding highest variance genes
  obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = num_most_variable_features)

  # Plotting highly variable genes
  var_features <- Seurat::VariableFeatures(obj)
  top10 <- head(var_features, 10)

  plot1 <- Seurat::VariableFeaturePlot(obj)
  plot2 <- Seurat::LabelPoints(plot = plot1, points = top10, repel = TRUE)

  # Save
  if (image_save_mode) {
    ggplot2::ggsave(filename = file.path(outdir, "highly_variable_genes.pdf"), plot = plot2, width = 6, height = 6)
  }

  # Checking composition of var features
  var_mt <- var_features[grepl("^MT-", var_features)]
  var_rb <- var_features[grepl("^RP[SL]", var_features)]
  var_ig <- var_features[grepl("^IG[HKL]", var_features)]

  message((length(var_mt) / length(var_features)) * 100, "%", " of the top 2000 most variably expressed genes are MT genes")
  message((length(var_rb) / length(var_features)) * 100, "%", " of the top 2000 most variably expressed genes are Ribo genes")
  message((length(var_ig) / length(var_features)) * 100, "%", " of the top 2000 most variably expressed genes are IG genes")

  # Saving list of variable features
  if (table_save_mode) {
    write.table(var_features, file = file.path(outdir, "top2000_variable_genes_list.txt"),
                col.names = F, row.names = F, quote = F, sep = "\t")
  }

  # Scaling
  obj <- Seurat::ScaleData(obj)

  # PCA
  obj <- Seurat::RunPCA(obj)

  # Saving PCA tables
  if (table_save_mode) {
    write.csv(obj@reductions$pca@cell.embeddings,
              file = file.path(outdir, "pca_cell_embeddings_top2000_genes.csv"), quote = F)
    write.csv(obj@reductions$pca@feature.loadings,
              file = file.path(outdir, "pca_feature_loadings_top2000_genes.csv"), quote = F)
  }

  # Saving the loading plots
  if (image_save_mode) {
    pdf(file = file.path(outdir, "pc1_2_loadings.pdf"), width = 8, height = 6)
    print(plot1 <- Seurat::VizDimLoadings(obj, dims = 1:2, reduction = "pca"))
    dev.off()
    pdf(file = file.path(outdir, "pc3_4_loadings.pdf"), width = 8, height = 6)
    print(plot1 <- Seurat::VizDimLoadings(obj, dims = 3:4, reduction = "pca"))
    dev.off()
  }

  # Getting perc_variance per PC
  # Calculating variance
  pc_stdev <- obj@reductions$pca@stdev
  pc_var <- pc_stdev^2
  percent_var <- (pc_var / sum(pc_var)) * 100

  # Putting in to a df
  pca_var_df <- data.frame(
    PC = 1:length(percent_var),
    Variance = pc_var,
    PercentVariance = percent_var
  )

  # Saving PCA related output
  if (image_save_mode) {
    pdf(file = file.path(outdir, "pca_elbow_plot.pdf"), width = 4, height = 4)
    print(Seurat::ElbowPlot(obj))
    dev.off()
    pdf(file = file.path(outdir, "pca_plot_pc1_2.pdf"), width = 4, height = 4)
    print(Seurat::DimPlot(obj, reduction = "pca", dims = c(1,2)) + Seurat::NoLegend())
    dev.off()

    pdf(file = file.path(outdir, "dim_heatmap_pc1.pdf"), width = 12, height = 11)
    print(Seurat::DimHeatmap(obj, dims = 1, cells = 500, balanced = TRUE, nfeatures = 40))
    dev.off()
    pdf(file = file.path(outdir, "dim_heatmap_pc2.pdf"), width = 12, height = 11)
    print(Seurat::DimHeatmap(obj, dims = 2, cells = 500, balanced = TRUE, nfeatures = 40))
    dev.off()
  }
  if (table_save_mode) {
    write.csv(pca_var_df, file = file.path(outdir, "pca_explained_variance.csv"), row.names = F,
              quote = F)
  }

  sink() # Bring output back to console

  # Getting user input for Clustering PCs
  logger::log_info("Please check the Elbow plot written to ", file.path(outdir, "pca_elbow_plot.pdf"),
           " to decide how many PCs you would like to use for clustering and further downstream analyses")

  user_input <- readline(prompt = "Enter your chosen number of PCs to proceed: ")
  user_input <- as.numeric(user_input)
  while (is.na(user_input)) {
    user_input <- readline(prompt = "Your entry was invalid. Please enter a valid number: ")
    user_input <- as.numeric(user_input)
  }

  sink(logfile) # Send output back to the logfile

  # Clustering
  obj <- Seurat::FindNeighbors(obj, dims = 1:user_input)
  obj <- Seurat::FindClusters(obj, resolution = 0.5)

  sink() # Bring output back to console

  # Getting user input for dims for UMAP
  user_input2 <- readline(prompt = "Enter your chosen number dims to use for UMAP: ")
  user_input2 <- as.numeric(user_input2)
  while (is.na(user_input2)) {
    user_input2 <- readline(prompt = "Your entry was invalid. Please enter a valid number: ")
    user_input2 <- as.numeric(user_input2)
  }

  sink(logfile) # Send output back to logfile

  # UMAP
  obj <- Seurat::RunUMAP(obj, dims = 1:user_input2)

  # UMAP plots
  if (image_save_mode) {
    seed_colors <- c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46",
                     "#008941", "#006FA6", "#A30059", "#FFDBE5", "#7A4900",
                     "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87")
    for (name in vars_to_color) {
      cluster_colors <- grDevices::colorRampPalette(seed_colors)(length(unique(obj@meta.data[[name]])))
      pdf(file = paste0(outdir, "/umap_plot_by", name, ".pdf"), width = 14, height = 11)
      print(Seurat::DimPlot(obj, reduction = "umap", cols = cluster_colors, group.by = name))
      dev.off()
    }
  }

  sink()

  logger::log_info("Processing Complete. Check: ", outdir, " for all saved outputs")

  return(obj)
}
