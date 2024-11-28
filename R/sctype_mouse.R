mouse_sctype <- function(path_to_file){
  
  cells.data <- Read10X(data.dir = path_to_file)
  cells <- CreateSeuratObject(counts = cells.data, project = "cells_mouse", min.cells = 3, min.features = 200)
  
  # normalize data
  cells[["percent.mt"]] <- PercentageFeatureSet(cells, pattern = "^MT-")
  cells <- subset(cells, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) # filter based on QC metrics  
  cells <- NormalizeData(cells, normalization.method = "LogNormalize", scale.factor = 10000)
  cells <- FindVariableFeatures(cells, selection.method = "vst", nfeatures = 2000)
  
  # scale and run PCA
  cells <- ScaleData(cells, features = rownames(cells))
  cells <- RunPCA(cells, features = VariableFeatures(object = cells))
  
  # cluster
  cells <- FindNeighbors(cells, dims = 1:10)
  cells <- FindClusters(cells, resolution = 0.8)
  cells <- RunUMAP(cells, dims = 1:10)
  
  ## Annotate cell types
  
  # load gene set preparation function
  source("https://raw.githubusercontent.com/natacow/sc-type/master/R/gene_sets_prepare.R")
  # load cell type annotation function
  source("https://raw.githubusercontent.com/natacow/sc-type/master/R/sctype_score_.R")
  
  # DB file
  db_ <- "https://raw.githubusercontent.com/natacow/sc-type/master/mouse_sctype_markers.xlsx"
  tissue <- "Lung" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
  
  # prepare gene sets
  gene_list <- gene_sets_prepare(db_, tissue)
  
  # check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
  seurat_package_v5 <- isFALSE('counts' %in% names(attributes(pbmc[["RNA"]])));
  
  # extract scaled scRNA-seq matrix
  scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(pbmc[["RNA"]]$scale.data) else as.matrix(pbmc[["RNA"]]@scale.data)
  
  # run ScType
  es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gene_list$gs_positive, gs2 = gene_list$gs_negative)
   
  es.max # sctype-sorted cells
  
}
