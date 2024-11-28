
singler_mouse <- function(path_to_file){
    #### annotating cells against built-in references ####
  library(SingleR) # accepts any log-normalized expression matrix
  library(celldex) # provides access to bulk RNA-seq datasets
  library(Seurat)
  mouseRNA <- MouseRNAseqData() # summarizedExperiment object
  
  
  # label cells
  library(scRNAseq) # test datasets 
  library(XML)
  mSCs <- Read10X(data.dir = path_to_file)
  
  pred.hesc <- SingleR(test = mSCs, ref = mouseRNA, assay.type.test=1,
                       labels = mouseRNA$label.main) # using mouseRNA to annotate each cell in mSCs
  pred.hesc # each row is a prediction result for a single cell
  
  pred.msc <- pred.hesc[,c(1, 4)]
  names(pred.msc)[1] <- 'scores'
  names(pred.msc)[2] <- 'Celltype'
  
  pred.msc # cells with assigned types
}
