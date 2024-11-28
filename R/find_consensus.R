find_consensus <- function(es.max, singler, sample_name){

  # rearrange sctype table
  sctype <- as.data.frame(t(es.max))
  Celltype <- colnames(sctype)[max.col(sctype)]
  sctype$scores <- apply(sctype, 1, max, na.rm=TRUE)
  sctype <- cbind(sctype, Celltype)
  sctype <- sctype[, -c(1:22)] # delete individual cell scores
  
  # find common ids
  ids_sctype <- rownames(sctype)
  ids_singler <- rownames(singler)
  common_ids <- intersect(ids_sctype, ids_singler)
  
  # subset the common cells
  sctype_subset <- sctype[common_ids, , drop = FALSE]  
  singler_subset <- singler[common_ids, , drop = FALSE]
  celltype_combined <- c()
  celltype_combined <- ifelse(sctype_subset$Celltype == singler_subset$pruned.labels, 
                              sctype_subset$Celltype, "none")
  
  # Create a new dataframe with the matched ids and the compared Celltype values
  consensus <- data.frame(
    CellID = common_ids,
    Celltype = celltype_combined,
    sctype = sctype_subset$Celltype,
    singler = singler_subset$pruned.labels
  )

  # replace NA values
  consensus <- consensus %>% mutate(Celltype = ifelse(!is.na(sctype) & is.na(singler), sctype, Celltype))
  consensus <- consensus %>% mutate(Celltype = ifelse(is.na(sctype) & !is.na(singler), singler, Celltype))
  
  # merge by cluster
  cL_resutls <- do.call("rbind", lapply(unique(cells@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(cells@meta.data[cells@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(cells@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  
  # cell clusters
  clusters <- cells@meta.data$seurat_clusters
  
  # list cluster types
  cluster_types <- sctype_scores$type
  
  # replace unknowns based on cluster type
  for (i in 1:length(cluster_types)) {
  
    cluster_index <- which(clusters == i, arr.ind = TRUE) # find the rows of cells that belong to cluster
    cluster <- sctype[cluster_index, ]
    cluster <- cluster[common_ids, , drop = FALSE] # cells that are not in the cluster are NA
    consensus <- cbind(consensus, cluster$Celltype) # add cluster types to consensus table
    names(consensus)[5] <- "cluster"
    consensus <- consensus %>% mutate(Celltype = ifelse(Celltype == "none" & !is.na(cluster), cluster_types[i], Celltype))
    consensus <- consensus[, -c(5)]
    
  }
  
  # Format table
  consensus <- subset(consensus, select = -sctype)
  consensus <- subset(consensus, select = -singler)
  consensus$SampleID <- sample_name
  consensus <- cbind(consensus %>% select(1, 3), consensus$Celltype)
  names(consensus)[3] <- "Celltype"
  
  consensus
  
}
