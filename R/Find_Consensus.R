Find_Consensus <- function(sctype, singler){

  ids_sctype <- rownames(sctype)
  ids_single2 <- rownames(singler)
  
  common_ids <- intersect(ids_sctype, ids_single2)
  
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

  consensus
}
