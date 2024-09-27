ids_sctype <- rownames(es.max)
ids_single2 <- rownames(pred.msc)

common_ids <- intersect(ids_sctype, ids_single2)

sctype_subset <- es.max[common_ids, , drop = FALSE]  
single2_subset <- pred.msc[common_ids, , drop = FALSE]

celltype_combined <- c()

celltype_combined <- ifelse(sctype_subset$Celltype == single2_subset$Celltype, 
                            sctype_subset$Celltype, "none")

# Create a new dataframe with the matched ids and the compared Celltype values
consensus <- data.frame(
  CellID = common_ids,
  Celltype = celltype_combined
)

consensus <- result[rowSums(is.na(result)) != ncol(result), ] # final consensus table

