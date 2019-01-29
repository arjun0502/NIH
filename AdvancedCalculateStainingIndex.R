#Arjun Jain 
#July 26, 2018
#AdvancedCalculateStainingIndex function 
#Purpose: Calculates the staining index for each antibody in a CITE-seq experiment
#Staining index measures the effectiveness of an antibody

#Precondition: Run FindCluster() on seurat objects
#Arguments: A Seurat object, an ADT matrix, output of FindClusters()
##Output: Data frame with calculated staining indexes for each antibody


AdvancedCalculateStainingIndex <- function(seurat_object, ADT_matrix, pos_neg_clusters) {
  #Vectors to contain staining indexes for the seurat object across the antibodies in the ADT matrix 
  staining_indexes <- vector(mode = "numeric")
  
  #Vectors to contain positive and negative cluster numbers (nos)
  positive_cluster_nos <- vector(mode = "integer")
  negative_cluster_nos <- vector(mode = "integer")
  m <- CreateCounter(0)
  surface_markers <- rownames(ADT_matrix)
  for(marker in surface_markers) {
    m$increment(1)
    #Access the positive and negative population cluster numbers to subset seurat object on 
    cluster_pos <- pos_neg_clusters$positive_indices[m$value()]
    positive_cluster_nos <- append(positive_cluster_nos, cluster_pos)
    
    cluster_neg <- pos_neg_clusters$negative_indices[m$value()]
    negative_cluster_nos <- append(negative_cluster_nos, cluster_neg)
    
    #Subset data for specific marker 
    specific_marker <- data.frame(FetchData(object = seurat_object, vars.all = marker))
    
    #Subset the positive and negative clusters 
    cells.to.include_pos <- WhichCells(object = seurat_object, ident = cluster_pos)
    pos_cluster_data <- specific_marker[cells.to.include_pos, , drop = FALSE]
    
    cells.to.include_neg <- WhichCells(object = seurat_object, ident = cluster_neg)
    neg_cluster_data <- specific_marker[cells.to.include_neg, , drop = FALSE]
    
    colnames(pos_cluster_data) <- "col_name"
    colnames(neg_cluster_data) <- "col_name"
    #Calculate staining index
    neg_median <- median(neg_cluster_data$col_name)
    pos_median <- median(pos_cluster_data$col_name)
    neg_sd <- sd(neg_cluster_data$col_name)
    staining_index = (pos_median - neg_median)/(2*neg_sd)
    staining_indexes <- append(staining_indexes, staining_index)
  }
  
  staining_indexResult <- data.frame(surface_markers, staining_indexes, negative_cluster_nos, positive_cluster_nos, stringsAsFactors = F)
  staining_indexResult
  
}

