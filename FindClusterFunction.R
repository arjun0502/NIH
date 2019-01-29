#Arjun Jain 
#July 26, 2018
#FindCluster function 
#Arguments: a list of Seurat objects and ADT matrix which contains marker names
#Ouput: Returns a dataframe with the cluster number that represents the positive and negative populations for each antibody

FindClusters <- function(seurat_objects, ADT_matrix) {
  #Make vectors to contain positive indices and negative indices
  positive_indices <- vector()
  negative_indices <- vector()
  listofdfs <- list()
  
  #Set the identity class of all seurat objects to cluster identity
  for(i in seurat_objects) {
    i <- SetIdent(object = i, ident.use = (i@meta.data)$cluster.ident)
  }
  
  #Make vector with levels of "cluster.ident" factor
  cluster_levels <- as.numeric(levels(seurat_objects[[1]]@ident))
  
  #Make vector for surface markers
  surface_markers <- rownames(ADT_used)
  
  #For each marker 
  for(marker in surface_markers) {
    
    #make a vector to contain all rownames/cell names across seurat objects
    all_rownames <- vector()
    #make a vector to contain all expression values across seurat objects
    all_exp_values <- vector()
    #make a vector to contain cluster medians 
    cluster_medians <- vector()
    #For each seurat object
    for(j in seurat_objects) {
      #Subset the marker for each seurat object
      marker_subset <- data.frame(FetchData(object = j, vars.all = marker))
      #Get cell names of that subset
      dilution_rownames <- rownames(marker_subset)
      #Add cell names to vector
      all_rownames <- append(all_rownames, dilution_rownames)
      colnames(marker_subset) <- "col_name"
      #Get expression values of subset and add to new vector
      dilution_exp_values <- marker_subset$col_name
      all_exp_values <- append(all_exp_values, dilution_exp_values)
    }
    #add columns "all_rownames" and "all_exp_values" to combined data frame 
    combined_dataframe <- data.frame(all_rownames, all_exp_values)
    #Find number of cells/rows
    length <- as.numeric(length(combined_dataframe$all_rownames))
    #Make a vector that is equal size to number of rows in combined dataframe 
    #Will fill to contain different cluster numbers at corresponding indices
    cluster_num <- vector(mode = "integer", length = length)
    #For each cluster.........
    for(cluster in cluster_levels) {
      #Have vector for cells that match 
      match_cluster_indices <- vector()
      #Vector containing cells for each cluster 
      cluster_n <- vector()
      #For each seurat object 
      for(k in seurat_objects) {
        #Subset a vector for just that cluster 
        cells.to.include <- WhichCells(object = k, ident = cluster)
        #Add to vector 
        cluster_n <- append(cluster_n, cells.to.include)
      }
      #Get rid of duplicates (cluster_n now contains cells from all seurat objects)
      cluster_n <- unique(cluster_n)
      #Find index where cells match up for each cluster and then add to a vector
      for(cell_name in combined_dataframe$all_rownames) {
        if(cell_name %in% cluster_n) {
          index <- which(combined_dataframe$all_rownames == cell_name)
          match_cluster_indices <- append(match_cluster_indices, index)
        }
      }
      #A vector called cluster_num used to insert cluster numbers at specified indices 
      for(num in match_cluster_indices) {
        cluster_num[num] <- cluster
      }
    }
    #Out of for loop, now make cluster number a new column 
    combined_dataframe <- cbind(combined_dataframe, cluster_num)
    
    #Find medians of each cluster and add to a vector "cluster_medians"
    #The length of cluster medians is equal to number of clusters 
    for(cluster in cluster_levels) {
      temp <- combined_dataframe[combined_dataframe$cluster_num == cluster,]
      clustern <- temp$all_exp_values
      clustern_median <- median(clustern)
      cluster_medians <- append(cluster_medians, clustern_median)
    }
    #Get cluster number for positive and negative population 
    positive_index <- as.numeric(which.max(cluster_medians) - 1)
    positive_indices <- append(positive_indices, positive_index)
    negative_index <- as.numeric(which.min(cluster_medians) -1 )
    negative_indices <- append(negative_indices, negative_index)
  }
  #Make a dataframe with surface markers and positive population and negative population cluster numbers 
  cluster_final_indices <- data.frame(surface_markers, positive_indices, negative_indices)
}