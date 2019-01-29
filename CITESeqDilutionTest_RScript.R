#Arjun Jain 
#July 16, 2018 - July 26, 2018

#I conducted an experiment to determine the lowest amount of various antibodies which can effectively separate cell populations
#In this experiment, CITE-seq staining was run with 8 million healthy PBMCs (peripheral blood mononuclear cells), with 4 pairs of
#hashtag-oligonucleotides (HTOs) to label 4 different antibody amounts.

#Each dilution of antibodies is represented by a pair of HTOs (hashtag-oligos)
  #HTO 1,2 -> 1ug
  #HTO 3,4 -> 0.5ug
  #HTO 5,6 -> 0.25ug
  #HTO 7,8 -> 0.125ug

#For this analysis, I used Seurat which is an R package designed for QC, analysis, and
#exploration of multi-modal data. 

#I also independently created two functions: "AdvancedCalculatingStainingIndex.R" and "FindClusterFunction.R" 

#The purpose of this R script: Analysis of CITE-seq antibody dilution data with Seurat package

#____________________________________________________________________________________________________________
#Load Seurat package
library(Seurat)

#_____________________________________________________________________________________________________________
#(Step 1) Creating a Seurat object with RNA, ADT (protein expression), and HTO data in "raw.data" slot

##Load in the RNA UMI expression matrix 
##Rownames: Genes
##Column names: Cells (cell barcodes, 18698 cells)
expression.matrix = Read10X("/hpcdata/sg/sg_data/users/jaina10/jaina10/CITESeqDilutionTest/CITEseqDilutionTestData/20180616_CITEseqdilutionTest/Lane1cDNA/outs/filtered_gene_bc_matrices/hg19")

##Load in ADT UMI matrix (protein expression )
##Rownames: Proteins (31)
##Column names: Cells (as cell barcodes)
ADT_result <- read.csv("/hpcdata/sg/sg_data/users/jaina10/jaina10/CITESeqDilutionTest/CITEseqDilutionTestData/20180616_CITEseqdilutionTest/CITEseqCountOut/ADT_Result_WL.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

##Cleaning up ADT UMI matrix 
###Rownames -> Proteins  
ADT_used = ADT_result[c(1:31),]
ADT_used = ADT_used[,order(colnames(ADT_used))]

##Load in HTO matrix (cell barcodes as columns, HTO as rows)
HTO_result <- read.csv("/hpcdata/sg/sg_data/users/jaina10/jaina10/CITESeqDilutionTest/CITEseqDilutionTestData/20180616_CITEseqdilutionTest/CITEseqCountOut/HTO_Result_WL.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

##Clean up HTO matrix  
HTO_used = HTO_result[c(1,4:11),]
HTO_used = HTO_used[,order(colnames(HTO_used))]

##Creating a Seurat object with raw expression data 
pbmc.dilution <- CreateSeuratObject(raw.data = expression.matrix)

##Performing log normalization on RNA data
pbmc.dilution <- NormalizeData(pbmc.dilution)

##Store CITE (ADT) data in raw.data slot of "CITE" assay slot of pbmc.dilution Seurat object
pbmc.dilution <- SetAssayData(pbmc.dilution, assay.type = "CITE", slot = "raw.data", new.data = ADT_used)

##Storing HTO data in raw.data slot of "HTO" assay slot of pbmc.dilution Seurat Object
pbmc.dilution <- SetAssayData(pbmc.dilution, assay.type = "HTO", slot = "raw.data", new.data = HTO_used)

#_________________________________________________________________________________________________________
#(Step 2) Demulitplexing the Seurat object based on HTOs and visualization.
##Also, repeating the preprocessing  (normalization and scaling) steps that we typically run with RNA, but modifying the 'assay.type' argument.  
#For CITE-seq data and HTO data, we use a centered log-ratio (CLR) normalization, computed independently for each gene. 

pbmc.dilution <- NormalizeData(pbmc.dilution, assay.type = "CITE", normalization.method = "genesCLR")
pbmc.dilution <- ScaleData(pbmc.dilution, assay.type = "CITE", display.progress = FALSE)
pbmc.dilution <- NormalizeData(pbmc.dilution, assay.type = "HTO", normalization.method = "genesCLR")
pbmc.dilution <- ScaleData(pbmc.dilution, assay.type = "HTO", display.progress = FALSE)

# Demultiplex based on HTOs
pbmc.dilution_demux = HTODemux(pbmc.dilution, assay.type = "HTO", positive_quantile = 0.99,init_centers = NULL, cluster_nstarts = 100, k_function = "clara", nsamples = 100, print.output = TRUE)

# Visualize HTO Di-multiplexing results using HTOs
HTOHeatmap(pbmc.dilution_demux, hto.classification = "hto_classification",global.classification = "hto_classification_global", assay.type = "HTO",num.cells = 5000, singlet.names = NULL)

# Subset dilutions
OneUG = SubsetData(pbmc.dilution_demux, ident.use = c("HTO_1", "HTO_2"))
p5UG = SubsetData(pbmc.dilution_demux, ident.use = c("HTO_3", "HTO_4"))
p25UG = SubsetData(pbmc.dilution_demux, ident.use = c("HTO_5", "HTO_6"))
p125UG = SubsetData(pbmc.dilution_demux, ident.use = c("HTO_7", "HTO_8"))

# Subset data without doublets, negative
#5044 doublets, 469 negatives, 13185 singlets 
demux_cleaned = SubsetData(pbmc.dilution_demux, ident.use = c("HTO_1", "HTO_2","HTO_3", "HTO_4","HTO_5", "HTO_6","HTO_7", "HTO_8"))

#Check Out Antibody results

# RidgePlot(pbmc.dilution_demux, features.plot = rownames(ADT_used)[1:9],nCol = 3)
# RidgePlot(pbmc.dilution_demux, features.plot = rownames(ADT_used)[10:18], nCol = 3)
# RidgePlot(pbmc.dilution_demux, features.plot = rownames(ADT_used)[19:27], nCol = 3)
# RidgePlot(pbmc.dilution_demux, features.plot = rownames(ADT_used)[28:31], nCol = 3)
# 
# GenePlot(OneUG, gene1 = "CCR4_PROT", gene2 = "CCR7_PROT", cex = 0.5)
# GenePlot(OneUG, gene1 = "CD4_PROT", gene2 = "CD8_PROT", cex = 0.5)
# GenePlot(OneUG, gene1 = "CD14_PROT", gene2 = "CD8_PROT", cex = 0.5)
# GenePlot(OneUG, gene1 = "ICOS_PROT", gene2 = "CD8_PROT", cex = 0.5)
# GenePlot(OneUG, gene1 = "IgM_PROT", gene2 = "IgD_PROT", cex = 0.5)

#_____________________________________________________________________________________________________
#(Step 3) Clustering cells based on RNA expression and visualization and finding corresponding cell subtypes 

# choose ~1k variable genes
demux_cleaned <- FindVariableGenes(demux_cleaned, do.plot = FALSE, y.cutoff = 0.5)

# standard scaling (no regression)
demux_cleaned <- ScaleData(demux_cleaned, display.progress = FALSE)

# Run PCA, select 13 PCs for tSNE visualization and graph-based clustering
demux_cleaned <- RunPCA(demux_cleaned, pcs.print = 0)
#PCElbowPlot(demux_cleaned)

demux_cleaned <- FindClusters(demux_cleaned, dims.use = 1:9, print.output = FALSE)

demux_cleaned <- RunTSNE(demux_cleaned, dims.use = 1:9)

# Find the markers that define each cluster, and use these to annotate the
# clusters, we use max.cells.per.ident to speed up the process (cluster based on expression)
demux_cleaned.rna.markers <- FindAllMarkers(demux_cleaned, max.cells.per.ident = 100, logfc.threshold = log(2), only.pos = TRUE, min.diff.pct = 0.3, do.print = F)

#Identify cell subsets by comparing Seurat Multimodal Vignette RNA clustering results to clusters (only numbers, not cell names) we got using RNA
#Looked at genes that were found in each cluster and saw in which cluster of Vignette clustering were they expressed
#Return a list with each component being a list of the clusters (on left) and their abundance (on right)
#Function GetCellTypes_RNA is also a separate RNA script in "CITESeqDilutionTest" folder 
GetCellTypes_RNA <- function(our_markers, seurat_markers) {
  list <- list()
  cluster_range <- as.numeric(as.character(unique(our_markers$cluster)))
  index <- CreateCounter(0)
  for(clust in cluster_range) {
    index$increment(1)
    #Subset our_markers based on specific cluster number 
    cluster_subset <- our_markers[our_markers$cluster == clust, ]
    cluster_vector <- vector()
    rownames_vector <- vector()
    for(gene in cluster_subset$gene) {
      for(row in rownames(seurat_markers)) {
        if(grepl(gene, row)) {
          row_num <- which(rownames(seurat_markers) == row)
          cluster_num <- as.numeric(as.character(seurat_markers[[row_num, 6]]))
          cluster_vector <- append(cluster_vector, cluster_num)
          #rownames_vector <- append(rownames_vector, row)
        }
      }
    }  
    ux <- unique(cluster_vector) 
    a <- table(cluster_vector)
    list[[index$value()]] <- a
  }
  
  return(list)
}  

shared_clusters_RNA <- GetCellTypes_RNA(demux_cleaned.rna.markers, cbmc.rna.markers)

#mRNA vs Protein plots 
par(mfrow = c(2,2))
GenePlot(pbmc.dilution_demux, gene1 = "CCR4_PROT", gene2 = "CCR4", cex.use = 0.5)
GenePlot(demux_cleaned, gene1 = "CCR7_PROT", gene2 = "CCR7", cex.use = 0.5)
GenePlot(demux_cleaned, gene1 = "CD10_PROT", gene2 = "MME", cex.use = 0.5)
GenePlot(demux_cleaned, gene1 = "CD11c_PROT", gene2 = "ITGAX", cex.use = 0.5)
par(mfrow = c(2,2))
GenePlot(demux_cleaned, gene1 = "CD138_PROT", gene2 = "SDC1", cex.use = 0.5)
GenePlot(demux_cleaned, gene1 = "CD14_PROT", gene2 = "CD14", cex.use = 0.5)
GenePlot(demux_cleaned, gene1 = "CD16_PROT", gene2 = "FCGR3A", cex.use = 0.5)
GenePlot(demux_cleaned, gene1 = "CD185_PROT", gene2 = "CXCR5", cex.use = 0.5)
par(mfrow = c(2,2))
GenePlot(demux_cleaned, gene1 = "CD195_PROT", gene2 = "CCR5", cex.use = 0.5)
GenePlot(demux_cleaned, gene1 = "CD196_PROT", gene2 = "CCR6", cex.use = 0.5)
GenePlot(demux_cleaned, gene1 = "CD20_PROT", gene2 = "MS4A1", cex.use = .5)
GenePlot(demux_cleaned, gene1 = "CD25_PROT", gene2 = "IL2RA", cex.use = .5)
par(mfrow = c(2,2))
GenePlot(demux_cleaned, gene1 = "CD27_PROT", gene2 = "CD27", cex.use = .5)
GenePlot(demux_cleaned, gene1 = "CD31_PROT", gene2 = "VEGFA", cex.use = .5)
GenePlot(demux_cleaned, gene1 = "CD34_PROT", gene2 = "CD34", cex.use = .5)
GenePlot(demux_cleaned, gene1 = "CD38_PROT", gene2 = "CD38", cex.use = .5)
par(mfrow = c(2,2))
GenePlot(demux_cleaned, gene1 = "CD45RA_PROT", gene2 = "PTPRC", cex.use = .5)
GenePlot(demux_cleaned, gene1 = "CD4_PROT", gene2 = "CD4", cex.use = 0.5)
GenePlot(demux_cleaned, gene1 = "CD56_PROT", gene2 = "NCAM1", use.scaled = TRUE, cex.use = 0.5)
GenePlot(demux_cleaned, gene1 = "CD8_PROT", gene2 = "CD8A", cex.use = 0.5)
par(mfrow = c(2,2))
GenePlot(demux_cleaned, gene1 = "CTLA4_PROT", gene2 = "CTLA4", cex.use = 0.5)
GenePlot(demux_cleaned, gene1 = "CXCR3_PROT", gene2 = "CXCR3", cex.use = 0.5)
GenePlot(demux_cleaned, gene1 = "ICOS_PROT", gene2 = "ICOS", cex.use = 0.5)
GenePlot(demux_cleaned, gene1 = "IgM_PROT", gene2 = "CD40LG", cex.use = 0.5)
par(mfrow = c(1,1))
GenePlot(demux_cleaned, gene1 = "CD8_PROT", gene2 = "CD8A", cex.use = 0.5)
GenePlot(demux_cleaned, gene1 = "NKG2D_PROT", gene2 = "KLRK1", cex.use = 0.5)
GenePlot(demux_cleaned, gene1 = "PD-1_PROT", gene2 = "PDCD1", cex.use = 0.5)

TSNEPlot(demux_cleaned, do.label = TRUE, pt.size = 0.5)
#using scaled data
# DoHeatmap(demux_cleaned, assay.type = "CITE", 
#           slim.col.label = TRUE, remove.key = TRUE, group.label.rot = TRUE)
# 
# #Not using scaled data
DoHeatmap(demux_cleaned, assay.type = "CITE", slim.col.label = TRUE, remove.key = TRUE, group.label.rot = TRUE, use.scaled = FALSE)

demux_cite = RunPCA(demux_cleaned, pc.genes = rownames(ADT_used), assay.type = "CITE", 
                    pcs.print = 0)
#PCAPlot(demux_cite, pt.size = 0.5)
#____________________________________________________________________________________________________________________
#(Step 4) Clustering based on surface protein markers and finding cell types 

#standard scaling (no regression)
demux_cleaned <- ScaleData(demux_cleaned, display.progress = FALSE)

#Actual clustering based on protein (CITE) data

pbmc_cite <- RunPCA(demux_cleaned, pc.genes = rownames(ADT_used), assay.type = "CITE", 
                    pcs.print = 0)

# Use a standard euclidean distance matrix here. 
adt.data <- GetAssayData(pbmc_cite, assay.type = "CITE", slot = "data")
adt.dist <- as.matrix(dist(t(adt.data)))


# Now, we rerun tSNE using our distance matrix defined only on ADT (protein)
# levels.
pbmc_cite <- RunTSNE(pbmc_cite, distance.matrix = adt.dist)

# We can also rerun clustering using the same distance matrix. We'll start
# with a very coarse clustering (resolution=0.2)
pbmc_cite <- FindClusters(pbmc_cite, distance.matrix = adt.dist, print.output = FALSE, resolution = 0.2)

tsne_adtClusters <- TSNEPlot(pbmc_cite, do.return = TRUE, pt.size = 0.5)
tsne_adtClusters <- tsne_adtClusters + ggtitle("Clustering based on ADT signal") + 
  theme(plot.title = element_text(hjust = 0.5))

#Plot of different cell clusters identified based on protein expression levels
#BETTER THAN RNA EXPRESSION CLUSTERING 
plot(tsne_adtClusters)

DoHeatmap(pbmc_cite, assay.type = "CITE", 
          slim.col.label = TRUE, remove.key = TRUE, group.label.rot = TRUE, use.scaled = FALSE)


pbmc_cite.prot.markers <- FindAllMarkers(pbmc_cite, max.cells.per.ident = 100, logfc.threshold = log(2), only.pos = TRUE, min.diff.pct = 0.3, do.print = F)


shared_clusters_PROT2 <- GetCellTypes_RNA(pbmc_cite.prot.markers, cbmc.rna.markers)

#_______________________________________________________________________________________________
#(Step 5) Separate cleaned and processed Seurat RNA object based on dilutions. Visualize separate dilutions with ridge plots. 

#StashIdent: Stashes the cluster identity in meta.data to be retrieved later
demux_cleaned <- StashIdent(demux_cleaned, save.name = 'cluster.ident')

#Set identity class of demux_cleaned object to hash_maxID: "HTO_1", "HTO_2", ........."HTO_8"
demux_cleaned <- SetIdent(object = demux_cleaned, ident.use = demux_cleaned@meta.data$hash_maxID)

#Dilution 1ug
demux_cleaned_dilution1 = SubsetData(object = demux_cleaned, ident.use = c("HTO_1", "HTO_2"))
demux_cleaned_dilution1 <- StashIdent(demux_cleaned_dilution1, save.name = "HTO1_HTO2")
demux_cleaned_dilution1 <- SetIdent(object = demux_cleaned_dilution1, ident.use = demux_cleaned_dilution1@meta.data$cluster.ident)
# RidgePlot(demux_cleaned_dilution1, features.plot = rownames(ADT_used)[1:9], nCol = 3)
# RidgePlot(demux_cleaned_dilution1, features.plot = rownames(ADT_used)[10:18], nCol = 3)
# RidgePlot(demux_cleaned_dilution1, features.plot = rownames(ADT_used)[19:27], nCol = 3)
# RidgePlot(demux_cleaned_dilution1, features.plot = rownames(ADT_used)[28:31], nCol = 3)

#Dilution .5ug
demux_cleaned_dilution2 = SubsetData(object = demux_cleaned, ident.use = c("HTO_3", "HTO_4"))
demux_cleaned_dilution2 <- StashIdent(demux_cleaned_dilution2, save.name = "HTO3_HTO4")
demux_cleaned_dilution2 <- SetIdent(object = demux_cleaned_dilution2, ident.use = demux_cleaned_dilution2@meta.data$cluster.ident)
# RidgePlot(demux_cleaned_dilution2, features.plot = rownames(ADT_used)[1:9], nCol = 3)
# RidgePlot(demux_cleaned_dilution2, features.plot = rownames(ADT_used)[10:18], nCol = 3)
# RidgePlot(demux_cleaned_dilution2, features.plot = rownames(ADT_used)[19:27], nCol = 3)
# RidgePlot(demux_cleaned_dilution2, features.plot = rownames(ADT_used)[28:31], nCol = 3)

#Dilution .25ug
demux_cleaned_dilution3 = SubsetData(object = demux_cleaned, ident.use = c("HTO_5", "HTO_6"))
demux_cleaned_dilution3 <- StashIdent(demux_cleaned_dilution3, save.name = "HTO5_HTO6")
demux_cleaned_dilution3 <- SetIdent(object = demux_cleaned_dilution3, ident.use = demux_cleaned_dilution3@meta.data$cluster.ident)
# RidgePlot(demux_cleaned_dilution3, features.plot = rownames(ADT_used)[1:9], nCol = 3)
# RidgePlot(demux_cleaned_dilution3, features.plot = rownames(ADT_used)[10:18], nCol = 3)
# RidgePlot(demux_cleaned_dilution3, features.plot = rownames(ADT_used)[19:27], nCol = 3)
# RidgePlot(demux_cleaned_dilution3, features.plot = rownames(ADT_used)[28:31], nCol = 3)

#Dilution .125ug
demux_cleaned_dilution4 = SubsetData(object = demux_cleaned, ident.use = c("HTO_7", "HTO_8"))
demux_cleaned_dilution4 <- StashIdent(demux_cleaned_dilution4, save.name = "HTO7_HTO8")
demux_cleaned_dilution4 <- SetIdent(object = demux_cleaned_dilution4, ident.use = demux_cleaned_dilution4@meta.data$cluster.ident)
# RidgePlot(demux_cleaned_dilution4, features.plot = rownames(ADT_used)[1:9], nCol = 3)
# RidgePlot(demux_cleaned_dilution4, features.plot = rownames(ADT_used)[10:18], nCol = 3)
# RidgePlot(demux_cleaned_dilution4, features.plot = rownames(ADT_used)[19:27], nCol = 3)
# RidgePlot(demux_cleaned_dilution4, features.plot = rownames(ADT_used)[28:31], n1Col = 3)

#_________________________________________________________________________________________________________________________
#(Step 6) Separate cleaned and processed Seurat object (Protein) based on dilutions. Visualize separate dilutions with ridge plots. 

#StashIdent: Stashes the cluster identity in meta.data to be retrieved later
pbmc_cite <- StashIdent(pbmc_cite, save.name = 'cluster.ident')

#Set identity class of pbmc_cite object to hash_maxID: "HTO_1", "HTO_2", ........."HTO_8"
pbmc_cite <- SetIdent(object = pbmc_cite, ident.use = pbmc_cite@meta.data$hash_maxID)

#Dilution 1ug
pbmc_cite_dilution1 = SubsetData(object = pbmc_cite, ident.use = c("HTO_1", "HTO_2"))
# RidgePlot(pbmc_cite_dilution1, features.plot = rownames(ADT_used)[1:9], nCol = 3)
# RidgePlot(pbmc_cite_dilution1, features.plot = rownames(ADT_used)[10:18], nCol = 3)
# RidgePlot(pbmc_cite_dilution1, features.plot = rownames(ADT_used)[19:27], nCol = 3)
# RidgePlot(pbmc_cite_dilution1, features.plot = rownames(ADT_used)[28:31], nCol = 3)

#Dilution .5ug
pbmc_cite_dilution2 = SubsetData(object = pbmc_cite, ident.use = c("HTO_3", "HTO_4"))
# RidgePlot(pbmc_cite_dilution2, features.plot = rownames(ADT_used)[1:9], nCol = 3)
# RidgePlot(pbmc_cite_dilution2, features.plot = rownames(ADT_used)[10:18], nCol = 3)
# RidgePlot(pbmc_cite_dilution2, features.plot = rownames(ADT_used)[19:27], nCol = 3)
# RidgePlot(pbmc_cite_dilution2, features.plot = rownames(ADT_used)[28:31], nCol = 3)

#Dilution .25ug
pbmc_cite_dilution3 = SubsetData(object = pbmc_cite, ident.use = c("HTO_5", "HTO_6"))
pbmc_cite_dilution3 <- SetIdent(object = pbmc_cite_dilution3, ident.use = pbmc_cite_dilution3@meta.data$cluster.ident)
# RidgePlot(pbmc_cite_dilution3, features.plot = rownames(ADT_used)[1:9], nCol = 3)
# RidgePlot(pbmc_cite_dilution3, features.plot = rownames(ADT_used)[10:18], nCol = 3)
# RidgePlot(pbmc_cite_dilution3, features.plot = rownames(ADT_used)[19:27], nCol = 3)
# RidgePlot(pbmc_cite_dilution3, features.plot = rownames(ADT_used)[28:31], nCol = 3)

#Dilution .125ug
pbmc_cite_dilution4 = SubsetData(object = pbmc_cite, ident.use = c("HTO_7", "HTO_8"))
pbmc_cite_dilution4 <- SetIdent(object = pbmc_cite_dilution4, ident.use = pbmc_cite_dilution4@meta.data$cluster.ident)
# RidgePlot(pbmc_cite_dilution4, features.plot = rownames(ADT_used)[1:9], nCol = 3)
# RidgePlot(pbmc_cite_dilution4, features.plot = rownames(ADT_used)[10:18], nCol = 3)
# RidgePlot(pbmc_cite_dilution4, features.plot = rownames(ADT_used)[19:27], nCol = 3)
# RidgePlot(pbmc_cite_dilution4, features.plot = rownames(ADT_used)[28:31], nCol = 3)


#_____________________________________________________________________________________________
#(Step 7) For each dilution and for both protein and RNA datasets, calculate the staining index for 31 antibodies 

#Created two functions: FindCluster() and AdvancedCalculateStainingIndex() functions to calculate stain index
#Functions stored in hpcdata/sg/sg_data/users/jaina10/jaina10/CITESeqDilutionTest/StainingIndexFunctions

#FindCluster() function
##Arguments: a list of Seurat objects and ADT matrix which contains marker names
##Ouput: Returns a dataframe with the cluster number for positive and negative populations for each antibody

#Applying the FindCluster() function with all 4 dilutions and both 
#Set identity to cluster identity (Factor with 6 levels)
pbmc_cite_dilution1 <- SetIdent(object = pbmc_cite_dilution1, ident.use = pbmc_cite_dilution1@meta.data$cluster.ident)
pbmc_cite_dilution2 <- SetIdent(object = pbmc_cite_dilution2, ident.use = pbmc_cite_dilution2@meta.data$cluster.ident)
pbmc_cite_dilution3 <- SetIdent(object = pbmc_cite_dilution3, ident.use = pbmc_cite_dilution3@meta.data$cluster.ident)
pbmc_cite_dilution4 <- SetIdent(object = pbmc_cite_dilution4, ident.use = pbmc_cite_dilution4@meta.data$cluster.ident)

#For both protein and RNA clustering, create a list containing dilution objects 
listofdilutiondata <- list(pbmc_cite_dilution1, pbmc_cite_dilution2, pbmc_cite_dilution3, pbmc_cite_dilution4)
dilution_datasets <- list(demux_cleaned_dilution1, demux_cleaned_dilution2, demux_cleaned_dilution3, demux_cleaned_dilution4)

protein_clusternos <- FindClusters(listofdilutiondata, ADT_used)

clusternums_acrossdilutions <- FindClusters(dilution_datasets, ADT_used)

#Run AdvanceCalculateStainingIndex for both RNA and protein clustering 
#Given the cluster numbers for positive and negative populations for each antibody,
#calculate the staining index across all antibodies for each seurat object
#Arguments: A seurat object, an ADT matrix, output of FindClusters() AKA dataframe with positive and negative pop clusters 
#of each antibody

##Antibody staining indexes based on RNA expression clustering 
RNAexp_dilution.125ug_stainingindices <- AdvancedCalculateStainingIndex(demux_cleaned_dilution4, ADT_used, clusternums_acrossdilutions) 
RNAexp_dilution.25ug_stainingindices <- AdvancedCalculateStainingIndex(demux_cleaned_dilution3, ADT_used, clusternums_acrossdilutions)
RNAexp_dilution.5ug_stainingindices <-  AdvancedCalculateStainingIndex(demux_cleaned_dilution2, ADT_used, clusternums_acrossdilutions)
RNAexp_dilution1ug_stainingindices <-  AdvancedCalculateStainingIndex(demux_cleaned_dilution1, ADT_used, clusternums_acrossdilutions)

Proteinexp_dilution1ug_stainingindices <- AdvancedCalculateStainingIndex(pbmc_cite_dilution1, ADT_used, protein_clusternos)
Proteinexp_dilution0.5ug_stainingindices <- AdvancedCalculateStainingIndex(pbmc_cite_dilution2, ADT_used, protein_clusternos)
Proteinexp_dilution0.25ug_stainingindices<- AdvancedCalculateStainingIndex(pbmc_cite_dilution3, ADT_used, protein_clusternos)
Proteinexp_dilution0.125ug_stainingindices <- AdvancedCalculateStainingIndex(pbmc_cite_dilution4, ADT_used, protein_clusternos)


#Data frames with staining indices for each marker for each dilution

x <- c(0.125, .25, .5, 1)

#Plots for staining index for each antibody (RNA)
CCR4_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[1,2], ADVANCED_dilution.25ug_stainingindices[1,2], ADVANCED_dilution.5ug_stainingindices[1,2], ADVANCED_dilution1ug_stainingindices[1,2]), main = "CCR4_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
CCR7_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[2,2], ADVANCED_dilution.25ug_stainingindices[2,2], ADVANCED_dilution.5ug_stainingindices[2,2], ADVANCED_dilution1ug_stainingindices[2,2]), main = "CCR7_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
CD10_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[3,2], ADVANCED_dilution.25ug_stainingindices[3,2], ADVANCED_dilution.5ug_stainingindices[3,2], ADVANCED_dilution1ug_stainingindices[3,2]), main = "CD10_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
CD11c_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[4,2], ADVANCED_dilution.25ug_stainingindices[4,2], ADVANCED_dilution.5ug_stainingindices[4,2], ADVANCED_dilution1ug_stainingindices[4,2]), main = "CD11c_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
CD138_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[5,2], ADVANCED_dilution.25ug_stainingindices[5,2], ADVANCED_dilution.5ug_stainingindices[5,2], ADVANCED_dilution1ug_stainingindices[5,2]), main = "CD138_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")

CD14_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[6,2], ADVANCED_dilution.25ug_stainingindices[6,2], ADVANCED_dilution.5ug_stainingindices[6,2], ADVANCED_dilution1ug_stainingindices[6,2]), main = "CD14_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
CD16_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[7,2], ADVANCED_dilution.25ug_stainingindices[7,2], ADVANCED_dilution.5ug_stainingindices[7,2], ADVANCED_dilution1ug_stainingindices[7,2]), main = "CD16_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
CD185_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[8,2], ADVANCED_dilution.25ug_stainingindices[8,2], ADVANCED_dilution.5ug_stainingindices[8,2], ADVANCED_dilution1ug_stainingindices[8,2]), main = "CD185_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
CD195_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[9,2], ADVANCED_dilution.25ug_stainingindices[9,2], ADVANCED_dilution.5ug_stainingindices[9,2], ADVANCED_dilution1ug_stainingindices[9,2]), main = "CD195_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
CD196_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[10,2], ADVANCED_dilution.25ug_stainingindices[10,2], ADVANCED_dilution.5ug_stainingindices[10,2], ADVANCED_dilution1ug_stainingindices[10,2]), main = "CD196_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")

CD20_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[11,2], ADVANCED_dilution.25ug_stainingindices[11,2], ADVANCED_dilution.5ug_stainingindices[11,2], ADVANCED_dilution1ug_stainingindices[11,2]), main = "CD20_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
CD25_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[12,2], ADVANCED_dilution.25ug_stainingindices[12,2], ADVANCED_dilution.5ug_stainingindices[12,2], ADVANCED_dilution1ug_stainingindices[12,2]), main = "CD25_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
CD27_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[13,2], ADVANCED_dilution.25ug_stainingindices[13,2], ADVANCED_dilution.5ug_stainingindices[13,2], ADVANCED_dilution1ug_stainingindices[13,2]), main = "CD27_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
CD31_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[14,2], ADVANCED_dilution.25ug_stainingindices[14,2], ADVANCED_dilution.5ug_stainingindices[14,2], ADVANCED_dilution1ug_stainingindices[14,2]), main = "CD31_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
CD34_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[15,2], ADVANCED_dilution.25ug_stainingindices[15,2], ADVANCED_dilution.5ug_stainingindices[15,2], ADVANCED_dilution1ug_stainingindices[15,2]), main = "CD34_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")

CD38_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[16,2], ADVANCED_dilution.25ug_stainingindices[16,2], ADVANCED_dilution.5ug_stainingindices[16,2], ADVANCED_dilution1ug_stainingindices[16,2]), main = "CD38_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
CD45RA_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[17,2], ADVANCED_dilution.25ug_stainingindices[17,2], ADVANCED_dilution.5ug_stainingindices[17,2], ADVANCED_dilution1ug_stainingindices[17,2]), main = "CD45RA_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
CD4_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[18,2], ADVANCED_dilution.25ug_stainingindices[18,2], ADVANCED_dilution.5ug_stainingindices[18,2], ADVANCED_dilution1ug_stainingindices[18,2]), main = "CD4_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
CD56_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[19,2], ADVANCED_dilution.25ug_stainingindices[19,2], ADVANCED_dilution.5ug_stainingindices[19,2], ADVANCED_dilution1ug_stainingindices[19,2]), main = "CD56_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
CD8_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[20,2], ADVANCED_dilution.25ug_stainingindices[20,2], ADVANCED_dilution.5ug_stainingindices[20,2], ADVANCED_dilution1ug_stainingindices[20,2]), main = "CD8_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")

CTLA4_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[21,2], ADVANCED_dilution.25ug_stainingindices[21,2], ADVANCED_dilution.5ug_stainingindices[21,2], ADVANCED_dilution1ug_stainingindices[21,2]), main = "CTLA4_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
CXCR3_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[22,2], ADVANCED_dilution.25ug_stainingindices[22,2], ADVANCED_dilution.5ug_stainingindices[22,2], ADVANCED_dilution1ug_stainingindices[22,2]), main = "CXCR3_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
ICOS_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[23,2], ADVANCED_dilution.25ug_stainingindices[23,2], ADVANCED_dilution.5ug_stainingindices[23,2], ADVANCED_dilution1ug_stainingindices[23,2]), main = "ICOS_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
IgD_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[24,2], ADVANCED_dilution.25ug_stainingindices[24,2], ADVANCED_dilution.5ug_stainingindices[24,2], ADVANCED_dilution1ug_stainingindices[24,2]), main = "IgD_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
IgM_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[25,2], ADVANCED_dilution.25ug_stainingindices[25,2], ADVANCED_dilution.5ug_stainingindices[25,2], ADVANCED_dilution1ug_stainingindices[25,2]), main = "IgM_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")

MIgG1kISO_ISO_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[26,2], ADVANCED_dilution.25ug_stainingindices[26,2], ADVANCED_dilution.5ug_stainingindices[26,2], ADVANCED_dilution1ug_stainingindices[26,2]), main = "MIgG1kISO_ISO_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
MIgG2ak_ISO_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[27,2], ADVANCED_dilution.25ug_stainingindices[27,2], ADVANCED_dilution.5ug_stainingindices[27,2], ADVANCED_dilution1ug_stainingindices[27,2]), main = "MIgG2ak_ISO_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
MIgG2bk_ISO_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[28,2], ADVANCED_dilution.25ug_stainingindices[28,2], ADVANCED_dilution.5ug_stainingindices[28,2], ADVANCED_dilution1ug_stainingindices[28,2]), main = "MIgG2bk_ISO_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
NKG2D_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[29,2], ADVANCED_dilution.25ug_stainingindices[29,2], ADVANCED_dilution.5ug_stainingindices[29,2], ADVANCED_dilution1ug_stainingindices[29,2]), main = "NKG2D_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
PD.1_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[30,2], ADVANCED_dilution.25ug_stainingindices[30,2], ADVANCED_dilution.5ug_stainingindices[30,2], ADVANCED_dilution1ug_stainingindices[30,2]), main = "PD.1_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")
RIgG2bk_ISO_stainingplot <- plot(x, c(ADVANCED_dilution.125ug_stainingindices[31,2], ADVANCED_dilution.25ug_stainingindices[31,2], ADVANCED_dilution.5ug_stainingindices[31,2], ADVANCED_dilution1ug_stainingindices[31,2]), main = "RIgG2bk_ISO_staining_indices", xlab = "Antibody amount (ug)", ylab = "Staining index")

