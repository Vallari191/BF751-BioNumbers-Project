# Script to get final UMAP 
# Annotation is till remaining 
# The research paper that are analysis is based on https://pubmed.ncbi.nlm.nih.gov/36195615/ 
# This Script uses Seurat version 5.1.0 and Seurat Object version 5.0.1
# Since we already performed QC (pre-processing) on our data , I started the analysis by Normalizing our data 


# Installing Seurat version 5.1.0
# This is where I installed Seurat v5 https://satijalab.org/seurat/ 
install.packages('Seurat')
setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))
# Install the remotes package
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
install.packages('Signac')
remotes::install_github("satijalab/seurat-data", quiet = TRUE)
remotes::install_github("satijalab/azimuth", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)
# install.packages("devtools")
devtools::install_github("immunogenomics/presto")


#Loading the required packages 
library(dplyr) 
library(Seurat) 
library(patchwork)  
library(ggplot2) 
#library(presto) #We'll need this package when we perform annotation 

# Reading into the Seurat Object which is saved in .rds format using the readRDS() function

seurat_bionumbers <- readRDS(file = '/usr4/bf527/vallarin/3sample_seurat_object.rds')

# Looking into the Metadata of the Seurat Object
View(seurat_bionumbers@meta.data)

#Checking for NA values in our Seurat Metadata 
is.na(seurat_bionumbers@meta.data) #All values come out as FALSE 
sum(is.na(seurat_bionumbers@meta.data)) # 0 NA values detected 

# Normalizing the data 
# Using the NormalizeData function and assigning it back to our Seurat Object 
seurat_bionumbers <- NormalizeData(seurat_bionumbers, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_bionumbers@assays$RNA@meta.data # Unsure of this 

#Identification of highly variable features (feature selection) using the FindVariableFeatures
seurat_bionumbers <- FindVariableFeatures(seurat_bionumbers, selection.method = "vst", nfeatures = 2000)
View(seurat_bionumbers@assays$RNA@meta.data) # The slot in the Seurat Object to view the variable feature information 
sum(is.na(seurat_bionumbers@assays$RNA@meta.data)) #Alot of NA values 176025
#seurat_bionumbers_filtered <- seurat_bionumbers[, !is.na(seurat_bionumbers@assays$RNA@meta.data$vf_vst_counts.2_variance.expected)] 


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_bionumbers), 10)
top10 # Top 10 genes are "SPRR1B"  "SPRR3"   "S100A8"  "S100A9"  "SPP1"    "WFDC21P" "TPSB2"   "IGHG4"   "IGKC"    "SLPI"  

# Plotting the Variable features with and without labels 
plot1 <- VariableFeaturePlot(seurat_bionumbers)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = FALSE, xnudge = 0, ynudge = 0)
plot1 + plot2 #There are a few NA values which need to be removed. Still working on that 

#Scaling the data 
all.genes <- rownames(seurat_bionumbers)
seurat_bionumbers <- ScaleData(seurat_bionumbers, features = all.genes)

#Determine the dimensionality of the dataset 
plot5 <- ElbowPlot(seurat_bionumbers)
plot5 #According to the elbow plot 5 seems like a good number of dimensions

# Performing Linear Dimensional reduction using the RunPCA function
seurat_bionumbers <- RunPCA(seurat_bionumbers, features = VariableFeatures(object = seurat_bionumbers))
# Examine and visualize PCA results a few different ways
# Using 5 as the number of dimensions according to the Elbow plot above
print(seurat_bionumbers[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize the PCA 
plot3 <- VizDimLoadings(seurat_bionumbers, dims = 1:2, reduction = "pca")
plot3

plot4 <- DimPlot(seurat_bionumbers, reduction = "pca", label = TRUE)
plot4

# Clustering the cells 
seurat_bionumbers <- FindNeighbors(seurat_bionumbers, dims = 1:5)
seurat_bionumbers <- FindClusters(seurat_bionumbers, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(seurat_bionumbers), 5) # 13 clusters found by the Louvain algorithm
seurat_bionumbers <- RunUMAP(seurat_bionumbers, dims = 1:10)
plot6 <- DimPlot(seurat_bionumbers, group.by='seurat_clusters', label=TRUE, label.size = 2.8) + ggtitle("Bionumbers UMAP by seurat clusters")
plot6 #The final UMAP with the 13 clusters 


# ~~~~~~~~~~~~~~~~~~~ We have the Final UMAP, Now we need to annotate it to find marker genes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#test script to find the Markers. Still working on it

# # Find Markers for each of the clusters 
# # find all markers of cluster 2
# seurat_bionumbers <- JoinLayers(seurat_bionumbers)
# cluster2.markers <- FindMarkers(seurat_bionumbers, ident.1 = 2)
# head(cluster2.markers, n = 5)
# 
# 
# # find all markers of cluster 2
# cluster2.markers <- FindMarkers(seurat_bionumbers, ident.1 = 2)
# head(cluster2.markers, n = 5)
# 
# # find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster5.markers <- FindMarkers(seurat_bionumbers, ident.1 = 5, ident.2 = c(0, 3))
# head(cluster5.markers, n = 5)
# 
# # find markers for every cluster compared to all remaining cells, report only the positive
# # ones
# pbmc.markers <- FindAllMarkers(seurat_bionumbers, only.pos = TRUE)
# pbmc.markers %>%
#   group_by(cluster) %>%
#   dplyr::filter(avg_log2FC > 1)
# 
# cluster0.markers <- FindMarkers(seurat_bionumbers, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
# plot6 <- VlnPlot(seurat_bionumbers, features = c("MS4A1", "CD79A"))
# # you can plot raw counts as well
# plot7 <- VlnPlot(seurat_bionumbers, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
# plot8 <- FeaturePlot(seurat_bionumbers, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
#                                "CD8A"))
# 
# new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
#                      "NK", "DC", "Platelet")
# names(new.cluster.ids) <- levels(seurat_bionumbers)
# seurat_bionumbers <- RenameIdents(seurat_bionumbers, new.cluster.ids)
# plot9 <-DimPlot(seurat_bionumbers, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# 
