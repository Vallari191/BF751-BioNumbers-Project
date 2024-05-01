
library(Seurat)
library(SeuratObject)
library(dplyr)
library(Seurat)
library(patchwork)
# Set the working directory to the path for storing data files
setwd("/Users/irisxu/Desktop/BU/Spring/BF751 molecular bio/project/GSE148071")

# Gets a list of all txt.gz files
file_list <- list.files(".", pattern = "\\.txt\\.gz$")

# Creates an empty list to store Seurat objects
seurat_list <- list()

# Loop through the data for each txt.gz file and create the Seurat object
for (file in file_list) {
  # Concatenated file path
  data.path <- paste0("./", file)
  seurat_data <- read.table(gzfile(data.path), row.names = 1,sep = '\t')
  
  # Create the Seurat object and specify the project name as the file name
  sample_name <- tools::file_path_sans_ext(basename(file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   project = "Seurat obj",
                                   min.features = 200,
                                   min.cells = 3)
  # Add the Seurat object to the list
  seurat_list <- append(seurat_list, seurat_obj)
}

# Extract the sample id (eg: P1, p2)
sample_names <- sub(".*_(P[0-9]+)_.*", "\\1",file_list)
# Merge all Seurat objects into one object
seurat_combined <- merge(seurat_list[[1]],
                         y = seurat_list[-1],
                         add.cell.ids = sample_names)

print(seurat_combined)

# QC
seurat_combined[["percent.mt"]] <- PercentageFeatureSet(seurat_combined, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(seurat_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# discarded cells with more than 30,000 UMIs
seurat_combined <- subset(seurat_combined, subset = nCount_RNA <= 30000)
# removed cells that had either lower than 200 or higher than 5000 expressed genes and  mitochondria content higher than 30%
seurat_combined <- subset(seurat_combined, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)

plot1 <- FeatureScatter(seurat_combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# save the comnbined files
saveRDS(seurat_combined, file = "/Users/irisxu/Desktop/BU/Spring/BF751 molecular bio/project/3sample_seurat_object.rds")

# load RDS file
seurat_combined <- readRDS(file = "/Users/irisxu/Desktop/BU/Spring/BF751 molecular bio/project/3sample_seurat_object.rds")

#Checking for NA values in our Seurat Metadata 
is.na(seurat_combined@meta.data) #All values come out as FALSE 
sum(is.na(seurat_combined@meta.data)) # 0 NA values detected 

# Normalizing the data 
# Using the NormalizeData function and assigning it back to our Seurat Object 
seurat_combined <- NormalizeData(seurat_combined, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_combined@assays$RNA@meta.data # Unsure of this 

#Identification of highly variable features (feature selection) using the FindVariableFeatures
seurat_combined <- FindVariableFeatures(seurat_combined, selection.method = "vst", nfeatures = 2000)
View(seurat_combined@assays$RNA@meta.data) # The slot in the Seurat Object to view the variable feature information 
sum(is.na(seurat_combined@assays$RNA@meta.data)) #Alot of NA values 176025


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_combined), 10)
top10 # Top 10 genes are "SPRR1B"  "SPRR3"   "S100A8"  "S100A9"  "SPP1"    "WFDC21P" "TPSB2"   "IGHG4"   "IGKC"    "SLPI"  

# Plotting the Variable features with and without labels 
plot1 <- VariableFeaturePlot(seurat_combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = FALSE, xnudge = 0, ynudge = 0)
plot1 + plot2 #There are a few NA values which need to be removed. Still working on that 

#Scaling the data 
all.genes <- rownames(seurat_combined)
seurat_combined <- ScaleData(seurat_combined, features = all.genes)

# Performing Linear Dimensional reduction using the RunPCA function
seurat_combined <- RunPCA(seurat_combined, features = VariableFeatures(object = seurat_combined))


#Determine the dimensionality of the dataset 
plot5 <- ElbowPlot(seurat_combined)
plot5 #According to the elbow plot 5 seems like a good number of dimensions

# Examine and visualize PCA results a few different ways
# Using 5 as the number of dimensions according to the Elbow plot above
print(seurat_combined[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize the PCA 
plot3 <- VizDimLoadings(seurat_combined, dims = 1:2, reduction = "pca")
plot3

plot4 <- DimPlot(seurat_combined, reduction = "pca", label = TRUE)
plot4

# Clustering the cells 
seurat_combined <- FindNeighbors(seurat_combined, dims = 1:5)
seurat_combined <- FindClusters(seurat_combined, resolution = 0.5)


# Look at cluster IDs of the first 5 cells
head(Idents(seurat_combined), 5) # 13 clusters found by the Louvain algorithm
seurat_combined <- RunUMAP(seurat_combined, dims = 1:10)
plot6 <- DimPlot(seurat_combined, group.by='seurat_clusters', label=TRUE, label.size = 2.8) + ggtitle("Bionumbers UMAP by seurat clusters")
plot6 #The final UMAP with the 13 clusters 

#number of cells in each cluster
table(Idents(seurat_combined))
#persentage of cells in each clusters
prop.table(table(Idents(seurat_combined)))

#DimPlot(seurat_combined,label = T,split.by = "orig.ident",ncol = 3)

DefaultAssay(seurat_combined) <- "RNA"
#Once integrative analysis is complete, you can rejoin the layers - which collapses the individual datasets together and 
#recreates the original counts and data layers. You will need to do this before performing any differential expression analysis. 
#However, you can always resplit the layers in case you would like to reperform integrative analysis.
seurat_combined <- JoinLayers(seurat_combined)
#FindAllMarkers takes time and CPU
all.markers <- FindAllMarkers(seurat_combined, #parameters can be changed 
                               only.pos = TRUE, 
                               min.pct = 0.01, 
                               logfc.threshold = 0.1)
head(x = all.markers)
#get the significant genes for each cluster
significant.markers  <- all.markers [all.markers $p_val_adj < 0.05, ]
write.csv(significant.markers, file = "significant.markers_padj<0.05.csv")






