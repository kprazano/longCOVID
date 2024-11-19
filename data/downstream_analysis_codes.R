#Convert h5ad object to a Seurat object
Convert("/path/to/file/adata_merged_raw.h5ad", dest = "h5seurat", overwrite = TRUE)
LongCovid_seu <- LoadH5Seurat("/path/to/file/adata_merged_raw.h5seurat", meta.data = FALSE, misc = FALSE)
LongCovid_seu

nCount = colSums(x = LongCovid_seu, slot = "counts")  # nCount_RNA
nFeature = colSums(x = GetAssayData(object = LongCovid_seu, slot = "counts") > 0)  # nFeatureRNA
LongCovid_seu$nCount_RNA <- nCount
LongCovid_seu$nFeature_RNA <- nFeature


#Quality Control and visualization
LongCovid_seu[["Percent_mt"]] <- PercentageFeatureSet(LongCovid_seu, pattern = "^MT-")
LongCovid_seu <- subset(LongCovid_seu, subset = nFeature_RNA > 900 & Percent_mt < 25) #& nCount_RNA > 500

#Violin plots showing QC metrics (per sample/ per sequencing set)
Idents(LongCovid_seu) <- LongCovid_seu$Sample #Set
VlnPlot(LongCovid_seu, features = c("nFeature_RNA", "nCount_RNA", "Percent_mt"), ncol = 1, raster = F, pt.size = 0)

#Pearson correlation heatmap
av.exp <- AverageExpression(LongCovid_seu)$RNA
cor.exp <- as.data.frame(cor(av.exp))
cor.exp$x <- rownames(cor.exp)
order_x <- c('set4', 'set5', 'set6', 'set7', 'set8', 'set9', 'set10')
cor.exp$x <- factor(cor.exp$x, levels = order_x)
cor.df <- tidyr::gather(data = cor.exp, key = "y", value = "correlation", -x)
cor.df$y <- factor(cor.df$y, levels = c('set4', 'set5', 'set6', 'set7', 'set8', 'set9', 'set10'))
ggplot(cor.df, aes(x, y, fill = correlation)) + geom_tile() + theme_bw()
ggplot(cor.df, aes(x, y, fill = correlation)) + 
  geom_tile(color = "black") +  # Add black borders
  scale_fill_gradient2(low = "steelblue3", mid = "white", high = "coral1", midpoint = 0.993) + theme_bw()

#Scatter plot showing relationship between the QC metrics
plot1 <- FeatureScatter(LongCovid_seu, feature1 = "nCount_RNA", feature2 = "Percent_mt", raster = FALSE)
plot2 <- FeatureScatter(LongCovid_seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = FALSE)
plot1 + plot2


#Normalization, identification of variable features
LongCovid_seu <- NormalizeData(LongCovid_seu, normalization.method = "LogNormalize", scale.factor = 10000)
LongCovid_seu <- FindVariableFeatures(LongCovid_seu, selection.method = "vst", nfeatures = 2000)

#Plot showing top 10 variable features 
top10 <- head(VariableFeatures(LongCovid_seu), 10)
plot1 <- VariableFeaturePlot(LongCovid_seu)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2


#Data scaling, PCA
all.genes <- rownames(LongCovid_seu)
LongCovid_seu <- ScaleData(LongCovid_seu, features = all.genes)
LongCovid_seu <- RunPCA(LongCovid_seu, features = VariableFeatures(object = LongCovid_seu))

#PCA plot grouped by sequencing set, heatmap of top 5 PCs
Idents(LongCovid_seu) <- LongCovid_seu$Set
DimPlot(LongCovid_seu, reduction = "pca", raster = FALSE) + theme_bw()
DimHeatmap(LongCovid_seu, dims = 1:5, cells = 500, balanced = TRUE, raster = FALSE)


#Determination of the number of PCs to include
ElbowPlot(LongCovid_seu)
LongCovid_seu <- JackStraw(LongCovid_seu)
LongCovid_seu <- ScoreJackStraw(LongCovid_seu, dims = 1:20)
JackStrawPlot(LongCovid_seu, dims = 1:20)


#Cell clustering and visualization
LongCovid_seu <- FindNeighbors(LongCovid_seu, dims = 1:20)
LongCovid_seu <- FindClusters(LongCovid_seu, resolution = 1.0) #resolution 1.0 for level 1 annotation, 0.2 for level 1
LongCovid_seu <- RunUMAP(LongCovid_seu, dims = 1:20)

Idents(LongCovid_seu) <- LongCovid_seu$seurat_clusters
DimPlot(LongCovid_seu, reduction = "umap", raster = FALSE, label = TRUE) + theme_bw() 

#UMAP plot showing distribution of cells from different sequencing sets
color_mapping <- c("set4"="darkgoldenrod2","set5"="springgreen3","set6"="turquoise3","set7" = "deepskyblue", "set8" = "mediumpurple1", "set9" = "violet", "set10" = "coral1")
DimPlot(LongCovid_seu, reduction = "umap", raster=FALSE, cols = color_mapping, order = c("set4", "set5", "set6","set7","set8","set9","set10")) + theme_bw()
DimPlot(LongCovid_seu, reduction = "umap", split.by = "Set", raster=FALSE) + theme_bw() 

#UMAP plot grouped by different metadata
DimPlot(LongCovid_seu, reduction = "umap", raster = F, label = F, group.by = "CIS_score") + theme_bw() #group.by = "Treatment"/"Time" etc.


#Marker identification
markers0.2 <- FindAllMarkers(LongCovid_seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #level 1 markers
markers1.0 <- FindAllMarkers(LongCovid_seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #level 2 markers

#Level 1 cell type annotation
Idents(LongCovid_seu) <- LongCovid_seu$seurat_clusters
new.cluster.ids <- c("0" = "T","1" = "Classic Mono","2" = "T","3" = "NK","4" = "B","5" = "Unassigned", 
                     "6" = "Non-classic Mono","7" = "DC", "8" = "Classic Mono","9" = "Mixed", "10" = "DC", 
                     "11" = "Proliferating","12" = "Classic Mono","13" = "Plasma", "14" = "Platelets", "15" = "HSPC")
names(new.cluster.ids) <- levels(LongCovid_seu)
LongCovid_seu <- RenameIdents(LongCovid_seu, new.cluster.ids)
DimPlot(LongCovid_seu, reduction = "umap", label = T, pt.size = 0.5, raster = F) + theme_bw() #+ NoLegend()
#Save the new annotation
LongCovid_seu[["new.cluster.ids"]] <- Idents(object = LongCovid_seu)
LongCovid_seu <- AddMetaData(LongCovid_seu, LongCovid_seu$new.cluster.ids, col.name = "Annotation_level1") 
Idents(LongCovid_seu) <- LongCovid_seu$Annotation_level1
#Dot plot of the gene markers
markers.to.plot <- c("CD3D","CD3E","CD3G","NCAM1","PRF1","KLRD1","TOP2A","MKI67","CENPF","ATG7","SLC8A1",
                     "FOXO1","IGHG1","IGHG2","IGHA1","MS4A1","BANK1","CD79B","FCN1","S100A8","CD14","FCGR3A",
                     "CSF1R","MS4A7","FCER1A","CD1C","IRF7","PF4","PPBP","TUBB1","CD34","SOX4","SPINK2")
Idents(LongCovid_seu) <- factor(LongCovid_seu$Annotation_level1, levels= c("HSPC","Platelets","DC","Non-classic Mono","Classic Mono","Mixed","B","Plasma","Unassigned","Proliferating","NK","T")) 
DotPlot(LongCovid_seu, features = markers.to.plot, cols = c("white", "blue"), dot.scale = 8) + RotatedAxis() #cols = c("white", "blue")

#Level 2 cell type annotation
Idents(LongCovid_seu) <- LongCovid_seu$seurat_clusters
new.cluster.ids <- c("0" = "CD4+ Memory T", "1" = "dimNK", "2" = "CD14+ TPPP3+ Monocytes", "3" = "CD4+ Naive T", "4" = "CD8+ Memory T", 
                     "5" = "CD8+ Teff", "6" = "CD14+ CYP1B1+ Monocytes", "7" = "CD14+ IFI+ Monocytes", 
                     "8" = "CD4+ Activated T", "9" = "Memory B", "10" = "Autophagic cells", 
                     "11" = "CD16+ Monocytes", "12" = "Mature Naive B", "13" = "CD14+ Early Monocytes", 
                     "14" = "CD8+ Naive T", "15" = "moDC", "16" = "NKT", "17" = "Immature B", 
                     "18" = "Exhausted T", "19" = "CD4+ Low quality T", "20" = "Mixed B/Myeloid", 
                     "21" = "Proliferating T/NK", "22" = "pDC", "23" = "Neutrophils", "24" = "Plasma", 
                     "25" = "Platelets", "26" = "HSPC")
names(new.cluster.ids) <- levels(LongCovid_seu)
LongCovid_seu <- RenameIdents(LongCovid_seu, new.cluster.ids)
DimPlot(LongCovid_seu, reduction = "umap", label = T, pt.size = 0.5, raster = F) + theme_bw() #+ NoLegend()
#Save the new annotation
LongCovid_seu[["new.cluster.ids"]] <- Idents(object = LongCovid_seu)
LongCovid_seu <- AddMetaData(LongCovid_seu, LongCovid_seu$new.cluster.ids, col.name = "Annotation_level2") 
Idents(LongCovid_seu) <- LongCovid_seu$Annotation_level2
#Dot plot of the gene markers
markers.to.plot <- c("CD3D","CD4","LEF1","CCR7","RORA","LRP1B","CD40LG","LTB","GATA3","CD8A","CD8B","GZMK",
                     "GZMH","TRGC2","TOX","STAT4","KLRF1","NCAM1","NCR1","MKI67","TOP2A","CENPF","ATG7",
                     "FOXO1","JCHAIN","IGKC","MS4A1","IL7R","FCER2","BTLA","AIM2","TNFRSF13B","FCN1",
                     "CD14","TPPP3","CYP1B1","IFIT2","OAS2","FCGR3A","MS4A7","CD1C","FCER1A","LILRA4",
                     "G0S2","AQP9","CXCL8","PPBP","PF4","CD34","SOX4")
Idents(LongCovid_seu) <- factor(LongCovid_seu$Annotation_level2, levels= c("CD4+ Naive T","CD4+ Activated T",
                       "CD4+ Low Quality T","CD4+ Memory T", "CD8+ Naive T", "CD8+ Memory T", "CD8+ Teff", 
                       "Exhausted T", "NKT", "NK", "Proliferating T/NK", "Autophagic cells", "Plasma", 
                       "Immature B", "Mature Naive B", "Memory B", "Mixed B/Myeloid","CD14+ TPPP3+ Monocytes", 
                       "CD14+ CYP1B1+ Monocytes", "CD14+ IFI+ Monocytes", "CD14+ Early Monocytes", 
                       "CD16+ Monocytes", "moDC", "pDC", "Neutrophils", "Platelets", "HSPC")) 
DotPlot(LongCovid_seu, features = markers.to.plot, cols = c("white", "blue"), dot.scale = 8) + RotatedAxis() #cols = c("white", "blue")


#Barplot showing proportion of cell types (example for level 2 annotation)
#Generate a dataframe with % of the cell types for each sample
metadata <- LongCovid_seu@meta.data
unique_cell_types <- unique(metadata$Annotation_level2)
unique_orig_idents <- unique(metadata$Group_type)

percentage_matrix <- matrix(NA, nrow = length(unique_cell_types), ncol = length(unique_orig_idents))
rownames(percentage_matrix) <- unique_cell_types
colnames(percentage_matrix) <- unique_orig_idents

for (j in seq_along(unique_orig_idents)) {
  orig_ident <- unique_orig_idents[j]
  subset_data <- metadata[metadata$Group_type == orig_ident, ]
  for (i in seq_along(unique_cell_types)) {
    cell_type <- unique_cell_types[i]
    cell_type_data <- subset_data[subset_data$Annotation_level2 == cell_type, ]
    total_cells <- nrow(cell_type_data)
    percentage_matrix[i, j] <- (total_cells / nrow(subset_data)) * 100}}

percentage_df <- as.data.frame(percentage_matrix)

#Additional settings (cell type order, color)
desired_order <- c("CD4+ Naive T","CD4+ Activated T","CD4+ Low Quality T","CD4+ Memory T", "CD8+ Naive T", 
                   "CD8+ Memory T", "CD8+ Teff", "Exhausted T", "NKT", "NK", "Proliferating T/NK", 
                   "Autophagic cells", "Plasma", "Immature B", "Mature Naive B", "Memory B", "Mixed B/Myeloid",
                   "CD14+ TPPP3+ Monocytes", "CD14+ CYP1B1+ Monocytes", "CD14+ IFI+ Monocytes", 
                   "CD14+ Early Monocytes", "CD16+ Monocytes", "moDC", "pDC", "Neutrophils", "Platelets", "HSPC")
color_mapping <- c("CD4+ Naive T"="palevioletred1","CD4+ Activated T"="lightslateblue","CD4+ Low Quality T"="skyblue2","CD4+ Memory T" = "thistle2", 
                   "CD8+ Naive T" = "orange2", "CD8+ Memory T" = "aquamarine3", "CD8+ Teff" = "plum", "Exhausted T" = "violet", "NKT" = "brown3", 
                   "NK" = "lightblue", "Proliferating T/NK" = "darkslategray3", "Autophagic cells" = "lightsalmon", "Plasma" = "gray", 
                   "Immature B" = "mediumaquamarine", "Mature Naive B" = "goldenrod1", "Memory B" = "lightpink1", "Mixed B/Myeloid" = "deepskyblue",
                   "CD14+ TPPP3+ Monocytes" = "paleturquoise3", "CD14+ CYP1B1+ Monocytes" = "salmon", "CD14+ IFI+ Monocytes" = "lightseagreen", 
                   "CD14+ Early Monocytes" = "hotpink", "CD16+ Monocytes" = "mediumpurple1", "moDC" = "plum3", "pDC" = "steelblue3", 
                   "Neutrophils" = "tan1", "Platelets" = "peru", "HSPC" = "red2")

percentage_df$Cell_type <- factor(percentage_df$Cell_type, levels = desired_order)

#Plot
ggplot(percentage_df, aes(x = "Cell_type", y = sc, fill = Cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_mapping) +  theme_bw()


#Knee-plot graphs using RSEC counts
#Calculate total UMI counts per cell
<samplename>_RSEC_MolsPerCell <- read_csv("<samplename>_RSEC_MolsPerCell.csv", comment = "#")
rownames(<samplename>_RSEC_MolsPerCell) <- <samplename>_RSEC_MolsPerCell$Cell_Index
data <- as.data.frame(<samplename>_RSEC_MolsPerCell)
umi_counts_per_cell <- data %>%
  rowwise() %>%
  mutate(Total_UMI = sum(c_across(-Cell_Index), na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(Cell_Index, Total_UMI)

#Sort UMI counts in descending order
sorted_umi_counts <- umi_counts_per_cell %>%
  arrange(desc(Total_UMI)) %>%
  mutate(Rank = row_number())

#Plot
ggplot(sorted_umi_counts, aes(x = Rank, y = Total_UMI)) +
  geom_line() +
  scale_y_log10() + 
  labs(title = "Knee Plot of UMI Counts per Cell",
       x = "Cell Rank",
       y = "Total UMI Counts (log10 scale)") +
  theme_minimal()
