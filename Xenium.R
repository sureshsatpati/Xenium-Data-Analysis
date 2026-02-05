library(dplyr)
library(Seurat)
library(future)
plan("multisession", workers = 10)
library(ggplot2)
options(future.seed = TRUE)
options(future.globals.maxSize = 18000 * 1024^2)

path.0024875_R1 <- "/rsrch3/home/genomic_med/ssatpati/1.LPS_Xenium/1.Raw/output-XETG00074__0024875__Region_1__20241009__193405"
xenium.obj.0024875_R1 <- LoadXenium(path.0024875_R1, fov = "fov")
xenium.obj.0024875_R1 <- subset(xenium.obj.0024875_R1, subset = nCount_Xenium > 0)



path.0024875_R3 <- "/rsrch3/home/genomic_med/ssatpati/1.LPS_Xenium/1.Raw/output-XETG00074__0024875__Region_3__20241009__193405"
xenium.obj.0024875_R3 <- LoadXenium(path.0024875_R3, fov = "fov")
xenium.obj.0024875_R3 <- subset(xenium.obj.0024875_R3, subset = nCount_Xenium > 0)


xenium.obj.0024875_R1$type = "LPS_10_DD"
xenium.obj.0024875_R3$type = "LPS_10_WD"


xenium.obj.0024875_R1 <- SCTransform(xenium.obj.0024875_R1, assay = "Xenium", verbose = FALSE)
xenium.obj.0024875_R3 <- SCTransform(xenium.obj.0024875_R3, assay = "Xenium", verbose = FALSE)


alldata <- merge(xenium.obj.0024875_R1,c(xenium.obj.0024875_R3), add.cell.ids=c("0024875_R1","0024875_R3"))


DefaultAssay(alldata) <- "Xenium"

pdf("QC_Plot.pdf")
VlnPlot(alldata, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
dev.off()

# Define the custom color palette (ensure there are 20 colors)
custom_palette <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33",
  "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
  "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
  "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072"
)

# Create a named vector for the genes and their corresponding colors
gene_colors <- setNames(custom_palette, c("CDK4", "HMGA2", "PPARG", "ADIPOQ", "CEBPA",
                                          "ATRX", "TP53", "RB1", "CDKN2A", "CDKN2B",
                                          "MKI67", "MTRNR2L12", "YAP1", "MRC1",
                                          "TIGIT", "ICOS", "HAVCR2", "PDCD1",
                                          "CD274", "LAG3"))

# Open the PDF device to save the plot
pdf("Gene_Dim_Plot.pdf", width = 10, height = 8)

# Generate the plot and save it into the 'p' object
p <- ImageDimPlot(
  alldata,
  fov = "fov",
  molecules = c("CDK4", "HMGA2", "PPARG", "ADIPOQ", "CEBPA", "ATRX", "TP53", "RB1",
                "CDKN2A", "CDKN2B", "MKI67", "MTRNR2L12", "YAP1", "MRC1", "TIGIT",
                "ICOS", "HAVCR2", "PDCD1", "CD274", "LAG3"),
  nmols = 20000
)

# Apply the custom color palette and ensure no conflicting scales
p +
  scale_color_manual(values = gene_colors) +
  theme(legend.position = "right")  # Optional: Adjust legend position

# Close the PDF device
dev.off()

pdf("Gene_Feature_Plot.pdf")
ImageFeaturePlot(alldata, features = c("CDK4", "HMGA2", "PPARG", "ADIPOQ", "CEBPA", "ATRX", "TP53", "RB1","CDKN2A", "CDKN2B", "MKI67", "MTRNR2L12", "YAP1", "MRC1", "TIGIT","ICOS", "HAVCR2", "PDCD1", "CD274", "LAG3"), max.cutoff = c(25, 35, 12, 10), size = 0.75, cols = c("white", "red"))
dev.off()

pdf("Gene_Feature_Plot1.pdf")                                                                                                                                                                                                                                               > ImageFeaturePlot(alldata, features = c("CDK4", "HMGA2", "PPARG", "ADIPOQ"), max.cutoff = c(25, 35, 12, 10), size = 0.75, cols = c("white", "red"))
dev.off()
  
pdf("Gene_Feature_Plot2.pdf")
ImageFeaturePlot(alldata, features = c("CEBPA", "ATRX", "TP53", "RB1"), max.cutoff = c(25, 35, 12, 10), size = 0.75, cols = c("white", "red"))
dev.off()

pdf("Gene_Feature_Plot3.pdf")
ImageFeaturePlot(alldata, features = c("CDKN2A", "CDKN2B", "MKI67", "MTRNR2L12"), max.cutoff = c(25, 35, 12, 10), size = 0.75, cols = c("white", "red"))
dev.off()

pdf("Gene_Feature_Plot4.pdf")
ImageFeaturePlot(alldata, features = c("YAP1", "MRC1", "TIGIT","ICOS"), max.cutoff = c(25, 35, 12, 10), size = 0.75, cols = c("white", "red"))
dev.off()

pdf("Gene_Feature_Plot5.pdf")
ImageFeaturePlot(alldata, features = c("HAVCR2", "PDCD1", "CD274", "LAG3"), max.cutoff = c(25, 35, 12, 10), size = 0.75, cols = c("white", "red"))
dev.off()


#options(future.seed = TRUE)
#options(future.globals.maxSize = 18000 * 1024^2)
#alldata <- SCTransform(alldata, assay = "Xenium")
#https://satijalab.org/seurat/reference/prepsctfindmarkers

alldata <- NormalizeData(alldata)
alldata <- ScaleData(alldata)
alldata <- RunPCA(alldata, npcs = 30, features = rownames(alldata))
alldata <- RunUMAP(alldata, dims = 1:30)
alldata <- FindNeighbors(alldata, reduction = "pca", dims = 1:30)
alldata <- FindClusters(alldata, resolution = 0.3)


alldata <- PrepSCTFindMarkers(object = alldata)
alldata <- JoinLayers(alldata)



pdf("UMAP.pdf")
DimPlot(alldata)
dev.off()

pdf("Cluster_Gene.pdf")
FeaturePlot(alldata, features = c("CDK4", "HMGA2", "PPARG", "ADIPOQ", "CEBPA", "ATRX", "TP53", "RB1","CDKN2A", "CDKN2B", "MKI67", "MTRNR2L12", "YAP1", "MRC1", "TIGIT","ICOS", "HAVCR2", "PDCD1", "CD274", "LAG3"))
dev.off()

pdf("Cluster_Gene1.pdf")
FeaturePlot(alldata, features = c("CDK4", "HMGA2", "PPARG", "ADIPOQ"))
dev.off()

pdf("Cluster_Gene2.pdf")
FeaturePlot(alldata, features = c("CEBPA", "ATRX", "TP53", "RB1"))
dev.off()

pdf("Cluster_Gene3.pdf")
FeaturePlot(alldata, features = c("CDKN2A", "CDKN2B", "MKI67", "MTRNR2L12"))
dev.off()

pdf("Cluster_Gene4.pdf")
FeaturePlot(alldata, features = c("YAP1", "MRC1", "TIGIT","ICOS"))
dev.off()

pdf("Cluster_Gene5.pdf")
FeaturePlot(alldata, features = c("HAVCR2", "PDCD1", "CD274", "LAG3"))
dev.off()


alldata.markers <- FindAllMarkers(alldata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
alldata.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
top100_markers <- alldata.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
write.table(top100_markers, file="Reg1-3_MP_LPS.markers.xls", quote=FALSE, sep="\t", row.names=FALSE)

alldata.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_markers <- alldata.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10_markers, file="Reg1-3_MP_LPS_10.markers.xls", quote=FALSE, sep="\t", row.names=FALSE)


saveRDS(alldata,'Reg1-3_Matched_LPS10_Pair_Xenium_final.rds')

library(Seurat)
library(ShinyCell)


scConf = createConfig(alldata)
DefaultAssay(alldata) = "SCT"

makeShinyApp(alldata, scConf,gex.assay = "SCT",gene.mapping = TRUE, shiny.dir = "Reg1-3_Xenium_LPS10_MP_together_Shiny", shiny.title = "Reg1-3_Xenium_LPS10_Matched_Pair_Together",shiny.footnotes = "by Rai Lab")





