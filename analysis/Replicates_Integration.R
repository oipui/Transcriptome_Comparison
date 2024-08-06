library(Matrix)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(patchwork)

setwd("/home/OiPui/share/OiPui/Ciona_data")
#this will create a df
ciona <- read_tsv(file="expression_matrix_10stage.tsv.gz")
ciona <- as.data.frame(ciona)
#This assigns rownames to the data.frame
#and then removes the first column which contains the gene names that we just moved to rownames
rownames(ciona) <- ciona$GENE
ciona <- ciona[, -1]

#transform into sparce matrix
ciona_matrix <- Matrix(as.matrix(ciona), sparse = TRUE)

ciona_meta1 <- read_tsv(file = "ciona10stage.cluster.upload.new.txt")
ciona_meta1 <- ciona_meta1 %>% 
  dplyr::filter(NAME != "TYPE") 

ciona_meta2 <- read_tsv(file= "ciona10stage.meta.upload.new.MSTRG.txt")
ciona_meta2 <- ciona_meta2 %>% 
  dplyr::filter(NAME != "TYPE")

#joins the two metadata columns by name 
ciona_metadata <- ciona_meta1 %>% inner_join(ciona_meta2, by = "NAME") 
unique(ciona_metadata$orig.ident)

ciona_metadata <- ciona_metadata %>% 
  dplyr::rename("Tissue_type" = `Tissue Type`) %>% 
  dplyr::mutate(X = as.numeric(X),
                Y = as.numeric(Y),
                nGene = as.numeric(nGene),
                nUMI = as.numeric(nUMI),
                percent.mito = as.numeric(percent.mito))

stages.list <- ciona_metadata %>% group_by(stage) %>% group_split()
names(stages.list) <- c("a.initial_gastrula",
                        "b.middle_gastrula",
                        "c.early_neurula",
                        "d.late_neurula",
                        "e.initial_tailbudI",
                        "f.early_tailbudI",   
                        "g.middle_tailbudII",
                        "h.late_tailbudI",    
                        "i.late_tailbudII",   
                        "j.larva")


SeuratObjects <- list()
for (name in names(stages.list)) {
  # print(name)
  stages.list[[name]] <- as.data.frame(stages.list[[name]])
  rownames(stages.list[[name]]) <- stages.list[[name]]$NAME
  matrix <- ciona_matrix[,stages.list[[name]]$NAME]
  meta = stages.list[[name]]
  # print(SeuratObjects[[name]])
  SeuratObjects[[name]] <- CreateSeuratObject(matrix, meta.data = meta)
}

list2env(SeuratObjects,envir=.GlobalEnv)

VlnPlot(a.initial_gastrula, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
a_plt1 <- FeatureScatter(a.initial_gastrula, feature1 = "nCount_RNA", feature2= "percent.mito")
a_plt2 <- FeatureScatter(a.initial_gastrula, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
a_plt1 + a_plt2

a.initial_gastrula <- subset(a.initial_gastrula, subset = percent.mito < 0.2)

VlnPlot(b.middle_gastrula , features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
b_plt1 <- FeatureScatter(b.middle_gastrula, feature1 = "nCount_RNA", feature2= "percent.mito")
b_plt2 <- FeatureScatter(b.middle_gastrula, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
b_plt1 + b_plt2

b.middle_gastrula <- subset(b.middle_gastrula, subset = percent.mito < 0.2)

VlnPlot(c.early_neurula , features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
c_plt1 <- FeatureScatter(c.early_neurula, feature1 = "nCount_RNA", feature2= "percent.mito")
c_plt2 <- FeatureScatter(c.early_neurula, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
c_plt1 + c_plt2

c.early_neurula <- subset(c.early_neurula, subset = percent.mito < 0.25 )

VlnPlot(d.late_neurula , features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
d_plt1 <- FeatureScatter(d.late_neurula, feature1 = "nCount_RNA", feature2= "percent.mito")
d_plt2 <- FeatureScatter(d.late_neurula, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
d_plt1 + d_plt2

d.late_neurula <- subset(d.late_neurula, subset = percent.mito < 0.25)

VlnPlot(e.initial_tailbudI , features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
e_plt1 <- FeatureScatter(e.initial_tailbudI, feature1 = "nCount_RNA", feature2= "percent.mito")
e_plt2 <- FeatureScatter(e.initial_tailbudI, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
e_plt1 + e_plt2

e.initial_tailbudI <- subset(e.initial_tailbudI, subset = percent.mito < 0.25)

VlnPlot(f.early_tailbudI , features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
f_plt1 <- FeatureScatter(f.early_tailbudI, feature1 = "nCount_RNA", feature2= "percent.mito")
f_plt2 <- FeatureScatter(f.early_tailbudI, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
f_plt1 + f_plt2

f.early_tailbudI <- subset(f.early_tailbudI, subset = percent.mito < 0.25)

VlnPlot(g.middle_tailbudII , features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
g_plt1 <- FeatureScatter(g.middle_tailbudII, feature1 = "nCount_RNA", feature2= "percent.mito")
g_plt2 <- FeatureScatter(g.middle_tailbudII, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
g_plt1 + g_plt2

g.middle_tailbudII <- subset(g.middle_tailbudII, subset = percent.mito < 0.25)

VlnPlot(h.late_tailbudI , features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
h_plt1 <- FeatureScatter(h.late_tailbudI, feature1 = "nCount_RNA", feature2= "percent.mito")
h_plt2 <- FeatureScatter(h.late_tailbudI, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
h_plt1 + h_plt2 

h.late_tailbudI <- subset(h.late_tailbudI, subset = percent.mito < 0.25)

VlnPlot(i.late_tailbudII , features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
i_plt1 <- FeatureScatter(i.late_tailbudII, feature1 = "nCount_RNA", feature2= "percent.mito")
i_plt2 <- FeatureScatter(i.late_tailbudII, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
i_plt1 + i_plt2

i.late_tailbudII <- subset(i.late_tailbudII, subset = percent.mito < 0.15)

VlnPlot(j.larva , features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
j_plt1 <- FeatureScatter(j.larva, feature1 = "nCount_RNA", feature2= "percent.mito")
j_plt2 <- FeatureScatter(j.larva, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
j_plt1 + j_plt2

j.larva <- subset(j.larva, subset = percent.mito < 0.2)

#split object into three replicates
init_gastrula.list <- SplitObject(a.initial_gastrula, split.by = "orig.ident")
middle_gastrula.list <- SplitObject(b.middle_gastrula , split.by = "orig.ident")
early_neurula.list <- SplitObject(c.early_neurula , split.by = "orig.ident")
late_neurula.list <- SplitObject(d.late_neurula, split.by = "orig.ident")
init_tail.list <- SplitObject(e.initial_tailbudI, split.by = "orig.ident")
early_tail.list <- SplitObject(f.early_tailbudI, split.by = "orig.ident")
mid_tail.list <- SplitObject(g.middle_tailbudII, split.by = "orig.ident")
late_tail.list <- SplitObject(h.late_tailbudI, split.by = "orig.ident")
late_tailII.list <- SplitObject(i.late_tailbudII, split.by = "orig.ident")
larval.list <- SplitObject(j.larva, split.by = "orig.ident")

for (i in 1:length(init_gastrula.list)) {
  init_gastrula.list[[i]] <- NormalizeData(init_gastrula.list[[i]], verbose = FALSE)
  init_gastrula.list[[i]] <- FindVariableFeatures(init_gastrula.list[[i]], selection.method = "vst", 
                                                  nfeatures = 2000, verbose = FALSE)
}

reference.list <- init_gastrula.list[c("C110.1","C110.2")]
init_gastrula.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
init_gastrula.integrated <- IntegrateData(anchorset = init_gastrula.anchors, dims = 1:30)

DefaultAssay(init_gastrula.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
init_gastrula.integrated <- ScaleData(init_gastrula.integrated, verbose = FALSE)
init_gastrula.integrated <- RunPCA(init_gastrula.integrated, npcs = 30, verbose = FALSE)
init_gastrula.integrated <- RunUMAP(init_gastrula.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(init_gastrula.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(init_gastrula.integrated, reduction = "umap", group.by = "Tissue_type", label = TRUE, 
              repel = TRUE)
init_gastrula.plots <- p1 + p2
filtered_init_gastrula.plots <- p1 + p2


####mid gastrula 
for (i in 1:length(middle_gastrula.list)) {
  middle_gastrula.list[[i]] <- NormalizeData(middle_gastrula.list[[i]], verbose = FALSE)
  middle_gastrula.list[[i]] <- FindVariableFeatures(middle_gastrula.list[[i]], selection.method = "vst", 
                                                    nfeatures = 2000, verbose = FALSE)
}

reference.list <- middle_gastrula.list[c("midG.1","midG.2")]
middle_gastrula.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
middle_gastrula.integrated <- IntegrateData(anchorset = middle_gastrula.anchors, dims = 1:30)

DefaultAssay(middle_gastrula.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
middle_gastrula.integrated <- ScaleData(middle_gastrula.integrated, verbose = FALSE)
middle_gastrula.integrated <- RunPCA(middle_gastrula.integrated, npcs = 30, verbose = FALSE)
middle_gastrula.integrated <- RunUMAP(middle_gastrula.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(middle_gastrula.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(middle_gastrula.integrated, reduction = "umap", group.by = "Tissue_type", label = TRUE, 
              repel = TRUE)
middle_gastrula.plots <- p1 + p2
filtered_middle_gastrula.plots <- p1 + p2

###early neurula
for (i in 1:length(early_neurula.list)) {
  early_neurula.list[[i]] <- NormalizeData(early_neurula.list[[i]], verbose = FALSE)
  early_neurula.list[[i]] <- FindVariableFeatures(early_neurula.list[[i]], selection.method = "vst", 
                                                  nfeatures = 2000, verbose = FALSE)
}

reference.list <- early_neurula.list[c("earlyN.1","earlyN.2")]
early_neurula.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
early_neurula.integrated <- IntegrateData(anchorset = early_neurula.anchors, dims = 1:30)

DefaultAssay(early_neurula.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
early_neurula.integrated <- ScaleData(early_neurula.integrated, verbose = FALSE)
early_neurula.integrated <- RunPCA(early_neurula.integrated, npcs = 30, verbose = FALSE)
early_neurula.integrated <- RunUMAP(early_neurula.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(early_neurula.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(early_neurula.integrated, reduction = "umap", group.by = "Tissue_type", label = TRUE, 
              repel = TRUE)
early_neurula.plots <- p1 + p2
filtered_early_neurula.plots <- p1 + p2

###late neurula
for (i in 1:length(late_neurula.list)) {
  late_neurula.list[[i]] <- NormalizeData(late_neurula.list[[i]], verbose = FALSE)
  late_neurula.list[[i]] <- FindVariableFeatures(late_neurula.list[[i]], selection.method = "vst", 
                                                 nfeatures = 2000, verbose = FALSE)
}

reference.list <- late_neurula.list[c("lateN.1","lateN.2")]
late_neurula.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
late_neurula.integrated <- IntegrateData(anchorset = late_neurula.anchors, dims = 1:30)

DefaultAssay(late_neurula.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
late_neurula.integrated <- ScaleData(late_neurula.integrated, verbose = FALSE)
late_neurula.integrated <- RunPCA(late_neurula.integrated, npcs = 30, verbose = FALSE)
late_neurula.integrated <- RunUMAP(late_neurula.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(late_neurula.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(late_neurula.integrated, reduction = "umap", group.by = "Tissue_type", label = TRUE, 
              repel = TRUE)
late_neurula.plots <- p1 + p2
filtered_late_neurula.plots <- p1 + p2

###initial tailbud
for (i in 1:length(init_tail.list)) {
  init_tail.list[[i]] <- NormalizeData(init_tail.list[[i]], verbose = FALSE)
  init_tail.list[[i]] <- FindVariableFeatures(init_tail.list[[i]], selection.method = "vst", 
                                              nfeatures = 2000, verbose = FALSE)
}

reference.list <- init_tail.list[c("ITB.1","ITB.3")]
init_tail.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
init_tail.integrated <- IntegrateData(anchorset = init_tail.anchors, dims = 1:30)

DefaultAssay(init_tail.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
init_tail.integrated <- ScaleData(init_tail.integrated, verbose = FALSE)
init_tail.integrated <- RunPCA(init_tail.integrated, npcs = 30, verbose = FALSE)
init_tail.integrated <- RunUMAP(init_tail.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(init_tail.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(init_tail.integrated, reduction = "umap", group.by = "Tissue_type", label = TRUE, 
              repel = TRUE)
init_tail.plots <- p1 + p2
filtered_init_tail.plots <- p1 + p2

###early tailbud
for (i in 1:length(early_tail.list)) {
  early_tail.list[[i]] <- NormalizeData(early_tail.list[[i]], verbose = FALSE)
  early_tail.list[[i]] <- FindVariableFeatures(early_tail.list[[i]], selection.method = "vst", 
                                               nfeatures = 2000, verbose = FALSE)
}

reference.list <- early_tail.list[c("ETB.1","ETB.2")]
early_tail.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
early_tail.integrated <- IntegrateData(anchorset = early_tail.anchors, dims = 1:30)

DefaultAssay(early_tail.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
early_tail.integrated <- ScaleData(early_tail.integrated, verbose = FALSE)
early_tail.integrated <- RunPCA(early_tail.integrated, npcs = 30, verbose = FALSE)
early_tail.integrated <- RunUMAP(early_tail.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(early_tail.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(early_tail.integrated, reduction = "umap", group.by = "Tissue_type", label = TRUE, 
              repel = TRUE)
early_tail.plots <- p1 + p2
filtered_early_tail.plots <- p1 + p2

###middle tailbud
for (i in 1:length(mid_tail.list)) {
  mid_tail.list[[i]] <- NormalizeData(mid_tail.list[[i]], verbose = FALSE)
  mid_tail.list[[i]] <- FindVariableFeatures(mid_tail.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}

reference.list <- mid_tail.list[c("MTB.1","MTB.2")]
mid_tail.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
mid_tail.integrated <- IntegrateData(anchorset = mid_tail.anchors, dims = 1:30)

DefaultAssay(mid_tail.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
mid_tail.integrated <- ScaleData(mid_tail.integrated, verbose = FALSE)
mid_tail.integrated <- RunPCA(mid_tail.integrated, npcs = 30, verbose = FALSE)
mid_tail.integrated <- RunUMAP(mid_tail.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(mid_tail.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(mid_tail.integrated, reduction = "umap", group.by = "Tissue_type", label = TRUE, 
              repel = TRUE)
mid_tail.plots <- p1 + p2
filtered_mid_tail.plots <- p1 + p2 

###late tailbud I
for (i in 1:length(late_tail.list)) {
  late_tail.list[[i]] <- NormalizeData(late_tail.list[[i]], verbose = FALSE)
  late_tail.list[[i]] <- FindVariableFeatures(late_tail.list[[i]], selection.method = "vst", 
                                              nfeatures = 2000, verbose = FALSE)
}

reference.list <- late_tail.list[c("LTB1.1","LTB1.2")]
late_tail.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
late_tail.integrated <- IntegrateData(anchorset = late_tail.anchors, dims = 1:30)

DefaultAssay(late_tail.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
late_tail.integrated <- ScaleData(late_tail.integrated, verbose = FALSE)
late_tail.integrated <- RunPCA(late_tail.integrated, npcs = 30, verbose = FALSE)
late_tail.integrated <- RunUMAP(late_tail.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(late_tail.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(late_tail.integrated, reduction = "umap", group.by = "Tissue_type", label = TRUE, 
              repel = TRUE)
late_tail.plots <- p1 + p2
filtered_late_tail.plots <- p1 + p2

###late tailbud II
for (i in 1:length(late_tailII.list)) {
  late_tailII.list[[i]] <- NormalizeData(late_tailII.list[[i]], verbose = FALSE)
  late_tailII.list[[i]] <- FindVariableFeatures(late_tailII.list[[i]], selection.method = "vst", 
                                                nfeatures = 2000, verbose = FALSE)
}

reference.list <- late_tailII.list[c("LTB2.1","LTB2.2")]
late_tailII.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
late_tailII.integrated <- IntegrateData(anchorset = late_tailII.anchors, dims = 1:30)

DefaultAssay(late_tailII.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
late_tailII.integrated <- ScaleData(late_tailII.integrated, verbose = FALSE)
late_tailII.integrated <- RunPCA(late_tailII.integrated, npcs = 30, verbose = FALSE)
late_tailII.integrated <- RunUMAP(late_tailII.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(late_tailII.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(late_tailII.integrated, reduction = "umap", group.by = "Tissue_type", label = TRUE, 
              repel = TRUE)
late_tailII.plots <- p1 + p2
filtered_late_tailII.plots <- p1 + p2 

###larval
for (i in 1:length(larval.list)) {
  larval.list[[i]] <- NormalizeData(larval.list[[i]], verbose = FALSE)
  larval.list[[i]] <- FindVariableFeatures(larval.list[[i]], selection.method = "vst", 
                                           nfeatures = 2000, verbose = FALSE)
}

reference.list <- larval.list[c("lv.3","lv.4","lv.1")]
larval.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
larval.integrated <- IntegrateData(anchorset = larval.anchors, dims = 1:30)

DefaultAssay(larval.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
larval.integrated <- ScaleData(larval.integrated, verbose = FALSE)
larval.integrated <- RunPCA(larval.integrated, npcs = 30, verbose = FALSE)
larval.integrated <- RunUMAP(larval.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(larval.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(larval.integrated, reduction = "umap", group.by = "Tissue_type", label = TRUE, 
              repel = TRUE)
larval.plots <- p1 + p2
filtered_larval.plots <- p1 + p2

saveRDS(init_gastrula.integrated, file = "IntegratedReplicates/init_gastrula.int")
saveRDS(middle_gastrula.integrated, file = "IntegratedReplicates/middle_gastrula.int")
saveRDS(early_neurula.integrated, file = "IntegratedReplicates/early_neurula.int")
saveRDS(late_neurula.integrated, file = "IntegratedReplicates/late_neurula.int")
saveRDS(init_tail.integrated, file = "IntegratedReplicates/init_tail.int")
saveRDS(early_tail.integrated, file = "IntegratedReplicates/early_tail.int")
saveRDS(mid_tail.integrated, file = "IntegratedReplicates/mid_tailII.int")
saveRDS(late_tail.integrated, file = "IntegratedReplicates/late_tail.int")
saveRDS(late_tailII.integrated, file = "IntegratedReplicates/lat_tailII.int")
saveRDS(larval.integrated, file = "IntegratedReplicates/larval.int")

larval <- readRDS("/home/OiPui/share/OiPui/Ciona_data/IntegratedReplicates/larval.int")
svg("/home/OiPui/share/OiPui/Ciona_data/IntegratedReplicates/UMAP/larval_all")

DimPlot(larval, reduction = "umap", group.by = "Tissue_type", label = TRUE, 
        repel = TRUE) + 
    NoLegend()

dev.off()


