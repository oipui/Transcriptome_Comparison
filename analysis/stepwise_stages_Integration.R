#Liger

library('rliger')

path <- "/home/share/OiPui/Ciona_data/data/"
out_path <- "~/share/OiPui/Ciona_data/NS_ReplicatesIntegration/subsequent_Integration/"

#this will create a df 
ciona <- read_tsv(file="expression_matrix_10stage.tsv.gz")
ciona <- as.data.frame(ciona)
#This assigns rownames to the data.frame
#and then removes the first column which contains the gene names that we just moved to rownames
rownames(ciona) <- ciona$GENE
ciona <- ciona[, -1]

#transform into sparse matrix
ciona_matrix <- Matrix(as.matrix(ciona), sparse = TRUE)

ciona_meta1 <- read_tsv(file = "ciona10stage.cluster.upload.new.txt")
ciona_meta1 <- ciona_meta1 %>% 
  dplyr::filter(NAME != "TYPE") %>% 
  dplyr::rename("Tissue_type" = 'Tissue Type') %>% 
  dplyr::filter(Tissue_type =="nervous system")

ciona_meta2 <- read_tsv(file= "ciona10stage.meta.upload.new.MSTRG.txt")
ciona_meta2 <- ciona_meta2 %>% 
  dplyr::filter(NAME != "TYPE") 

ciona_metadata <- ciona_meta1 %>% inner_join(ciona_meta2, by = "NAME")

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
  print(name)
  stages.list[[name]] <- as.data.frame(stages.list[[name]])
  rownames(stages.list[[name]]) <- stages.list[[name]]$NAME
  matrix <- ciona_matrix[,stages.list[[name]]$NAME]
  meta = stages.list[[name]]
  print(SeuratObjects[[name]])
  SeuratObjects[[name]] <- CreateSeuratObject(matrix, meta.data = meta)
  #SeuratObjects[[name]] <- subset(SeuratObjects[[name]], subset = 'Tissue Type' == "nervous system")
}

list2env(SeuratObjects,envir=.GlobalEnv)

#CNS <- readRDS(file ="/home/OiPui/share/OiPui/Ciona_data/SAMap/data/Ci_CNS_larval.rds")
CNS.lv <- readRDS("~/share/OiPui/Ciona_data/SAMap/data/CNS_larval_modified.rds")

####################################################
C110 <- readRDS(file = "~/share/OiPui/Ciona_data/NS_ReplicatesIntegration/subsequent_Integration/C110_NS.integrated.rds" )
DefaultAssay(C110) <- "RNA"
a_b <- list(C110, b.middle_gastrula)
names(a_b) <- c("C110", "midGastrula")

#create liger object using the list of matrices
ciona_liger <- rliger::seuratToLiger(a_b, names = c('initGastrula','midGastrula'),remove.missing = FALSE)

#preprocessing and normalization
#normalization of the expression to account for differences in sequencing depth and efficiency between cells
#identify variably expressed genes 
#scale date so each gene has the same variance
#NOTE: because of nonnegative matrix factorization requires positive values, no centering of data by subtracting the mean 
# also no log transformation of the data 
ciona_liger <- rliger::normalize(ciona_liger)
ciona_liger <- rliger::selectGenes(ciona_liger)
ciona_liger <- rliger::scaleNotCenter(ciona_liger, remove.missing = FALSE)


#Joint matrix factorization
ciona_liger_k17 <- optimizeALS(ciona_liger ,k = 17)

#Error with cell data not being updated when removing the genes in previous step 
#solution:https://github.com/welch-lab/liger/issues/217
names<-rbind(as.matrix(row.names(ciona_liger_k17@scale.data[[1]])), as.matrix(row.names(ciona_liger_k17@scale.data[[2]])))
tmp<- ciona_liger_k17@cell.data[names,]
cell.data.orig <- ciona_liger_k17@cell.data
ciona_liger_k17@cell.data <- tmp

ciona_liger_k17 <- louvainCluster(ciona_liger_k17, resolution = 0.3)

ciona_liger_k17 <- runUMAP(ciona_liger_k17, use.raw = TRUE, dims.use = 1:17)


plots <- plotByDatasetAndCluster(ciona_liger_k17, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = T)
t <- plots[[1]] + plots[[2]]


#########################################################

library(SeuratWrappers)

gastrula_neurula <- merge(gastrula_earNeurula, y= d.late_neurula)
gastrula_neurula <- NormalizeData(gastrula_neurula)
gastrula_neurula <- FindVariableFeatures(gastrula_neurula)
gastrula_neurula <- ScaleData(gastrula_neurula, split.by = "stage", do.center = FALSE)
gastrula_neurula <- RunOptimizeALS(gastrula_neurula, k =20, lambda = 5, split.by = "stage")
gastrula_neurula <- FindNeighbors(gastrula_neurula, reduction = "iNMF_raw", dims = 1:20)
gastrula_neurula <- FindClusters(gastrula_neurula, resolution = 0.3)
# Dimensional reduction and plotting
gastrula_neurula <- RunUMAP(gastrula_neurula, dims = 1:ncol(test[["iNMF_raw"]]), reduction = "iNMF_raw")
DimPlot(gastrula_neurula, group.by = c("stage", "ident", "Cell_type"), ncol = 3, label = TRUE)
DimPlot(gastrula_neurula, group.by = c("liger_cluster","Cell_type2", "Cell_type"), ncol = 3, label = TRUE) +
  NoLegend()

gastrula_neurula$liger_cluster <- Idents(gastrula_neurula)
gastrula_neurula$Cell_type2 <- Idents(gastrula_neurula)
Idents(gastrula_neurula) <- "Cell_type2"

cells.use <- WhichCells(object = gastrula_neurula, idents = '12')
gastrula_neurula <- SetIdent(object = gastrula_neurula, cells = cells.use, value = 'b.1.7')

gastrula_neurula$Cell_type2 <- Idents(gastrula_neurula)

saveRDS(init_midGastrula, file = "~/share/OiPui/Ciona_data/NS_ReplicatesIntegration/subsequent_Integration/init_midGastrula_ligerIntegration.rds")

########################################################

gastrula_to_lateTailII$Tissue_type <- NA

gastrula_to_larval <- merge(gastrula_to_lateTailII, y= CNS.lv)
gastrula_to_larval <- NormalizeData(gastrula_to_larval)
gastrula_to_larval <- FindVariableFeatures(gastrula_to_larval)
gastrula_to_larval <- ScaleData(gastrula_to_larval, split.by = "stage", do.center = FALSE)
gastrula_to_larval <- RunOptimizeALS(gastrula_to_larval, k =20, lambda = 5, split.by = "stage")
gastrula_to_larval <- FindNeighbors(gastrula_to_larval, reduction = "iNMF_raw", dims = 1:20)
gastrula_to_larval <- FindClusters(gastrula_to_larval, resolution = 0.3)
# Dimensional reduction and plotting
gastrula_to_larval <- RunUMAP(gastrula_to_larval, dims = 1:ncol(gastrula_to_larval[["iNMF_raw"]]), reduction = "iNMF_raw")
DimPlot(gastrula_to_larval, group.by = c("stage", "ident", "Cell_type"), ncol = 3, label = TRUE)

gastrula_to_larval$liger_cluster<- Idents(gastrula_to_larval)
gastrula_to_larval@meta.data$Cell_type2 <- ifelse(is.na(gastrula_to_larval@meta.data$Cell_type2), gastrula_to_larval@meta.data$liger_cluster, gastrula_to_larval@meta.data$Cell_type2)
gastrula_to_larval@meta.data$Cell_type3 <- ifelse(str_detect(gastrula_to_larval$Cell_type2, "^[:digit:]+$"), gastrula_to_larval@meta.data$Tissue_type, test@meta.data$Cell_type2)
test@meta.data$Cell_type2 <- ifelse(is.na(gastrula_to_larval$Cell_type2), test@meta.data$liger_cluster, test@meta.data$Cell_type2)

DimPlot(gastrula_to_lateTailII, group.by = c("Cell_type","liger_cluster","Cell_type2"), ncol = 3, label = TRUE)
DimPlot(test, group.by = "liger_cluster",label = TRUE) +
  NoLegend()
DimPlot(gastrula_to_larval, group.by = "stage",label = FALSE)
DimPlot(gastrula_to_larval, group.by = "Tissue_type",label = TRUE) +
  NoLegend()
DimPlot(gastrula_to_larval, group.by = "Cell_type3", label = TRUE)+
  NoLegend()


Idents(gastrula_to_larval) <- "Cell_type2"

DimPlot(gastrula_to_larval, group.by = c("ident", "Tissue_type"), ncol = 2, label = TRUE)

cells.use <- WhichCells(object = gastru, idents = '18')
gastrula_to_larval <- SetIdent(object = gastrula_to_larval, cells = cells.use, value = 'a.8.1')
DimPlot(gastrula_to_larval, group.by = "stage", label = TRUE)

gastrula_to_larval$Cell_type2 <- Idents(gastrula_to_larval)

test <- readRDS(file = "~/share/OiPui/Ciona_data/NS_ReplicatesIntegration/subsequent_Integration/gastrula_to_lateTailII_ligerIntegration.rds")

saveRDS(gastrula_to_lateTailII, file = "~/share/OiPui/Ciona_data/NS_ReplicatesIntegration/subsequent_Integration/gastrula_to_lateTailII_ligerIntegration.rds")


##########################
#changing cell types again 

#a.1.2, b.1.3, b.1.2 separate in initial gastrula stage from a and b lineage prospectively to PNS 

gastrula_to_larval$Cell_type3.1 <- gastrula_to_larval$Cell_type3
Idents(gastrula_to_larval) <- "Cell_type3.1"
cells.use <- WhichCells(object = gastrula_to_larval, idents = 'a.7.0')
gastrula_to_larval <- SetIdent(object = gastrula_to_larval, cells = cells.use, value = "PNS.7.0")
DimPlot(gastrula_to_larval, group.by = "ident", label = TRUE) + 
  NoLegend()

gastrula_to_larval$Cell_type3.1 <- Idents(gastrula_to_larval)
gastrula_to_larval$Cell_type3.1 <- as.character(gastrula_to_larval$Cell_type3.1)

DimPlot(gast_neur_initTail, group.by = c("Cell_type2","liger_cluster","Cell_type"),ncol = 3, label = TRUE) 
+
  NoLegend()

DimPlot(gastrula_to_larval, group.by = "Cell_lineages",label = TRUE) + 
  NoLegend()

###############
#somethings went wrong with the adding annotations will add them here again 
cell <- as.character(gastrula_to_lateTailII$Cell_type2[gastrula_to_lateTailII$stage == "i.late tailbud II"])
gastrula_to_larval$Cell_type2 <- ifelse(str_detect(gastrula_to_larval$Cell_type2, "^[:digit:]+$"), cell, gastrula_to_larval$Cell_type2)
gastrula_to_larval$Cell_type3.1 <-ifelse(is.na(gastrula_to_larval$Cell_type3.1), gastrula_to_larval$Cell_type2, gastrula_to_larval$Cell_type3.1)

WhichCells(object = gastrula_to_larval, idents = "ependymal cells")

DimPlot(gastrula_to_larval, group.by = "ident",label = TRUE, cells.highlight = WhichCells(object = gastrula_to_larval, idents = c("KCNB1+ motor ganglion","AMD+ motor ganglion", "GLRA1/2/3+ motor ganglion"))) +
  NoLegend()

DimPlot(gastrula_to_larval, group.by = "Cell_type3.1",label = TRUE, cells.highlight = WhichCells(object = gastrula_to_larval, idents = "b.2.4")) +
  NoLegend()

DimPlot(gastrula_to_larval, group.by = "Cell_type3.1",label = TRUE, cells.highlight = WhichCells(object = gastrula_to_larval, idents = c("KCNB1+ motor ganglion","AMD+ motor ganglion", "GLRA1/2/3+ motor ganglion"))) +
  NoLegend()

strings <- as.vector(gastrula_to_larval$Cell_type3.1[startsWith(gastrula_to_larval$Cell_type3.1, "PNS.")])
PNS <- subset(gastrula_to_larval, subset= Cell_type3.1 == unique(strings)| Cell_lineages == "PNS")
DimPlot(PNS, group.by = c("Cell_type3.1","stage"), ncol = 2, label = TRUE)

##############
#rename Cell lineages

gastrula_to_larval$Cell_lineages <- ifelse(is.na(gastrula_to_larval$Cell_lineages), gastrula_to_larval$Cell_type3.1, gastrula_to_larval$Cell_lineages)


Idents(gastrula_to_larval) <- "Cell_type3.1"
gastrula_to_larval$Cell_type3.1 <- as.character(gastrula_to_larval$Cell_type3.1)
strings <- as.vector(gastrula_to_larval$Cell_type3.1[startsWith(gastrula_to_larval$Cell_type3.1, "PNS.")])
cells.use <- WhichCells(gastrula_to_larval, idents = unique(strings))
gastrula_to_larval <- SetIdent(object = gastrula_to_larval, cells = cells.use, value = "PNS")

strings <- unique(grep("^a\\.[A-z]", gastrula_to_larval$Cell_type3.1, value = TRUE))
cells.use <- WhichCells(gastrula_to_larval, idents = strings)
gastrula_to_larval <- SetIdent(object = gastrula_to_larval, cells = cells.use, value = "a-lineage")

strings <- unique(grep("^A\\.[A-z]", gastrula_to_larval$Cell_lineages, value = TRUE))
cells.use <- WhichCells(gastrula_to_larval, idents = strings)
gastrula_to_larval <- SetIdent(object = gastrula_to_larval, cells = cells.use, value = "A-lineage")
gastrula_to_larval$Cell_lineages1 <- Idents(gastrula_to_larval)

Idents(gastrula_to_larval) <- "Cell_lineages"
strings <- unique(grep("[A-z]-lineage|PNS", gastrula_to_larval$Cell_lineages, value = TRUE, invert = TRUE))
strings <- unique(grep("PNS\\.", gastrula_to_larval$Cell_lineages, value = TRUE))

cells.use <- WhichCells(gastrula_to_larval, idents = strings)
gastrula_to_larval <- SetIdent(object = gastrula_to_larval, cells = cells.use, value = NA)

gastrula_to_larval$Cell_lineages <- Idents(gastrula_to_larval)

DimPlot(gastrula_to_larval,group.by = "Cell_lineages1", label = TRUE) 




saveRDS(gastrula_to_larval, file = "~/share/OiPui/Ciona_data/NS_ReplicatesIntegration/subsequent_Integration/gastrula_to_larval_ligerIntegration.rds")

#################
gastrula_to_larval$Cell_lineages1.1 <- as.character(gastrula_to_larval$Cell_lineages1)

Idents(gastrula_to_larval) <- "Cell_lineages1.1"

strings <- unique(grep("^A\\.\\d+", gastrula_to_larval$Cell_lineages1.1, value = TRUE))
cells.use <-WhichCells(gastrula_to_larval, idents = strings)
gastrula_to_larval <- SetIdent(object = gastrula_to_larval, cells = cells.use, value = "maybe.A-lineage")

strings <- unique(grep("^a\\.\\d+", gastrula_to_larval$Cell_lineages1.1, value = TRUE))
cells.use <-WhichCells(gastrula_to_larval, idents = strings)
gastrula_to_larval <- SetIdent(object = gastrula_to_larval, cells = cells.use, value = "maybe.a-lineage")

gastrula_to_larval$Cell_lineages1.1 <- Idents(gastrula_to_larval)

gastrula_to_larval$lineages_final <- gastrula_to_larval$Level3

gastrula_to_larval$Cell_lineages1.1 <- as.character(gastrula_to_larval$Cell_lineages1.1)
gastrula_to_larval$lineages_final <- ifelse(is.na(gastrula_to_larval$lineages_final), gastrula_to_larval$Cell_lineages1.1 ,gastrula_to_larval$lineages_final)

DimPlot(gastrula_to_larval,group.by = "Cell_lineages1.1", label = TRUE, cells.highlight = WhichCells(object = gastrula_to_larval, idents = "PNS")) 





