rm(list = ls())

#IMPORTANT KEEP THE ORDER OF LOADING LIBRARY
library(sceasy)
library(reticulate)
library(stringr)
library(Seurat)

path <- "~/share/OiPui/Ciona_data/SAMap/data/"
project_name <- "CNS_larval_modified"
# Load Seurat object
chromium_seurat <- readRDS(file = str_c(path, project_name, ".rds"))
# DefaultAssay(chromium_seurat) <- "RNA"
# gene_trans <- read_tsv("~/share/Francesco/Bassi/SAMap/gene_trans/mouse.tsv")
# chromium_seurat <- subset(chromium_seurat, features = gene_trans$external_gene_name)
setwd(path)

sceasy::convertFormat(chromium_seurat, 
                      from = "seurat", 
                      to = "anndata",
                      outFile = str_c(project_name, ".h5ad"))

gastrula_to_larval <- readRDS(file = "~/share/OiPui/Ciona_data/NS_ReplicatesIntegration/subsequent_Integration/data/gastrula_to_larval_ligerIntegration.rds")


sceasy::convertFormat(gastrula_to_larval, 
                      from = "seurat", 
                      to = "anndata",
                      outFile = "~/share/OiPui/Ciona_data/NS_ReplicatesIntegration/subsequent_Integration/data/Ci_allStages_ligerInt.h5ad")