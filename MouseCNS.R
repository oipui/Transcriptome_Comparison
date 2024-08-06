library("tidyverse")

df <- read_csv("/home/OiPui/share/OiPui/Ciona_data/SAMap/maps/CiMm/CiMm_GenePairs_n100.csv", col_names = FALSE)

df <- df[, -1]
df_sep <- names(df) %>%
  map(
    function(x) 
      df %>% 
      dplyr::select(x) %>% 
      separate(x, 
               into = paste0(x, c("s1", "s2")), 
               sep = ";")) %>%
  bind_cols()
colnames(df_sep) <- df_sep[1, ]
df_sep <- df_sep[-1,]

#this was to check whether there were differences with new SAMap version -> no 
#library(arsenal)
#summary(comparedf(df_sep,df2_sep))

pairs <- lapply(1:(ncol(df_sep)/2), function(x) df_sep[, c(2 * x - 1, 2 * x)])

IDs <- read_tsv("/home/share/OiPui/Ciona_data/CionaBlast/ENSEMBL_KH_ids.list")
IDs$KH_id <- sub("KH", "Ci_KH2012:KH",IDs$KH_id)

for(i in 1:length(pairs)) { 
  names(pairs)[[i]] = colnames(pairs[[i]][1]) #assigns column names to names of each df in list
  pairs[[i]] = dplyr::rename(pairs[[i]], "ciona_gene_ids" = colnames(pairs[[i]][1])) #renaming
}

for(i in 1:length(pairs)) { 
  pairs[[i]] = inner_join(pairs[[i]], IDs, by= c("ciona_gene_ids" = "KH_id"))
}


for(i in 1:length(pairs)) { 
  colnames(pairs[[i]])[colnames(pairs[[i]]) == "ciona_gene_ids"] <- names(pairs)[[i]]
}

mouseCNS <- readRDS("/home/OiPui/share/OiPui/Ciona_data/SAMap/data/mouse_CNS.rds")


library(ggplot2)
ggplot(mouseCNS, aes(.X, .Y, colour = Species)) +
  geom_point()

DimPlot(mouseCNS, reduction = ".tSNE", group.by = "ClusterNames", label = TRUE) + 
  NoLegend()

test <- ScaleData(mouseCNS, verbose = FALSE)
test <- RunPCA(test, npcs = 30, verbose = FALSE)


