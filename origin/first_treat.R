pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", "factoextra", "cluster", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "Seurat", "ComplexHeatmap")
for (pkg in pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "healthman"
dataset <- "chenlei"
species <- "human"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/") %>% checkdir()
setwd(workdir)

paths <- list.files("~/projects/healthman/data/chenlei/human/clinical/raw", full.names = T)
names(paths) <- paste0(str_extract(paths, pattern = "[0-9]+(?=月)"), "_month")
excellist <- lapply(paths, function(path){
  readxl::read_excel(path, 
                     col_types = c("text", "text", "text", "guess", "numeric","text","text","numeric","numeric", 
                                   rep("guess", times = 53)))
}) %>% list_rbind(names_to = "month") %>% dplyr::select(-2)


