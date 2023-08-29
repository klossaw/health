pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", "factoextra", "cluster",
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "Seurat", "ComplexHeatmap",
          "ggpointdensity", "ggpubr", "randomForest")
for (pkg in pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "healthman"
dataset <- "zhanglei"
species <- "human"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/") %>% checkdir()
setwd(workdir)
core_frame <- read_csv("updata_core.csv") %>% dplyr::select(-1)
library(rpart)
rddt <- function(modelf, df, item, table = T){
  v <- base::setdiff(c("left_level", "right_level", "max_level"), item)
  df <- df %>% dplyr::select(- any_of(v))
  pd <- predict(modelf, newdata = df)
  if(table){
    crosstable <- tibble(
      predict = apply(pd,1,function(x) c("level1", "level2", "normal")[which.max(x)]),
      read_data = df[[item]]
    ) %>% group_by(predict, read_data) %>% 
      dplyr::summarise(number = n()) %>% ungroup()
    return(crosstable)
  }else{
    return(pd)
  }
}
split_set <- function(df, item, pct = 0.7){
  idt <- df %>% group_by("{item}":=df[[item]]) %>% 
    slice_sample(prop = pct)
  
  train_data <- idt
  test_data <- anti_join(df, idt)
  list(
    train_data = train_data,
    test_data = test_data
  )
}
rdt <- function(df, item, imp = T){
  v <- base::setdiff(c("left_level", "right_level", "max_level"), item)
  df <- df %>% dplyr::select(- any_of(v))
  total_level_train.forest <- rpart(formula = glue::glue("{item} ~ .") %>% as.formula(),
                                    data = df, 
                                    model = TRUE)
  if(imp){
    f <- total_level_train.forest$variable.importance %>% as.data.frame() 
    colnames(f) <- "importance"
    f <- f %>% arrange(desc(importance))
    return(f)
  }else{
    return(total_level_train.forest)
  }
}

# 2400
tiny_core_frame <- core_frame %>% dplyr::select(tidyselect::contains("level")) %>% 
  dplyr::select(-c("xingbie_level", "total_level")) %>% na.omit() %>% 
  dplyr::select("left_level", "right_level", "max_level", everything())
tiny_core_frame$left_level <- factor(tiny_core_frame$left_level, levels = c("normal", "level1", "level2"))


bigtreelist <- lapply(1:20, function(x){
  set.seed(x)
  impf <- lapply(c("left_level", "right_level", "max_level"), function(item){
    df <- tiny_core_frame
    df <- df %>% dplyr::select(- "NL_level")
    df[[item]] <- as.factor(df[[item]])
    split_list <- split_set(df, item)
    train_modle <- rdt(split_list$train_data, item, F) 
    prd <- rddt(train_modle, split_list$test_data, item, table = T) 
    true_n <- prd %>% dplyr::filter(predict == read_data) %>% 
      dplyr::summarise(true = sum(number)) %>% pull
    zhunquedu <- true_n / nrow(split_list$test_data)
    full_modle <- rdt(df = df, item =  item, imp =  F)
    ipm <- full_modle$variable.importance %>% as.data.frame() 
    colnames(ipm) <- "importance"
    ipm <- ipm %>% arrange(desc(importance))
    tibble(accurate = zhunquedu, importance = list(ipm), full_modle = list(full_modle))
  })
  names(impf) <- c("left_level", "right_level", "max_level")
  impf
})
saveRDS(bigtreelist, file = "./treelist.rds")

bigframetree <- list_transpose(bigtreelist) %>% list_rbind(names_to = "id") %>% group_by(id) %>% 
  dplyr::slice_max(order_by = accurate, n = 1)
bigframetree
library(rpart.plot) 
rpart.plot(bigframetree[5,"full_modle"][[1]][[1]])















