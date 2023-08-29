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

rdf <- function(df, item, imp = T){
  v <- base::setdiff(c("left_level", "right_level", "max_level"), item)
  df <- df %>% dplyr::select(- any_of(v))
  total_level_train.forest <- randomForest(
    x = df %>% dplyr::select(- all_of(item)), 
    y = df[[item]], importance = TRUE, 
    ntree = 100, keep.forest = T,
    maxnodes = 10, mtry = 2)
  if(imp){
    f <- importance(total_level_train.forest) %>% as.data.frame() %>% 
      arrange(desc(MeanDecreaseGini))
    return(f)
  }else{
    return(total_level_train.forest)
  }
}


split_set <- function(df, item, pct = 0.7){
  idt <- sample(1:nrow(df), floor(pct * nrow(df)))
  train_data <- df[idt,]
  test_data <- df[-idt,]
  list(
    train_data = train_data,
    test_data = test_data
  )
}

rdp <- function(modelf, df, item, table = T){
  v <- base::setdiff(c("left_level", "right_level", "max_level"), item)
  df <- df %>% dplyr::select(- any_of(v))
  pd <- predict(modelf, newdata = df)
  if(table){
    crosstable <- tibble(
      predict = pd,
      read_data = df[[item]]
    ) %>% group_by(predict, read_data) %>% 
      dplyr::summarise(number = n()) %>% ungroup()
    return(crosstable)
  }else{
    return(pd)
  }
}




# 2400
tiny_core_frame <- core_frame %>% dplyr::select(tidyselect::contains("level")) %>% 
  dplyr::select(-c("xingbie_level", "total_level")) %>% na.omit() %>% 
  dplyr::select("left_level", "right_level", "max_level", everything())

tiny_core_frame$left_level <- factor(tiny_core_frame$left_level, levels = c("normal", "level1", "level2"))
rdf(tiny_core_frame, "left_level")

impfl <- parallel::mclapply(1:40, function(seed){
  set.seed(seed)
  df <- tiny_core_frame
  impf <- lapply(c("left_level", "right_level", "max_level"), function(item){
    df[[item]] <- as.factor(df[[item]])
    
    split_list <- split_set(df, item)
    train_modle <- rdf(split_list$train_data, item, F)
    prd <- rdp(train_modle, split_list$test_data, item, table = T) 
    true_n <- prd %>% dplyr::filter(predict == read_data) %>% 
      dplyr::summarise(true = sum(number)) %>% pull
    zhunquedu <- true_n / nrow(split_list$test_data)
    
    full_modle <- rdf(df = df, item =  item, imp =  F)
    ipm <- importance(full_modle) %>% as.data.frame() %>% 
      arrange(desc(MeanDecreaseGini))
    tibble(accurate = zhunquedu, importance = list(ipm), full_modle = list(full_modle))
  })
  names(impf) <- c("left_level", "right_level", "max_level")
  impf
}, mc.cores = 4L)
write_rds(impfl, file = "randomf.rds")
bigframe <- list_transpose(impfl) %>% list_rbind(names_to = "id") %>% group_by(id) %>% 
  dplyr::slice_max(order_by = accurate, n = 1)

# =========================== 左侧 =================================
bigframe[1,4] %>% pull %>% .[[1]] -> md
tree <- getTree(md, k = 1, labelVar=TRUE)
tree %>% class
#rename the name of the column
colnames(tree) <- sapply(colnames(tree),collapse)
rules <- getConds(tree)
print(rules)



addf <- function(x){
  paste0("{", x, "}")
}
for(i in 1:ncol(tiny_core_frame)){
  assign(colnames(tiny_core_frame)[i], tiny_core_frame[[i]])
}
all_tree <- parallel::mclapply(1:100, function(nt){
  tree <- getTree(md, k = nt, labelVar=TRUE)
  tree %>% class
  #rename the name of the column
  colnames(tree) <- sapply(colnames(tree),collapse)
  rules <- getConds(tree)
  # print(rules)
  
  tree1 <- lapply(1:length(rules), function(i){
    xstr <- rules[[i]]
    text <- xstr %>% str_remove(pattern = " => .+$") %>% 
      str_replace_all(pattern = paste0("\\b", colnames(tiny_core_frame), "\\b") %>% 
                        paste0(collapse = "|"), replacement = addf)
    ex <- glue::glue(text) %>% as.character()
    df <- lapply(ex, function(exi){
      branch <- eval(parse(text = exi))
      branchy <- xstr %>% str_extract(pattern = "(?<= =>  ).+$") 
      tibble(
        "{branchy}_{i}" := branch
      )
    }) %>% list_rbind()
  }) %>% list_cbind()
  tree1$id <- 1:nrow(tree1)
  tree1_ans <- tree1 %>% pivot_longer(-"id") %>% dplyr::filter(value) %>% 
    tidyr::separate_wider_delim(cols = name , delim = "_", names = c("ans", "branch")) %>% 
    dplyr::select(- "value")
  
  tree1_ans
}, mc.cores = 10L)

trees_frame <- lapply(all_tree, function(x){
  x[2]
}) %>% list_cbind()
trees_frame$id <- 1:nrow(trees_frame)
trees_frame <- trees_frame %>% pivot_longer(-"id")

strees_frame <- trees_frame %>% group_by(id, value) %>% summarise(number = n()) %>% 
  mutate(order = row_number(desc(number))) %>% dplyr::filter(order == 1)
  


trees_error_frame <- lapply(all_tree, function(x){
  data.frame(
    error = as.numeric(x[[2]] == strees_frame$value)
    )
}) %>% list_cbind()
mx <- trees_error_frame %>% as.matrix()
colnames(mx) <- 1:ncol(mx)
pdf("error_heatmap.pdf")
ComplexHeatmap::Heatmap(mx, use_raster = F, show_row_names = F,
                        clustering_method_rows = "ward.D2",
                        clustering_method_columns = "ward.D2")
dev.off()

tibble(
  tree_id = 1:ncol(mx),
  csm = colSums(mx)
) %>% arrange(desc(csm)) %>% 
  as.data.frame() %>% head(20)



trees_true_error_frame <- lapply(all_tree, function(x){
  data.frame(
    error = as.numeric(x[[2]] == tiny_core_frame$left_level)
  )
}) %>% list_cbind()
mx <- trees_true_error_frame %>% as.matrix()
colnames(mx) <- 1:ncol(mx)
pdf("true_error_heatmap.pdf")
ComplexHeatmap::Heatmap(mx, use_raster = F, show_row_names = F,
                        clustering_method_rows = "ward.D2",
                        clustering_method_columns = "ward.D2")
dev.off()


colnames(trees_true_error_frame) <- 1:ncol(trees_true_error_frame)
trees_true_error_framex <- trees_true_error_frame %>% mutate(lb = tiny_core_frame$left_level) %>% 
  pivot_longer(-"lb") %>% group_by(lb, name) %>% summarise(true_number = sum(value)) %>% 
  arrange(name, lb)
fullan <- tibble(
  lb = tiny_core_frame$left_level
) %>% group_by(lb) %>% summarise(total_nub = n())
truepct <- trees_true_error_framex %>% left_join(fullan) %>% mutate(pct = true_number/total_nub)
truepct %>% dplyr::filter(lb == "level1") %>% arrange(desc(pct))
truepct %>% dplyr::filter(lb == "level2") %>% arrange(desc(pct))


tibble(
  tree_id = 1:ncol(mx),
  csm = colSums(mx)
) %>% arrange(desc(csm)) %>% 
  as.data.frame() %>% head(20)
# 53
# 7
# 15
# 34
# 35
# 43
# 46


tree <- getTree(md, k = 18, labelVar=TRUE)
tree %>% class
#rename the name of the column
colnames(tree) <- sapply(colnames(tree),collapse)
rules <- getConds(tree)
print(rules)

tree <- getTree(md, k = 4, labelVar=TRUE)
tree %>% class
#rename the name of the column
colnames(tree) <- sapply(colnames(tree),collapse)
rules <- getConds(tree)
print(rules)

tree <- getTree(md, k = 39, labelVar=TRUE)
tree %>% class
#rename the name of the column
colnames(tree) <- sapply(colnames(tree),collapse)
rules <- getConds(tree)
print(rules)

# =========================== 二侧 =================================
bigframe[2,4] %>% pull %>% .[[1]] -> md


addf <- function(x){
  paste0("{", x, "}")
}
for(i in 1:ncol(tiny_core_frame)){
  assign(colnames(tiny_core_frame)[i], tiny_core_frame[[i]])
}
all_tree <- parallel::mclapply(1:100, function(nt){
  tree <- getTree(md, k = nt, labelVar=TRUE)
  tree %>% class
  #rename the name of the column
  colnames(tree) <- sapply(colnames(tree),collapse)
  rules <- getConds(tree)
  # print(rules)
  
  tree1 <- lapply(1:length(rules), function(i){
    xstr <- rules[[i]]
    text <- xstr %>% str_remove(pattern = " => .+$") %>% 
      str_replace_all(pattern = paste0("\\b", colnames(tiny_core_frame), "\\b") %>% 
                        paste0(collapse = "|"), replacement = addf)
    ex <- glue::glue(text) %>% as.character()
    df <- lapply(ex, function(exi){
      branch <- eval(parse(text = exi))
      branchy <- xstr %>% str_extract(pattern = "(?<= =>  ).+$") 
      tibble(
        "{branchy}_{i}" := branch
      )
    }) %>% list_rbind()
  }) %>% list_cbind()
  tree1$id <- 1:nrow(tree1)
  tree1_ans <- tree1 %>% pivot_longer(-"id") %>% dplyr::filter(value) %>% 
    tidyr::separate_wider_delim(cols = name , delim = "_", names = c("ans", "branch")) %>% 
    dplyr::select(- "value")
  
  tree1_ans
}, mc.cores = 10L)

trees_frame <- lapply(all_tree, function(x){
  x[2]
}) %>% list_cbind()
trees_frame$id <- 1:nrow(trees_frame)
trees_frame <- trees_frame %>% pivot_longer(-"id")

strees_frame <- trees_frame %>% group_by(id, value) %>% summarise(number = n()) %>% 
  mutate(order = row_number(desc(number))) %>% dplyr::filter(order == 1)



trees_error_frame <- lapply(all_tree, function(x){
  data.frame(
    error = as.numeric(x[[2]] == strees_frame$value)
  )
}) %>% list_cbind()
mx <- trees_error_frame %>% as.matrix()
colnames(mx) <- 1:ncol(mx)
pdf("error_heatmap.pdf")
ComplexHeatmap::Heatmap(mx, use_raster = F, show_row_names = F,
                        clustering_method_rows = "ward.D2",
                        clustering_method_columns = "ward.D2")
dev.off()

tibble(
  tree_id = 1:ncol(mx),
  csm = colSums(mx)
) %>% arrange(desc(csm)) %>% 
  as.data.frame() %>% head(20)



trees_true_error_frame <- lapply(all_tree, function(x){
  data.frame(
    error = as.numeric(x[[2]] == tiny_core_frame$left_level)
  )
}) %>% list_cbind()
mx <- trees_true_error_frame %>% as.matrix()
colnames(mx) <- 1:ncol(mx)
pdf("true_error_heatmap.pdf")
ComplexHeatmap::Heatmap(mx, use_raster = F, show_row_names = F,
                        clustering_method_rows = "ward.D2",
                        clustering_method_columns = "ward.D2")
dev.off()


colnames(trees_true_error_frame) <- 1:ncol(trees_true_error_frame)
trees_true_error_framex <- trees_true_error_frame %>% mutate(lb = tiny_core_frame$left_level) %>% 
  pivot_longer(-"lb") %>% group_by(lb, name) %>% summarise(true_number = sum(value)) %>% 
  arrange(name, lb)
fullan <- tibble(
  lb = tiny_core_frame$left_level
) %>% group_by(lb) %>% summarise(total_nub = n())
truepct <- trees_true_error_framex %>% left_join(fullan) %>% mutate(pct = true_number/total_nub)
truepct %>% dplyr::filter(lb == "level1") %>% arrange(desc(pct))
truepct %>% dplyr::filter(lb == "level2") %>% arrange(desc(pct))


tibble(
  tree_id = 1:ncol(mx),
  csm = colSums(mx)
) %>% arrange(desc(csm)) %>% 
  as.data.frame() %>% head(20)






tree <- getTree(md, k = 61, labelVar=TRUE)
tree %>% class
#rename the name of the column
colnames(tree) <- sapply(colnames(tree),collapse)
rules <- getConds(tree)
print(rules)

tree <- getTree(md, k = 39, labelVar=TRUE)
tree %>% class
#rename the name of the column
colnames(tree) <- sapply(colnames(tree),collapse)
rules <- getConds(tree)
print(rules)

tree <- getTree(md, k = 74, labelVar=TRUE)
tree %>% class
#rename the name of the column
colnames(tree) <- sapply(colnames(tree),collapse)
rules <- getConds(tree)
print(rules)



