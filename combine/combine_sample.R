pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", "factoextra", "cluster",
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "Seurat", "ComplexHeatmap",
          "ggpointdensity", "ggpubr")
for (pkg in pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "healthman"
dataset <- "zhanglei"
species <- "human"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/") %>% checkdir()
setwd(workdir)
# ================================ 第一部分 过滤离散化 ====================================
# 9133 个
fixsheet <- readxl::read_excel("~/projects/healthman/analysis/zhanglei/human/fix/fix.xlsx", 
                                col_types = c("text", "text", "text", "guess", "numeric","text","text","numeric","numeric", 
                                              rep("guess", times = 75)))
rawobj <-  "~/projects/healthman/data/zhanglei/human/clinical/raw.rds" %>%
  read_rds
rawobj <- rawobj[rawobj$"是否采纳" == "采纳",]
for(i in fixsheet$"匹配"){
  rawobj[rawobj$"匹配" == i,c("右侧颈动脉", "左侧颈动脉" )] <- fixsheet[fixsheet$"匹配" == i,c("右侧颈动脉", "左侧颈动脉")]
}
nac <- apply(rawobj, 2, function(x) sum(is.na(x)))
rawobj <- rawobj[,nac != nrow(rawobj)]
c("脂肪肝",
  "尿素","同型半胱氨酸","甘油三酯",
  "总胆固醇","高密度脂蛋白C","低密度脂蛋白C",
  "载脂蛋白A1","载脂蛋白B","尿酸",
  "空腹血糖","糖化血红蛋白A1","糖化血红蛋白A1C",
  "空腹C肽","空腹胰岛素","超敏C反应蛋白",
  "唾液酸","游离脂肪酸","载脂蛋白E","羟基维生素D3")
fixsheet <- function(sheet){
  sheet <- sheet %>% dplyr::mutate(
    right = `右侧颈动脉` %>% as.numeric(),
    left = `左侧颈动脉` %>% as.numeric(),
    NHDL = as.numeric(`总胆固醇`) - as.numeric(`高密度脂蛋白C`),
    NL =  `NL` %>% str_remove("岁") %>% as.numeric(),
    LDL =  `低密度脂蛋白C` %>% as.numeric(),
    NHDL_NL = NHDL * NL,
    LDL_NL = LDL * NL,
    smoking =  `吸烟`,
    shousuoya =  `收缩压` %>% as.numeric(),
    shuzhangya =  `舒张压` %>% as.numeric(),
    bmi =  `体重指数` %>% as.numeric(),
    fuwei =  `腹围` %>% as.numeric(),
    xingbie = `XB`,
    left_level = case_when(
      is.na(left) ~ NA,
      left < 1 ~ "normal",
      left < 1.5 ~ "level1",
      T  ~ "level2"
    ),
    right_level = case_when(
      is.na(right) ~ NA,
      right < 1 ~ "normal",
      right < 1.5 ~ "level1",
      T ~ "level2"
    ),
    smoking_level = case_when(
      is.na(smoking) ~ NA,
      smoking == "不吸" ~ 0,
      T ~ 1),
    NL_level = floor((NL - 20)/10) + 1,
    fuwei_level = case_when(
      is.na(fuwei) | is.na(xingbie) ~ NA,
      (fuwei < 90 & xingbie == "男") | (fuwei < 85 & xingbie == "女") ~ 0,
      T ~ 1
    ),
    bmi_level = case_when(
      is.na(bmi) ~ NA,
      bmi < 24 ~ 0,
      bmi < 28 ~ 1,
      T ~ 2
    ),
    LDL_level = case_when(
      is.na(LDL) ~ NA,
      LDL < 1.8 ~ 0,
      LDL < 2.6 ~ 1,
      LDL < 3.4 ~ 2,
      LDL < 4.9 ~ 3,
      T ~ 4
    ),
    xueya_level = case_when(
      is.na(shousuoya) | is.na(shuzhangya) ~ NA,
      (shousuoya < 130 & shuzhangya < 80) ~ 0,
      (shousuoya < 140 & shuzhangya < 90) ~ 1,
      (shousuoya < 160 & shuzhangya < 100) ~ 2,
      T ~ 3
    ),
    xingbie_level = case_when(
      is.na(xingbie) ~ NA,
      xingbie == "女" ~ 0,
      xingbie == "男" ~ 1
    )
  )
  sheet <- sheet %>% rowwise() %>% mutate(max_data = max(left, right))
  sheet <- sheet %>% dplyr::mutate(max_level = case_when(
    is.na(max_data) ~ NA,
    max_data < 1 ~ "normal",
    max_data < 1.5 ~ "level1",
    T ~ "level2"
  ))
  sheet$total_level <- rowSums(sheet[,c("smoking_level", "NL_level",
                                          "fuwei_level", "bmi_level",
                                          "LDL_level" , "xueya_level" )])
  sheet
}
rawobj <- fixsheet(rawobj)
rawobj[is.na(rawobj$right),]$"右侧颈动脉" <- c("0.4", "0.5")
rawobj <- fixsheet(rawobj)
for(col_i in colnames(rawobj[,-c(1, 2, 3, 4)])){
  if(sum(as.numeric(rawobj[[col_i]]) %>% is.na()) != nrow(rawobj)){
    rawobj[[col_i]] <- as.numeric(rawobj[[col_i]])
  }
}
rawobj <- rawobj[!is.na(rawobj$right) & !is.na(rawobj$left),]

rawobj$"脂肪肝"[rawobj$"脂肪肝" == 5] <- NA
rawobj$NL_level[rawobj$NL_level >= 6] <- 6 
core_frame <- rawobj[,c(1,2,3,4, 13, 61:84)]


dir.create("./combine")
write_rds(rawobj, file = "./combine/updata_raw.rds")
write_rds(core_frame, file = "./combine/updata_core.rds")










