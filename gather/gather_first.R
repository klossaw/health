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
rawobj <-  "/cluster/home/jhuang/projects/healthman/data/chenlei/human/clinical/raw.rds" %>%
  read_rds
rawobj <- rawobj[rawobj$"是否采纳" == "采纳",]
nac <- apply(rawobj, 2, function(x) sum(is.na(x)))
rawobj <- rawobj[,nac != nrow(rawobj)]
c("尿素","同型半胱氨酸","甘油三酯",
  "总胆固醇","高密度脂蛋白C","低密度脂蛋白C",
  "载脂蛋白A1","载脂蛋白B","尿酸",
  "空腹血糖","糖化血红蛋白A1","糖化血红蛋白A1C",
  "空腹C肽","空腹胰岛素","超敏C反应蛋白",
  "唾液酸","游离脂肪酸","载脂蛋白E","羟基维生素D3")
{
  rawobj <- rawobj %>% dplyr::mutate(
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
  rawobj <- rawobj %>% rowwise() %>% mutate(max_data = max(left, right))
  rawobj <- rawobj %>% dplyr::mutate(max_level = case_when(
    is.na(max_data) ~ NA,
    max_data < 1 ~ "normal",
    max_data < 1.5 ~ "level1",
    T ~ "level2"
  ))
  rawobj$total_level <- rowSums(rawobj[,c("smoking_level", "NL_level",
                                          "fuwei_level", "bmi_level",
                                          "LDL_level" , "xueya_level" )])
}

rawobj_numeric <- apply(rawobj, 2, as.numeric)
cna_number <- apply(rawobj_numeric, 2, function(x) sum(is.na(x)))
rawobj_numeric <- rawobj_numeric[,cna_number != max(cna_number)]
rawobj_numeric <- rawobj_numeric[,-c(1, 2)] %>% as_tibble()

for(i in colnames(rawobj_numeric)){
  rawobj[[i]] <- rawobj_numeric[[i]]
}
rawobj <- rawobj[!is.na(rawobj$right) & !is.na(rawobj$left),]
core_frame <- rawobj[,c(13, 61:84)]
write.csv(rawobj, file = "updata_raw.csv")
write.csv(core_frame, file = "updata_core.csv")

# pdf("score_boxplot.pdf")
# my_comparisons <- list(c("0", "1"))
# ggboxplot(leftight_level, x = "smoking_level", y = "left", fill = "smoking_level",
#           palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
#   stat_compare_means(comparisons = my_comparisons)
# ggboxplot(leftight_level, x = "smoking_level", y = "right", fill = "smoking_level",
#           palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
#   stat_compare_means(comparisons = my_comparisons)
#
# ggboxplot(leftight_level, x = "NL_level", y = "left", fill = "NL_level")+
#   stat_compare_means(method = "anova")
# ggboxplot(leftight_level, x = "NL_level", y = "right", fill = "NL_level")+
#   stat_compare_means(method = "anova")
#
# ggboxplot(leftight_level, x = "fuwei_level", y = "left", fill = "fuwei_level",
#           palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
#   stat_compare_means(comparisons = my_comparisons)
# ggboxplot(leftight_level, x = "fuwei_level", y = "right", fill = "fuwei_level",
#           palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
#   stat_compare_means(comparisons = my_comparisons)
#
# my_comparisons <- list(c("0", "1"),c("1","2"),c("0","2"))
# ggboxplot(leftight_level, x = "bmi_level", y = "left", fill = "bmi_level",
#           palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
#   stat_compare_means(comparisons = my_comparisons)
# ggboxplot(leftight_level, x = "bmi_level", y = "right", fill = "bmi_level",
#           palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
#   stat_compare_means(comparisons = my_comparisons)
#
# ggboxplot(leftight_level, x = "LDL_level", y = "left", fill = "LDL_level")+
#   stat_compare_means(method = "anova")
# ggboxplot(leftight_level, x = "LDL_level", y = "right", fill = "LDL_level")+
#   stat_compare_means(method = "anova")
#
# ggboxplot(leftight_level, x = "xueya_level", y = "left", fill = "xueya_level")+
#   stat_compare_means(method = "anova")
# ggboxplot(leftight_level, x = "xueya_level", y = "right", fill = "xueya_level")+
#   stat_compare_means(method = "anova")
# dev.off()
# sum(is.na(rawobj$left))
# sum(is.na(rawobj$right))
# sum(is.na(rawobj$smoking))
# sum(is.na(rawobj$NL))
# sum(is.na(rawobj$fuwei))
# sum(is.na(rawobj$bmi))
# sum(is.na(rawobj$LDL))
# sum(is.na(rawobj$shuzhangya))
# sum(is.na(rawobj$shousuoya))
#
# sum(is.na(rawobj$left_level))
# sum(is.na(rawobj$right_level))
# sum(is.na(rawobj$smoking_level))
# sum(is.na(rawobj$NL_level))
# sum(is.na(rawobj$fuwei_level))
# sum(is.na(rawobj$bmi_level))
# sum(is.na(rawobj$LDL_level))






