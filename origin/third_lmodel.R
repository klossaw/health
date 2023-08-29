pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", "factoextra", "cluster", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "Seurat", "ComplexHeatmap",
          "ggpointdensity", "ggpubr")
for (pkg in pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "healthman"
dataset <- "chenlei"
species <- "human"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/") %>% checkdir()
setwd(workdir)
# ================================ 第一部分 点散图数值回归 ====================================
rawobj <- "/cluster/home/jhuang/projects/healthman/data/chenlei/human/clinical/raw.rds" %>% read_rds
rawobj <- rawobj[rawobj$"是否采纳" == "采纳",]
rawobj$"非高密度脂蛋白C" <- as.numeric(rawobj$"总胆固醇") - as.numeric(rawobj$"高密度脂蛋白C")

leftight_level <- tibble(
  right = rawobj$"右侧颈动脉" %>% as.numeric(),
  left = rawobj$"左侧颈动脉" %>% as.numeric(),
  NL =  rawobj$NL %>% str_remove("岁") %>% as.numeric(),
  LDL =  rawobj$"低密度脂蛋白C" %>% as.numeric(),
  smoking =  rawobj$"吸烟",
  shousuoya =  rawobj$"收缩压" %>% as.numeric(),
  shuzhangya =  rawobj$"舒张压" %>% as.numeric(),
  bmi =  rawobj$"体重指数" %>% as.numeric(),
  fuwei =  rawobj$"腹围" %>% as.numeric(),
  xingbie = rawobj$XB,
  NHDL = rawobj$"非高密度脂蛋白C",
  NHDL_NL = NHDL * NL,
  LDL_NL = LDL * NL
) %>% na.omit()

# 2390 
leftight_level <- leftight_level %>% 
  dplyr::mutate(
    smoking_level = case_when(
      smoking == "不吸" ~ 0,
      T ~ 1),
    NL_level = floor((NL - 20)/10) + 1,
    fuwei_level = case_when(
      (fuwei < 90 & xingbie == "男") | (fuwei < 85 & xingbie == "女") ~ 0,
      T ~ 1
    ),
    bmi_level = case_when(
      bmi < 24 ~ 0,
      bmi <= 28 ~ 1,
      T ~ 2
    ),
    LDL_level = case_when(
      LDL < 1.8 ~ 0,
      LDL <= 2.6 ~ 1,
      LDL <= 3.4 ~ 2,
      LDL <= 4.9 ~ 3,
      T ~ 4
    ),
    xueya_level = case_when(
      (shousuoya < 130 & shuzhangya < 80) ~ 0,
      (shousuoya < 140 & shuzhangya < 90) ~ 1,
      (shousuoya < 160 & shuzhangya < 100) ~ 2,
      T ~ 3
    )
)
leftight_level$max_data <- apply(leftight_level, 1, function(df){
  max(df[[1]] %>% as.numeric(), df[[2]] %>% as.numeric())
}) 


v <- c("NL","LDL","shousuoya","shuzhangya","bmi","fuwei","NHDL","NHDL_NL","LDL_NL")
pdf("data_lm.pdf")
for(i in v){
  print(ggscatter(leftight_level, x = i, y = "right", add = "reg.line") +
    stat_cor(label.y = 7, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
    stat_regline_equation(label.y = 6))
  print(ggscatter(leftight_level, x = i, y = "left", add = "reg.line") +
    stat_cor(label.y = 7, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
    stat_regline_equation(label.y = 6))
  print(ggscatter(leftight_level, x = i, y = "max_data", add = "reg.line") +
    stat_cor(label.y = 7, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
    stat_regline_equation(label.y = 6))
}
dev.off()

pdf("high_NL_data_lm.pdf")
for(i in v){
  print(ggscatter(leftight_level %>% dplyr::filter(NL >= 40), x = i, y = "right", add = "reg.line", color = "NL") +
          stat_cor(label.y = 7, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
          stat_regline_equation(label.y = 6))
  print(ggscatter(leftight_level %>% dplyr::filter(NL >= 40), x = i, y = "left", add = "reg.line", color = "NL") +
          stat_cor(label.y = 7, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
          stat_regline_equation(label.y = 6))
  print(ggscatter(leftight_level %>% dplyr::filter(NL >= 40), x = i, y = "max_data", add = "reg.line", color = "NL") +
          stat_cor(label.y = 7, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
          stat_regline_equation(label.y = 6))
}
dev.off()

pdf("low_NL_data_lm.pdf")
for(i in v){
  print(ggscatter(leftight_level %>% dplyr::filter(NL < 40), x = i, y = "right", add = "reg.line", color = "NL") +
          stat_cor(label.y = 7, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
          stat_regline_equation(label.y = 6))
  print(ggscatter(leftight_level %>% dplyr::filter(NL < 40), x = i, y = "left", add = "reg.line", color = "NL") +
          stat_cor(label.y = 7, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
          stat_regline_equation(label.y = 6))
  print(ggscatter(leftight_level %>% dplyr::filter(NL < 40), x = i, y = "max_data", add = "reg.line", color = "NL") +
          stat_cor(label.y = 7, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
          stat_regline_equation(label.y = 6))
}
dev.off()


pdf("mid_NL_data_lm.pdf")
for(i in v){
  print(ggscatter(leftight_level %>% dplyr::filter(NL <= 50 & NL > 40), x = i, y = "right", add = "reg.line", color = "NL") +
          stat_cor(label.y = 7, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
          stat_regline_equation(label.y = 6))
  print(ggscatter(leftight_level %>% dplyr::filter(NL <= 50 & NL > 40), x = i, y = "left", add = "reg.line", color = "NL") +
          stat_cor(label.y = 7, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
          stat_regline_equation(label.y = 6))
  print(ggscatter(leftight_level %>% dplyr::filter(NL <= 50 & NL > 40), x = i, y = "max_data", add = "reg.line", color = "NL") +
          stat_cor(label.y = 7, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
          stat_regline_equation(label.y = 6))
}
dev.off()
# ================================ 第二部分 分等级堆积图 ====================================
leftight_level <- leftight_level %>% dplyr::mutate(
  right_level = case_when(
    right < 1 ~ 0,
    right < 1.5 ~ 1,
    T ~ 2
  ),
  left_level = case_when(
    left < 1 ~ 0,
    left < 1.5 ~ 1,
    T ~ 2
  ),
  max_data_level = case_when(
    max_data < 1 ~ 0,
    max_data < 1.5 ~ 1,
    T ~ 2
  ),
  xingbie_level = case_when(
    xingbie == "女" ~ 0,
    xingbie == "男" ~ 1
  )
)

v2 <- c("xingbie_level", "smoking_level", "NL_level","fuwei_level",
        "bmi_level","LDL_level","xueya_level")
pdf("level_data_jitter.pdf",width = 12)
for(i in v2){
  print(
    ggplot(leftight_level, aes(x = .data[[i]], fill = .data[["right_level"]] %>% as.character())) + 
      geom_bar(position = "fill") + theme_classic())
  print(
    ggplot(leftight_level, aes( x = .data[[i]], fill = .data[["left_level"]] %>% as.character())) + 
      geom_bar(position = "fill") + theme_classic())
  print(
    ggplot(leftight_level, aes(x = .data[[i]], fill = .data[["max_data_level"]] %>% as.character())) + 
      geom_bar(position = "fill") + theme_classic())
}
dev.off()
# ================================ 第三部分 总分的点图堆积图 ====================================
leftight_level$total_level <- rowSums(leftight_level[,c("smoking_level", "NL_level", 
                                                        "fuwei_level", "bmi_level", 
                                                        "LDL_level" , "xueya_level" )])

pdf("total_point_bar.pdf",width = 12)
print(
  ggplot(leftight_level, aes(x = .data[["total_level"]], fill = .data[["right_level"]] %>% as.character())) + 
    geom_bar(position = "fill") + theme_classic())
print(
  ggplot(leftight_level, aes(x = .data[["total_level"]], fill = .data[["left_level"]] %>% as.character())) + 
    geom_bar(position = "fill") + theme_classic())
print(
  ggplot(leftight_level, aes(x = .data[["total_level"]], fill = .data[["max_data_level"]] %>% as.character())) + 
    geom_bar(position = "fill") + theme_classic())
print(ggscatter(leftight_level, x = "total_level", y = "right", add = "reg.line", color = "NL") +
        stat_cor(label.y = 7, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
        stat_regline_equation(label.y = 8))
print(ggscatter(leftight_level, x = "total_level", y = "left", add = "reg.line", color = "NL") +
        stat_cor(label.y = 7, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
        stat_regline_equation(label.y = 8))
print(ggscatter(leftight_level, x = "total_level", y = "max_data", add = "reg.line", color = "NL") +
        stat_cor(label.y = 7, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
        stat_regline_equation(label.y = 8))
dev.off()

# ================================ 第四部分 随机森林 ====================================
leftight_levelx <- leftight_level
leftight_levelx$right_level <- factor(leftight_levelx$right_level, levels = 0:2)
leftight_levelx$left_level <- factor(leftight_levelx$left_level, levels = 0:2)
leftight_levelx$max_data_level <- factor(leftight_levelx$max_data_level, levels = 0:2)

train <- sample(1:nrow(leftight_levelx),floor(nrow(leftight_levelx) * 0.7))
data_train <- leftight_levelx[train, ]
data_test <- leftight_levelx[-train, ]

library(randomForest)
total_level_train.forest <- randomForest(total_level ~ ., data = data_train[,c(14:19, 24, 25)], importance = TRUE, ntree=1000)
right_level.forest <- randomForest(right_level ~ ., data = data_train[,c(14:19, 24, 21)], importance = TRUE, ntree=1000)
left_level.forest <- randomForest(left_level ~ ., data = data_train[,c(14:19, 24, 22)], importance = TRUE, ntree=1000)
max_data_level.forest <- randomForest(max_data_level ~ ., data = data_train[,c(14:19, 24, 23)], importance = TRUE, ntree=1000)



doxezt <- function(data_train, data_test, model, realtrain, realtest, fn){
  train_total_level_predict <- predict(model, data_train)
  test_total_level_predict <- predict(model, data_test)
  dftrain <- data.frame(
    predict = train_total_level_predict,
    realdata = realtrain
  )
  dftest <- data.frame(
    predict = test_total_level_predict,
    realdata = realtest
  )
  pdf(fn)
  print(
    ggplot(dftrain, aes(x = predict, y = realdata)) + geom_bin_2d() + 
      geom_abline(mapping = aes(intercept = 0, slope = 1)) + 
      ggtitle("train") + theme_classic() + scale_fill_viridis_c()
    )
  print(
    ggplot(dftest, aes(x = predict, y = realdata)) + geom_bin_2d() + 
      geom_abline(mapping = aes(intercept = 0, slope = 1)) + 
      ggtitle("test") + theme_classic() + scale_fill_viridis_c()
  )
  dev.off()
}

doxezt(data_train, data_test, model = total_level_train.forest, realtrain = data_train$total_level, realtest = data_test$total_level, fn = "total_level_train.forest.pdf")
doxezt(data_train, data_test, model = right_level.forest, realtrain = data_train$right_level, realtest = data_test$right_level, fn = "right_level.forest.pdf")
doxezt(data_train, data_test, model = left_level.forest, realtrain = data_train$left_level, realtest = data_test$left_level, fn = "left_level.forest.pdf")
doxezt(data_train, data_test, model = max_data_level.forest, realtrain = data_train$max_data_level, realtest = data_test$max_data_level, fn = "max_data_level.forest.forest.pdf")


importance(total_level_train.forest)
importance(right_level.forest)
importance(left_level.forest)
importance(max_data_level.forest)
