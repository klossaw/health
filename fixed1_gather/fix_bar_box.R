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
rawobj <- read_rds("./fixdata/updata_raw.rds")
core_frame <- read_rds("./fixdata/updata_core.rds")

# ==========================  厚度做分组  ======================================
boxsub1_core_frame <- core_frame[, c("NHDL", "LDL", "NHDL_NL", "LDL_NL", "shousuoya", "shuzhangya", "bmi", "fuwei", "NL",
                                     "left_level", "right_level", "max_level", "total_level")] %>% na.omit()
nrow(boxsub1_core_frame)
# 2391
boxsub1_core_frame$left_level <- factor(boxsub1_core_frame$left_level, levels = c("normal", "level1", "level2"))
boxsub1_core_frame$right_level <- factor(boxsub1_core_frame$right_level, levels = c("normal", "level1", "level2"))
boxsub1_core_frame$max_level <- factor(boxsub1_core_frame$max_level, levels = c("normal", "level1", "level2"))
dir.create("./fix_first/")
ic <- c("NHDL", "LDL", "NHDL_NL", "LDL_NL", "shousuoya", "shuzhangya", "bmi", "fuwei")
pdf("./fix_first/updata_score1_boxplot.pdf")
my_comparisons <- list(c("normal", "level1"),c("level1","level2"),c("normal","level2"))
for(i in ic){
  p1 <- ggboxplot(boxsub1_core_frame, x = "left_level", y = i, fill = "left_level",
                  palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
    stat_compare_means(comparisons = my_comparisons) + NoLegend()
  p2 <- ggboxplot(boxsub1_core_frame, x = "right_level", y = i, fill = "right_level",
                  palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
    stat_compare_means(comparisons = my_comparisons) + NoLegend()
  p3 <- ggboxplot(boxsub1_core_frame, x = "max_level", y = i, fill = "max_level",
                  palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
    stat_compare_means(comparisons = my_comparisons) + NoLegend()
  p4 <- ggboxplot(boxsub1_core_frame, x = "total_level", y = i, fill = "total_level") +
    stat_compare_means(method = "anova") + NoLegend()
  pt <- (p1 + p2)/(p3 + p4)
  print(pt)
}
dev.off()

# ============================  厚度做分组的其他数据  ==========================
vm <- c("脂肪肝", "心率", "身高", "体重",
        "白细胞计数", "中性粒细胞","中性粒细胞百分比","淋巴细胞百分比","单核细胞百分比","血红蛋白","血小板计数","血小板压积",
        "尿微量白蛋白","尿微量白蛋白尿肌酐比值",
        "总蛋白","白蛋白","谷丙转氨酶","谷草转氨酶","肌酐","尿素","肾小球滤过率",
        "同型半胱氨酸","甘油三酯","总胆固醇","高密度脂蛋白C","载脂蛋白A1","载脂蛋白B","尿酸",
        "空腹血糖","糖化血红蛋白A1","糖化血红蛋白A1C","空腹C肽","空腹胰岛素","超敏C反应蛋白",
        "唾液酸","游离脂肪酸","载脂蛋白E","羟基维生素D3")
names(vm) <- c(
  "zhifanggan", "xinlv","shengao","tizhong",
  "baijisuh","zhongjishu","zhongbili","linbabili","danhebili","xuehong","xiaobanjishu","xiaobanjiya",
  "niaoweidanbai", "niaoweidanbaijigan",
  "zongdanbai","baidanbai","gubing","gucao","jigan","niaosu","lvguolv",
  "tongxing","ganyousanzhi","zongdanguchun","gaomidu","zaizhia1","zaizhib","niaosuan",
  "kongfuxuetang","tanghuaA1","tanghuaA1C","kongfuC","kongfuyidaosu","chaomingC",
  "tuoyesuan","youlizhifangsuan","zaizhiE","qiangjiVD3"
)
boxsub1x_core_frame <- rawobj[, c(vm, "left_level", "right_level", "max_level", "total_level")]
nrow(boxsub1x_core_frame)

boxsub1x_core_frame$left_level <- factor(boxsub1x_core_frame$left_level, levels = c("normal", "level1", "level2"))
boxsub1x_core_frame$right_level <- factor(boxsub1x_core_frame$right_level, levels = c("normal", "level1", "level2"))
boxsub1x_core_frame$max_level <- factor(boxsub1x_core_frame$max_level, levels = c("normal", "level1", "level2"))
ic <- vm
pdf("./fix_first/updata_score1x_boxplot.pdf")
my_comparisons <- list(c("normal", "level1"),c("level1","level2"),c("normal","level2"))
for(ix in 1:length(ic)){
  i = ic[[ix]]
  yn <- names(ic)[[ix]]
  sub_frame <- boxsub1x_core_frame[,c("left_level", "right_level","max_level", "total_level", i )] %>% na.omit()
  nsp <- nrow(sub_frame)
  p1 <- ggboxplot(sub_frame, x = "left_level", y = i, fill = "left_level",
                  palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
    stat_compare_means(comparisons = my_comparisons) + NoLegend() + ylab(yn) +
    ggtitle(glue("sample:{nsp}"))
  p2 <- ggboxplot(sub_frame, x = "right_level", y = i, fill = "right_level",
                  palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
    stat_compare_means(comparisons = my_comparisons) + NoLegend() + ylab(yn) +
    ggtitle(glue("sample:{nsp}"))
  p3 <- ggboxplot(sub_frame, x = "max_level", y = i, fill = "max_level",
                  palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
    stat_compare_means(comparisons = my_comparisons) + NoLegend()+ ylab(yn) +
    ggtitle(glue("sample:{nsp}"))
  p4 <- ggboxplot(sub_frame, x = "total_level",
                  y = i, fill = "total_level") +
    stat_compare_means(method = "anova") + NoLegend() + ylab(yn) +
    ggtitle(glue("sample:{nsp}"))
  pt <- (p1 + p2)/(p3 + p4)
  print(pt)
}
dev.off()
# ==========================  协变量做分组  ====================================
boxsub2_core_frame <- core_frame[, c("smoking_level","NL_level","fuwei_level","bmi_level","LDL_level","xueya_level","xingbie_level",
                                     "max_data","right","left","total_level")] %>% na.omit()
# 2401
ic <- c("smoking_level","NL_level","fuwei_level","bmi_level","LDL_level","xueya_level","xingbie_level")
pdf("./fix_first/updata_score2_boxplot.pdf",width = 10, height = 10)
for(i in ic){
  mx <- 3
  if(max(boxsub2_core_frame[[i]]) == 2){
    my_comparisons <- list(c("normal", "level1"),c("level1","level2"),c("normal","level2"))
    boxsub2_core_frame[[i]] <- case_when(
      boxsub2_core_frame[[i]] == 0 ~ "normal",
      boxsub2_core_frame[[i]] == 1 ~ "level1",
      boxsub2_core_frame[[i]] == 2 ~ "level2"
    )
    boxsub2_core_frame[[i]] <- factor(boxsub2_core_frame[[i]], levels = c("normal", "level1", "level2"))
    mx <- 2
  }else if(max(boxsub2_core_frame[[i]]) == 1){
    my_comparisons <- list(c("normal", "level1"))
    boxsub2_core_frame[[i]] <- case_when(
      boxsub2_core_frame[[i]] == 0 ~ "normal",
      boxsub2_core_frame[[i]] == 1 ~ "level1"
    )
    boxsub2_core_frame[[i]] <- factor(boxsub2_core_frame[[i]], levels = c("normal", "level1"))
    mx <- 1
  }
  if(mx <= 2){
    p1 <- ggboxplot(boxsub2_core_frame, x = i, y = "left", fill = i,
                    palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
      stat_compare_means(comparisons = my_comparisons) + NoLegend()
    p2 <- ggboxplot(boxsub2_core_frame, x = i, y = "right", fill = i,
                    palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
      stat_compare_means(comparisons = my_comparisons) + NoLegend()
    p3 <- ggboxplot(boxsub2_core_frame, x = i, y = "max_data", fill = i,
                    palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
      stat_compare_means(comparisons = my_comparisons) + NoLegend()
    p4 <- ggboxplot(boxsub2_core_frame, x = i, y = "total_level", fill = i,
                    palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
      stat_compare_means(comparisons = my_comparisons) + NoLegend()
  }else{
    p1 <- ggboxplot(boxsub2_core_frame, x = i, y = "left", fill = i) +
      stat_compare_means(method = "anova") + NoLegend()
    p2 <- ggboxplot(boxsub2_core_frame, x = i, y = "right", fill = i) +
      stat_compare_means(method = "anova") + NoLegend()
    p3 <- ggboxplot(boxsub2_core_frame, x = i, y = "max_data", fill = i) +
      stat_compare_means(method = "anova") + NoLegend()
    p4 <- ggboxplot(boxsub2_core_frame, x = i, y = "total_level", fill = i) +
      stat_compare_means(method = "anova") + NoLegend()
  }
  pt <- (p1 + p2)/(p3 + p4)
  print(pt)
}
dev.off()
# ============================  堆积柱装图卡方检验  ======================================
boxsub3_core_frame <- core_frame[, c("smoking_level","NL_level","fuwei_level","bmi_level","LDL_level","xueya_level","xingbie_level",
                                     "left_level", "right_level", "max_level", "total_level")] %>%
  na.omit()
v2 <- c("smoking_level","NL_level","fuwei_level","bmi_level","LDL_level","xueya_level","xingbie_level")
myggtitle <- function(x, y){
  fisher.test.out <- table(x, y) %>% fisher.test(simulate.p.value=TRUE)
  ggtitle(glue::glue("fisher.test : {signif(fisher.test.out$p.value, 3)}"))
}
boxsub3_core_frame$right_level <- factor(boxsub3_core_frame$right_level, levels = c("normal", "level1", "level2" ))
boxsub3_core_frame$left_level <- factor(boxsub3_core_frame$left_level, levels = c("normal", "level1", "level2" ))
boxsub3_core_frame$max_level <- factor(boxsub3_core_frame$max_level, levels = c("normal", "level1", "level2" ))
# 2401
pdf("./fix_first/level_barplot.pdf",width = 12)
for(i in v2){
  p1 <- ggplot(boxsub3_core_frame, aes(x = .data[[i]], fill = .data[["right_level"]])) +
    geom_bar(position = "fill") + theme_classic() + myggtitle(boxsub3_core_frame[[i]], boxsub3_core_frame[["right_level"]])
  p2 <- ggplot(boxsub3_core_frame, aes( x = .data[[i]], fill = .data[["left_level"]])) +
    geom_bar(position = "fill") + theme_classic()+ myggtitle(boxsub3_core_frame[[i]], boxsub3_core_frame[["left_level"]])
  p3 <- ggplot(boxsub3_core_frame, aes(x = .data[[i]], fill = .data[["max_level"]])) +
    geom_bar(position = "fill") + theme_classic()+ myggtitle(boxsub3_core_frame[[i]], boxsub3_core_frame[["max_level"]])
  p4 <- ggplot(boxsub3_core_frame, aes(x = .data[["total_level"]], fill = .data[[i]] %>% as.character())) +
    geom_bar(position = "fill") + theme_classic() +
    labs(fill = i)+ myggtitle(boxsub3_core_frame[[i]], boxsub3_core_frame[["total_level"]])
  print((p1 + p2)/(p3 + p4))
}
dev.off()

pdf("./fix_first/level_barplot_count.pdf",width = 12)
for(i in v2){
  p1 <- ggplot(boxsub3_core_frame, aes(x = .data[[i]], fill = .data[["right_level"]])) +
    geom_bar() + theme_classic() + myggtitle(boxsub3_core_frame[[i]], boxsub3_core_frame[["right_level"]])
  p2 <- ggplot(boxsub3_core_frame, aes( x = .data[[i]], fill = .data[["left_level"]])) +
    geom_bar() + theme_classic()+ myggtitle(boxsub3_core_frame[[i]], boxsub3_core_frame[["left_level"]])
  p3 <- ggplot(boxsub3_core_frame, aes(x = .data[[i]], fill = .data[["max_level"]])) +
    geom_bar() + theme_classic()+ myggtitle(boxsub3_core_frame[[i]], boxsub3_core_frame[["max_level"]])
  p4 <- ggplot(boxsub3_core_frame, aes(x = .data[["total_level"]], fill = .data[[i]] %>% as.character())) +
    geom_bar() + theme_classic() +
    labs(fill = i)+ myggtitle(boxsub3_core_frame[[i]], boxsub3_core_frame[["total_level"]])
  print((p1 + p2)/(p3 + p4))
}
dev.off()

# ================================  脂肪肝 =====================================
boxsub3x_core_frame <- rawobj[, c("脂肪肝", "left_level", "right_level", "max_level", "total_level")] %>%
  na.omit()
v2 <- c("脂肪肝")
myggtitle <- function(x, y){
  fisher.test.out <- table(x, y) %>% fisher.test(simulate.p.value=TRUE)
  ggtitle(glue::glue("fisher.test : {signif(fisher.test.out$p.value, 3)}"))
}
boxsub3x_core_frame$right_level <- factor(boxsub3x_core_frame$right_level, levels = c("normal", "level1", "level2" ))
boxsub3x_core_frame$left_level <- factor(boxsub3x_core_frame$left_level, levels = c("normal", "level1", "level2" ))
boxsub3x_core_frame$max_level <- factor(boxsub3x_core_frame$max_level, levels = c("normal", "level1", "level2" ))
# 2382
pdf("./fix_first/zhifanggan_level_barplot.pdf",width = 12)
i = v2
p1 <- ggplot(boxsub3x_core_frame, aes(x = .data[[i]], fill = .data[["right_level"]])) +
  geom_bar(position = "fill") + theme_classic() + myggtitle(boxsub3x_core_frame[[i]], boxsub3x_core_frame[["right_level"]]) + 
  xlab("zhifanggan")
p2 <- ggplot(boxsub3x_core_frame, aes( x = .data[[i]], fill = .data[["left_level"]])) +
  geom_bar(position = "fill") + theme_classic()+ myggtitle(boxsub3x_core_frame[[i]], boxsub3x_core_frame[["left_level"]]) + 
  xlab("zhifanggan")
p3 <- ggplot(boxsub3x_core_frame, aes(x = .data[[i]], fill = .data[["max_level"]])) +
  geom_bar(position = "fill") + theme_classic()+ myggtitle(boxsub3x_core_frame[[i]], boxsub3x_core_frame[["max_level"]]) + 
  xlab("zhifanggan")
p4 <- ggplot(boxsub3x_core_frame, aes(x = .data[["total_level"]], fill = .data[[i]] %>% as.character())) +
  geom_bar(position = "fill") + theme_classic() +
  labs(fill = i)+ myggtitle(boxsub3x_core_frame[[i]], boxsub3x_core_frame[["total_level"]]) 
print((p1 + p2)/(p3 + p4))
dev.off()

pdf("./fix_first/zhifanggan_level_barplot_count.pdf",width = 12)
i = v2
p1 <- ggplot(boxsub3x_core_frame, aes(x = .data[[i]], fill = .data[["right_level"]])) +
  geom_bar() + theme_classic() + myggtitle(boxsub3x_core_frame[[i]], boxsub3x_core_frame[["right_level"]]) + 
  xlab("zhifanggan")
p2 <- ggplot(boxsub3x_core_frame, aes( x = .data[[i]], fill = .data[["left_level"]])) +
  geom_bar() + theme_classic()+ myggtitle(boxsub3x_core_frame[[i]], boxsub3x_core_frame[["left_level"]]) + 
  xlab("zhifanggan")
p3 <- ggplot(boxsub3x_core_frame, aes(x = .data[[i]], fill = .data[["max_level"]])) +
  geom_bar() + theme_classic()+ myggtitle(boxsub3x_core_frame[[i]], boxsub3x_core_frame[["max_level"]]) + 
  xlab("zhifanggan")
p4 <- ggplot(boxsub3x_core_frame, aes(x = .data[["total_level"]], fill = .data[[i]] %>% as.character())) +
  geom_bar() + theme_classic() +
  labs(fill = i)+ myggtitle(boxsub3x_core_frame[[i]], boxsub3x_core_frame[["total_level"]]) 
print((p1 + p2)/(p3 + p4))
dev.off()
# ================================  厚度为y的核心点图  =========================
point_core_frame <- core_frame[, c("NHDL", "LDL", "NHDL_NL", "LDL_NL", "shousuoya", "shuzhangya", "bmi", "fuwei", "NL",
                                   "max_data","right","left","total_level")] %>% na.omit()
v2 <- c("NHDL", "LDL", "NHDL_NL", "LDL_NL", "shousuoya", "shuzhangya", "bmi", "fuwei", "NL")
dim(point_core_frame)
# 2391
pdf("./fix_first/core_pointplot.pdf")
for(i in v2){
  p1 <- ggscatter(point_core_frame, x = i, y = "right", add = "reg.line") +
    stat_cor(label.y = 7, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
    stat_regline_equation(label.y = 6)
  p2 <- ggscatter(point_core_frame, x = i, y = "left", add = "reg.line") +
    stat_cor(label.y = 7, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
    stat_regline_equation(label.y = 6)
  p3 <- ggscatter(point_core_frame, x = i, y = "max_data", add = "reg.line") +
    stat_cor(label.y = 7, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
    stat_regline_equation(label.y = 6)
  p4 <- ggscatter(point_core_frame, x = i, y = "total_level", add = "reg.line") +
    stat_cor(label.y = 18, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
    stat_regline_equation(label.y = 17)
  print((p1 + p2)/(p3 + p4))
}
dev.off()
# =============================  厚度为y的其他数据的点图  ======================
vm <- c("脂肪肝", "心率", "身高", "体重",
        "白细胞计数", "中性粒细胞","中性粒细胞百分比","淋巴细胞百分比","单核细胞百分比","血红蛋白","血小板计数","血小板压积",
        "尿微量白蛋白","尿微量白蛋白尿肌酐比值",
        "总蛋白","白蛋白","谷丙转氨酶","谷草转氨酶","肌酐","尿素","肾小球滤过率",
        "同型半胱氨酸","甘油三酯","总胆固醇","高密度脂蛋白C","载脂蛋白A1","载脂蛋白B","尿酸",
        "空腹血糖","糖化血红蛋白A1","糖化血红蛋白A1C","空腹C肽","空腹胰岛素","超敏C反应蛋白",
        "唾液酸","游离脂肪酸","载脂蛋白E","羟基维生素D3")
names(vm) <- c(
  "zhifanggan", "xinlv","shengao","tizhong",
  "baijisuh","zhongjishu","zhongbili","linbabili","danhebili","xuehong","xiaobanjishu","xiaobanjiya",
  "niaoweidanbai", "niaoweidanbaijigan",
  "zongdanbai","baidanbai","gubing","gucao","jigan","niaosu","lvguolv",
  "tongxing","ganyousanzhi","zongdanguchun","gaomidu","zaizhia1","zaizhib","niaosuan",
  "kongfuxuetang","tanghuaA1","tanghuaA1C","kongfuC","kongfuyidaosu","chaomingC",
  "tuoyesuan","youlizhifangsuan","zaizhiE","qiangjiVD3"
)
point_other_frame <- rawobj[, c(vm, "max_data","right","left","total_level")]
ic <- vm
pdf("./fix_first/other_pointplot.pdf")
for(ix in 1:length(ic)){
  i = ic[[ix]]
  yn <- names(ic)[[ix]]
  sub_frame <- point_other_frame[,c("left", "right","max_data", "total_level", i )] %>% na.omit()
  nsp <- nrow(sub_frame)
  p1 <- ggscatter(sub_frame, x = i, y = "right", add = "reg.line") +
    stat_cor(label.y = 7, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
    stat_regline_equation(label.y = 6) + xlab(yn) +
    ggtitle(glue("sample:{nsp}"))
  p2 <- ggscatter(sub_frame, x = i, y = "left", add = "reg.line") +
    stat_cor(label.y = 7, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
    stat_regline_equation(label.y = 6) + xlab(yn) +
    ggtitle(glue("sample:{nsp}"))
  p3 <- ggscatter(sub_frame, x = i, y = "max_data", add = "reg.line") +
    stat_cor(label.y = 7, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
    stat_regline_equation(label.y = 6) + xlab(yn) +
    ggtitle(glue("sample:{nsp}"))
  p4 <- ggscatter(sub_frame, x = i, y = "total_level", add = "reg.line") +
    stat_cor(label.y = 18, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
    stat_regline_equation(label.y = 17) + xlab(yn)
  print((p1 + p2)/(p3 + p4))
}
dev.off()







