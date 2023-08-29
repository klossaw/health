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
rawobj <- read_rds("./combine/updata_raw.rds")
core_frame <- read_rds("./combine/updata_core.rds")
cross_chisquare <- function(year, gender, x, y){
  black <- tibble(
    year = year,
    gender = gender,
    x = x,
    y = y
  ) %>% group_by(year, gender, x, y) %>% 
    summarise(count = n()) %>% 
    ungroup() %>% group_by(year, gender) %>% 
    summarise(lx = length(unique(x)), ly = length(unique(y))) %>% ungroup() %>% 
    dplyr::filter(lx == 1 | ly == 1)
  
  tb <- tibble(
    year = year,
    gender = gender,
    x = x,
    y = y
  ) %>% anti_join(black) %>% group_by(year, gender) %>% 
    summarise(number = n(), fisher.test = list(table(pick(x, y)) %>% fisher.test(simulate.p.value=TRUE))) %>% 
    rowwise() %>% mutate(p.value = fisher.test$p.value) %>% dplyr::select(-"fisher.test") %>% ungroup()
  tb
}
cal_fsheet <- function(df, v){
  fftb <- lapply(v, function(i){
    lapply(c("right_level", "left_level", "max_level"), function(side){
      ftb <- cross_chisquare(df$NL_level, df$xingbie_level, df[[side]], df[[i]])
      colnames(ftb)[[4]] <- side
      ftb
    }) %>% purrr::reduce(., left_join)
  }) 
  names(fftb) <- v
  fftb <- fftb %>% list_rbind(names_to = "feature_level")
  fftb
}
# ==================================== 核心指标 按照指标分柱子 =================
boxsub3_core_frame <- core_frame[, c("smoking_level","NL_level","fuwei_level","bmi_level","LDL_level","xueya_level","xingbie_level",
                                     "left_level", "right_level", "max_level", "total_level")] %>%
  na.omit()
v2 <- c("smoking_level","NL_level","fuwei_level","bmi_level","LDL_level","xueya_level","xingbie_level")
boxsub3_core_frame$right_level <- factor(boxsub3_core_frame$right_level, levels = c("normal", "level1", "level2" ))
boxsub3_core_frame$left_level <- factor(boxsub3_core_frame$left_level, levels = c("normal", "level1", "level2" ))
boxsub3_core_frame$max_level <- factor(boxsub3_core_frame$max_level, levels = c("normal", "level1", "level2" ))
# 2401
cal_fsheet(boxsub3_core_frame, base::setdiff(v2, c("xingbie_level","NL_level"))) %>% write_csv("./combine/face_level_barplot.csv")
pdf("./combine/face_level_barplot.pdf",width = 12)
for(i in v2){
  p1 <- ggplot(boxsub3_core_frame, aes(x = .data[[i]], fill = .data[["right_level"]])) +
    geom_bar(position = "fill") + theme_classic() + 
    facet_grid(xingbie_level ~  NL_level)
    
  p2 <- ggplot(boxsub3_core_frame, aes( x = .data[[i]], fill = .data[["left_level"]])) +
    geom_bar(position = "fill") + theme_classic() + 
    facet_grid(xingbie_level ~  NL_level)
  
  p3 <- ggplot(boxsub3_core_frame, aes(x = .data[[i]], fill = .data[["max_level"]])) +
    geom_bar(position = "fill") + theme_classic() + 
    facet_grid(xingbie_level ~  NL_level)
  print(p1 / p2 / p3)
}
dev.off()

pdf("./combine/face_count_level_barplot.pdf",width = 12)
for(i in v2){
  p1 <- ggplot(boxsub3_core_frame, aes(x = .data[[i]], fill = .data[["right_level"]])) +
    geom_bar() + theme_classic() + 
    facet_grid(xingbie_level ~  NL_level)
  p2 <- ggplot(boxsub3_core_frame, aes( x = .data[[i]], fill = .data[["left_level"]])) +
    geom_bar() + theme_classic() + 
    facet_grid(xingbie_level ~  NL_level)
  p3 <- ggplot(boxsub3_core_frame, aes(x = .data[[i]], fill = .data[["max_level"]])) +
    geom_bar() + theme_classic() + 
    facet_grid(xingbie_level ~  NL_level)
  print(p1 / p2 / p3)
}
dev.off()
# ================================  脂肪肝 =====================================
boxsub3x_core_frame <- rawobj[, c("NL_level", "xingbie_level", "脂肪肝", "left_level", "right_level", "max_level", "total_level")] %>%
  na.omit()
v2 <- c("脂肪肝")
myggtitle <- function(x, y){
  fisher.test.out <- table(x, y) %>% fisher.test(simulate.p.value=TRUE)
  ggtitle(glue::glue("fisher.test : {signif(fisher.test.out$p.value, 3)}"))
}
boxsub3x_core_frame$right_level <- factor(boxsub3x_core_frame$right_level, levels = c("normal", "level1", "level2" ))
boxsub3x_core_frame$left_level <- factor(boxsub3x_core_frame$left_level, levels = c("normal", "level1", "level2" ))
boxsub3x_core_frame$max_level <- factor(boxsub3x_core_frame$max_level, levels = c("normal", "level1", "level2" ))
cal_fsheet(boxsub3x_core_frame, v2) %>% write_csv("./combine/face_zhifanggan_level_barplot.csv")

# 2382
pdf("./combine/face_zhifanggan_level_barplot.pdf",width = 12)
i = v2
p1 <- ggplot(boxsub3x_core_frame, aes(x = .data[[i]], fill = .data[["right_level"]])) +
  geom_bar(position = "fill") + theme_classic() + 
  xlab("zhifanggan") +
  facet_grid(xingbie_level ~  NL_level)
p2 <- ggplot(boxsub3x_core_frame, aes( x = .data[[i]], fill = .data[["left_level"]])) +
  geom_bar(position = "fill") + theme_classic()+ 
  xlab("zhifanggan") +
  facet_grid(xingbie_level ~  NL_level)
p3 <- ggplot(boxsub3x_core_frame, aes(x = .data[[i]], fill = .data[["max_level"]])) +
  geom_bar(position = "fill") + theme_classic()+ 
  xlab("zhifanggan") +
  facet_grid(xingbie_level ~  NL_level)
print(p1 / p2 / p3)
dev.off()

pdf("./combine/face_zhifanggan_level_barplot_count.pdf",width = 12)
i = v2
p1 <- ggplot(boxsub3x_core_frame, aes(x = .data[[i]], fill = .data[["right_level"]])) +
  geom_bar() + theme_classic() + 
  xlab("zhifanggan") +
  facet_grid(xingbie_level ~  NL_level)
p2 <- ggplot(boxsub3x_core_frame, aes( x = .data[[i]], fill = .data[["left_level"]])) +
  geom_bar() + theme_classic()+ 
  xlab("zhifanggan") +
  facet_grid(xingbie_level ~  NL_level)
p3 <- ggplot(boxsub3x_core_frame, aes(x = .data[[i]], fill = .data[["max_level"]])) +
  geom_bar() + theme_classic()+
  xlab("zhifanggan") +
  facet_grid(xingbie_level ~  NL_level)
print(p1 / p2 / p3)
dev.off()

# ==================================== 核心指标 按照厚度分柱子 =================
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

cal_fsheet(boxsub3_core_frame, v2) %>% write_csv("./combine/face_thin_level_barplot.csv")

pdf("./combine/face_thin_level_barplot.pdf",width = 12)
for(i in v2){
  p1 <- ggplot(boxsub3_core_frame, aes(x = .data[["right_level"]], fill = .data[[i]] %>% as.character())) +
    geom_bar(position = "fill") + theme_classic() + 
    facet_grid(xingbie_level ~  NL_level)+ 
    scale_fill_viridis_d()+ 
    guides(fill = guide_legend(title=i))
  
  p2 <- ggplot(boxsub3_core_frame, aes( x = .data[["left_level"]], fill = .data[[i]] %>% as.character())) +
    geom_bar(position = "fill") + theme_classic() + 
    facet_grid(xingbie_level ~  NL_level)+ 
    scale_fill_viridis_d()+ 
    guides(fill = guide_legend(title=i))
  
  p3 <- ggplot(boxsub3_core_frame, aes(x = .data[["max_level"]], fill = .data[[i]] %>% as.character())) +
    geom_bar(position = "fill") + theme_classic() + 
    facet_grid(xingbie_level ~  NL_level)+ 
    scale_fill_viridis_d()+ 
    guides(fill = guide_legend(title=i))
  print(p1 / p2 / p3)
}
dev.off()

pdf("./combine/face_count_thin_level_barplot.pdf",width = 12)
for(i in v2){
  p1 <- ggplot(boxsub3_core_frame, aes(x = .data[["right_level"]], fill = .data[[i]] %>% as.character())) +
    geom_bar() + theme_classic() + 
    facet_grid(xingbie_level ~  NL_level)+ 
    scale_fill_viridis_d()+ 
    guides(fill = guide_legend(title=i))
  p2 <- ggplot(boxsub3_core_frame, aes( x = .data[["left_level"]], fill = .data[[i]] %>% as.character())) +
    geom_bar() + theme_classic() + 
    facet_grid(xingbie_level ~  NL_level)+ 
    scale_fill_viridis_d()+ 
    guides(fill = guide_legend(title=i))
  p3 <- ggplot(boxsub3_core_frame, aes(x = .data[["max_level"]], fill = .data[[i]] %>% as.character())) +
    geom_bar() + theme_classic() + 
    facet_grid(xingbie_level ~  NL_level)+ 
    scale_fill_viridis_d() + 
    guides(fill = guide_legend(title=i))
  print(p1 / p2 / p3)
}
dev.off()
# ================================  厚度脂肪肝 =====================================
boxsub3x_core_frame <- rawobj[, c("NL_level", "xingbie_level", "脂肪肝", "left_level", "right_level", "max_level", "total_level")] %>%
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

cal_fsheet(boxsub3x_core_frame, v2) %>% write_csv("./combine/face_zhifanggan_thin_level_barplot.csv")

pdf("./combine/face_zhifanggan_thin_level_barplot.pdf",width = 12)
i = v2
p1 <- ggplot(boxsub3x_core_frame, aes(x = .data[["right_level"]], fill = .data[[i]] %>% as.character())) +
  geom_bar(position = "fill") + theme_classic() + 
  facet_grid(xingbie_level ~  NL_level)+ 
  scale_fill_viridis_d() + 
  guides(fill = guide_legend(title="zhifanggan"))
p2 <- ggplot(boxsub3x_core_frame, aes( x = .data[["left_level"]], fill = .data[[i]] %>% as.character())) +
  geom_bar(position = "fill") + theme_classic()+ 
  facet_grid(xingbie_level ~  NL_level)+ 
  scale_fill_viridis_d() + 
  guides(fill = guide_legend(title="zhifanggan"))
p3 <- ggplot(boxsub3x_core_frame, aes(x = .data[["max_level"]], fill = .data[[i]] %>% as.character())) +
  geom_bar(position = "fill") + theme_classic()+
  facet_grid(xingbie_level ~  NL_level)+ 
  scale_fill_viridis_d() + 
  guides(fill = guide_legend(title="zhifanggan"))
print(p1 / p2 / p3)
dev.off()

pdf("./combine/face_zhifanggan_thin_level_barplot_count.pdf",width = 12)
i = v2
p1 <- ggplot(boxsub3x_core_frame, aes(x = .data[["right_level"]], fill = .data[[i]] %>% as.character())) +
  geom_bar() + theme_classic() + 
  facet_grid(xingbie_level ~  NL_level) + 
  scale_fill_viridis_d() + 
  guides(fill = guide_legend(title="zhifanggan"))
p2 <- ggplot(boxsub3x_core_frame, aes( x = .data[["left_level"]], fill = .data[[i]] %>% as.character())) +
  geom_bar() + theme_classic()+ 
  facet_grid(xingbie_level ~  NL_level)+ 
  scale_fill_viridis_d()+ 
  guides(fill = guide_legend(title="zhifanggan"))
p3 <- ggplot(boxsub3x_core_frame, aes(x = .data[["max_level"]], fill = .data[[i]] %>% as.character())) +
  geom_bar() + theme_classic()+ 
  facet_grid(xingbie_level ~  NL_level) + 
  scale_fill_viridis_d()+ 
  guides(fill = guide_legend(title="zhifanggan"))
print(p1 / p2 / p3)
dev.off()

# =============================== 核心指标 厚度分组 盒子 =======================
boxsub1_core_frame <- core_frame[, c("NL_level", "xingbie_level",
                                     "NHDL", "LDL", "NHDL_NL", "LDL_NL", 
                                     "shousuoya", "shuzhangya", "bmi", "fuwei", "NL",
                                     "left_level", "right_level", "max_level", "total_level")] %>% na.omit()
nrow(boxsub1_core_frame)
# 2391
boxsub1_core_frame$left_level <- factor(boxsub1_core_frame$left_level, levels = c("normal", "level1", "level2"))
boxsub1_core_frame$right_level <- factor(boxsub1_core_frame$right_level, levels = c("normal", "level1", "level2"))
boxsub1_core_frame$max_level <- factor(boxsub1_core_frame$max_level, levels = c("normal", "level1", "level2"))
dir.create("./combine/")
ic <- c("NHDL", "LDL", "NHDL_NL", "LDL_NL", "shousuoya", "shuzhangya", "bmi", "fuwei")
my_comparisons <- list(c("normal", "level1"),c("level1","level2"),c("normal","level2"))
pdf("./combine/face_boxplot.pdf", height = 12, width = 16)
for(i in ic){
  p1 <- ggboxplot(boxsub1_core_frame, x = "left_level", y = i, fill = "left_level",
                  palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
    stat_compare_means(comparisons = my_comparisons) + NoLegend() +
    facet_grid(xingbie_level ~ NL_level)
  p2 <- ggboxplot(boxsub1_core_frame, x = "right_level", y = i, fill = "right_level",
                  palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
    stat_compare_means(comparisons = my_comparisons) + NoLegend() +
    facet_grid(xingbie_level ~ NL_level)
  p3 <- ggboxplot(boxsub1_core_frame, x = "max_level", y = i, fill = "max_level",
                  palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
    stat_compare_means(comparisons = my_comparisons) + NoLegend() +
    facet_grid(xingbie_level ~ NL_level)
  pt <- p1 / p2 / p3
  print(pt)
}
dev.off()

# ==================================== 其他指标 厚度分组盒子 ===================
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
boxsub1x_core_frame <- rawobj[, c(vm,"NL_level", "xingbie_level", "left_level", "right_level", "max_level", "total_level")]
nrow(boxsub1x_core_frame)

boxsub1x_core_frame$left_level <- factor(boxsub1x_core_frame$left_level, levels = c("normal", "level1", "level2"))
boxsub1x_core_frame$right_level <- factor(boxsub1x_core_frame$right_level, levels = c("normal", "level1", "level2"))
boxsub1x_core_frame$max_level <- factor(boxsub1x_core_frame$max_level, levels = c("normal", "level1", "level2"))
ic <- vm
my_comparisons <- list(c("normal", "level1"),c("level1","level2"),c("normal","level2"))

pdf("./combine/face_boxplotx.pdf", height = 12, width = 16)
for(i in 1:length(ic)){
  yn <- names(ic)[[i]]
  i = ic[[i]]
  sub_frame <- boxsub1x_core_frame[,c("NL_level", "xingbie_level", "left_level", "right_level","max_level", "total_level", i )] %>% na.omit()
  nsp <- nrow(sub_frame)
  p1 <- ggboxplot(sub_frame, x = "left_level", y = i, fill = "left_level",
                  palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
    stat_compare_means(comparisons = my_comparisons) + NoLegend() +
    facet_grid(xingbie_level ~ NL_level)  + ylab(yn) +
    ggtitle(glue("sample:{nsp}"))
  p2 <- ggboxplot(sub_frame, x = "right_level", y = i, fill = "right_level",
                  palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
    stat_compare_means(comparisons = my_comparisons) + NoLegend() +
    facet_grid(xingbie_level ~ NL_level) + ylab(yn) +
    ggtitle(glue("sample:{nsp}"))
  p3 <- ggboxplot(sub_frame, x = "max_level", y = i, fill = "max_level",
                  palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
    stat_compare_means(comparisons = my_comparisons) + NoLegend() +
    facet_grid(xingbie_level ~ NL_level) + ylab(yn) +
    ggtitle(glue("sample:{nsp}"))
  pt <- p1 / p2 / p3
  print(pt)
}
dev.off()











