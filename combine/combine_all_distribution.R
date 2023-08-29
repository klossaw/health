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

# ==================================== 核心指标 按照指标分柱子 =================
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
pdf("./combine/core_distribution_fill_barplot.pdf",width = 7, height = 8)
for(i in v2){
  p1 <- ggplot(boxsub3_core_frame, aes(x =  NL_level, fill = .data[[i]] %>% as.character())) +
    geom_bar(position = "fill") + theme_classic() + 
    facet_grid( ~  xingbie_level) + 
    scale_fill_viridis_d() + 
    guides(fill = guide_legend(title=i))+
    theme(legend.position = "top")
  
  p2 <- ggplot(boxsub3_core_frame, aes(x =  NL_level, fill = .data[[i]] %>% as.character())) +
    geom_bar() + theme_classic() + 
    facet_grid( ~  xingbie_level) + 
    scale_fill_viridis_d() + 
    guides(fill = guide_legend(title=i))+
    theme(legend.position = "top")
  print(p1/p2)
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
# 2382
pdf("./combine/core_distribution_zhifanggan_fill_barplot.pdf",width = 7, height = 8)
i = v2
p1 <- ggplot(boxsub3x_core_frame, aes(x =  NL_level, fill = .data[[i]] %>% as.character())) +
  geom_bar(position = "fill") + theme_classic() + 
  xlab("zhifanggan") +
  facet_grid( ~  xingbie_level) + 
  scale_fill_viridis_d() + 
  guides(fill = guide_legend(title = "zhifanggan"))+
  theme(legend.position = "top")

p2 <- ggplot(boxsub3x_core_frame, aes(x =  NL_level, fill = .data[[i]] %>% as.character())) +
  geom_bar() + theme_classic() + 
  xlab("zhifanggan") +
  facet_grid( ~  xingbie_level) + 
  scale_fill_viridis_d() + 
  guides(fill = guide_legend(title = "zhifanggan"))+
  theme(legend.position = "top")
print(p1/p2)
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

ic <- c("NHDL", "LDL", "NHDL_NL", "LDL_NL", "shousuoya", "shuzhangya", "bmi", "fuwei")

pdf("./combine/face_distribution_boxplot.pdf", height = 6, width = 8)
for(i in ic){
  p1 <- ggboxplot(boxsub1_core_frame, x = "NL_level", y = i, fill = "NL_level") +
    NoLegend() +
    facet_grid(rows = vars(xingbie_level))
  pt <- p1
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
pdf("./combine/face_distribution_boxplotx.pdf", height = 6, width = 8)
for(i in 1:length(ic)){
  yn <- names(ic)[[i]]
  i = ic[[i]]
  sub_frame <- boxsub1x_core_frame[,c("NL_level", "xingbie_level", "left_level", "right_level","max_level", "total_level", i )] %>% na.omit()
  nsp <- nrow(sub_frame)
  p1 <- ggboxplot(sub_frame, x = "NL_level", y = i, fill = "NL_level") +
     NoLegend()  +
    facet_grid(rows = vars(xingbie_level)) + 
    ylab(yn) +
    ggtitle(glue("sample:{nsp}"))
  pt <- p1
  print(pt)
}
dev.off()



