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

abnormalsheet <- rawobj %>% dplyr::mutate(abnormal = case_when(
  LDL < 1.8 & max_data >= 1 ~ "loLDL_hiMax",
  LDL < 1.8 & max_data < 1 ~ "loLDL_loMax",
  LDL >= 1.8 & max_data >= 1 ~ "hiLDL_hiMax",
  LDL >= 1.8 & max_data < 1 ~ "hiLDL_loMax",
  T ~ "other"))
core_abnormalsheet <- abnormalsheet[,c(1,2,3,4, 13, 61:85)]
write_rds(abnormalsheet, file = "./combine/abnormalsheet.rds")
write_rds(core_abnormalsheet, file = "./combine/core_abnormalsheet.rds")

# ====================== 其他指标 ==============================================
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
abnormalsheetx <- abnormalsheet[, c(vm,"abnormal")] %>% dplyr::filter(abnormal != "other")
abnormalsheetx$abnormal <- factor(abnormalsheetx$abnormal, levels = c("loLDL_hiMax", "hiLDL_hiMax", "loLDL_loMax", "hiLDL_loMax"))

ic <- vm
my_comparisons <- list(c("loLDL_hiMax","hiLDL_loMax"),c("loLDL_hiMax","loLDL_loMax"), c("loLDL_hiMax", "hiLDL_hiMax"))


pl <- list()
for(i in 1:length(ic)){
  yn <- names(ic)[[i]]
  i = ic[[i]]
  sub_frame <- abnormalsheetx[,c("abnormal", i )] %>% na.omit()
  nsp <- nrow(sub_frame)
  pl[[i]] <- ggboxplot(sub_frame, x = "abnormal", y = i, fill = "abnormal") +
    stat_compare_means(comparisons = my_comparisons) + NoLegend() +
    ylab(yn) +
    ggtitle(glue("sample:{nsp}")) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2))
}
multi_plot(pl, fig_fn = "./combine/other_loLDL_hiMax_boxplotx.pdf", width = 12, height = 10)

# ====================== 核心指标 ==============================================
boxsub1_core_frame <- core_abnormalsheet[, c("NL_level", "xingbie_level",
                                     "NHDL", "LDL", "NHDL_NL", "LDL_NL", 
                                     "shousuoya", "shuzhangya", "bmi", "fuwei", "NL",
                                     "abnormal")] %>% na.omit()
nrow(boxsub1_core_frame)
# 2391
boxsub1_core_frame$abnormal <- factor(boxsub1_core_frame$abnormal, levels = c("loLDL_hiMax", "hiLDL_hiMax", "loLDL_loMax", "hiLDL_loMax"))

ic <- c("NHDL", "LDL", "NHDL_NL", "LDL_NL", "shousuoya", "shuzhangya", "bmi", "fuwei","NL")
my_comparisons <- list(c("loLDL_hiMax","hiLDL_loMax"),c("loLDL_hiMax","loLDL_loMax"), c("loLDL_hiMax", "hiLDL_hiMax"))

pl <- list()
for(i in ic){
  pl[[i]] <- ggboxplot(boxsub1_core_frame, x = "abnormal", y = i, fill = "abnormal") +
    stat_compare_means(comparisons = my_comparisons) + NoLegend() +
    scale_x_discrete(guide = guide_axis(n.dodge = 2))
}
multi_plot(pl, fig_fn = "./combine/core_loLDL_hiMax_boxplotx.pdf", width = 12, height = 10)



