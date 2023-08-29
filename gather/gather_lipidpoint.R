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
rawobj <- read_csv( "updata_raw.csv") %>% dplyr::select(-1)
core_frame <- read_csv("updata_core.csv") %>% dplyr::select(-1)
dir.create("./lipidpoint/")

vm <- c("LDL", "脂肪肝",
        "心率", "身高", "体重",
        "同型半胱氨酸","甘油三酯","总胆固醇","高密度脂蛋白C","载脂蛋白A1","载脂蛋白B","尿酸",
        "空腹血糖","糖化血红蛋白A1","糖化血红蛋白A1C","空腹C肽","空腹胰岛素","超敏C反应蛋白",
        "唾液酸","游离脂肪酸","载脂蛋白E","羟基维生素D3")
names(vm) <- c(
  "LDL", "zhifanggan",
  "xinlv","shengao","tizhong",
  "tongxing","ganyousanzhi","zongdanguchun","gaomidu","zaizhia1","zaizhib","niaosuan",
  "kongfuxuetang","tanghuaA1","tanghuaA1C","kongfuC","kongfuyidaosu","chaomingC",
  "tuoyesuan","youlizhifangsuan","zaizhiE","qiangjiVD3"
)
point_core_frame <- rawobj[, c(vm, "max_level")]
for(i in vm){
  for(j in vm){
    pdf(glue::glue("./lipidpoint/core_{i}_{j}_pointplot.pdf"))
    p1 <- ggplot(point_core_frame %>% arrange(point_core_frame$max_data), aes(x = .data[[i]], y = .data[[j]], color = max_level)) +
      geom_point() + scale_color_viridis_d()
    print(p1)
    dev.off()
  }
}


vm0 <- c("LDL", "脂肪肝",
        "心率", "身高", "体重",
        "同型半胱氨酸","甘油三酯","总胆固醇","高密度脂蛋白C","载脂蛋白A1","载脂蛋白B","尿酸",
        "空腹血糖","糖化血红蛋白A1","糖化血红蛋白A1C","空腹C肽","空腹胰岛素","超敏C反应蛋白",
        "唾液酸","游离脂肪酸","载脂蛋白E","羟基维生素D3",
        "smoking","shousuoya","shuzhangya", "bmi","fuwei","NL",
        "max_data")
point_core_frame <- rawobj[, c(vm0)]
apply(point_core_frame, 2, function(x) sum(!is.na(x)))



vm0 <- c("LDL", "脂肪肝",
         "心率", "身高", "体重",
         "同型半胱氨酸","甘油三酯","总胆固醇","高密度脂蛋白C","载脂蛋白A1","载脂蛋白B","尿酸",
         "smoking","shousuoya","shuzhangya", "bmi","fuwei","NL",
         "max_data")
point_core_frame <- rawobj[, c(vm0)]%>% na.omit
dim(point_core_frame)
apply(point_core_frame[,c("smoking","shousuoya","shuzhangya", "bmi","fuwei","NL")] %>% na.omit,
      2, function(x) sum(!is.na(x)))

apply(point_core_frame[!is.na(point_core_frame$smoking) & 
                         !is.na(point_core_frame$shousuoya) & 
                         !is.na(point_core_frame$shuzhangya) & 
                         !is.na(point_core_frame$bmi) & 
                         !is.na(point_core_frame$fuwei) & 
                         !is.na(point_core_frame$NL) & 
                         !is.na(point_core_frame$LDL),],
      2, function(x) sum(!is.na(x)))
pcaobj <- prcomp(point_core_frame[,c("LDL", "脂肪肝",
                           "身高", "体重",
                           "同型半胱氨酸","甘油三酯","总胆固醇","高密度脂蛋白C",
                           "载脂蛋白A1","载脂蛋白B","NL")])

pcaf <- pcaobj$x[,c(1,2,3)] %>% as.data.frame()
pcaf$max_data <- point_core_frame$max_data
lm(max_data ~ ., data = pcaf) %>% summary()
pcaobj$rotation

pcalm <- lm(max_data ~ ., data = pcaf)
pdf("pca.pdf")
data.frame(
  data = pcaf$max_data,
  predict = predict(pcalm, pcaf)
) %>% ggplot(aes(x = data, y = predict)) + geom_point()
dev.off()




cor(point_core_frame$NL, pcaobj$x, method = "spearman")
cor(pcaobj$x)



