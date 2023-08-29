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

regdata <- core_frame[,c(22, 23, 1:7, 9:12, 16)] %>% dplyr::select(right, left, max_data, xingbie_level, everything())
lm(max_data ~ ., data = regdata %>% 
     dplyr::select(-c("right", "left","NHDL_NL","LDL_NL","xingbie_level")) %>% 
     na.omit() %>% scale %>% as.data.frame()) %>% 
  summary()

regdata %>% 
  dplyr::select(-c("right", "left","NHDL_NL","LDL_NL","xingbie_level", "max_data")) %>% 
  na.omit() %>%  prcomp() -> pcaobj
pcaobj$x
pcaobj$rotation
sdev <- pcaobj$sdev
(sdev * sdev) /sum(sdev * sdev) * 100

df <- pcaobj$x[,1:2] %>% as.data.frame()
df$max_data <- regdata %>% 
  dplyr::select(-c("right", "left","NHDL_NL","LDL_NL","xingbie_level")) %>% 
  na.omit() %>% pull("max_data")
pdf("pca.pdf")
ggplot(df %>% arrange(max_data), aes(x = PC1, y = PC2, color = log(max_data))) + 
  geom_point() + scale_color_viridis_c()
dev.off()


regdata %>% 
  dplyr::select(-c("right", "left","NHDL_NL","LDL_NL","xingbie_level", "max_data")) %>% 
  na.omit() %>% cor(method = "spearman")

lm(LDL ~ NHDL, data = regdata %>% 
     dplyr::select(-c("right", "left","NHDL_NL","LDL_NL","xingbie_level")) %>% 
     na.omit() %>% scale %>% as.data.frame()) %>% 
  summary()
# 正相关
lm(bmi ~ NHDL*LDL, data = regdata %>% 
     dplyr::select(-c("right", "left","NHDL_NL","LDL_NL","xingbie_level")) %>% 
     na.omit() %>% scale %>% as.data.frame()) %>% 
  summary()
# 正相关
lm(fuwei ~ NHDL*LDL, data = regdata %>% 
     dplyr::select(-c("right", "left","NHDL_NL","LDL_NL","xingbie_level")) %>% 
     na.omit() %>% scale %>% as.data.frame()) %>% 
  summary()
# 正相关

checkdata <- regdata %>% dplyr::filter(LDL > 1.8) %>% 
  dplyr::select(-c("right", "left","NHDL_NL","LDL_NL","xingbie_level")) %>% 
  na.omit() %>% as.data.frame()
# lm(fuwei ~ bmi, data = checkdata) %>% 
#   summary()
# checkdata$max_data[checkdata$max_data > 2.5] <- 2.5
# checkdata$LDL[checkdata$LDL > 2.5] <- 2.5
model_nl <- lm(max_data ~ LDL + NL + shousuoya + smoking_level, data = checkdata)
summary(model_nl)
checkdata$delta <- checkdata$max_data - predict(model_nl, checkdata)
pdf("test.pdf")
ggscatter(checkdata %>% arrange(LDL), x = "LDL", y = "delta", add = "reg.line") +
  geom_smooth(method = "lm") +
  stat_cor(label.y = 5, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  stat_regline_equation(label.y = 4)
ggscatter(checkdata %>% arrange(LDL), x = "LDL", y = "max_data", add = "reg.line") +
  geom_smooth(method = "lm") +
  stat_cor(label.y = 5, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  stat_regline_equation(label.y = 4)
ggscatter(checkdata %>% arrange(LDL), x = "LDL", y = "bmi", add = "reg.line") +
  geom_smooth(method = "lm") +
  stat_cor(label.y = 5, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  stat_regline_equation(label.y = 4)
ggscatter(checkdata %>% arrange(LDL), x = "LDL", y = "fuwei", add = "reg.line") +
  geom_smooth(method = "lm") +
  stat_cor(label.y = 5, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  stat_regline_equation(label.y = 4)
ggscatter(checkdata %>% arrange(LDL), x = "LDL", y = "NL", add = "reg.line") +
  geom_smooth(method = "lm") +
  stat_cor(label.y = 5, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  stat_regline_equation(label.y = 4)
dev.off()


checkdata <- regdata %>% dplyr::filter(LDL > 0) %>% 
  dplyr::select(-c("right", "left","NHDL_NL","LDL_NL","xingbie_level")) %>% 
  na.omit() %>% as.data.frame()
pdf("test0.pdf")
checkdata$max_data_level <- case_when(
  is.na(checkdata$max_data) ~ NA,
  checkdata$max_data < 1 ~ 0,
  checkdata$max_data < 1.5 ~ 1,
  T ~ 2
)
ggplot(checkdata %>% arrange(max_data), aes(x = LDL, y = fuwei, color = max_data_level)) + geom_point() + scale_color_viridis_c()
ggplot(checkdata %>% arrange(max_data), aes(x = LDL, y = bmi, color = max_data_level)) + geom_point() + scale_color_viridis_c()
ggplot(checkdata %>% arrange(max_data), aes(x = LDL, y = NHDL, color = max_data_level)) + geom_point() + scale_color_viridis_c()
dev.off()


pdf("testx.pdf")
checkdata <- regdata %>% dplyr::filter(LDL > 0) %>% 
  dplyr::select(-c("right", "left","NHDL_NL","LDL_NL","xingbie_level")) %>% 
  na.omit() %>% as.data.frame()
ggscatter(checkdata %>% arrange(LDL), x = "LDL", y = "max_data", add = "reg.line") +
  geom_smooth(method = "lm") +
  stat_cor(label.y = 5, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  stat_regline_equation(label.y = 4) + scale_color_viridis_c()
checkdata <- regdata %>% dplyr::filter(LDL > 1.8) %>% 
  dplyr::select(-c("right", "left","NHDL_NL","LDL_NL","xingbie_level")) %>% 
  na.omit() %>% as.data.frame()
ggscatter(checkdata %>% arrange(LDL), x = "LDL", y = "max_data", add = "reg.line") +
  geom_smooth(method = "lm") +
  stat_cor(label.y = 5, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  stat_regline_equation(label.y = 4)+ ggtitle("> 1.8")
checkdata <- regdata %>% dplyr::filter(LDL < 1.8) %>% 
  dplyr::select(-c("right", "left","NHDL_NL","LDL_NL","xingbie_level")) %>% 
  na.omit() %>% as.data.frame()
ggscatter(checkdata %>% arrange(LDL), x = "LDL", y = "max_data", add = "reg.line") +
  geom_smooth(method = "lm") +
  stat_cor(label.y = 5, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  stat_regline_equation(label.y = 4)+ ggtitle("< 1.8")
dev.off()





rawobj[,c("脂肪肝", "NL_level", "max_data","LDL")] %>% group_by(NL_level, `脂肪肝`) %>% 
  summarise(max_data = mean(max_data), mean_LDL = mean(LDL %>% na.omit())) %>% 
  dplyr::filter( `脂肪肝` == 3)

rawobj[,c("脂肪肝", "NL_level", "max_data","LDL")] %>% group_by(NL_level, `脂肪肝`) %>% 
  summarise(max_data = mean(max_data), mean_LDL = mean(LDL %>% na.omit()), count = n()) %>% 
  dplyr::filter( NL_level == 3)

(rawobj[rawobj$LDL < 1.8 & rawobj$max_data > 1.5,"脂肪肝"] %>% table)/(rawobj[rawobj$LDL < 1.8 & rawobj$max_data > 1.5,"脂肪肝"] %>% table %>% sum)
(rawobj[rawobj$LDL > 1.8 & rawobj$max_data > 1.5,"脂肪肝"] %>% table)/(rawobj[rawobj$LDL > 1.8 & rawobj$max_data > 1.5,"脂肪肝"] %>% table %>% sum)
(rawobj[rawobj$max_data > 1.5,"脂肪肝"] %>% table)/(rawobj[rawobj$max_data > 1.5,"脂肪肝"] %>% table %>% sum)

write_csv(regdata %>% dplyr::filter(regdata$LDL < 1.8 & regdata$max_data > 1.5), file = "low_LDL.csv")

checkdata <- regdata %>% dplyr::filter(LDL > 0) %>% 
  dplyr::select(-c("right", "left","NHDL_NL","LDL_NL","xingbie_level")) %>% 
  na.omit() %>% as.data.frame()
checkdata$LDL < 1.8 & checkdata$max_data > 1.5
