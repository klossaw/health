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

boxsub3_core_frame <- core_frame[, c("smoking_level","NL_level","fuwei_level","bmi_level","LDL_level","xueya_level","xingbie_level",
                                     "left_level", "right_level", "max_level", "total_level","NL")] %>%
  na.omit()
v2 <- c("smoking_level","NL_level","fuwei_level","bmi_level","LDL_level","xueya_level","xingbie_level")
boxsub3_core_frame$right_level <- factor(boxsub3_core_frame$right_level, levels = c("normal", "level1", "level2" ))
boxsub3_core_frame$left_level <- factor(boxsub3_core_frame$left_level, levels = c("normal", "level1", "level2" ))
boxsub3_core_frame$max_level <- factor(boxsub3_core_frame$max_level, levels = c("normal", "level1", "level2" ))

# ======== 分面：年龄x性别，分组：厚度等级，x：特征等级，y：厚度百分比 =========
pdf("./combine/face_yeargender.pdf",width = 12,height = 4)
for(i in v2){
  subf1 <- boxsub3_core_frame %>% dplyr::select("xingbie_level","NL_level", "max_level", {{i}}) %>% 
    group_by(xingbie_level, NL_level, max_level, "{i}" := get(i)) %>% 
    dplyr::summarise(count = n()) %>% ungroup()
  subf2 <- subf1 %>% group_by(xingbie_level, NL_level, "{i}" := get(i)) %>% dplyr::summarise(total_number = sum(count)) %>% ungroup()
  full_sub <- subf1 %>% left_join(subf2) %>% dplyr::mutate(pct = count / total_number * 100)
  print(full_sub %>% 
          ggplot(aes(x = .data[[i]], y = pct, group = max_level, color = max_level)) +
          geom_point(mapping = aes(size = count))+
          geom_smooth(se = F, span = 1) + 
          theme_classic() + 
          facet_grid(xingbie_level ~ NL_level))
}
dev.off()
pdf("./combine/face_yeargender_combine.pdf",width = 12,height = 4)
for(i in v2){
  subf1 <- boxsub3_core_frame %>% 
    mutate(max_level = str_extract(max_level, pattern = "normal|level")) %>% 
    dplyr::select("xingbie_level","NL_level", "max_level", {{i}}) %>% 
    group_by(xingbie_level, NL_level, max_level, "{i}" := get(i)) %>% 
    dplyr::summarise(count = n()) %>% ungroup()
  subf2 <- subf1 %>% group_by(xingbie_level, NL_level, "{i}" := get(i)) %>% dplyr::summarise(total_number = sum(count)) %>% ungroup()
  full_sub <- subf1 %>% left_join(subf2) %>% dplyr::mutate(pct = count / total_number * 100)
  print(full_sub %>% 
          ggplot(aes(x = .data[[i]], y = pct, group = max_level, color = max_level)) +
          geom_point(mapping = aes(size = count))+
          geom_smooth(se = F, span = 1) + 
          theme_classic() + 
          facet_grid(xingbie_level ~ NL_level))
}
dev.off()
# ======== 分面：性别x厚度等级，分组：特征，x：年龄，y：厚度百分比 =============
pdf("./combine/face_genderthin.pdf",width = 12,height = 4)
for(i in v2){
  subf1 <- boxsub3_core_frame %>% dplyr::select("xingbie_level","NL_level", "max_level", {{i}}) %>% 
    group_by(xingbie_level, NL_level, max_level, "{i}" := get(i)) %>% 
    dplyr::summarise(count = n()) %>% ungroup()
  subf2 <- subf1 %>% group_by(xingbie_level, NL_level, "{i}" := get(i)) %>% dplyr::summarise(total_number = sum(count)) %>% ungroup()
  full_sub <- subf1 %>% left_join(subf2) %>% dplyr::mutate(pct = count / total_number * 100)
  full_sub[[i]] <- as.character(full_sub[[i]] )
  print(full_sub %>% 
          ggplot(aes(x = NL_level, y = pct, group = .data[[i]], color = .data[[i]])) +
          geom_point(mapping = aes(size = count))+
          geom_smooth(se = F, span = 1) + 
          theme_classic() + 
          facet_grid(xingbie_level ~ max_level))
}
dev.off()

pdf("./combine/face_genderthincombine.pdf",width = 12,height = 4)
for(i in v2){
  subf1 <- boxsub3_core_frame %>% mutate(max_level = str_extract(max_level, pattern = "normal|level")) %>% 
    dplyr::select("xingbie_level","NL_level", "max_level", {{i}}) %>% 
    group_by(xingbie_level, NL_level, max_level, "{i}" := get(i)) %>% 
    dplyr::summarise(count = n()) %>% ungroup()
  subf2 <- subf1 %>% group_by(xingbie_level, NL_level, "{i}" := get(i)) %>% dplyr::summarise(total_number = sum(count)) %>% ungroup()
  full_sub <- subf1 %>% left_join(subf2) %>% dplyr::mutate(pct = count / total_number * 100)
  full_sub[[i]] <- as.character(full_sub[[i]] )
  print(full_sub %>% 
          ggplot(aes(x = NL_level, y = pct, group = .data[[i]], color = .data[[i]])) +
          geom_point(mapping = aes(size = count))+
          geom_smooth(se = F, span = 1) + 
          theme_classic() + 
          facet_grid(xingbie_level ~ max_level))
}
dev.off()
# ======================= 分层显著性 ===========================================
boxsub3_core_frame <- core_frame[, c("smoking_level","NL_level","fuwei_level","bmi_level","LDL_level","xueya_level","xingbie_level",
                                     "left_level", "right_level", "max_level", "total_level","NL")] %>%
  na.omit()
v2 <- c("smoking_level","NL_level","fuwei_level","bmi_level","LDL_level","xueya_level","xingbie_level")
boxsub3_core_frame$right_level <- factor(boxsub3_core_frame$right_level, levels = c("normal", "level1", "level2" ))
boxsub3_core_frame$left_level <- factor(boxsub3_core_frame$left_level, levels = c("normal", "level1", "level2" ))
boxsub3_core_frame$max_level <- factor(boxsub3_core_frame$max_level, levels = c("normal", "level1", "level2" ))

i <- "smoking_level"
tmpdt <- boxsub3_core_frame %>% 
  dplyr::select(NL_level, xingbie_level, max_level, i) %>% 
  dplyr::mutate(max_level = factor(str_extract(max_level, pattern = "normal|level"), levels = c("level", "normal"))) %>% 
  group_by(NL_level, xingbie_level) %>% 
  dplyr::summarise(glmn = list(glm(formula = glue("max_level ~ {i}") %>% as.formula() , data = tmpdt, family = binomial()))) %>% 
  rowwise() %>% 
  dplyr::mutate()




# ======================= 脂肪肝单做 ===========================================
boxsub3_core_frame <- rawobj[, c("NL_level", "xingbie_level", "脂肪肝", "left_level", "right_level", "max_level", "total_level")] %>%
  na.omit()
v2 <- c("脂肪肝")
boxsub3_core_frame$right_level <- factor(boxsub3_core_frame$right_level, levels = c("normal", "level1", "level2" ))
boxsub3_core_frame$left_level <- factor(boxsub3_core_frame$left_level, levels = c("normal", "level1", "level2" ))
boxsub3_core_frame$max_level <- factor(boxsub3_core_frame$max_level, levels = c("normal", "level1", "level2" ))

# ======== 分面：年龄x性别，分组：厚度等级，x：特征等级，y：厚度百分比 =========
pdf("./combine/face_yeargender_zhifanggan.pdf",width = 12,height = 4)
for(i in v2){
  subf1 <- boxsub3_core_frame %>% dplyr::select("xingbie_level","NL_level", "max_level", {{i}}) %>% 
    group_by(xingbie_level, NL_level, max_level, "{i}" := get(i)) %>% 
    dplyr::summarise(count = n()) %>% ungroup()
  subf2 <- subf1 %>% group_by(xingbie_level, NL_level, "{i}" := get(i)) %>% dplyr::summarise(total_number = sum(count)) %>% ungroup()
  full_sub <- subf1 %>% left_join(subf2) %>% dplyr::mutate(pct = count / total_number * 100)
  print(full_sub %>% 
          ggplot(aes(x = .data[[i]], y = pct, group = max_level, color = max_level)) +
          geom_point(mapping = aes(size = count))+
          geom_smooth(se = F, span = 1) + 
          theme_classic() + 
          facet_grid(xingbie_level ~ NL_level))
}
dev.off()
pdf("./combine/face_yeargender_zhifanggan_combine.pdf",width = 12,height = 4)
for(i in v2){
  subf1 <- boxsub3_core_frame %>% 
    mutate(max_level = str_extract(max_level, pattern = "normal|level")) %>% 
    dplyr::select("xingbie_level","NL_level", "max_level", {{i}}) %>% 
    group_by(xingbie_level, NL_level, max_level, "{i}" := get(i)) %>% 
    dplyr::summarise(count = n()) %>% ungroup()
  subf2 <- subf1 %>% group_by(xingbie_level, NL_level, "{i}" := get(i)) %>% dplyr::summarise(total_number = sum(count)) %>% ungroup()
  full_sub <- subf1 %>% left_join(subf2) %>% dplyr::mutate(pct = count / total_number * 100)
  print(full_sub %>% 
          ggplot(aes(x = .data[[i]], y = pct, group = max_level, color = max_level)) +
          geom_point(mapping = aes(size = count))+
          geom_smooth(se = F, span = 1) + 
          theme_classic() + 
          facet_grid(xingbie_level ~ NL_level))
}
dev.off()
# ======== 分面：性别x厚度等级，分组：特征，x：年龄，y：厚度百分比 =============
pdf("./combine/face_genderthin_zhifanggan.pdf",width = 12,height = 4)
for(i in v2){
  subf1 <- boxsub3_core_frame %>% dplyr::select("xingbie_level","NL_level", "max_level", {{i}}) %>% 
    group_by(xingbie_level, NL_level, max_level, "{i}" := get(i)) %>% 
    dplyr::summarise(count = n()) %>% ungroup()
  subf2 <- subf1 %>% group_by(xingbie_level, NL_level, "{i}" := get(i)) %>% dplyr::summarise(total_number = sum(count)) %>% ungroup()
  full_sub <- subf1 %>% left_join(subf2) %>% dplyr::mutate(pct = count / total_number * 100)
  full_sub[[i]] <- as.character(full_sub[[i]] )
  print(full_sub %>% 
          ggplot(aes(x = NL_level, y = pct, group = .data[[i]], color = .data[[i]])) +
          geom_point(mapping = aes(size = count))+
          geom_smooth(se = F, span = 1) + 
          theme_classic() + 
          facet_grid(xingbie_level ~ max_level))
}
dev.off()

pdf("./combine/face_genderthincombine_zhifanggan.pdf",width = 12,height = 4)
for(i in v2){
  subf1 <- boxsub3_core_frame %>% mutate(max_level = str_extract(max_level, pattern = "normal|level")) %>% 
    dplyr::select("xingbie_level","NL_level", "max_level", {{i}}) %>% 
    group_by(xingbie_level, NL_level, max_level, "{i}" := get(i)) %>% 
    dplyr::summarise(count = n()) %>% ungroup()
  subf2 <- subf1 %>% group_by(xingbie_level, NL_level, "{i}" := get(i)) %>% dplyr::summarise(total_number = sum(count)) %>% ungroup()
  full_sub <- subf1 %>% left_join(subf2) %>% dplyr::mutate(pct = count / total_number * 100)
  full_sub[[i]] <- as.character(full_sub[[i]] )
  print(full_sub %>% 
          ggplot(aes(x = NL_level, y = pct, group = .data[[i]], color = .data[[i]])) +
          geom_point(mapping = aes(size = count))+
          geom_smooth(se = F, span = 1) + 
          theme_classic() + 
          facet_grid(xingbie_level ~ max_level))
}
dev.off()


