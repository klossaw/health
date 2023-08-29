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
rawobj <- read_rds("./fourtime/updata_raw.rds") %>% 
  mutate(
    NL = `NL` %>% str_remove("岁") %>% as.numeric(),
    age_bin = str_extract(age_level, pattern = "\\d+") %>% as.numeric()
  ) %>% dplyr::filter(NL < 70)

level1 <- c(
  "max_level",
  "gender_level","age_level", "age",# 性别 年龄
  "drinking_level","smoking_level", # 抽烟 喝酒
  "systolic_blood_pressure_level", "diastolic_pressure_level","blood_pressure_level", # 血压
  "abdominal_circumference_level", "bmi_level", # 腹围bmi
  "white_blood_cell_count_level","neutrophils_level","neutrophil_percentage_level", # 血
  "lymphocyte_percentage_level","monocyte_percentage_level",
  "hemoglobin_level","platelet_count_level","platelet_volume_level","total_protein_level",
  "albumin_level","alanine_aminotransferase_level", "aspartate_aminotransferase_level", # 总蛋白 转氨酶
  "creatinine_level", "urea_level","uric_acid_level", # 肌酐尿素尿酸
  "glomerular_filtration_rate_level", "triglycerides_level", # 肾小球滤过率  甘油三酯
  "TC_level","HDL_level", # 总胆固醇 HDL
  "LDL_level","fasting_blood_sugar_level", # LDL 空腹血糖
  "fatty_liver" # 脂肪肝
)
boxsub3_core_frame <- rawobj[, level1] 
v2 <- level1[-c(1, 2, 3, 4)]
boxsub3_core_frame$max_level <- factor(boxsub3_core_frame$max_level, levels = c("normal", "level1", "level2" ))


# ======== 分面：性别x厚度等级，分组：特征，x：年龄，y：厚度百分比 =============
pdf("./fourtime/fourtime/fig2/age_smooth/face_genderthin.pdf",width = 12,height = 4)
for(i in v2){
  subf1 <- boxsub3_core_frame %>% 
    dplyr::select("gender_level","age_level", "max_level", {{i}}) %>% na.omit() %>% 
    group_by(gender_level, age_level, max_level, "{i}" := get(i)) %>% 
    dplyr::summarise(count = n()) %>% ungroup()
  subf2 <- subf1 %>% group_by(gender_level, age_level, "{i}" := get(i)) %>% dplyr::summarise(total_number = sum(count)) %>% ungroup()
  full_sub <- subf1 %>% left_join(subf2) %>% dplyr::mutate(pct = count / total_number * 100)
  full_sub[[i]] <- as.character(full_sub[[i]] )
  print(full_sub %>% 
          ggplot(aes(x = age_level, y = pct, group = .data[[i]], color = .data[[i]])) +
          geom_point(mapping = aes(size = count))+
          geom_smooth(se = F, method = "lm") + 
          theme_classic() + 
          facet_grid(gender_level ~ max_level))
}
dev.off()
pdf("./fourtime/fourtime/fig2/age_smooth/smooth_face_genderthin.pdf",width = 12,height = 4)
for(i in v2){
  subf1 <- boxsub3_core_frame %>% 
    dplyr::select("gender_level","age_level", "max_level", {{i}}) %>% na.omit() %>% 
    group_by(gender_level, age_level, max_level, "{i}" := get(i)) %>% 
    dplyr::summarise(count = n()) %>% ungroup()
  subf2 <- subf1 %>% group_by(gender_level, age_level, "{i}" := get(i)) %>% dplyr::summarise(total_number = sum(count)) %>% ungroup()
  full_sub <- subf1 %>% left_join(subf2) %>% dplyr::mutate(pct = count / total_number * 100)
  full_sub[[i]] <- as.character(full_sub[[i]] )
  print(full_sub %>% 
          ggplot(aes(x = age_level, y = pct, group = .data[[i]], color = .data[[i]])) +
          geom_point(mapping = aes(size = count))+
          geom_smooth(se = F, span = 1) + 
          theme_classic() + 
          facet_grid(gender_level ~ max_level))
}
dev.off()
pdf("./fourtime/fourtime/fig2/age_smooth/age_face_genderthin.pdf",width = 12,height = 4)
for(i in v2){
  subf1 <- boxsub3_core_frame %>% 
    dplyr::select("gender_level","age", "max_level", {{i}}) %>% na.omit() %>% 
    group_by(gender_level, age, max_level, "{i}" := get(i)) %>% 
    dplyr::summarise(count = n()) %>% ungroup()
  subf2 <- subf1 %>% group_by(gender_level, age, "{i}" := get(i)) %>% dplyr::summarise(total_number = sum(count)) %>% ungroup()
  full_sub <- subf1 %>% left_join(subf2) %>% dplyr::mutate(pct = count / total_number * 100)
  full_sub[[i]] <- as.character(full_sub[[i]] )
  print(full_sub %>% 
          ggplot(aes(x = age, y = pct, group = .data[[i]], color = .data[[i]])) +
          geom_point(mapping = aes(size = count))+
          geom_smooth(se = F, span = 1) + 
          theme_classic() + 
          facet_grid(gender_level ~ max_level))
}
dev.off()




pdf("./fourtime/fourtime/fig2/age_smooth/face_genderthincombine.pdf",width = 12,height = 4)
for(i in v2){
  subf1 <- boxsub3_core_frame %>% mutate(max_level = str_extract(max_level, pattern = "normal|level")) %>% 
    dplyr::select("gender_level","age_level", "max_level", {{i}}) %>% na.omit() %>% 
    group_by(gender_level, age_level, max_level, "{i}" := get(i)) %>% 
    dplyr::summarise(count = n()) %>% ungroup()
  subf2 <- subf1 %>% group_by(gender_level, age_level, "{i}" := get(i)) %>% 
    dplyr::summarise(total_number = sum(count)) %>% ungroup()
  full_sub <- subf1 %>% left_join(subf2) %>% dplyr::mutate(pct = count / total_number * 100)
  full_sub[[i]] <- as.character(full_sub[[i]] )
  print(full_sub %>% 
          ggplot(aes(x = age_level, y = pct, group = .data[[i]], color = .data[[i]])) +
          geom_point(mapping = aes(size = count))+
          geom_smooth(se = F, method = "lm") + 
          theme_classic() + 
          facet_grid(gender_level ~ max_level))
}
dev.off()
pdf("./fourtime/fourtime/fig2/age_smooth/age_face_genderthincombine.pdf",width = 12,height = 4)
for(i in v2){
  subf1 <- boxsub3_core_frame %>% mutate(max_level = str_extract(max_level, pattern = "normal|level")) %>% 
    dplyr::select("gender_level","age", "max_level", {{i}}) %>% na.omit() %>% 
    group_by(gender_level, age, max_level, "{i}" := get(i)) %>% 
    dplyr::summarise(count = n()) %>% ungroup()
  subf2 <- subf1 %>% group_by(gender_level, age, "{i}" := get(i)) %>% 
    dplyr::summarise(total_number = sum(count)) %>% ungroup()
  full_sub <- subf1 %>% left_join(subf2) %>% dplyr::mutate(pct = count / total_number * 100)
  full_sub[[i]] <- as.character(full_sub[[i]] )
  print(full_sub %>% 
          ggplot(aes(x = age, y = pct, group = .data[[i]], color = .data[[i]])) +
          geom_point(mapping = aes(size = count))+
          geom_smooth(se = F, span = 1) + 
          theme_classic() + 
          facet_grid(gender_level ~ max_level))
}
dev.off()
pdf("./fourtime/fourtime/fig2/age_smooth/smooth_face_genderthincombine.pdf",width = 12,height = 4)
for(i in v2){
  subf1 <- boxsub3_core_frame %>% mutate(max_level = str_extract(max_level, pattern = "normal|level")) %>% 
    dplyr::select("gender_level","age_level", "max_level", {{i}}) %>% na.omit() %>% 
    group_by(gender_level, age_level, max_level, "{i}" := get(i)) %>% 
    dplyr::summarise(count = n()) %>% ungroup()
  subf2 <- subf1 %>% group_by(gender_level, age_level, "{i}" := get(i)) %>% 
    dplyr::summarise(total_number = sum(count)) %>% ungroup()
  full_sub <- subf1 %>% left_join(subf2) %>% dplyr::mutate(pct = count / total_number * 100)
  full_sub[[i]] <- as.character(full_sub[[i]] )
  print(full_sub %>% 
          ggplot(aes(x = age_level, y = pct, group = .data[[i]], color = .data[[i]])) +
          geom_point(mapping = aes(size = count))+
          geom_smooth(se = F, span = 1) + 
          theme_classic() + 
          facet_grid(gender_level ~ max_level))
}
dev.off()



# ======== 分面：性别，分组：特征，x：年龄，y：厚度分位数 =============
level1 <- c(
  "max_data",
  "gender_level","age_level", "age",# 性别 年龄
  "drinking_level","smoking_level", # 抽烟 喝酒
  "systolic_blood_pressure_level", "diastolic_pressure_level","blood_pressure_level", # 血压
  "abdominal_circumference_level", "bmi_level", # 腹围bmi
  "white_blood_cell_count_level","neutrophils_level","neutrophil_percentage_level", # 血
  "lymphocyte_percentage_level","monocyte_percentage_level",
  "hemoglobin_level","platelet_count_level","platelet_volume_level","total_protein_level",
  "albumin_level","alanine_aminotransferase_level", "aspartate_aminotransferase_level", # 总蛋白 转氨酶
  "creatinine_level", "urea_level","uric_acid_level", # 肌酐尿素尿酸
  "glomerular_filtration_rate_level", "triglycerides_level", # 肾小球滤过率  甘油三酯
  "TC_level","HDL_level", # 总胆固醇 HDL
  "LDL_level","fasting_blood_sugar_level", # LDL 空腹血糖
  "fatty_liver" # 脂肪肝
)
boxsub3_core_frame <- rawobj[, level1] 
v2 <- level1[-c(1, 2, 3, 4)]

pdf("./fourtime/fourtime/fig2/age_smooth/face_genderthin_qt.pdf",width = 12,height = 4)
for(i in v2){
  framelist_gender <- lapply(seq(0.5, 0.95, 0.05), function(qt){
    boxsub3_core_frame %>% 
      group_by(gender_level, age_level, "{i}" := get(i)) %>% 
      dplyr::summarise(
        "qt_max_dataht_{qt}" := quantile(max_data, {{ qt }})
        )
  })
  framesheet_gender <- purrr::reduce(framelist_gender, left_join)
  maxframe_gender <- framesheet_gender %>% 
    dplyr::select(1, 2, 3, contains("max")) %>% na.omit() %>% 
    pivot_longer(-c("age_level", "gender_level", i))
  maxframe_gender[[i]] <- as.character(maxframe_gender[[i]])
  print(maxframe_gender %>% 
          ggplot(aes(x = age_level, y = value, group = .data[[i]], color = .data[[i]])) +
          geom_point()+
          geom_smooth(se = F, method = "lm") + 
          theme_classic() + 
          facet_grid(gender_level ~ name))
}
dev.off()


pdf("./fourtime/fourtime/fig2/age_smooth/age_face_genderthin_qt.pdf",width = 12,height = 4)
for(i in v2){
  framelist_gender <- lapply(seq(0.5, 0.95, 0.05), function(qt){
    boxsub3_core_frame %>% 
      group_by(gender_level, age, "{i}" := get(i)) %>% 
      dplyr::summarise(
        "qt_max_dataht_{qt}" := quantile(max_data, {{ qt }})
      )
  })
  framesheet_gender <- purrr::reduce(framelist_gender, left_join)
  maxframe_gender <- framesheet_gender %>% 
    dplyr::select(1, 2, 3, contains("max")) %>% na.omit() %>% 
    pivot_longer(-c("age", "gender_level", i))
  maxframe_gender[[i]] <- as.character(maxframe_gender[[i]])
  print(maxframe_gender %>% 
          ggplot(aes(x = age, y = value, group = .data[[i]], color = .data[[i]])) +
          geom_point()+
          geom_smooth(se = F, span = 1) + 
          theme_classic() + 
          facet_grid(gender_level ~ name))
}
dev.off()

pdf("./fourtime/fourtime/fig2/age_smooth/smooth_face_genderthin_qt.pdf",width = 12,height = 4)
for(i in v2){
  framelist_gender <- lapply(seq(0.5, 0.95, 0.05), function(qt){
    boxsub3_core_frame %>% 
      group_by(gender_level, age_level, "{i}" := get(i)) %>% 
      dplyr::summarise(
        "qt_max_dataht_{qt}" := quantile(max_data, {{ qt }})
      )
  })
  framesheet_gender <- purrr::reduce(framelist_gender, left_join)
  maxframe_gender <- framesheet_gender %>% 
    dplyr::select(1, 2, 3, contains("max")) %>% na.omit() %>% 
    pivot_longer(-c("age_level", "gender_level", i))
  maxframe_gender[[i]] <- as.character(maxframe_gender[[i]])
  print(maxframe_gender %>% 
          ggplot(aes(x = age_level, y = value, group = .data[[i]], color = .data[[i]])) +
          geom_point()+
          geom_smooth(se = F, span = 1) + 
          theme_classic() + 
          facet_grid(gender_level ~ name))
}
dev.off()





