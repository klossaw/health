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
rawobj <- read_rds("./fourtime/updata_raw.rds")
# core_frame <- read_rds("./fourtime/updata_core.rds")


# level1 <- c(
#   "left_level","right_level","max_level",
#   "gender_level","age_level", # 性别 年龄
#   "drinking_level","smoking_level", # 抽烟 喝酒
#   "systolic_blood_pressure_level", "diastolic_pressure_level","blood_pressure_level", # 血压
#   "abdominal_circumference_level", "bmi_level", # 腹围bmi
#   "white_blood_cell_count_level","neutrophils_level","neutrophil_percentage_level", # 血
#   "lymphocyte_percentage_level","monocyte_percentage_level",
#   "hemoglobin_level","platelet_count_level","platelet_volume_level","total_protein_level",
#   "albumin_level","alanine_aminotransferase_level", "aspartate_aminotransferase_level", # 总蛋白 转氨酶
#   "creatinine_level", "urea_level","uric_acid_level", # 肌酐尿素尿酸
#   "glomerular_filtration_rate_level", "triglycerides_level", # 肾小球滤过率  甘油三酯
#   "TC_level","HDL_level", # 总胆固醇 HDL
#   "LDL_level","fasting_blood_sugar_level", # LDL 空腹血糖
#   "fatty_liver" # 脂肪肝
# )

# nolevel1 <- c(
#   "right", "left", "max_data",
#   "gender_level","age_level", "age", # 性别 年龄
#   "systolic_blood_pressure","diastolic_pressure", # 收缩压 舒张压
#   "heart_rate","height","weight", # 心率 升高 体重
#   "abdominal_circumference","bmi", # 腹围 BMI
#   "white_blood_cell_count","neutrophils","neutrophil_percentage", # 血
#   "lymphocyte_percentage","monocyte_percentage","hemoglobin","platelet_count","platelet_volume",
#   "urinary_microalbumin", "urine_microalbuminuria_creatinine_ratio", "total_protein", "albumin", #尿微量蛋白  白蛋白 总蛋白 
#   "alanine_aminotransferase", # 谷丙转氨酶
#   "aspartate_aminotransferase", # 谷草转氨酶
#   "creatinine", # 肌酐
#   "urea","uric_acid", # 尿素 尿酸
#   "glomerular_filtration_rate","homocysteine", # 肾小球滤过率 同型半胱氨酸
#   "triglycerides","fasting_blood_sugar", # 甘油三酯  空腹血糖
#   "apolipoproteinA1","apolipoproteinB", # 载脂蛋白A1 载脂蛋白B
#   "apolipoproteinE","glycated_hemoglobinA1", # 载脂蛋白E 糖化血红蛋白A1
#   "glycated_hemoglobinA1C","fasting_C_peptide", # 糖化血红蛋白A1C  空腹C肽
#   "fasting_insulin","hsc_reactive_protein", # 空腹胰岛素  超敏C反应蛋白
#   "sialic_acid","free_fatty_acid",# 唾液酸 游离脂肪酸
#   "Hydroxyvitamin_D3", # 羟基维生素D3 
#   "TC","HDL","LDL","NHDL", # TC HDL LDL NHDL
#   "NHDL_age", "LDL_age"  # NHDL*age LDL*age
# )
# =========================== 总体盒子图 =======================================
nolevel <- c(
  "max_data",
  "gender_level","age_level", # 性别 年龄
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

boxsub1_nolevel_frame <- rawobj[, nolevel] 
boxsub1_nolevel_frame$max_data <- log10(boxsub1_nolevel_frame$max_data)
ic <- colnames(boxsub1_nolevel_frame)[-1]
# my_comparisons <- list(c("normal", "level1"),c("level1","level2"),c("normal","level2"))

pt <- list()
for(i in ic){
  boxsub1_nolevel_framex <- boxsub1_nolevel_frame %>% 
    dplyr::select("max_data", i) %>% 
    na.omit()
  boxsub1_nolevel_framex[[i]] <- as.factor(boxsub1_nolevel_framex[[i]])
  pt[[i]] <- ggboxplot(boxsub1_nolevel_framex, x = i, y = "max_data", fill = i) +
    theme_base() + ggtitle(glue("sample_size:{nrow(boxsub1_nolevel_framex)}")) +
    theme(legend.position = "top")
}
multi_plot(pt, fig_fn = "./fourtime/fourtime/fig1/total_box/updata_score1_boxplot2.pdf", width = 12, height = 10)


# ============================ 盒子分页 ========================================
boxsub1_nolevel_frame <- rawobj[, nolevel] 
boxsub1_nolevel_frame$max_data <- log10(boxsub1_nolevel_frame$max_data)
ic <- colnames(boxsub1_nolevel_frame)[-c(1, 2, 3)]
# my_comparisons <- list(c("normal", "level1"),c("level1","level2"),c("normal","level2"))
pdf(glue("./fourtime/fourtime/fig1/total_box/updata_score1_boxplot_face2.pdf"), width = 12, height = 4)
for(i in ic){
  boxsub1_nolevel_framex <- boxsub1_nolevel_frame %>% 
    dplyr::select("max_data", i,"gender_level","age_level") %>% 
    na.omit()
  boxsub1_nolevel_framex[[i]] <- as.factor(boxsub1_nolevel_framex[[i]])
  print(
    ggboxplot(boxsub1_nolevel_framex, x = i, y = "max_data", fill = i) +
      # stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
      theme_base() + 
      ggtitle(glue("sample_size:{nrow(boxsub1_nolevel_framex)}")) +
      facet_grid(gender_level ~ age_level)
  )
}
dev.off()

# =========================== 总体分位数图 =======================================
nolevel <- c(
  "max_data",
  "gender_level","age_level",
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

boxsub1_pct_frame <- rawobj[, nolevel] 
ic <- colnames(boxsub1_pct_frame)[-1]
pt <- list()
for(i in ic){
  boxsub1_pctl_framex <- boxsub1_pct_frame %>% 
    dplyr::select("max_data", i) %>% 
    na.omit()

  framelist <- lapply(c(0.25, 0.5, 0.75, 0.90), function(qt){
    boxsub1_pctl_framex %>% group_by(!!sym(i)) %>% dplyr::summarise(
      "max_data_{qt}_{i}" := quantile(max_data, {{ qt }})
    )
  })
  framesheet <- purrr::reduce(framelist, left_join)
  maxframe <- framesheet %>% pivot_longer(-i) 
  pt[[i]] <- ggplot(maxframe, aes(x = .data[[i]], y = value, color = name, group = name)) +
    geom_point() + 
    theme_base() +
    theme(legend.position = "top") +
    scale_color_viridis_d() + 
    geom_line() +
    ggtitle(glue("sample_size:{nrow(boxsub1_pctl_framex)}"))
}
multi_plot(pt, fig_fn = "./fourtime/fourtime/fig1/total_box/updata_pct_plot2.pdf",
           width = 12, height = 10,
           ncol = 2, nrow = 2)


# ============================ 分位数分页 ========================================
boxsub1_pct_frame <- rawobj[, nolevel] 
ic <- colnames(boxsub1_pct_frame)[-c(1, 2, 3)]
pdf(glue("./fourtime/fourtime/fig1/total_box/updata_pct_plot_face2.pdf"), width = 12, height = 4)
for(i in ic){
  boxsub1_nolevel_framex <- boxsub1_pct_frame %>% 
    dplyr::select("max_data", i,"gender_level","age_level") %>% 
    na.omit()
  framelist <- lapply(c(0.25, 0.5, 0.75, 0.90), function(qt){
    boxsub1_nolevel_framex %>% group_by(!!sym(i), gender_level, age_level) %>% dplyr::summarise(
      "max_data_{qt}_{i}" := quantile(max_data, {{ qt }})
    )
  })
  framesheet <- purrr::reduce(framelist, left_join)
  maxframe <- framesheet %>% 
    pivot_longer(-c(i,"gender_level","age_level")) 
  
  print(
    ggplot(maxframe, aes(x = .data[[i]], y = value, color = name, group = name)) +
      geom_point() + 
      theme_base() +
      theme(legend.position = "top") +
      scale_color_viridis_d() + 
      geom_line() +
      ggtitle(glue("sample_size:{nrow(boxsub1_pctl_framex)}")) + 
      facet_grid(gender_level ~ age_level)
  )
}
dev.off()
# =========================== 总体柱子图 =======================================
level1 <- c(
  "max_level",
  "gender_level","age_level", # 性别 年龄
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
boxsub1_level_frame <- rawobj[, level1] 
boxsub1_level_frame$max_level <- factor(boxsub1_level_frame$max_level, levels = c("normal", "level1", "level2"))

ic <- colnames(boxsub1_level_frame)[-1]
pt <- list()
for(i in ic){
  boxsub1_level_framex <- boxsub1_level_frame %>% 
    dplyr::select("max_level", i) %>% 
    na.omit() 
  boxsub1_level_framex[[i]] <- as.character(boxsub1_level_framex[[i]])
  pt[[i]] <- ggplot(boxsub1_level_framex, aes(x = .data[[i]], fill = max_level)) +
    geom_bar(position = "fill") + 
    theme_base()  + theme (legend.position= "top") + 
    ggtitle(glue("sample_size:{nrow(boxsub1_level_framex)}"))
}
multi_plot(pt, 
           fig_fn = "./fourtime/fourtime/fig1/total_box/updata_score1_barplot2.pdf", 
           width = 12, height = 10,
           ncol = 2, nrow = 2)

# ============================ 柱子分页 ========================================
level1 <- c(
  "max_level",
  "gender_level","age_level", # 性别 年龄
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
boxsub1_nolevel_frame <- rawobj[, level1] 
boxsub1_nolevel_frame$max_level <- factor(boxsub1_nolevel_frame$max_level, 
                                          levels = c("normal", "level1", "level2"))
ic <- colnames(boxsub1_nolevel_frame)[-c(1, 2, 3)]
pdf(glue("./fourtime/fourtime/fig1/total_box/updata_score1_barplot_face2.pdf"), width = 12, height = 4)
for(i in ic){
  boxsub1_level_framex <- boxsub1_level_frame %>% 
    dplyr::select("max_level", i,"gender_level","age_level") %>% 
    na.omit() 
  boxsub1_level_framex[[i]] <- as.character(boxsub1_level_framex[[i]])
  print(
    ggplot(boxsub1_level_framex, aes(x = .data[[i]], fill = max_level)) +
      geom_bar(position = "fill") + 
      theme_base()  + theme (legend.position= "top") + 
      ggtitle(glue("sample_size:{nrow(boxsub1_level_framex)}")) +
      facet_grid(gender_level ~ age_level)
  )
}
dev.off()
# ============================= 百分比总体 =====================================
level1 <- c(
  "max_level",
  "gender_level","age_level", # 性别 年龄
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
boxsub1_nolevel_frame <- rawobj[, level1] 
boxsub1_nolevel_frame$max_level <- factor(boxsub1_nolevel_frame$max_level, 
                                          levels = c("normal", "level1", "level2"))
ic <- colnames(boxsub1_nolevel_frame)[-c(1, 2, 3)]
pdf("./fourtime/fourtime/fig1/total_box/updata_score1_pct2.pdf",width = 5, height = 4)
for(i in ic){
  subf1 <- boxsub1_nolevel_frame %>% dplyr::select("max_level", {{i}}) %>% 
    group_by(max_level, "{i}" := get(i)) %>% 
    dplyr::summarise(count = n()) %>% ungroup()
  subf2 <- subf1 %>% group_by("{i}" := get(i)) %>% 
    dplyr::summarise(total_number = sum(count)) %>% ungroup()
  full_sub <- subf1 %>% left_join(subf2) %>% dplyr::mutate(pct = count / total_number * 100)
  print(full_sub %>% na.omit() %>% 
          ggplot(aes(x = .data[[i]], y = pct, group = max_level, color = max_level)) +
          geom_point(mapping = aes(size = count))+
          geom_line() + 
          theme_classic() )
}
dev.off()

# ======================== 百分比分面 ==========================================
level1 <- c(
  "max_level",
  "gender_level","age_level", # 性别 年龄
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
boxsub1_nolevel_frame <- rawobj[, level1] 
boxsub1_nolevel_frame$max_level <- factor(boxsub1_nolevel_frame$max_level, 
                                          levels = c("normal", "level1", "level2"))
ic <- colnames(boxsub1_nolevel_frame)[-c(1, 2, 3)]
pdf("./fourtime/fourtime/fig1/total_box/updata_score1_pct_face2.pdf",width = 12,height = 4)
for(i in ic){
  subf1 <- boxsub1_nolevel_frame %>% dplyr::select("gender_level","age_level", "max_level", {{i}}) %>% 
    group_by(gender_level, age_level, max_level, "{i}" := get(i)) %>% 
    dplyr::summarise(count = n()) %>% ungroup()
  subf2 <- subf1 %>% group_by(gender_level, age_level, "{i}" := get(i)) %>% dplyr::summarise(total_number = sum(count)) %>% ungroup()
  full_sub <- subf1 %>% left_join(subf2) %>% dplyr::mutate(pct = count / total_number * 100)
  print(full_sub %>% na.omit() %>% 
          ggplot(aes(x = .data[[i]], y = pct, group = max_level, color = max_level)) +
          geom_point(mapping = aes(size = count))+
          geom_line() + 
          theme_classic() + 
          facet_grid(gender_level ~ age_level))
}
dev.off()


