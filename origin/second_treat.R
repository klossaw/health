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
# ================================ 第一部分 ====================================
# 9133 个
rawobj <- "/cluster/home/jhuang/projects/healthman/data/chenlei/human/clinical/raw.rds" %>% read_rds
rawobj <- rawobj[rawobj$"是否采纳" == "采纳",]

leftight <- tibble(
  right = rawobj$"右侧颈动脉" %>% as.numeric(),
  left = rawobj$"左侧颈动脉" %>% as.numeric()
)

leftight <- leftight %>% dplyr::mutate(
  left_level = case_when(
    left < 1 ~ "normal",
    left < 1.5 ~ "level1",
    T ~ "level2"
  ),
  right_level = case_when(
    right < 1 ~ "normal",
    right < 1.5 ~ "level1",
    T ~ "level2"
  ),
  total_level = glue("{left_level}_{right_level}")
)
leftight$max_data <- apply(leftight, 1, function(df){
  max(df[[1]] %>% as.numeric(), df[[2]] %>% as.numeric())
}) 
leftight <- leftight %>% dplyr::mutate(max_level = case_when(
  max_data < 1 ~ "normal",
  max_data < 1.5 ~ "level1",
  T ~ "level2"
))

pdf("patients_density.pdf")
ggplot(leftight, aes(x = left, y = right)) +
  geom_pointdensity(size = 0.8) + scale_color_viridis_c() +
  theme_base()
ggplot(leftight, aes(x = left, y = right)) +
  geom_pointdensity(size = 0.8) + scale_color_viridis_c() +
  theme_base() +
  geom_hline(yintercept = c(1, 1.5)) +
  geom_vline(xintercept = c(1, 1.5))
ggplot(leftight, aes(x = left, fill = after_stat(density))) +  stat_density(geom = "bar") +
  theme_base() + scale_fill_viridis_c() + ggtitle("left")
ggplot(leftight, aes(x = left, fill = after_stat(x))) + stat_density(geom = "bar") +
  geom_vline(xintercept = c(1, 1.5)) + ggtitle("left") +
  theme_base() + scale_fill_viridis_c()
ggplot(leftight, aes(x = right, fill = after_stat(density))) + stat_density(geom = "bar") +
  theme_base() + scale_fill_viridis_c() + ggtitle("right")
ggplot(leftight, aes(x = right, fill = after_stat(x))) +  stat_density(geom = "bar")+
  geom_vline(xintercept = c(1, 1.5)) + ggtitle("right") + scale_fill_viridis_c() +
  theme_base() + scale_fill_viridis_c()
ggplot(leftight, aes(x = max_data, fill = after_stat(density))) +  stat_density(geom = "bar") +
  theme_base() + scale_fill_viridis_c() + ggtitle("max_data")
ggplot(leftight, aes(x = max_data, fill = after_stat(x))) +  stat_density(geom = "bar") +
  geom_vline(xintercept = c(1, 1.5)) + ggtitle("max_data") +
  theme_base() + scale_fill_viridis_c()
dev.off()

pdf("patients_bar.pdf")
leftight %>% dplyr::select(c("left_level", "right_level", "max_level")) %>% 
  pivot_longer(cols = c("left_level", "right_level", "max_level")) %>% 
  ggplot(aes(x = name, fill = value)) + geom_bar() +
  theme_base()
dev.off()

leftight$left_level %>% table()
leftight$right_level %>% table()
table(left_level = leftight$left_level, right_level = leftight$right_level)
leftight$max_level %>% table()

# ================================ 第二部分 ====================================
# 9029 个
rawobj$"非高密度脂蛋白C" <- as.numeric(rawobj$"总胆固醇") - as.numeric(rawobj$"高密度脂蛋白C")
# rawobj$NL
leftight_year_NHDC <- tibble(
  right = rawobj$"右侧颈动脉" %>% as.numeric(),
  left = rawobj$"左侧颈动脉" %>% as.numeric(),
  NHDC = rawobj$"非高密度脂蛋白C",
  NL =  rawobj$NL %>% str_remove("岁") %>% as.numeric(),
  LDL =  rawobj$"低密度脂蛋白C" %>% as.numeric()
) %>% na.omit()

leftight_year_NHDC <- leftight_year_NHDC %>% dplyr::mutate(
  left_level = case_when(
    left < 1 ~ "normal",
    left < 1.5 ~ "level1",
    T ~ "level2"
  ),
  right_level = case_when(
    right < 1 ~ "normal",
    right < 1.5 ~ "level1",
    T ~ "level2"
  ),
  total_level = glue("{left_level}_{right_level}")
)
leftight_year_NHDC$max_data <- apply(leftight_year_NHDC, 1, function(df){
  max(df[[1]] %>% as.numeric(), df[[2]] %>% as.numeric())
}) 
leftight_year_NHDC <- leftight_year_NHDC %>% dplyr::mutate(max_level = case_when(
  max_data < 1 ~ "normal",
  max_data < 1.5 ~ "level1",
  T ~ "level2"
))

leftight_year_NHDC$left_level <- factor(leftight_year_NHDC$left_level, levels = c("level2", "level1", "normal"))
leftight_year_NHDC$right_level <- factor(leftight_year_NHDC$right_level, levels = c("level2", "level1", "normal"))
leftight_year_NHDC$max_level <- factor(leftight_year_NHDC$max_level, levels = c("level2", "level1", "normal"))
pdf("patients_year_NHDC.pdf")
ggplot(leftight_year_NHDC %>% arrange(desc(left_level)), aes(x = log(LDL), y = NL, color = left_level)) +
  geom_point(size = 0.8) + scale_color_viridis_d(direction = -1) +
  theme_base() +
  ggtitle("left_level")
ggplot(leftight_year_NHDC %>% arrange(desc(right_level)), aes(x = log(LDL), y = NL, color = right_level)) +
  geom_point(size = 0.8) + scale_color_viridis_d(direction = -1) +
  theme_base() +
  ggtitle("right_level")
ggplot(leftight_year_NHDC %>% arrange(desc(max_level)), aes(x = log(LDL), y = NL, color = max_level)) +
  geom_point(size = 0.8) + scale_color_viridis_d(direction = -1) +
  theme_base() +
  ggtitle("max_level")
ggplot(leftight_year_NHDC %>% arrange(desc(left_level)), aes(x = log(NHDC), y = NL, color = left_level)) +
  geom_point(size = 0.8) + scale_color_viridis_d(direction = -1) +
  theme_base() +
  ggtitle("left_level")
ggplot(leftight_year_NHDC %>% arrange(desc(right_level)), aes(x = log(NHDC), y = NL, color = right_level)) +
  geom_point(size = 0.8) + scale_color_viridis_d(direction = -1) +
  theme_base() +
  ggtitle("right_level")
ggplot(leftight_year_NHDC %>% arrange(desc(max_level)), aes(x = log(NHDC), y = NL, color = max_level)) +
  geom_point(size = 0.8) + scale_color_viridis_d(direction = -1) +
  theme_base() +
  ggtitle("max_level")
dev.off()


lm(left ~ NHDC * NL + LDL * NL,data = leftight_year_NHDC) %>% summary()
lm(right ~ NHDC * NL + LDL * NL,data = leftight_year_NHDC) %>% summary()
lm(max_data ~ NHDC * NL + LDL * NL,data = leftight_year_NHDC) %>% summary()

lm(left ~ NHDC * NL,data = leftight_year_NHDC) %>% summary()
lm(right ~ NHDC * NL,data = leftight_year_NHDC) %>% summary()
lm(max_data ~ NHDC * NL,data = leftight_year_NHDC) %>% summary()

lm(left ~ LDL * NL,data = leftight_year_NHDC) %>% summary()
lm(right ~ LDL * NL,data = leftight_year_NHDC) %>% summary()
lm(max_data ~ LDL * NL,data = leftight_year_NHDC) %>% summary()



# ================================ 第三部分 ====================================
# 2400 个
rawobj <- "/cluster/home/jhuang/projects/healthman/data/chenlei/human/clinical/raw.rds" %>% read_rds
rawobj <- rawobj[rawobj$"是否采纳" == "采纳",]
colnames(rawobj)
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
  xingbie = rawobj$XB
) %>% na.omit()


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

pdf("score.pdf")
leftight_level %>% arrange(smoking_level) %>% 
  ggplot(., aes(x = right, y = left, color = smoking_level)) +
  geom_point() + scale_color_viridis_c() + theme_base()
leftight_level %>% arrange(NL_level) %>% 
  ggplot(., aes(x = right, y = left, color = NL_level)) +
  geom_point() + scale_color_viridis_c() + theme_base()
leftight_level %>% arrange(fuwei_level) %>% 
  ggplot(., aes(x = right, y = left, color = fuwei_level)) +
  geom_point() + scale_color_viridis_c() + theme_base()
leftight_level %>% arrange(bmi_level) %>% 
  ggplot(., aes(x = right, y = left, color = bmi_level)) +
  geom_point() + scale_color_viridis_c() + theme_base()
leftight_level %>% arrange(LDL_level) %>% 
  ggplot(., aes(x = right, y = left, color = LDL_level)) +
  geom_point() + scale_color_viridis_c() + theme_base()
leftight_level %>% arrange(xueya_level) %>% 
  ggplot(., aes(x = right, y = left, color = xueya_level)) +
  geom_point() + scale_color_viridis_c() + theme_base()
dev.off()

leftight_level[,c("left", "right", "smoking_level", "NL_level", "fuwei_level",
                  "bmi_level", "LDL_level" , "xueya_level" )]
lm1 <- lm(left ~ smoking_level + NL_level + fuwei_level + bmi_level + LDL_level + xueya_level, data = leftight_level)
lm1 %>% summary

lm2 <- lm(right ~ smoking_level + NL_level + fuwei_level + bmi_level + LDL_level + xueya_level, data = leftight_level)
lm2 %>% summary

cor(rowSums(leftight_level[,c("smoking_level", "NL_level", "fuwei_level",
                          "bmi_level", "LDL_level" , "xueya_level" )]), 
    leftight_level$right, method = "spearman")
cor(rowSums(leftight_level[,c("smoking_level", "NL_level", "fuwei_level",
                              "bmi_level", "LDL_level" , "xueya_level" )]), 
    leftight_level$left, method = "spearman")
level_dt <- data.frame(
  total_level = rowSums(leftight_level[,c("smoking_level", "NL_level", "fuwei_level",
                                          "bmi_level", "LDL_level" , "xueya_level" )]),
  left = leftight_level$left,
  right = leftight_level$right
)
lm(left ~ total_level,data = level_dt) %>% summary
lm(right ~ total_level,data = level_dt) %>% summary

pdf("score_boxplot.pdf")
my_comparisons <- list(c("0", "1"))
ggboxplot(leftight_level, x = "smoking_level", y = "left", fill = "smoking_level", 
          palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
  stat_compare_means(comparisons = my_comparisons)
ggboxplot(leftight_level, x = "smoking_level", y = "right", fill = "smoking_level", 
          palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
  stat_compare_means(comparisons = my_comparisons)

ggboxplot(leftight_level, x = "NL_level", y = "left", fill = "NL_level")+
  stat_compare_means(method = "anova")
ggboxplot(leftight_level, x = "NL_level", y = "right", fill = "NL_level")+
  stat_compare_means(method = "anova")

ggboxplot(leftight_level, x = "fuwei_level", y = "left", fill = "fuwei_level", 
          palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
  stat_compare_means(comparisons = my_comparisons)
ggboxplot(leftight_level, x = "fuwei_level", y = "right", fill = "fuwei_level", 
          palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
  stat_compare_means(comparisons = my_comparisons)

my_comparisons <- list(c("0", "1"),c("1","2"),c("0","2"))
ggboxplot(leftight_level, x = "bmi_level", y = "left", fill = "bmi_level", 
          palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
  stat_compare_means(comparisons = my_comparisons)
ggboxplot(leftight_level, x = "bmi_level", y = "right", fill = "bmi_level", 
          palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
  stat_compare_means(comparisons = my_comparisons)

ggboxplot(leftight_level, x = "LDL_level", y = "left", fill = "LDL_level")+
  stat_compare_means(method = "anova")
ggboxplot(leftight_level, x = "LDL_level", y = "right", fill = "LDL_level")+
  stat_compare_means(method = "anova")

ggboxplot(leftight_level, x = "xueya_level", y = "left", fill = "xueya_level")+
  stat_compare_means(method = "anova")
ggboxplot(leftight_level, x = "xueya_level", y = "right", fill = "xueya_level")+
  stat_compare_means(method = "anova")
dev.off()

# ================================ 第四部分 ====================================
rawobj <- "/cluster/home/jhuang/projects/healthman/data/chenlei/human/clinical/raw.rds" %>% read_rds
right = rawobj$"右侧颈动脉" %>% as.numeric()
left = rawobj$"左侧颈动脉" %>% as.numeric()
ft <- is.na(right) | is.na(left)
rawobj <- rawobj %>%  
  dplyr::filter(!ft) %>% 
  group_by(TIJIANKABM) %>% 
  dplyr::mutate(test_number = n())
rawobj_manytimes <- rawobj %>% dplyr::filter(test_number > 1)

rawobj_manytimes$TIJIANKABM %>% unique() %>% length()
# 388

rawobj_manytimes <- rawobj_manytimes %>% mutate(check_times = dense_rank(desc(`TIJIANRQ`)))

pdf("rawobj_manytimes.pdf")
ggplot(rawobj_manytimes, aes(x = check_times %>% as.numeric(), y = `左侧颈动脉` %>% as.numeric(), group = TIJIANKABM)) +  
  geom_path() +
  geom_point() + 
  geom_hline(yintercept = c(1, 1.5)) + 
  ggtitle("manytimes") +
  xlab("check_times") +
  ylab("left") +
  theme_base() 
ggplot(rawobj_manytimes, aes(x = check_times %>% as.numeric(), y = `右侧颈动脉` %>% as.numeric(), group = TIJIANKABM)) +  
  geom_path() +
  geom_point() + 
  geom_hline(yintercept = c(1, 1.5)) + 
  ggtitle("manytimes") +
  xlab("check_times") +
  ylab("right") +
  theme_base() 
dev.off()
















