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
getrootframe <- function(df, y = 1){
  df %>% group_by(group, PANEL) %>% 
    mutate(delta = abs(y - {{y}}))  %>% 
    slice_min(order_by = delta, n = 2) %>% 
    summarise( root = ifelse(x[1] * x[2] < 0, mean(x), x[1])) %>% 
    ungroup() 
}

rawobj <- read_rds("./fourtime/updata_raw.rds") %>% 
  mutate(
    NL = `NL` %>% str_remove("岁") %>% as.numeric()
) %>% dplyr::filter(NL < 70)

framelist <- lapply(seq(0.5, 0.95, 0.05), function(qt){
  rawobj %>% group_by(NL) %>% dplyr::summarise(
    "qt_max_dataht_{qt}" := quantile(max_data, {{ qt }})
  )
})
framesheet <- purrr::reduce(framelist, left_join)
maxframe <- framesheet %>% dplyr::select(1, contains("max")) %>% 
  pivot_longer(-"NL")
p2 <- ggplot() +
  geom_point(rawobj, mapping = aes(x = NL, y = max_data)) +
  geom_smooth(maxframe,
              mapping = aes(x = NL, y = value, color = name),
              linewidth = 1, se = F) + theme_base() +
  geom_hline(yintercept = c(1,1.5), color = "orange", size = 1) +
  geom_vline(xintercept = c(30, 40), color = "orange", size = 1) +
  theme(legend.position = "top")+
  scale_color_viridis_d()
pdf("./fourtime/fourtime/fig2/age_smooth/NL_lm3.pdf")
print(p2)
dev.off()
bp <- ggplot_build(p2)
data.frame(
  lable = glue::glue("qt_max_{seq(0.5, 0.95, 0.05)}"),
  level1 = getrootframe(bp$data[[2]], 1)$root,
  level2 = getrootframe(bp$data[[2]], 1.5)$root
) %>% arrange(desc(lable)) %>% 
  write_csv(file = "./fourtime/fourtime/fig2/age_smooth/NL_lm3.csv")





framelist_gender <- lapply(seq(0.5, 0.95, 0.05), function(qt){
  rawobj %>% group_by(gender_level, NL) %>% dplyr::summarise(
    "qt_max_dataht_{qt}" := quantile(max_data, {{ qt }})
  )
})
framesheet_gender <- purrr::reduce(framelist_gender, left_join)
maxframe_gender <- framesheet_gender %>% 
  dplyr::select(1, 2, contains("max")) %>% pivot_longer(-c("NL", "gender_level"))
p3 <- ggplot() +
  geom_point(rawobj, mapping = aes(x = NL, y = max_data)) +
  geom_smooth(maxframe_gender,
              mapping = aes(x = NL, y = value, color = name),
              linewidth = 1, se = F) + theme_base() +
  geom_hline(yintercept = c(1, 1.5), color = "orange", size = 1) +
  geom_vline(xintercept = c(30, 40), color = "orange", size = 1) +
  theme(legend.position = "top") +
  scale_color_viridis_d() + 
  facet_wrap(~ gender_level)
pdf("./fourtime/fourtime/fig2/age_smooth/NL_lm4.pdf", width = 7, height = 4)
print(p3)
dev.off()

bp <- ggplot_build(p3)
df <- getrootframe(bp$data[[2]], 1)
df$level2 <- getrootframe(bp$data[[2]], 1.5)$root
colnames(df) <- c("lable", "panel", "level1","level2")
lable_frame <- tibble(
  lable = 1:10,
  new_lable = glue::glue("qt_max_{seq(0.5, 0.95, 0.05)}")
)
df <- df %>% left_join(lable_frame, by = "lable")
df$panel <- df$panel %>% as.character() %>% as.numeric()
gender_frame <- tibble(
  panel = 1:2,
  gender = c("female", "male")
)
df <- df %>% left_join(gender_frame, by = "panel")
df %>% arrange(desc(lable)) %>% dplyr::select(gender, new_lable, level1, level2) %>% 
  write_csv(file = "./fourtime/fourtime/fig2/age_smooth/NL_lm4.csv")



max_level_frame <- rawobj %>% 
  group_by(gender_level, age_level, max_level) %>% 
  dplyr::summarise(
    max_level_number = n()
  ) %>% 
  dplyr::mutate(max_level_total_count = sum(max_level_number), 
                max_level_pct = max_level_number / max_level_total_count * 100,
                gender_level = as.character(gender_level))
max_level_frame$age_bin <- str_extract(max_level_frame$age_level, pattern = "\\d+") %>% as.numeric()
p4 <- ggplot(max_level_frame, aes(x = age_bin, y = max_level_pct, color = gender_level)) + 
  geom_smooth(se = F, span = 1) +
  geom_point(mapping = aes(size = max_level_number)) +
  facet_grid(~ max_level) +
  theme_base()
pdf("./fourtime/fourtime/fig2/age_smooth/NL_lm5.pdf", width = 10, height = 4)
print(p4)
dev.off()

max_level_frame0 <- rawobj %>% 
  mutate(max_level0 = str_extract(max_level, pattern = "level|normal")) %>%  
  group_by(gender_level, age_level, max_level0) %>% dplyr::summarise(
    max_level_number = n()
  ) %>% 
  dplyr::mutate(max_level_total_count = sum(max_level_number), 
                max_level_pct = max_level_number / max_level_total_count * 100,
                gender_level = as.character(gender_level))
max_level_frame0$age_bin <- str_extract(max_level_frame0$age_level, pattern = "\\d+") %>% as.numeric()
p5 <- ggplot(max_level_frame0, aes(x = age_bin, y = max_level_pct, color = gender_level)) + 
  geom_point(mapping = aes(size = max_level_number)) +
  geom_smooth(se = F, span = 0.8) +
  facet_grid(~ max_level0) +
  theme_base()
pdf("./fourtime/fourtime/fig2/age_smooth/NL_lm6.pdf", width = 10, height = 4)
print(p5)
dev.off()



