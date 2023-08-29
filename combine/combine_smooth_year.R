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

rawobj <- read_rds("./combine/updata_raw.rds") %>% dplyr::filter(NL <= 90)
core_frame <- read_rds("./combine/updata_core.rds") %>% dplyr::filter(NL <= 90)
framelist <- lapply(seq(0.5, 0.95, 0.05), function(qt){
  core_frame %>% group_by(NL) %>% dplyr::summarise(
    "qt_max_dataht_{qt}" := quantile(max_data, {{ qt }})
  )
})
framesheet <- purrr::reduce(framelist, left_join)
maxframe <- framesheet %>% dplyr::select(1, contains("max")) %>% pivot_longer(-"NL")
p1 <- ggplot(maxframe, aes(x = NL, y = value, color = name)) +
  geom_smooth(se = F) + 
  theme_base() +
  geom_hline(yintercept = c(1, 1.5)) +
  geom_vline(xintercept = c(30, 40)) +
  theme(legend.position = "top") +
  scale_color_viridis_d()
pdf("./combine/NL_lm1.pdf")
print(p1)
dev.off()

bp <- ggplot_build(p1)
data.frame(
  lable = glue::glue("qt_max_{seq(0.5, 0.95, 0.05)}"),
  level1 = getrootframe(bp$data[[1]], 1)$root,
  level2 = getrootframe(bp$data[[1]], 1.5)$root
) %>% arrange(desc(lable)) %>% write_csv(file = "./combine/NL_lm1.csv")

p2 <- ggplot() +
  geom_point(core_frame, mapping = aes(x = NL, y = max_data)) +
  geom_smooth(maxframe %>%
                dplyr::filter(name %in% c("qt_max_dataht_0.95", "qt_max_dataht_0.9", "qt_max_dataht_0.85", "qt_max_dataht_0.8", "qt_max_dataht_0.5")),
              mapping = aes(x = NL, y = value, color = name),
              linewidth = 1, se = F) + theme_base() +
  geom_hline(yintercept = c(1,1.5), color = "orange", size = 1) +
  geom_vline(xintercept = c(30, 40), color = "orange", size = 1) +
  theme(legend.position = "top")+
  scale_color_viridis_d()
pdf("./combine/NL_lm3.pdf")
print(p2)
dev.off()
bp <- ggplot_build(p2)
data.frame(
  lable = rev(c("qt_max_dataht_0.95", "qt_max_dataht_0.9", "qt_max_dataht_0.85", "qt_max_dataht_0.8", "qt_max_dataht_0.5")),
  level1 = getrootframe(bp$data[[2]], 1)$root,
  level2 = getrootframe(bp$data[[2]], 1.5)$root
) %>% arrange(desc(lable)) %>% write_csv(file = "./combine/NL_lm3.csv")

framelist_gender <- lapply(seq(0.5, 0.95, 0.05), function(qt){
  core_frame %>% group_by(xingbie_level, NL) %>% dplyr::summarise(
    "qt_max_dataht_{qt}" := quantile(max_data, {{ qt }})
  )
})
framesheet_gender <- purrr::reduce(framelist_gender, left_join)
maxframe_gender <- framesheet_gender %>% dplyr::select(1, 2, contains("max")) %>% pivot_longer(-c("NL", "xingbie_level"))
p3 <- ggplot() +
  geom_point(core_frame, mapping = aes(x = NL, y = max_data)) +
  geom_smooth(maxframe_gender %>%
                dplyr::filter(name %in% c("qt_max_dataht_0.95", "qt_max_dataht_0.9","qt_max_dataht_0.85", "qt_max_dataht_0.8","qt_max_dataht_0.5")),
              mapping = aes(x = NL, y = value, color = name),
              linewidth = 1, se = F) + theme_base() +
  geom_hline(yintercept = c(1, 1.5), color = "orange", size = 1) +
  geom_vline(xintercept = c(30, 40), color = "orange", size = 1) +
  theme(legend.position = "top") +
  scale_color_viridis_d() + 
  facet_wrap(~ xingbie_level)
pdf("./combine/NL_lm4.pdf", width = 7, height = 4)
print(p3)
dev.off()
bp <- ggplot_build(p3)
df <- getrootframe(bp$data[[2]], 1)
df$level2 <- getrootframe(bp$data[[2]], 1.5)$root
colnames(df)[1:3] <- c("lable", "panel", "level1")
df$lable <- rev(c("qt_max_dataht_0.95", "qt_max_dataht_0.95",
                  "qt_max_dataht_0.9", "qt_max_dataht_0.9",
                  "qt_max_dataht_0.85", "qt_max_dataht_0.85",
                  "qt_max_dataht_0.8", "qt_max_dataht_0.8",
                  "qt_max_dataht_0.5", "qt_max_dataht_0.5"))
df$panel <- c(0,1,0,1,0,1,0,1,0,1)
df %>% arrange(desc(lable)) %>% write_csv(file = "./combine/NL_lm4.csv")

max_level_frame <- core_frame %>% 
  group_by(xingbie_level, NL_level, max_level) %>% 
  dplyr::summarise(
    max_level_number = n()
    ) %>% 
  dplyr::mutate(max_level_total_count = sum(max_level_number), 
                max_level_pct = max_level_number / max_level_total_count * 100,
                xingbie_level = as.character(xingbie_level))
p4 <- ggplot(max_level_frame, aes(x = NL_level, y = max_level_pct, color = xingbie_level)) + 
  geom_smooth(se = F, span = 0.8) +
  geom_point(mapping = aes(size = max_level_number)) +
  facet_grid(~ max_level) +
  theme_base()
pdf("./combine/NL_lm5.pdf", width = 10, height = 4)
print(p4)
dev.off()
max_level_frame0 <- core_frame %>% mutate(max_level0 = str_extract(max_level, pattern = "level|normal")) %>%  
  group_by(xingbie_level, NL_level, max_level0) %>% dplyr::summarise(
  max_level_number = n()
) %>% 
  dplyr::mutate(max_level_total_count = sum(max_level_number), 
                max_level_pct = max_level_number / max_level_total_count * 100,
                xingbie_level = as.character(xingbie_level))
p5 <- ggplot(max_level_frame0, aes(x = NL_level, y = max_level_pct, color = xingbie_level)) + 
  geom_point(mapping = aes(size = max_level_number)) +
  geom_smooth(se = F) +
  facet_grid(~ max_level0) +
  theme_base()
pdf("./combine/NL_lm6.pdf", width = 10, height = 4)
print(p5)
dev.off()



