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
rawobj <- read_rds("./fixdata/updata_raw.rds")
core_frame <- read_rds("./fixdata/updata_core.rds")


framelist <- lapply(c(seq(0.1,0.9,0.1), 0.95), function(qt){
  core_frame %>% group_by(NL) %>% dplyr::summarise(
    mean_right = mean(right),
    mean_left = mean(left),
    mean_max_data = mean(max_data),
    "qt_right_{qt}" := quantile(right, {{ qt }}),
    "qt_left_{qt}" := quantile(left, {{ qt }}),
    "qt_max_dataht_{qt}" := quantile(max_data, {{ qt }})
  )
})

framesheet <- purrr::reduce(framelist, left_join)
rightframe <- framesheet %>% dplyr::select(1, contains("right")) %>% pivot_longer(-"NL")
leftframe <- framesheet %>% dplyr::select(1, contains("left")) %>% pivot_longer(-"NL")
maxframe <- framesheet %>% dplyr::select(1, contains("max")) %>% pivot_longer(-"NL")

pdf("./fix_first/NL_lm1.pdf")
ggplot(rightframe, aes(x = NL, y = value, color = name)) +
  geom_smooth(se = F) + theme_base() +
  geom_hline(yintercept = c(1,1.5)) +
  geom_vline(xintercept = c(30, 40)) +
  theme(legend.position = "top")+
  scale_color_viridis_d()
ggplot(leftframe, aes(x = NL, y = value, color = name)) +
  geom_smooth(se = F) + theme_base() +
  geom_hline(yintercept = c(1,1.5)) +
  geom_vline(xintercept = c(30, 40)) +
  theme(legend.position = "top")+
  scale_color_viridis_d()
ggplot(maxframe, aes(x = NL, y = value, color = name)) +
  geom_smooth(se = F) + theme_base() +
  geom_hline(yintercept = c(1,1.5)) +
  geom_vline(xintercept = c(30, 40)) +
  theme(legend.position = "top")+
  scale_color_viridis_d()
dev.off()
pdf("./fix_first/NL_lm2.pdf")
ggplot(rightframe, aes(x = NL, y = value, color = name)) +
  geom_point(size = 1) + geom_smooth(se = F) + theme_base() +
  geom_hline(yintercept = c(1,1.5)) + facet_wrap(~ name) +
  geom_vline(xintercept = c(30, 40)) +
  theme(legend.position = "top")+
  scale_color_viridis_d()
ggplot(leftframe, aes(x = NL, y = value, color = name)) +
  geom_point(size = 1) + geom_smooth(se = F) + theme_base() +
  geom_hline(yintercept = c(1,1.5)) + facet_wrap(~ name) +
  geom_vline(xintercept = c(30, 40)) +
  theme(legend.position = "top")+
  scale_color_viridis_d()
ggplot(maxframe, aes(x = NL, y = value, color = name)) +
  geom_point(size = 1) + geom_smooth(se = F) + theme_base() +
  geom_hline(yintercept = c(1,1.5)) + facet_wrap(~ name) +
  geom_vline(xintercept = c(30, 40)) +
  theme(legend.position = "top")+
  scale_color_viridis_d()
dev.off()
pdf("./fix_first/NL_lm3.pdf")
ggplot() +
  geom_point(core_frame, mapping = aes(x = NL, y = right)) +
  geom_smooth(rightframe %>%
                dplyr::filter(name %in% c("mean_right","qt_right_0.95","qt_right_0.8","qt_right_0.5")),
              mapping = aes(x = NL, y = value, color = name),
              linewidth = 1, se = F) + theme_base() +
  geom_hline(yintercept = c(1,1.5), color = "orange", size = 1) +
  geom_vline(xintercept = c(30, 40), color = "orange", size = 1)+
  theme(legend.position = "top") +
  scale_color_viridis_d()
ggplot() +
  geom_point(core_frame, mapping = aes(x = NL, y = left)) +
  geom_smooth(leftframe %>%
                dplyr::filter(name %in% c("mean_left","qt_left_0.95","qt_left_0.8","qt_left_0.5")),
              mapping = aes(x = NL, y = value, color = name),
              linewidth = 1, se = F) + theme_base() +
  geom_hline(yintercept = c(1,1.5), color = "orange", size = 1) +
  geom_vline(xintercept = c(30, 40), color = "orange", size = 1)+
  theme(legend.position = "top")+
  scale_color_viridis_d()
ggplot() +
  geom_point(core_frame, mapping = aes(x = NL, y = max_data)) +
  geom_smooth(maxframe %>%
                dplyr::filter(name %in% c("mean_max_data","qt_max_dataht_0.95","qt_max_dataht_0.8","qt_max_dataht_0.5")),
              mapping = aes(x = NL, y = value, color = name),
              linewidth = 1, se = F) + theme_base() +
  geom_hline(yintercept = c(1,1.5), color = "orange", size = 1) +
  geom_vline(xintercept = c(30, 40), color = "orange", size = 1) +
  theme(legend.position = "top")+
  scale_color_viridis_d()
dev.off()



framelist_gender <- lapply(c(seq(0.1,0.9,0.1), 0.95), function(qt){
  core_frame %>% group_by(xingbie_level, NL) %>% dplyr::summarise(
    mean_right = mean(right),
    mean_left = mean(left),
    mean_max_data = mean(max_data),
    "qt_right_{qt}" := quantile(right, {{ qt }}),
    "qt_left_{qt}" := quantile(left, {{ qt }}),
    "qt_max_dataht_{qt}" := quantile(max_data, {{ qt }})
  )
})

framesheet_gender <- purrr::reduce(framelist_gender, left_join)
rightframe_gender <- framesheet_gender %>% dplyr::select(1, 2, contains("right")) %>% pivot_longer(-c("NL", "xingbie_level"))
leftframe_gender <- framesheet_gender %>% dplyr::select(1, 2, contains("left")) %>% pivot_longer(-c("NL", "xingbie_level"))
maxframe_gender <- framesheet_gender %>% dplyr::select(1, 2, contains("max")) %>% pivot_longer(-c("NL", "xingbie_level"))

pdf("./fix_first/NL_lm4.pdf", width = 7, height = 4)
ggplot() +
  geom_point(core_frame, mapping = aes(x = NL, y = right)) +
  geom_smooth(rightframe_gender %>%
                dplyr::filter(name %in% c("mean_right","qt_right_0.95","qt_right_0.8","qt_right_0.5")),
              mapping = aes(x = NL, y = value, color = name),
              linewidth = 1, se = F) + theme_base() +
  geom_hline(yintercept = c(1,1.5), color = "orange", size = 1) +
  geom_vline(xintercept = c(30, 40), color = "orange", size = 1) +
  theme(legend.position = "top") +
  scale_color_viridis_d() + facet_wrap(~ xingbie_level)
ggplot() +
  geom_point(core_frame, mapping = aes(x = NL, y = left)) +
  geom_smooth(leftframe_gender %>%
                dplyr::filter(name %in% c("mean_left","qt_left_0.95","qt_left_0.8","qt_left_0.5")),
              mapping = aes(x = NL, y = value, color = name),
              linewidth = 1, se = F) + theme_base() +
  geom_hline(yintercept = c(1,1.5), color = "orange", size = 1) +
  geom_vline(xintercept = c(30, 40), color = "orange", size = 1)+
  theme(legend.position = "top")+
  scale_color_viridis_d() + facet_wrap(~ xingbie_level)
ggplot() +
  geom_point(core_frame, mapping = aes(x = NL, y = max_data)) +
  geom_smooth(maxframe_gender %>%
                dplyr::filter(name %in% c("mean_max_data","qt_max_dataht_0.95","qt_max_dataht_0.8","qt_max_dataht_0.5")),
              mapping = aes(x = NL, y = value, color = name),
              linewidth = 1, se = F) + theme_base() +
  geom_hline(yintercept = c(1,1.5), color = "orange", size = 1) +
  geom_vline(xintercept = c(30, 40), color = "orange", size = 1) +
  theme(legend.position = "top") +
  scale_color_viridis_d() + facet_wrap(~ xingbie_level)
dev.off()


max_level_frame <- core_frame %>% group_by(xingbie_level, NL_level, max_level) %>% dplyr::summarise(
  max_level_number = n()
) %>% 
  dplyr::mutate(max_level_total_count = sum(max_level_number), 
                max_level_pct = max_level_number / max_level_total_count * 100,
                xingbie_level = as.character(xingbie_level))
max_level_frame$NL_level[max_level_frame$NL_level > 7] <- 7

pdf("./fix_first/NL_lm5.pdf", width = 10, height = 4)
ggplot(max_level_frame, aes(x = NL_level, y = max_level_pct, color = xingbie_level)) + 
  geom_smooth(se = F) +
  geom_point(mapping = aes(size = max_level_number)) +
  facet_grid(~ max_level) +
  theme_base()
dev.off()

max_level_frame0 <- core_frame %>% mutate(max_level0 = str_extract(max_level, pattern = "level|normal")) %>%  
  group_by(xingbie_level, NL_level, max_level0) %>% dplyr::summarise(
  max_level_number = n()
) %>% 
  dplyr::mutate(max_level_total_count = sum(max_level_number), 
                max_level_pct = max_level_number / max_level_total_count * 100,
                xingbie_level = as.character(xingbie_level))
max_level_frame0$NL_level[max_level_frame0$NL_level > 7] <- 7

pdf("./fix_first/NL_lm6.pdf", width = 10, height = 4)
ggplot(max_level_frame0, aes(x = NL_level, y = max_level_pct, color = xingbie_level)) + 
  geom_point(mapping = aes(size = max_level_number)) +
  geom_smooth(se = F) +
  facet_grid(~ max_level0) +
  theme_base()
dev.off()



