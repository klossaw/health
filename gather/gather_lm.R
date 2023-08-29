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
core_frame <- read_csv("updata_core.csv") %>% dplyr::select(NL, right, left, max_data)

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

pdf("./secondtime/NL_lm1.pdf")
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

pdf("./secondtime/NL_lm2.pdf")
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

pdf("./secondtime/NL_lm3.pdf")
ggplot() +
  geom_point(core_frame, mapping = aes(x = NL, y = right)) +
  geom_smooth(rightframe %>%
                dplyr::filter(name %in% c("mean_right","qt_right_0.95")),
              mapping = aes(x = NL, y = value, color = name),
              linewidth = 1, se = F) + theme_base() +
  geom_hline(yintercept = c(1,1.5), color = "orange", size = 1) +
  geom_vline(xintercept = c(30, 40), color = "orange", size = 1)+
  theme(legend.position = "top") +
  scale_color_viridis_d()

ggplot() +
  geom_point(core_frame, mapping = aes(x = NL, y = left)) +
  geom_smooth(leftframe %>%
                dplyr::filter(name %in% c("mean_left","qt_left_0.95")),
              mapping = aes(x = NL, y = value, color = name),
              linewidth = 1, se = F) + theme_base() +
  geom_hline(yintercept = c(1,1.5), color = "orange", size = 1) +
  geom_vline(xintercept = c(30, 40), color = "orange", size = 1)+
  theme(legend.position = "top")+
  scale_color_viridis_d()

ggplot() +
  geom_point(core_frame, mapping = aes(x = NL, y = max_data)) +
  geom_smooth(maxframe %>%
                dplyr::filter(name %in% c("mean_max_data","qt_max_dataht_0.95")),
              mapping = aes(x = NL, y = value, color = name),
              linewidth = 1, se = F) + theme_base() +
  geom_hline(yintercept = c(1,1.5), color = "orange", size = 1) +
  geom_vline(xintercept = c(30, 40), color = "orange", size = 1) +
  theme(legend.position = "top")+
  scale_color_viridis_d()
dev.off()


















