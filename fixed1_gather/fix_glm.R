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


core_frame$group <- ifelse(str_detect(core_frame$max_level, pattern = "level"), 1, 0)

tmpobj <- glm(group ~ LDL, family = binomial(link = "logit"), data = core_frame)
summary(tmpobj)
confint(tmpobj)
coef(tmpobj)




