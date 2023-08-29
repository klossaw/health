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

boxsub3_core_frame <- core_frame[, c("smoking_level","NL_level","fuwei_level","bmi_level","LDL_level","xueya_level","xingbie_level",
                                     "left_level", "right_level", "max_level", "total_level","NL")] %>%
  na.omit()
v2 <- c("smoking_level","NL_level","fuwei_level","bmi_level","LDL_level","xueya_level","xingbie_level")
boxsub3_core_frame$right_level <- factor(boxsub3_core_frame$right_level, levels = c("normal", "level1", "level2" ))
boxsub3_core_frame$left_level <- factor(boxsub3_core_frame$left_level, levels = c("normal", "level1", "level2" ))
boxsub3_core_frame$max_level <- factor(boxsub3_core_frame$max_level, levels = c("normal", "level1", "level2" ))

i <- "smoking_level"

glmndf <- function(df, v){
  tmpl <- lapply(base::setdiff(v, c("xingbie_level","NL_level")), function(i){
    black <- tibble(
      NL_level = df$NL_level,
      xingbie_level = df$xingbie_level,
      x = df$max_level,
      y = df[[i]]
    ) %>% group_by(NL_level, xingbie_level, x, y) %>% 
      summarise(count = n()) %>% 
      ungroup() %>% group_by(NL_level, xingbie_level) %>% 
      summarise(lx = length(unique(x)), ly = length(unique(y))) %>% ungroup() %>% 
      dplyr::filter(lx == 1 | ly == 1)
    
    tmpdt <- df %>% 
      dplyr::select(NL_level, xingbie_level, max_level, i) %>% anti_join(black) %>% 
      dplyr::mutate(max_level = factor(str_extract(max_level, pattern = "normal|level"), levels = c("normal", "level"))) %>% 
      group_by(NL_level, xingbie_level) %>% 
      dplyr::summarise(sample_number = n(), glmn = list(glm(formula = glue("max_level ~ {i}") %>% as.formula() , data = pick("max_level", i), family = binomial()))) %>% 
      rowwise() %>% 
      dplyr::mutate(
        coef = list(glmn %>% summary() %>% .$coefficients %>% as.data.frame() %>% .[2,]),
      ) %>% unnest(coef) %>% dplyr::select(-"glmn")
    tmpdt
  })
  names(tmpl) <- base::setdiff(v, c("xingbie_level","NL_level"))
  tmpl <- tmpl %>% list_rbind(names_to = "feature")
}
glmndfx <- glmndf(boxsub3_core_frame, v2) %>% as.data.frame()



