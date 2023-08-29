pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "factoextra", "cluster","jhtools", "glue", "ggsci", "patchwork", 
          "tidyverse", "dplyr", "ComplexHeatmap", "ggpointdensity", "ggpubr",
          "tidygraph", "ggraph")
for (pkg in pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "healthman"
dataset <- "zhanglei"
species <- "human"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/") %>% checkdir()
setwd(workdir)
rawobj <- read_rds("./fourtime/updata_raw.rds")

nolevel1 <- c(
  "right", "left", "max_data",
  "gender_level","age_level", "age", # 性别 年龄
  "drinking","smoking", # 抽烟 喝酒
  "systolic_blood_pressure","diastolic_pressure", # 收缩压 舒张压
  "heart_rate","height","weight", # 心率 升高 体重
  "abdominal_circumference","bmi", # 腹围 BMI
  "white_blood_cell_count","neutrophils","neutrophil_percentage", # 血
  "lymphocyte_percentage","monocyte_percentage","hemoglobin","platelet_count","platelet_volume",
  "urinary_microalbumin", "urine_microalbuminuria_creatinine_ratio", "total_protein", "albumin", #尿微量蛋白  白蛋白 总蛋白 
  "alanine_aminotransferase", # 谷丙转氨酶
  "aspartate_aminotransferase", # 谷草转氨酶
  "creatinine", # 肌酐
  "urea","uric_acid", # 尿素 尿酸
  "glomerular_filtration_rate","homocysteine", # 肾小球滤过率 同型半胱氨酸
  "triglycerides","fasting_blood_sugar", # 甘油三酯  空腹血糖
  "apolipoproteinA1","apolipoproteinB", # 载脂蛋白A1 载脂蛋白B
  "apolipoproteinE","glycated_hemoglobinA1", # 载脂蛋白E 糖化血红蛋白A1
  "glycated_hemoglobinA1C","fasting_C_peptide", # 糖化血红蛋白A1C  空腹C肽
  "fasting_insulin","hsc_reactive_protein", # 空腹胰岛素  超敏C反应蛋白
  "sialic_acid","free_fatty_acid",# 唾液酸 游离脂肪酸
  "Hydroxyvitamin_D3", # 羟基维生素D3 
  "TC","HDL","LDL","NHDL", # TC HDL LDL NHDL
  "NHDL_age", "LDL_age"  # NHDL*age LDL*age
)
tw2cor <- rawobj[,nolevel1[-c(1,2,5,7,8)]] %>% as.data.frame()

cormf <- function(df, group_col = "gender_level"){
  lbs <- unique(df[[group_col]])
  bingl <- list()
  for(k in lbs){
    gft <- df[[group_col]] == k
    subdf <- df %>% dplyr::select(-group_col)
    n <- 0
    clmale <- list()
    for(i in 1:ncol(subdf)){
      for(j in 1:ncol(subdf)){
        n <- n + 1
        cdf <- subdf[gft,c(i, j)] %>% na.omit()
        if(nrow(cdf) <= 2){
          clmale[[n]] <- tibble( i = colnames(subdf)[i], j = colnames(subdf)[j], r = NA, p = NA)
        }else{
          cdf[[1]] <- as.numeric(cdf[[1]])
          cdf[[2]] <- as.numeric(cdf[[2]])
          tmp <- cor.test(cdf[[1]], cdf[[2]], method = "spearman")
          clmale[[n]] <- tibble( i = colnames(subdf)[i], j = colnames(subdf)[j], 
                                 r = tmp$estimate, p = tmp$p.value)
        }
      }
    }
    tmp <- list_rbind(clmale)
    tmp$padj <- p.adjust(tmp$p, method = "BH")
    colnames(tmp)[c(3, 4, 5)] <- glue::glue("{k}_{colnames(tmp)[c(3, 4, 5)]}")
    bingl[[k]] <- tmp
  }
  bingl
}
joinlist <- function(cl, main_f = "max_data", fig = T, filter_sig = 0.05){
  allframe <- purrr::reduce(cl, left_join, by = c("i", "j")) %>% 
    dplyr::filter(i == {{main_f}} & i != j) 
  Rf <- allframe %>%
    dplyr::select(tidyselect::matches("_r$|^j$")) %>% 
    pivot_longer(cols = -"j", values_to = "Rvalue") %>% 
    mutate(name = str_remove(name, pattern = "_.+"))
  pf <- allframe %>%
    dplyr::select(tidyselect::matches("_padj$|^j$")) %>% 
    pivot_longer(cols = -"j", values_to = "padjvalue")%>% 
    mutate(name = str_remove(name, pattern = "_.+"))
  f <- Rf %>% left_join(pf, by = c("j", "name"))
  
  f$padjvalue[f$padjvalue == 0] <- min(f$padjvalue[f$padjvalue != 0]) * 0.001
  if(fig){
    if(is.null(filter_sig)){
      f <- f %>% 
        ggplot(aes(x = j, y = Rvalue, color = name, size = -log(padjvalue))) + 
        geom_point(alpha = 0.8)+ coord_flip() + theme_bw()
    }else{
      f <- f %>% dplyr::filter(padjvalue < {{filter_sig}}) %>% 
        ggplot(aes(x = j, y = Rvalue, color = name, size = -log(padjvalue))) + 
        geom_point(alpha = 0.8)+ coord_flip() + theme_bw()
    }
  }
  return(f)
}
draw_matrix <- function(cl, pth = "./fourtime/fourtime/fig3/cor_matrix/", filter_sig = 0.05, clustermatrix = T){
  fns <- names(cl)
  for(fn in fns){
    cn <- grep(colnames(cl[[fn]]), pattern = "_padj", value = T)
    if(is.null(filter_sig)){
      clmalem <- cl[[fn]] %>% 
        pivot_wider(id_cols = 'i', names_from = "j", values_from = tidyselect::matches("_r")) %>% 
        column_to_rownames("i") %>% as.matrix()
    }else{
      clmalem <- cl[[fn]] %>% dplyr::filter(!!sym(cn) < {{filter_sig}}) %>% 
        pivot_wider(id_cols = 'i', names_from = "j", values_from = tidyselect::matches("_r")) %>% 
        column_to_rownames("i") %>% as.matrix()
    }
    pdf(glue("{pth}/heatmap_{fn}.pdf"), width = 10, height = 10)
    print(Heatmap(clmalem, cluster_columns = clustermatrix, cluster_rows = clustermatrix))
    dev.off()
  }
}
draw_graph <- function(cl, pth = "./fourtime/fourtime/fig3/cor_graph/", filter_sig = 0.05){
  fns <- names(cl)
  for(fn in fns){
    cn <- grep(colnames(cl[[fn]]), pattern = "_padj", value = T)
    if(is.null(filter_sig)){
      clmalem <- cl[[fn]] 
    }else{
      clmalem <- cl[[fn]] %>% dplyr::filter(!!sym(cn) < {{filter_sig}} & i != j) 
    }
    nodes <- tibble(node = unique(c(clmalem$i, clmalem$j)))
    net <- tbl_graph(edges = clmalem, nodes  = nodes)
    cnr <- grep(colnames(cl[[fn]]), pattern = "_r", value = T)
    fg <- ggraph(net)+
      geom_edge_link2(aes(color = .data[[cnr]])) + 
      scale_edge_color_gradient2() +
      geom_node_point() +
      geom_node_text(aes(label = node)) 
    
    pdf(glue("{pth}/graph_{fn}.pdf"), width = 10, height = 10)
    print(fg)
    dev.off()
  }
}

dir.create("./fourtime/fourtime/fig3/cor_matrix/", recursive = T)
tw2cor <- rawobj[,nolevel1[-c(1,2,4,7,8)]] %>% as.data.frame()
mmm <- cormf(tw2cor %>% dplyr::filter(age_level  != ">= 70"), group_col = "age_level")
p <- joinlist(cl = mmm, main_f = "max_data",fig = T,filter_sig = 0.05)
pdf("./fourtime/fourtime/fig3/cor_matrix/point_age_level.pdf")
print(p)
dev.off()
draw_matrix(mmm, filter_sig = NULL)
dir.create("./fourtime/fourtime/fig3/cor_graph/")
draw_graph(mmm, filter_sig = 0.05)

tw2cor <- rawobj[,nolevel1[-c(1,2,5,7,8)]] %>% as.data.frame()
mmm <- cormf(tw2cor, group_col = "gender_level")
p <- joinlist(cl = mmm, main_f = "max_data",fig = T,filter_sig = NULL)
pdf("./fourtime/fourtime/fig3/cor_matrix/point_gender_level.pdf")
print(p)
dev.off()

draw_matrix(mmm, filter_sig = NULL)
draw_graph(mmm, filter_sig = 0.05)






