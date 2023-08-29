#!/cluster/apps/anaconda3/2020.02/envs/R-4.2.1/bin/Rscript
load_libraries <- function() {
  pkgs <- c(
    "fs", "tidyverse", "futile.logger", "configr", "stringr", "optparse", "glue"
  )
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}
# dependencies
load_libraries()

project <- "healthman"
species <- "human"
dataset <- "liuzhong"

workdir <- glue("~/projects/{project}/analysis/")
rds_dir <- glue("{workdir}/{dataset}/{species}/clinical/rds")

df01 <- readRDS(glue("{rds_dir}/detail_2020.rds"))
df02 <- readRDS(glue("{rds_dir}/main_2020.rds"))
df11 <- readRDS(glue("{rds_dir}/detail_2021.rds"))
df12 <- readRDS(glue("{rds_dir}/main_2021.rds"))
df21 <- readRDS(glue("{rds_dir}/detail_2022.rds"))
df22 <- readRDS(glue("{rds_dir}/main_2022.rds"))
df31 <- readRDS(glue("{rds_dir}/detail_2023.rds"))
df32 <- readRDS(glue("{rds_dir}/main_2023.rds"))

write.csv(df01, file = glue("{workdir}/detail_2020.csv", row.names = FALSE))
write.csv(df02, file = glue("{workdir}/main_2020.csv", row.names = FALSE))
write.csv(df11, file = glue("{workdir}/detail_2021.csv", row.names = FALSE))
write.csv(df12, file = glue("{workdir}/main_2021.csv", row.names = FALSE))
write.csv(df21, file = glue("{workdir}/detail_2022.csv", row.names = FALSE))
write.csv(df22, file = glue("{workdir}/main_2022.csv", row.names = FALSE))
write.csv(df31, file = glue("{workdir}/detail_2023.csv", row.names = FALSE))
write.csv(df32, file = glue("{workdir}/main_2023.csv", row.names = FALSE))
