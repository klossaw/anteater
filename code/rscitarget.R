pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes",
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr","scMetabolism",
          "ggplot2","rsvd", "RcisTarget")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "anteater"
dataset <- "anteater"
species <- "mouse"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/scmetabolism")
setwd(workdir)



