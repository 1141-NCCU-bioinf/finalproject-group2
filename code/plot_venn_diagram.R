library(ggvenn)
library(readr)
library(dplyr)
library(tidyr)

file_name <- c(
  "DEG_table_Dsec_M_TW_vs_Dmel_M_TW.csv",
  "DEG_table_Dsec_M_JP_vs_Dsim_M_JP.csv",
  "DEG_table_Dsec_F_TW_vs_Dmel_F_TW.csv",
  "DEG_table_Dsec_F_JP_vs_Dsim_F_JP.csv"
)

for (i in 1:4){
  temp_file <- read_delim(file_name[i])
  temp_file <- temp_file |> select(Gene, M, prob)|> filter(prob >= .95)|> select(-prob)
  if (i == 1){
    dta <- temp_file
  }else{
    dta <- dta |> full_join(temp_file, by = "Gene")
  }
}

rm(temp_file)

colnames(dta) <- c("Gene", "Dsec_v._Dmel_M", "Dsec_v._Dsim_M", "Dsec_v._Dmel_F", "Dsec_v._Dsim_F")

dta <- dta |> mutate(
  Dsec_v._Dmel_M = Dsec_v._Dmel_M |> replace_na(0),
  Dsec_v._Dsim_M = Dsec_v._Dsim_M |> replace_na(0),
  Dsec_v._Dmel_F = Dsec_v._Dmel_F |> replace_na(0),
  Dsec_v._Dsim_F = Dsec_v._Dsim_F |> replace_na(0)
)

up_count <- dta > 0
up_count <- data.frame(up_count[,c(2:5)])

up <- ggvenn(up_count) 
ggsave("up.jpeg")

down_count <- dta < 0
down_count <- data.frame(down_count[,c(2:5)])
down <- ggvenn(down_count)
ggsave("down.jpeg")