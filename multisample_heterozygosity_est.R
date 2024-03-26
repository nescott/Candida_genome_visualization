## ---------------------------
## Purpose: Calculate genome heterozygosity % and plot distribution
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999)

# Input vars
spreadsheet_list <- "~/umn/data/metadata/Calbicans_snp-depth_paths.txt"
sample_list <- "~/umn/data/metadata/Calbicans_MEC_raxml_midpoint_tips.csv"
save_dir <- "~/umn/images/Calbicans/"
genome_size <- 14324315

# Load packages
library(readxl)
library(tidyverse)
library(ggplot2)
library(writexl)

# Read in SNP counts for all samples
snp_files <- scan(spreadsheet_list, what=character())

genome_snp <- read_xlsx(snp_files[1]) %>%
  select(index, pos, snp_count)
names(genome_snp)[names(genome_snp)=="snp_count"] <-str_extract(snp_files[1], "AMS[:digit:]+|MEC[:digit:]+")

for(i in 2:100){
  new_snp <- read_xlsx(snp_files[i]) %>%
    select(index, pos, snp_count)
  names(new_snp)[names(new_snp)=="snp_count"] <-str_extract(snp_files[i], "AMS[:digit:]+|MEC[:digit:]+")

  genome_snp <- genome_snp %>%
    left_join(new_snp, by = join_by(index,pos))
}

snp_again <- genome_snp %>%
  pivot_longer(names_to = "sample", values_to = "snp_count", cols=-c(index,pos))

snp_total <- snp_again %>% group_by(sample) %>%
  summarize(all_snps = sum(snp_count, na.rm = TRUE)) %>%
  mutate(het_percentage = all_snps/genome_size*100)

h <- ggplot(snp_total, aes(x = het_percentage)) +
  geom_histogram(bins=50, fill = "lightgrey", color = "black") +
  theme_bw() +
  xlab("Genome heterozygosity, %") +
  ylab("Count")

ggsave(paste0(save_dir,"Calbicans_genome_heterozygosity_estimate.png"),
       h,
       device = png,
       dpi = 300,
       bg = "white")
