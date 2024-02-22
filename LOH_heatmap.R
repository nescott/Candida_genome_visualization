## ---------------------------
## Script name: LOH_heatmap.R
##
## Purpose: Combine het SNP counts and generate heatmap across chrs
##
## Author: Nancy Scott
##
## Date Created: 2024-02-15
##
## Email: scot0854@umn.edu
## ---------------------------

# Input vars
spreadsheet_list <- "~/umn/data/metadata/Calbicans_snp-depth_paths.txt"
sample_list <- "~/umn/data/metadata/Calbicans_MEC_raxml_midpoint_tips.csv"
save_dir <- "~/umn/images/Calbicans/"

## ---------------------------
# Load packages
library(readxl)
library(tidyverse)
library(ggplot2)
library(paletteer)
library(writexl)

## ---------------------------
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

# Plotting: regular intervals for x-axis
genome_snp$x_num <- seq.int(nrow(genome_snp))

# Finish tidying data
snp_again <- genome_snp %>%
  pivot_longer(names_to = "sample", values_to = "snp_count", cols=-c(index,pos, x_num))

# Plotting: chr outlines and tick marks
chrs <- genome_snp %>%
  group_by(index) %>%
  summarise(border_start=min(x_num),
            border_stop=max(x_num),
            tick=min(x_num) + (max(x_num)-min(x_num))/2
            )

# Sort samples by tree order
sample_order <- read_csv(sample_list, show_col_types = FALSE)

# Plot
p <- snp_again %>%
  mutate(clustered_samples = fct_relevel(sample, rev(sample_order$sample))) %>%
    ggplot(aes(x=x_num, y=clustered_samples, fill=snp_count))+
    #scale_fill_gradient2(low = "white", mid="darkgoldenrod2", high = "dodgerblue3",
     #                    midpoint= 80, na.value = "white",
       #                  name="Heterozygous SNPs\n per 5kb window") +
    scale_fill_gradient(na.value = "white", high = "#800000", low = "white",
                      name="Heterozygous SNPs\n per 5kb window") +
    geom_tile() +
    theme_minimal() +
    geom_rect(data=chrs, aes(group=index, xmin=border_start, xmax=border_stop, ymin=0.5, ymax=Inf),
                             fill=NA, inherit.aes=FALSE,colour = "grey26", linejoin = "round") +
    scale_x_continuous(expand = c(0,0),
                       breaks = chrs$tick,
                       labels=c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "ChrR"))+
    theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=6.5),
        legend.title = element_text(size=10)
        )

# What are the axis positions per sample?
#details <- ggplot_build(p)

ggsave(paste0(save_dir,Sys.Date(),"_MEC_Calbicans_LOH_heatmap.png"), p, device=png, dpi=300, bg="white",
       width = 10, height = 7.5, units="in")

