## ---------------------------
## Script name: gene_cov.R
##
## Purpose: Summarize deleted and truncated genes (determined by read depth)
##
## Author: Nancy Scott
##
## Date Created: 2024-02-15
##
## Email: scot0854@umn.edu
## ---------------------------
# Input vars
truncation_file <- "~/umn/data/metadata/Calbicans_MEC_all_truncations.tab"
deletion_file <- "~/umn/data/metadata/Calbicans_MEC_all_deletions.tab"
include_file <- "~/Calbicans_persvade/scripts/Calbicans_genes_to_check_coverage.txt"

## ---------------------------
# Load packages
library(readxl)
library(tidyverse)
library(ggplot2)
library(paletteer)
library(writexl)

## ---------------------------
truncations <- read.delim(truncation_file)

deletions <- read.delim(deletion_file)

genes_to_keep <- scan(include_file, what = character())

del_summary <- deletions %>%
  filter(ORF %in% genes_to_keep) %>%
  group_by(ORF) %>%
  summarize(total_samples=n())

truncate_summary <- truncations %>%
  filter(ORF %in% genes_to_keep) %>%
  group_by(ORF) %>%
  summarize(total_samples=n())

# Quick bar plots
del_plot <- ggplot(del_summary, aes(x=ORF, y=total_samples, fill = ORF)) +
  geom_col(colour = "grey26") +
  scale_fill_manual(values = paletteer_c("grDevices::Purple-Yellow", 11) ) +
  theme_bw() +
  ylab("Number of samples") +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")

ggsave("images/Calbicans/2022_Calbicans_ORF_deletions.png", del_plot, device = png, dpi = 300, bg="white")


trunc_plot <- ggplot(truncate_summary, aes(x=ORF, y=total_samples, fill = ORF)) +
  geom_col(colour = "grey26") +
  scale_fill_manual(values = paletteer_c("grDevices::Purple-Yellow", 60) ) +
  theme_bw() +
  ylab("Number of samples") +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")

ggsave("images/Calbicans/2022_Calbicans_ORF_truncations.png", trunc_plot, device = png, dpi = 300, bg="white")
