## ---------------------------
## Purpose: Calculate binned homozygous SNP counts, merge with het SNP file
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
# Input files
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

snp_file <- ""
summary_file <- ""

# Load packages
library(tidyverse)
library(ggplot2)
library(writexl)
library(readxl)

# Set variables
window <- 5000 # size of window used for rolling mean and snp density
mito <- "Ca19-mtDNA"
species <- "Calbicans"
save_dir <- "~/umn/data/genome_plots/"

sample_id <- str_extract(snp_file, "AMS[:digit:]+|MEC[:digit:]+")
################################################################################
# Homozygous SNP freq calcs
genome_snp <- read.table(snp_file, header = TRUE)
genome_snp <- genome_snp %>%
  filter(chr != mito) %>%
  mutate(reads=rowSums(pick(A,T,G,C))) %>%  # sum read count per row
  mutate(across(c(A,T,G,C), ~ .x /reads, .names = "{.col}_freq")) # get allele freq

genome_snp <- genome_snp %>%
  filter(reads >= 10) %>%
  filter(!if_all(ends_with("_freq"), is.na)) %>%
  mutate(snp_bin=(pos %/% window) * window +1)

gaf <- genome_snp %>%
  mutate(non_A = case_when((ref=="A" & (T_freq >= 0.95 | C_freq >=0.95 | G_freq >= 0.95)) ~ 1, .default = 0),
         non_T = case_when((ref=="T" & (A_freq >= 0.95 | C_freq >=0.95 | G_freq >= 0.95)) ~ 1, .default = 0),
         non_G = case_when((ref=="G" & (A_freq >= 0.95 | C_freq >=0.95 | T_freq >= 0.95)) ~ 1, .default = 0),
         non_C = case_when((ref=="C" & (T_freq >= 0.95 | A_freq >=0.95 | G_freq >= 0.95)) ~ 1, .default = 0))

gaf <- gaf %>%
  group_by(chr, snp_bin) %>%
  summarise(hom_count=sum(non_A, non_T, non_G, non_C)) %>%
  mutate(concat_pos=paste(chr,snp_bin) )

het_freqs <- read_xlsx(summary_file) %>%
  mutate(concat_pos=paste(chr,pos))

het_freqs <- het_freqs %>%
  left_join(gaf, by=join_by(concat_pos)) %>%
  mutate(all_snps=rowSums(pick(snp_count, hom_count), na.rm=TRUE)) %>%
  select(chr.x, pos, depth, rolling_mean, index, relative_depth, copy_number, chr_sums, plot_pos, snp_count, hom_count, all_snps)

het_freqs <- het_freqs %>%
  rename(chr=chr.x, het_count=snp_count)
################################################################################
# Save dataframes as excel file

write_xlsx(het_freqs, path = sprintf("%s%s/%s_%s_all_snps.xlsx", save_dir, species, Sys.Date(), sample_id))
