## ---------------------------
## Purpose: Replot genome-view of LOH and CNV from saved excel file (see "genome_vis.R)
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# Input file variables
genome_df_file <- args[1]
sample_id <- args[2] # or "YourID"
feature_file <- args[3] # or "path/to/features.txt"
label_file <- args[4] # or "path/to/chr_labels.txt"

# Load packages
library(readxl)
library(tidyverse)
library(ggplot2)

# Set variables
window <- 5000 # size of window used for rolling mean and snp density
ploidy <- 2

y_axis_labels <- c(1,2)  # manual y-axis labels, adjust as needed
inter_chr_spacing <- 150000 # size of space between chrs

save_dir <- "images/Calbicans/genome_plots/" # path with trailing slash, or  "" to save locally
ref <- "SC5314" # short label for file name or "" to leave out

snp_low <- "white"  # snp LOH colors, plot function uses 2-color gradient scale
snp_high <- "black"  # snp LOH colors, plot function uses 2-color gradient scale
cnv_color <- "dodgerblue4"  # copy number color

feature_colors <- c("white", "grey26", "deepskyblue")
feature_shapes <- c(24,21,22)

ploidy_multiplier <- 2  # this multiplied by ploidy sets the max-y scale

chrom_outline_color <- "gray15"  # color of chromosome outlines

chrom_line_width <- 0.2  # line width of chromosome outlines

# X-axis labels overwrite input scaffold names in final plot
chr_ids <- scan(label_file, what = character())

# Dataframe of joined copy number, snps, and plotting positions per window
genome_depth <- read_xlsx(genome_df_file, sheet=1)

# Small dataframes for chrom. outlines and features
chroms <- genome_depth %>%
  group_by(index) %>%
  summarise(xmin=min(plot_pos), xmax=max(plot_pos), ymin=0, ymax=Inf)

features <- read_tsv(feature_file, show_col_types = FALSE)

features <- features %>%
     group_by(chr, index=consecutive_id(chr)) %>%
     left_join(chroms, by=join_by(index))

features <- features %>%
     mutate(plot_start = start + xmin, plot_end = end + xmin)

# Tick marks to center chromosome ID label
ticks <- tapply(genome_depth$plot_pos, genome_depth$index, quantile, probs =
                0.5, na.remove = TRUE)

# Plot linear genome
p <- ggplot(genome_depth) +
  scale_color_gradient(low=snp_low,high=snp_high, na.value = "white", guide = "none") +
  geom_segment(aes(x = plot_pos, y = 0, color = snp_count, xend = plot_pos, yend = Inf)) +
  geom_segment(aes(x = plot_pos,
                   y = ifelse(copy_number <= ploidy*ploidy_multiplier, copy_number, Inf),
                   xend = plot_pos, yend = ploidy), alpha = 0.9, color = cnv_color) +
  geom_rect(data=chroms, aes(group=index, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            linewidth = chrom_line_width, fill = NA,
            colour = chrom_outline_color, linejoin = "round", inherit.aes = FALSE) +
  geom_point(data = features, size = 2,
               aes(group=index, x=plot_start, y=ymin, shape = Feature, fill = Feature),
               position = position_nudge(y=0.07)) +
    scale_fill_manual(values = feature_colors) +
    scale_shape_manual(values = feature_shapes) +
  ylab(sample_id) +
  scale_x_continuous(name = NULL, expand = c(0, 0), breaks = ticks, labels=chr_ids) +
  scale_y_continuous(limits = c(0, ploidy*ploidy_multiplier), breaks = y_axis_labels) +
  theme_classic() +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        axis.ticks = element_line(color = NA),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size=11))

# Save plot
ggsave(sprintf("%s%s_%s_%s_%sbp.png", save_dir, Sys.Date(), sample_id, ref, window),
       p, width = 18, height = 1.7, units = "in", device = png, dpi = 300, bg = "white")
