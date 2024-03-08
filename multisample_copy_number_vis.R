### ---------------------------
## Purpose: plot all relative copy number data for a group of samples on a single linear ideogram
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
# Input file variables
spreadsheet_list <- "~/umn/data/metadata/Calbicans_snp-depth_paths.txt"
feature_file <- "ref_genome_files/Calbicans_SC5314_A21_plotting_features.txt"
label_file <-  "ref_genome_files/Calbicans_SC5314_A21_chr_labels.txt"
patient_data <- "~/umn/data/metadata/2022_Calbicans_sorted_patient_pop.csv"

# Load packages
library(readxl)
library(tidyverse)
library(ggplot2)
library(paletteer)

# Plotting variables
save_dir <-"images/Calbicans/"

ploidy <- 2

y_axis_labels <- c(1,2,3,4,5,6)  # manual y-axis labels, adjust as needed
inter_chr_spacing <- 150000 # size of space between chrs
patient_colors <- c(paletteer_d("colorBlindness::paletteMartin"))

ploidy_multiplier <- 3  # this multiplied by ploidy sets the max-y scale

chrom_outline_color <- "gray15"  # color of chromosome outlines

chrom_line_width <- 0.2  # line width of chromosome outlines

################################################################################
# X-axis labels overwrite input scaffold names in final plot
chr_ids <- scan(label_file, what = character())

# Read in patient metadata
pop.data <- read.table(patient_data,
                       sep = ",",
                       header = TRUE)

# Read in and combine depth data for all isolates
depth_files <- scan(spreadsheet_list, what=character())

genome_depth <- read_xlsx(depth_files[1]) %>%
  select(chr, index, pos, plot_pos, copy_number)
names(genome_depth)[names(genome_depth)=="copy_number"] <-str_extract(depth_files[1], "AMS[:digit:]+|MEC[:digit:]+")

for(i in 2:100){
  new_depth <- read_xlsx(depth_files[i]) %>%
    select(chr, index, pos, plot_pos, copy_number)
  names(new_depth)[names(new_depth)=="copy_number"] <-str_extract(depth_files[i], "AMS[:digit:]+|MEC[:digit:]+")

  genome_depth <- genome_depth %>%
    left_join(new_depth, by = join_by(chr,index,pos, plot_pos))
}

genome_depth <- rename(genome_depth, c("MEC103_2"="MEC103"))

genome_depth <- genome_depth %>%
  pivot_longer(names_to = "sample", values_to = "copy_number", cols=-c(chr,index,pos,plot_pos))

genome_depth <- genome_depth %>% inner_join(pop.data, by=join_by(sample))

genome_depth <- genome_depth%>%
  mutate(copy_number = ifelse(sample %in% c("MEC185", "MEC324"), copy_number +1, copy_number)) %>%
  arrange(ploidy_change)

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

################################################################################
# Plot linear genome
p <- genome_depth %>%
  group_by(sample) %>%
  ggplot() +
  geom_point(aes(x = plot_pos,
                 y = copy_number,
                 color=as.factor(ploidy_change)),alpha = 0.8, size=0.3, shape=20,
             show.legend = FALSE) +
  scale_color_manual(values=patient_colors) +
  geom_point(data = features, size = 2,
             aes(group=index, x=plot_start, y=ymin, shape = Feature),# fill = Feature),
             position = position_nudge(y=0.07)) +
  #scale_fill_manual(values = c("white", "grey26", "deepskyblue")) +
  scale_shape_manual(values = c(17,19,15), guide = guide_legend(title = NULL)) +
  ylab('Relative copy number') +
  geom_rect(data=chroms, aes(group=index, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            linewidth = chrom_line_width, fill = NA,
            colour = chrom_outline_color, linejoin = "round", inherit.aes = FALSE) +
  scale_x_continuous(name = NULL, expand = c(0, 0), breaks = ticks, labels=chr_ids) +
  scale_y_continuous(limits = c(0, ploidy*ploidy_multiplier), breaks = y_axis_labels) +
  theme_classic() +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        axis.ticks = element_line(color = NA),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size=12),
        legend.position = "bottom")

# Save plot
ggsave(paste0(save_dir,Sys.Date(),"_MEC_Calbicans_all_copy_number.png"),
       p, width = 12, height = 2.7, units = "in", device = png, dpi = 300, bg = "white")
