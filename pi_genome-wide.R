## ---------------------------
## Script name: pi_genome-wide.R
##
## Purpose: Plot calculated nucleotide diversity across chrs
##
## Author: Nancy Scott
##
## Date Created: 2024-02-21
##
## Email: scot0854@umn.edu
## ---------------------------
library(tidyverse)

feature_file <- "ref_genome_files/Cglabrata_cgd_s05m03r02_features.txt"
label_file <- "ref_genome_files/Cglabrata_chr_labels.txt"
mito <- "mito_C_glabrata_CBS138"  #"Ca19-mtDNA"
pi_text <- "/home/selmecki/shared/2022_Cglabrata_Erayil/align_variants_cgd_cbs138/bwa/pop_statistics/out.windowed.pi"

window <- "10kb" # size of calculation window in kb

# Chrom plotting and labeling
chrom_outline_color <- "gray15"
chr_ids <- scan(label_file, what = character())
chrom_line_width <- 0.2

# Saving
save_dir <- "images/Cglabrata/"# path with trailing slash or "" to save locally
species <- "Cglabrata"
ref <- "CBS138"

pi_file <- read_delim(pi_text, show_col_types = FALSE) %>%
  filter(CHROM != mito) %>%
  group_by(CHROM, index=consecutive_id(CHROM)) %>%
  group_by(index) %>%
  arrange(index)


pi_file$x_num <- seq.int(nrow(pi_file))

chrs <- pi_file %>%
  group_by(index) %>%
  summarise(border_start=min(x_num),
            border_stop=max(x_num),
            tick=min(x_num) + (max(x_num)-min(x_num))/2
  )

#features <- read_tsv(feature_file, show_col_types = FALSE)

#features <- features %>%
 # group_by(chr, index=consecutive_id(chr)) %>%
  #left_join(chroms, by=join_by(index))

#features <- features %>%
 # mutate(plot_start = start + xmin,
  #       plot_end = end + xmin)

p <- ggplot(pi_file,
       aes(x = x_num, y = PI)) +
  geom_point(alpha = 0.9, size = 0.5) +
  theme_minimal() +
  geom_rect(data=chrs,
            aes(group=index, xmin=border_start, xmax=border_stop, ymin=0, ymax=0.02),
            fill=NA,
            inherit.aes=FALSE,
            colour = "grey26",
            linejoin = "round") +
  scale_x_continuous(expand = c(0,0),
                    breaks = chrs$tick,
                   labels=chr_ids)+
  ylab("Pi, 10kb window") +
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title.y = element_text(size=14)
  )

ggsave(paste0(save_dir,Sys.Date(),"_",species,"_",ref,"_",window,".png"),
       p,
       device = png,
       dpi=300,
       bg = "white",
       width = 10,
       height = 3,
       units = "in")
