# script to visualize heterozygosity between two haplotype assemblies 

library(VariantAnnotation)
library(GenomicRanges)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

#Read VCF
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
vcf <- readVcf("../syri_snps_only.vcf", genome = "unknown")
gr <- rowRanges(vcf)

# Calculate chromosome lengths and cumulative starts
fai <- read.delim("../mat_filtered.fasta.fai",header=F)[,c(1,2)]
colnames(fai)[1] <- "chr_name"
colnames(fai)[2] <- "length"
# set length
chrom_lengths <- setNames(fai$length, fai$chr_name)
common_chroms <- intersect(seqlevels(gr), names(chrom_lengths))
seqlengths(gr)[common_chroms] <- chrom_lengths[common_chroms]

# sliding window
window_size <- 1000000   # 1 Mb
step_size <- 100000      # 100 kb

# Function to create sliding windows for a chromosome
makeSlidingWindows <- function(chr, chr_len, win_size, step) {
  starts <- seq(1, chr_len - win_size + 1, by = step)
  ends <- starts + win_size - 1
  GRanges(seqnames = chr, ranges = IRanges(start = starts, end = ends))
}

# Apply to all chromosomes in your GRanges object
sliding_windows <- do.call(
  c,
  lapply(names(seqlengths(gr)), function(chr) {
    chr_len <- seqlengths(gr)[chr]
    if (!is.na(chr_len) && chr_len >= window_size) {
      makeSlidingWindows(chr, chr_len, window_size, step_size)
    } else {
      NULL
    }
  })
)

# count SNPs
counts <- countOverlaps(sliding_windows, gr)
df <- data.frame(
  chr = as.character(seqnames(sliding_windows)),
  start = start(sliding_windows),
  end = end(sliding_windows),
  n_het = counts
)

# plot all chromosomes in one x axis
chr_lengths <- seqlengths(gr)[unique(df$chr)]
chr_lengths <- chr_lengths[!is.na(chr_lengths)]


# Compute cumulative start positions for each chromosome
library(gtools)
chr_order <- mixedsort(names(chr_lengths))
# Reorder chr_lengths
chr_lengths_ordered <- chr_lengths[chr_order]
chr_offsets <- cumsum(c(0, head(chr_lengths_ordered, -1)))
names(chr_offsets) <- chr_order

# Add cumulative position to df
df$cum_start <- df$start + chr_offsets[df$chr]
# Compute chromosome midpoints for x-axis labels
label_df <- data.frame(
  chr = chr_order,
  midpoint = chr_offsets + chr_lengths_ordered / 2
)


# make two colors
df$chr_color <- as.factor(as.integer(factor(df$chr)) %% 2)


# add centromere position
centro <- read.delim("../../data/bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.gff",header = F)
centro <- centro %>% filter(str_detect(V1,"mat")) %>% select(V1,V4,V5)
colnames(centro) <- c("chr","start","end")
centro$chr <- gsub("_mat","",centro$chr)
centro <- centro %>% filter(!chr %in% "chrW")
centro$pos <- (centro$start + centro$end )/2
# need to calculate the cumulative start
centro$cum_pos <- centro$pos + chr_offsets[centro$chr]


# Flip the chromosomes so that centromeres are on the left of each chromosome
flip_chr <- c("chr1A","chr3","chr4","chr4A","chr8","chr11","chr12","chr13","chr14","chr19","chr21","chr23",
              "chr26","chr27","chr28","chr30","chr31","chr32","chr33","chr35","chr36","chr37")
# flip window position
df$cum_start_flip <- NA
for (i in 1:length(df$chr)) {
  if (!df$chr[i] %in% flip_chr){
    df$cum_start_flip[i] <- df$cum_start[i]
  } else {
    # abs(length-start)+1+offset
    df$cum_start_flip[i] <- abs(chr_lengths[df$chr[i]] - df$end[i])+1+chr_offsets[df$chr[i]]
  }
}
# flip centromere
centro$cum_pos_flip <- NA
for (i in 1:length(centro$chr)) {
  if (!centro$chr[i] %in% flip_chr){
    centro$cum_pos_flip[i] <- centro$cum_pos[i]
  } else {
    # abs(length-start)+1+offset
    centro$cum_pos_flip[i] <- abs(chr_lengths[centro$chr[i]] - centro$pos[i])+1+chr_offsets[centro$chr[i]]
  }
}








#-------------------
# plot actual chr size
#png("minimap2_syri_SNP.png", units = "in", width = 9, height = 5, res=300)
ggplot(df, aes(x = cum_start_flip / 1e6, y = n_het/1000, fill=chr_color)) +
  geom_col(width = 1, position = "identity") +
  # Custom y-axis line
  geom_segment(x = 0, xend = 0, y = 0, yend = 25, color = "black", linewidth = 0.5) +
  # Custom x-axis line
  geom_segment(x = 0, xend = max(df$cum_start) / 1e6, y = 0, yend = 0, color = "black", linewidth = 0.5) +
  # add centromere
  annotate("point", x = centro$cum_pos_flip/1e6, y = 0, shape = 24, size =2, fill = "red", color = "black")+
  scale_fill_manual(values = c("lightblue","steelblue")) +
  scale_x_continuous(
    breaks = label_df$midpoint / 1e6,
    labels = label_df$chr,
    #limits=c(0, max(df$global_start)/1e6),
    expand = c(0, 0)
  ) +
  ylim(0, 25)+
  labs(y = "Het / 1kbp") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 0.5), # hjust = 0.5
    panel.grid = element_blank(),
    axis.line = element_blank(),
    legend.position = "none"
  )
#dev.off()




#-------------------
# plot actual chr size, split into two rows
df_up <- df %>% filter(chr %in% c("chr1","chr1A","chr2","chr3","chr4"))
centro_up <- centro %>% filter(chr %in% c("chr1","chr1A","chr2","chr3","chr4"))
df_down <- df %>% filter(!chr %in% c("chr1","chr1A","chr2","chr3","chr4"))
centro_down <- centro %>% filter(!chr %in% c("chr1","chr1A","chr2","chr3","chr4"))

# re-calculate the mid-point
up_order <- mixedsort(unique(df_up$chr))
up_lengths_ordered <- chr_lengths[up_order]
up_offsets <- cumsum(c(0, head(up_lengths_ordered, -1)))
names(up_offsets) <- up_order
label_df_up <- data.frame(
  chr = up_order,
  midpoint = up_offsets + up_lengths_ordered / 2
)
label_df_up$chr <- sub("^chr", "", label_df_up$chr)
# re-calculate start point
df_up$cum_start <- df_up$start + up_offsets[df_up$chr]
df_up$chr_color <- as.factor(as.integer(factor(df_up$chr)) %% 2)
# make 0 visible
df_up$plot_het <- ifelse(df_up$n_het < 300, 300, df_up$n_het)
# flip window position
df_up$cum_start_flip <- NA
for (i in 1:length(df_up$chr)) {
  if (!df_up$chr[i] %in% flip_chr){
    df_up$cum_start_flip[i] <- df_up$cum_start[i]
  } else {
    # abs(length-start)+1+offset
    df_up$cum_start_flip[i] <- abs(chr_lengths[df_up$chr[i]] - df_up$end[i])+1+up_offsets[df_up$chr[i]]
  }
}
# flip centromere
centro_up$cum_pos_flip <- NA
for (i in 1:length(centro_up$chr)) {
  if (!centro_up$chr[i] %in% flip_chr){
    centro_up$cum_pos_flip[i] <- centro_up$cum_pos[i]
  } else {
    # abs(length-start)+1+offset
    centro_up$cum_pos_flip[i] <- abs(chr_lengths[centro_up$chr[i]] - centro_up$pos[i])+1+up_offsets[centro_up$chr[i]]
  }
}
# plot
(p1 <- ggplot(df_up, aes(x = cum_start_flip / 1e6, y = n_het/1000, fill=chr_color)) +
    geom_col(width = 1, position = "identity") +
    scale_fill_manual(values = c("lightblue","steelblue")) +
    scale_x_continuous(breaks = label_df_up$midpoint / 1e6,labels = label_df_up$chr,expand = c(0, 0),
                       limits= c(-3, max(df_up$cum_start_flip) / 1e6)) +
    scale_y_continuous(expand = c(0, 0), limits = c(-1, 22))+
    geom_hline(yintercept = 0, color = "black", linewidth = 1)+
    geom_vline(xintercept = 0, color = "black", linewidth = 1)+
    labs(y = "") +
    annotate("point", x = centro_up$cum_pos_flip/1e6, y = 0, shape = 24, size =4, fill = "red", color = "black")+
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          #axis.text.x = element_text(angle = 45, vjust = 0.5), # hjust = 0.5
          panel.grid = element_blank(),
          axis.line = element_blank(),
          legend.position = "none",
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          base_family = "Arial"
          )
)




# re-calculate the mid-point
down_order <- mixedsort(unique(df_down$chr))
down_lengths_ordered <- chr_lengths[down_order]
down_offsets <- cumsum(c(0, head(down_lengths_ordered, -1)))
names(down_offsets) <- down_order
label_df_down <- data.frame(
  chr = down_order,
  midpoint = down_offsets + down_lengths_ordered / 2
)
label_df_down$chr <- sub("^chr", "", label_df_down$chr)
# re-calculate start point
df_down$cum_start <- df_down$start + down_offsets[df_down$chr]
df_down$chr_color <- as.factor(as.integer(factor(df_down$chr)) %% 2)
centro_down$cum_pos <- centro_down$pos + down_offsets[centro_down$chr]
# make 0 visible
#df_down$plot_het <- ifelse(df_down$n_het < 300, 300, df_down$n_het)
# flip window position
df_down$cum_start_flip <- NA
for (i in 1:length(df_down$chr)) {
  if (!df_down$chr[i] %in% flip_chr){
    df_down$cum_start_flip[i] <- df_down$cum_start[i]
  } else {
    # abs(length-start)+1+offset
    df_down$cum_start_flip[i] <- abs(chr_lengths[df_down$chr[i]] - df_down$end[i])+1+down_offsets[df_down$chr[i]]
  }
}
# flip centromere
centro_down$cum_pos_flip <- NA
for (i in 1:length(centro_down$chr)) {
  if (!centro_down$chr[i] %in% flip_chr){
    centro_down$cum_pos_flip[i] <- centro_down$cum_pos[i]
  } else {
    # abs(length-start)+1+offset
    centro_down$cum_pos_flip[i] <- abs(chr_lengths[centro_down$chr[i]] - centro_down$pos[i])+1+down_offsets[centro_down$chr[i]]
  }
}
# to make the x axis text prettier
#label_df_down$chr[which(label_df_down$chr %in% c("19"))] <- "\n19"
label_df_down$chr[which(label_df_down$chr %in% c("21"))] <- "\n21"
label_df_down$chr[which(label_df_down$chr %in% c("23"))] <- "\n23"
label_df_down$chr[which(label_df_down$chr %in% c("25"))] <- "\n25"
label_df_down$chr[which(label_df_down$chr %in% c("27"))] <- "\n27"
label_df_down$chr[which(label_df_down$chr %in% c("29"))] <- "\n29"
label_df_down$chr[which(label_df_down$chr %in% c("31"))] <- "\n31"
label_df_down$chr[which(label_df_down$chr %in% c("33"))] <- "\n33"
label_df_down$chr[which(label_df_down$chr %in% c("35"))] <- "\n35"
label_df_down$chr[which(label_df_down$chr %in% c("37"))] <- "\n37"
# plot
(p2 <- ggplot(df_down, aes(x = cum_start_flip / 1e6, y = n_het/1000, fill=chr_color)) +
    geom_col(width = 1, position = "identity") +
    scale_fill_manual(values = c("lightblue","steelblue")) +
    scale_x_continuous(breaks = label_df_down$midpoint / 1e6,expand = c(0, 0),limits= c(-3, max(df_down$cum_start_flip) / 1e6),
                       labels = label_df_down$chr) + #,guide=guide_axis(n.dodge=2)
    scale_y_continuous(expand = c(0, 0), limits = c(-1, 22))+
    geom_hline(yintercept = 0, color = "black", linewidth = 1)+
    geom_vline(xintercept = 0, color = "black", linewidth = 1)+
    labs(y = "") +
    annotate("point", x = centro_down$cum_pos_flip/1e6, y = 0, shape = 24, size =4, fill = "red", color = "black")+
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          legend.position = "none",
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          base_family = "Arial"
          )
)


# plot together
library(grid)
library(gridExtra)

# Convert to grobs
g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)

# Automatically calculate x-range widths
range1 <- diff(range(df_up$cum_start))  # 5000
range2 <- diff(range(df_down$cum_start))  # 4000

# Compute relative widths
max_width <- max(range1, range2)
width1 <- range1 / max_width
width2 <- range2 / max_width

# Wrap each plot in a viewport with the right width
g1_wrapped <- gTree(children = gList(g1), vp = viewport(width = unit(width1, "npc")))
g2_wrapped <- gTree(children = gList(g2), vp = viewport(width = unit(width2, "npc")))


png("minimap2_syri_SNP_v8.png", units = "in", width = 25, height = 4, res=300)
library(svglite)
svglite("minimap2_syri_SNP_v8.svg", width = 25, height = 4)

# Draw them stacked vertically
grid.newpage()
# Stack plots vertically with fixed relative widths
#pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1)))
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 2,widths = unit.c(unit(1, "lines"), unit(1, "null")))))

# Add y-axis label (rotated vertical text) in left margin (column 1, spans both rows)
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1:2))
grid.text("      Heterozygosity / 1kbp", rot = 90,gp = gpar(fontsize = 20))
popViewport()

# Plot 1, left aligned
pushViewport(viewport(layout.pos.row = 1))
pushViewport(viewport(width = unit(width1, "npc"), x = 0, just = "left"))
grid.draw(g1)
popViewport(2)

# Plot 2, left aligned
pushViewport(viewport(layout.pos.row = 2))
pushViewport(viewport(width = unit(width2, "npc"), x = 0, just = "left"))
grid.draw(g2)
popViewport(2)

dev.off()












#-------------------
#-------------------
#-------------------
# plot actual chr size, split macro and micro
df_macro <- df %>% filter(chr %in% c("chr1","chr1A","chr2","chr3","chr4","chr5","chr6","chr7","chr8"))

# re-calculate the mid-point
macro_order <- c("chr1","chr1A","chr2","chr3","chr4","chr5","chr6","chr7","chr8")
macro_lengths_ordered <- chr_lengths[macro_order]
macro_offsets <- cumsum(c(0, head(macro_lengths_ordered, -1)))
names(macro_offsets) <- macro_order
label_df_macro <- data.frame(
  chr = macro_order,
  midpoint = macro_offsets + macro_lengths_ordered / 2
)
label_df_macro$chr <- sub("^chr", "", label_df_macro$chr)
# re-calculate start point
df_macro$cum_start <- df_macro$start + macro_offsets[df_macro$chr]
df_macro$chr_color <- as.factor(as.integer(factor(df_macro$chr)) %% 2)
# make 0 visible
df_macro$plot_het <- ifelse(df_macro$n_het < 300, 300, df_macro$n_het)
# plot
(p1 <- ggplot(df_macro, aes(x = cum_start / 1e6, y = plot_het/1000, fill=chr_color)) +
  geom_col(width = 1) +
  scale_fill_manual(values = c("lightblue","steelblue")) +
  scale_x_continuous(breaks = label_df_macro$midpoint / 1e6,labels = label_df_macro$chr,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 22))+
  geom_hline(yintercept = 0, color = "black", linewidth = 1)+
  geom_vline(xintercept = 0, color = "black", linewidth = 1)+
  labs(y = "SNPs / 1kbp") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_text(angle = 45, vjust = 0.5), # hjust = 0.5
        panel.grid = element_blank(),
        axis.line = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))
)


# microchromosomes
df_micro <- df %>% filter(!chr %in% c("chr1","chr1A","chr2","chr3","chr4","chr5","chr6","chr7","chr8"))

# re-calculate the mid-point
micro_order <- mixedsort(unique(df_micro$chr))
micro_lengths_ordered <- chr_lengths[micro_order]
micro_offsets <- cumsum(c(0, head(micro_lengths_ordered, -1)))
names(micro_offsets) <- micro_order
label_df_micro <- data.frame(
  chr = micro_order,
  midpoint = micro_offsets + micro_lengths_ordered / 2
)
label_df_micro$chr <- sub("^chr", "", label_df_micro$chr)
# re-calculate start point
df_micro$cum_start <- df_micro$start + micro_offsets[df_micro$chr]
df_micro$chr_color <- as.factor(as.integer(factor(df_micro$chr)) %% 2)
# make 0 visible
df_micro$plot_het <- ifelse(df_micro$n_het < 300, 300, df_micro$n_het)
# plot
(p2 <- ggplot(df_micro, aes(x = cum_start / 1e6, y = plot_het/1000, fill=chr_color)) +
  geom_col(width = 1) +
  scale_fill_manual(values = c("steelblue","lightblue")) +
  scale_x_continuous(breaks = label_df_micro$midpoint / 1e6,labels = label_df_micro$chr,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 22))+
  geom_hline(yintercept = 0, color = "black", linewidth = 1)+
  labs(x = NULL, y = NULL)+
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_text(angle = 45, vjust = 0.5), # hjust = 0.5
        panel.grid = element_blank(),
        axis.line = element_blank(),
        legend.position = "none",
        axis.line.y = element_blank(),   # hide y axis line
        axis.ticks.y = element_blank(),  # hide y axis ticks
        axis.text.y = element_blank(),   # hide y axis tick labels
        axis.title.y = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))   # hide y axis title
)

# plot together
library(patchwork)
png("minimap2_syri_SNP_v3.png", units = "in", width = 32, height = 4, res=300)
p1 + plot_spacer() + p2 + plot_layout(widths = c(0.7, 0.03, 1.35))
p1 + plot_spacer() + p2 + plot_layout(widths = c(0.7, 0.01, 1.45)) 
dev.off()









#-------------------
# plot equal chr length
fixed_chr_width <- 1
n_chr <- length(chr_order)
# Offsets: chromosomes stacked one after another with fixed width
chr_offsets_fixed <- setNames((0:(n_chr - 1)) * fixed_chr_width, chr_order)

# Compute relative position inside each chromosome scaled to [0,1]
df$rel_pos <- df$start / seqlengths(gr)[df$chr]
# Then calculate cumulative position on fixed scale:
df$fixed_cum_pos <- df$rel_pos * fixed_chr_width + chr_offsets_fixed[df$chr]

# create a mid-position for x axis text
label_df_fixed <- data.frame(chr = chr_order,midpoint = chr_offsets_fixed + fixed_chr_width / 2)

# plot
# to make 0 visible
df$plot_het <- ifelse(df$n_het == 0, 100, df$n_het)

ggplot(df, aes(x = fixed_cum_pos, y = plot_het / 1000, fill = chr_color)) +
  geom_col(width=0.05) +  # width just under fixed width for gaps
  scale_x_continuous(breaks = label_df_fixed$midpoint,labels = label_df_fixed$chr,expand = c(0, 0)) +
  scale_fill_manual(values = c("lightblue","steelblue")) +
  labs(x = NULL,y = "Het / 1kbp") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    panel.grid = element_blank(),
    legend.position = "none"
  )


#-------------------
# plot log scale chr
scale_factor <- 1e7  # adjust as needed
df$cum_start_plog <- log1p(df$cum_start / scale_factor)  # log(1 + x)
label_df <- data.frame(
  chr = chr_order,
  midpoint = chr_offsets + chr_lengths_ordered / 2
)
label_df$midpoint_plog <- log1p(label_df$midpoint / scale_factor)

ggplot(df, aes(x = cum_start_plog, y = n_het / 1000, fill = chr_color)) +
  geom_col(width = 0.02) +  # width on transformed scale, adjust if needed
  scale_x_continuous(
    breaks = label_df$midpoint_plog,
    labels = label_df$chr,
    expand = c(0, 0)
  ) +
  labs(
    title = "Genome-wide Heterozygosity Scan (Pseudo-log scale on chromosomes)",
    x = "Genomic Position (pseudo-log scale)",
    y = "Het / 1kbp"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    panel.grid = element_blank(),
    legend.position = "none"
  )















