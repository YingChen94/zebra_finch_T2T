# Make a plot of the chromosome length difference between maternal and paternal haplotype

library(ggplot2)
library(scales)
library(ggrepel)
library(dplyr)

mat <- read.table("../data/mat.fai")
pat <- read.table("../data/pat.fai")
pat$V1 <- gsub("_pat", "", pat$V1)
mat$V1 <- gsub("_mat", "", mat$V1)
pat <- pat[,c(1:2)]
mat <- mat[,c(1:2)]
mat$V1[mat$V1=="chrW"] <- "ZW"
pat$V1[pat$V1=="chrZ"] <- "ZW"

# merge paternal and maternal
dat <- merge(pat,mat,by="V1")
colnames(dat) <- c("chromosome","paternal","maternal")
dat$larger <- pmax(dat$paternal, dat$maternal)
dat$diff <- abs(dat$paternal - dat$maternal)
dat$diff_per <- abs(dat$paternal - dat$maternal)/dat$larger

# define chr order
chr_order <- c("chr1","chr1A",paste("chr", 2:4, sep = ""),"chr4A",paste("chr", 5:37, sep = ""), "ZW")
dat$chromosome <- factor(dat$chromosome, levels = chr_order)

# remove sex chromosome because it's a different evolutionary force
dat_noZW <- dat[c(1:39),]


# label the point
library(ggrepel)
label_data <- subset(dat_noZW, diff_per > 0.1)
dat_noZW$label <- ifelse(dat_noZW$diff_per > 0,
                         as.character(dat_noZW$chromosome), NA)
chr_classifi <- read.csv("../data/chr_classification.csv")
dat_noZW <- left_join(dat_noZW, chr_classifi, by = "chromosome")
#png("diff_vs_hap_size_col.png", units = "in", width = 4, height = 3.5, res=300)
(p3 <- ggplot(dat_noZW,aes(larger,diff_per,color=classification))+
    geom_point(alpha=0.5,size=2.6)+
    scale_color_manual(name=NULL,values = c("dot" = "#D55E00", "macro" = "#0072B2", "micro" = "#F0E442")) +
    geom_text_repel(aes(label = label),segment.color = "gray40",segment.size = 0.5,nudge_x=16)+
    #geom_text(aes(label = label, vjust = -1, color = "#D55E00")) +
    scale_y_continuous(limits=c(0,0.3), breaks = seq(0, 1, by = 0.1),labels = scales::percent) +
    scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6))+
    labs(x="Haplotype size",y="Size difference")+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor.x = element_blank(),
          #panel.border = element_rect(color = "black", fill = NA),
          axis.line = element_line(), 
          axis.line.y.right = element_blank(),
          axis.line.x.top = element_blank(),
          legend.position = c(0.85, 0.8),
          legend.background = element_rect(fill = NA, color = NA),
          legend.text = element_text(size = 10),
          legend.key.height = unit(0.7, "lines"),
          axis.title.y = element_text(margin = margin(r = 10)), # Increase right margin to push title away from tick labels
          text = element_text(family = "Arial",size=10),
          plot.margin = margin(0, 0, 0, 0)))
#dev.off()