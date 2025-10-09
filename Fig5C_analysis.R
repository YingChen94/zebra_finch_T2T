# Load necessary libraries
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# TOTAL = SYN + NOTAL + INV + INVTR + DUP + TRANS + INVDP
# SNPs are called in the SYN, INV, INVTR and TRANS, not in DUP and INVDP

# SYN = "Syntenic region"  [Fig.S1] = SYNAL + INS + DEL + HDR
# SYNAL = "Syntenic alignment"   
# INV = "Inversion"   = INVAL + INS + DEL + HDR
# INVAL = "Inversion alignment"
# TRANS = "Translocation" = TRANSAL + INS + DEL
# TRANSAL = "Translocation alignment"
# INVTR = "Inverted Translocation"  
# INVTRAL = "Inverted Translocation alignment"
# DUP = "Duplication" 
# DUPAL = "Duplication alignment" = DUPAL +
# INVDP = "Inverted Duplication" 
# INVDPAL = "Inverted Duplication alignment"
# HDR = "Highly diverged regions"
# INS = "Insertion in non-reference genome"
# DEL = "Deletion in non-reference genome"
# CPG = "Copy gain in non-reference genome"
# CPL = "Copy loss in non-reference genome" 
# SNP = "Single nucleotide polymorphism"
# TDM = "Tandem repeat" 
# NOTAL = "Not Aligned region" 

# Read the syri.out file
syri <- read_delim("syri.out", col_names = FALSE)
# take a bit of time in reading it
colnames(syri) <- c("ref_chr", "ref_start", "ref_end", "ref_seq", "alt_seq",
                    "qry_chr", "qry_start", "qry_end", "ID", "parent", "type", "dup_copy_status")
# Note it's the length of maternal haplotype (reference)
syri <- syri %>%
  mutate(ref_start = as.numeric(ref_start), # will have NAs
         ref_end = as.numeric(ref_end), # will have NAs
         length_ref = abs(ref_end - ref_start) + 1)


# chromosome classification
chr_classifi <- read.csv("../../data/chr_classification.csv")
# chromosome length
mat <- read.delim("../mat_filtered.fasta.fai",header=F)[,c(1,2)]
colnames(mat) <- c("chr_name","length_ref")
#pat <- read.delim("../../data/pat.fai",header=F)[,c(1,2)]
#colnames(pat) <- c("chr_name","length_qry")
#pat$chr_name <- gsub("_pat", "", pat$chr_name)

# Total length of variants per type per chromosome
length_matrix <- syri %>%
  group_by(ref_chr, type) %>%
  summarise(total_length_ref = sum(length_ref), 
            .groups = "drop") %>%
  pivot_wider(names_from = type, values_from = c(total_length_ref), values_fill = 0)
length_matrix <- left_join(length_matrix, chr_classifi, by = c("ref_chr" = "chromosome"))
length_matrix <- left_join(length_matrix, mat, by = c("ref_chr" = "chr_name"))
#length_matrix <- left_join(length_matrix, pat, by = c("ref_chr" = "chr_name"))
# Save to file
#write_tsv(length_matrix, "syri_variant_lengths.tsv")

# percentage of length
length_matrix_pct <- length_matrix
length_matrix_pct[2:21] <- length_matrix[2:21] / length_matrix$length_ref


#-------------------------------------------
# plot length percentage
# TOTAL length is divded by (SYN + NOTAL + INV + INVTR + DUP + TRANS + INVDP) category
length_long_pct <- length_matrix_pct %>%
  pivot_longer(
    cols = 2:21,
    names_to = "variant_type",
    values_to = "value"
  )
# filter variant types you don't wanna show
length_long_pct <- length_long_pct %>% filter(ref_chr != "-") %>%
  filter(variant_type %in% c("SYN","NOTAL","INV","INVTR","TRANS","DUP","INVDP"))
# Now plot
# re-order
length_long_pct$variant_type <- factor(length_long_pct$variant_type, 
                                       levels = c("SYN","NOTAL","INV","INVTR","TRANS","DUP","INVDP"))
png("SV_v2.png", units = "in", width = 7, height = 3.5, res=300)
ggplot(length_long_pct, aes(x = variant_type, y = value)) +
  geom_point(position = position_jitter(width = 0.2, height = 0),aes(colour = classification),size=3,alpha=0.6) + 
  #geom_text(aes(label=ref_chr),position = position_jitter(width = 0.2, height = 0),size = 3)+
  #ggrepel::geom_text_repel(aes(label=ref_chr),max.overlaps=10)+
  scale_color_manual(name=NULL,values = c("dot" = "#D55E00", "macro" = "#0072B2", "micro" = "#F0E442")) +
  theme_bw() +
  labs(x = "Structural Variant Type",
       y = "% of haplotype length",
       color = "Chromosome Type") +
  theme(#axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        panel.grid = element_blank())
dev.off()






#-------------------------------------------
# create the stats summary for supplementary table
# Read the syri.out file
syri <- read_delim("syri.out", col_names = FALSE)
# take a bit of time in reading it
colnames(syri) <- c("ref_chr", "ref_start", "ref_end", "ref_seq", "alt_seq",
                    "qry_chr", "qry_start", "qry_end", "ID", "parent", "type", "dup_copy_status")

# calculate the stats
dat <- syri %>%
  mutate(ref_start = as.numeric(ref_start), # will have NAs
         ref_end = as.numeric(ref_end), # will have NAs
         length_mat = abs(ref_end - ref_start) + 1,
         length_pat = abs(as.numeric(qry_end) - as.numeric(qry_start)) + 1)

# Total length of variants per type per chromosome
length_matrix <- dat %>%
  group_by(ref_chr, type) %>%
  summarise(total_length_mat = sum(length_mat), 
            total_length_pat = sum(length_pat), 
            .groups = "drop") %>%
  pivot_wider(names_from = type, values_from = c(total_length_mat, total_length_pat), values_fill = 0)
length_matrix <- length_matrix[gtools::mixedorder(length_matrix$ref_chr), ]
write.csv(length_matrix,"length_matrix_per_chr.csv",quote=F,row.names = F)










