# Generate plots for zebra finch T2T
# investigate the associations between inversions and segmental duplications (mechanism non-allelic homologous recombination)

# Notes on BISER parameter:
# --max-error (total alignment error): This counts everything: mismatches + all indels (including large gaps)
# --max-edit-error (edit error): This counts mismatches and small indels but excludes large gaps.



library(tidyverse)
library(GenomicRanges)

# setting working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# -----------------------------------------------------------------------------
# input and output data
# -----------------------------------------------------------------------------

syri_file <- "data/syri.out"
category_file <- "data/chr_classification.csv"
fai_files <- c("data/mat.fai", "data/pat.fai")
biser_file <- "data/bTaeGut7v04.BISER.sd"


# -----------------------------------------------------------------------------
# 1. Load and prepare data
# -----------------------------------------------------------------------------

#------------------
# inversion dataset
#------------------
# function to read SyRI output
# 1) keep inversions only, 2) calculate length (average of maternal and paternal) from
read_syri_inversions <- function(syri_file, annotation_filter = "INV") {
  # SyRI output: 12 columns
  syri_data <- read_tsv(syri_file, 
                        col_names = c("ref_chrom", "ref_start", "ref_end", 
                                      "ref_seq", "ref_extra",
                                      "query_chrom", "query_start", "query_end",
                                      "ID", "parent_ID", "annotation", "extra"),
                        col_types = cols(
                          ref_chrom = col_character(),
                          ref_start = col_integer(),
                          ref_end = col_integer(),
                          ref_seq = col_character(),
                          ref_extra = col_character(),
                          query_chrom = col_character(),
                          query_start = col_integer(),
                          query_end = col_integer(),
                          ID = col_character(),
                          parent_ID = col_character(),
                          annotation = col_character(),
                          extra = col_character()
                        ),
                        comment = "#")
  
  # Filter for specified annotation category
  filtered_data <- syri_data %>%
    filter(annotation %in% annotation_filter) %>%
    mutate(
      # Lengths on each haplotype
      mat_length = ref_end - ref_start,
      pat_length = query_end - query_start,
      # Use average length for general stats
      length = (mat_length + pat_length) / 2,
      # Extract chromosome ID (ref and query should match)
      chrom_id = str_extract(ref_chrom, "(?<=chr)[^_]+"),
      chrom_base = ref_chrom
    ) %>%
    filter(!is.na(chrom_id))
  
  return(filtered_data)
}
read_syri_inversions_chrZ <- function(syri_file, annotation_filter = "INV") {
  # SyRI output: 12 columns
  syri_data <- read_tsv(syri_file, 
                        col_names = c("ref_chrom", "ref_start", "ref_end", 
                                      "ref_seq", "ref_extra",
                                      "query_chrom", "query_start", "query_end",
                                      "ID", "parent_ID", "annotation", "extra"),
                        col_types = cols(
                          ref_chrom = col_character(),
                          ref_start = col_integer(),
                          ref_end = col_integer(),
                          ref_seq = col_character(),
                          ref_extra = col_character(),
                          query_chrom = col_character(),
                          query_start = col_integer(),
                          query_end = col_integer(),
                          ID = col_character(),
                          parent_ID = col_character(),
                          annotation = col_character(),
                          extra = col_character()
                        ),
                        comment = "#")
  
  # Filter for specified annotation category
  filtered_data <- syri_data %>%
    filter(annotation %in% annotation_filter) %>%
    mutate(
      # Lengths on each haplotype
      mat_length = ref_end - ref_start,
      pat_length = query_end - query_start,
      # Use average length for general stats
      length = (mat_length + pat_length) / 2,
      # Extract chromosome ID (ref and query should match)
      chrom_id = "Z",
      chrom_base = "chrZ"
    ) %>%
    filter(!is.na(chrom_id))
  
  return(filtered_data)
}
# run
inv_data_autosome <- read_syri_inversions(syri_file, annotation_filter = c("INV"))
inv_data_1.4 <- read_syri_inversions_chrZ("data/syri_bTaeGut7_bTaeGut1.4.out", annotation_filter = c("INV"))
inv_data_3.2.4 <- read_syri_inversions_chrZ("data/syri_bTaeGut7_bTaeGut3.2.4.out", annotation_filter = c("INV"))
inv_data <- rbind(inv_data_autosome,inv_data_1.4,inv_data_3.2.4)


# function to read chromosome category: macro, micro, dot
read_chrom_categories <- function(category_file) {
  read_csv(category_file, col_types = cols(.default = col_character())) %>%
    set_names(c("chrom_base", "category")) %>%
    mutate(
      chrom_base = str_trim(chrom_base),
      category = str_trim(category)
    )
}
# run 
categories <- read_chrom_categories(category_file)

# function to read chromosome length from fai files
read_fai_files <- function(fai_files) {
  fai_data <- map_df(fai_files, function(f) {
    read_tsv(f,
             col_names = c("chrom", "chrom_length", "offset", "line_bases", "line_width"),
             col_types = cols(
               chrom = col_character(),
               chrom_length = col_integer(),
               offset = col_double(),
               line_bases = col_integer(),
               line_width = col_integer()
             ))
  }) %>%
    select(chrom, chrom_length) %>%
    mutate(
      chrom_id = str_extract(chrom, "(?<=chr)[^_]+"),
      haplotype = str_extract(chrom, "(?<=_)[a-z]+$")
    ) %>%
    filter(!is.na(chrom_id))
  
  return(fai_data)
}
# run
chrom_lengths <- read_fai_files(fai_files)


# function add chromosome category and length to inversion dataset
add_categories_and_lengths <- function(inv_data, categories, chrom_lengths) {
  # Get maternal and paternal lengths separately
  mat_lengths <- chrom_lengths %>%
    filter(haplotype == "mat") %>%
    select(chrom_id, mat_chrom_length = chrom_length)
  
  pat_lengths <- chrom_lengths %>%
    filter(haplotype == "pat") %>%
    select(chrom_id, pat_chrom_length = chrom_length)
  
  inv_data <- inv_data %>%
    left_join(categories, by = "chrom_base") %>%
    left_join(mat_lengths, by = "chrom_id") %>%
    left_join(pat_lengths, by = "chrom_id") %>%
    mutate(
      # Combined length for normalization
      total_chrom_length = mat_chrom_length + pat_chrom_length
    )
  
  inv_data <- inv_data %>%
    mutate(category = factor(category, levels = c("macro", "micro", "dot", "sex")))
  
  # Natural sort order
  chrom_order <- inv_data %>%
    distinct(chrom_id, category) %>%
    mutate(
      num_part = as.numeric(str_extract(chrom_id, "^[0-9]+")),
      letter_part = str_extract(chrom_id, "[A-Za-z]+$")
    ) %>%
    arrange(category, num_part, letter_part) %>%
    pull(chrom_id)
  
  inv_data <- inv_data %>%
    mutate(chrom_id = factor(chrom_id, levels = chrom_order))
  
  return(inv_data)
}
# run
inv_data <- add_categories_and_lengths(inv_data, categories, chrom_lengths)


# Note that three paternal chromosomes are flipped for syri analyses
# chr16, chr26, chr27
# flip the coordinates
# Function to flip query coordinates for specified chromosomes
# Function to flip query coordinates for specified chromosomes
# Run AFTER add_categories_and_lengths()
flip_query_coordinates <- function(inv_data, chroms_to_flip = c("chr16", "chr26", "chr27")) {
  inv_data <- inv_data %>%
    mutate(
      # Flip coordinates for specified chromosomes
      query_start_new = if_else(
        query_chrom %in% chroms_to_flip,
        pat_chrom_length - query_end + 1,
        query_start
      ),
      query_end_new = if_else(
        query_chrom %in% chroms_to_flip,
        pat_chrom_length - query_start + 1,
        query_end
      ),
      # Replace original coordinates
      query_start = query_start_new,
      query_end = query_end_new
    ) %>%
    select(-query_end_new,-query_start_new)
  
  return(inv_data)
}
# run
inv_data_flipped <- flip_query_coordinates(inv_data, chroms_to_flip = c("chr16", "chr26", "chr27"))

# Prepare inversions with matching chromosome format
inv_data_column_select <- bind_rows(
  # Maternal
  inv_data_flipped %>%
    transmute(
      inv_id = ID,
      chr = paste0(ref_chrom, "_mat"),
      start = ref_start,
      end = ref_end,
      inv_length = mat_length,
      chrom_id,
      category,
      haplotype = "mat"
    ),
  # Paternal
  inv_data_flipped %>%
    transmute(
      inv_id = ID,
      chr = paste0(query_chrom, "_pat"),
      start = query_start,
      end = query_end,
      inv_length = pat_length,
      chrom_id,
      category,
      haplotype = "pat"
    )
)

# Remove the chrZ coordinates for bTaeGut1.4 and bTaeGut3.2.4
inversions_for_analysis <- inv_data_column_select %>% filter(!chr %in% c("NC_044241.2_pat","NC_011493.1_pat"))
inversions_for_analysis[inversions_for_analysis=="NC_133063.1_mat"] <- "chrZ_pat" 
inversions_for_analysis$haplotype[which(inversions_for_analysis$chr=="chrZ_pat")] <- "pat"
inversions_for_analysis$category <- as.character(inversions_for_analysis$category)
inversions_for_analysis$category[inversions_for_analysis$chr == "chrZ_pat"] <- "chrZ"






#----------
# SD dataset
#----------
# get segmental duplication file from biser output
load_sds <- function(file_path) {
  #' Load BISER output file
  #'
  #' @param file_path Path to BISER .sd output file
  #' @return Tibble of all segmental duplications with calculated fields
  
  biser_cols <- c("chr1", "start1", "end1", "chr2", "start2", "end2", 
                  "reference", "error_rate", "strand1", "strand2",
                  "max_len", "aln_len", "cigar", "optional",
                  "mismatch_rate", "gap_rate")
  
  sds <- read_tsv(file_path, col_names = biser_cols, show_col_types = FALSE) %>%
    mutate(
      orientation = if_else(strand1 == strand2, "direct", "inverted"),
      identity = 100 - error_rate,
      length = end1 - start1
    )
  
  return(sds)
}
# run
sds <- load_sds(biser_file)


# - Filter SDs using: 
filter_sds <- function(sds, 
                       min_identity = 90, # identity = 100 - mismatch rate - gap rate
                       min_length = 1000, # length = end1 - start1
                       same_chromosome = TRUE, # same chromosome or not
                       orientation = "inverted", # keep inverted, directed, or both
                       remove_overlapping = FALSE, # whether SD1 and SD2 are overlapping
                       remove_identical = TRUE, # SD pairs are same coordinates, likely palindromes
                       deduplicate = TRUE, # remove duplicate rows, i.e. the SD1 and SD2 are swapped 
                       # make sure SD1 coordinates are always smaller than SD2
                       normalize = TRUE) {
  #' Filter and process segmental duplications
  #'
  #' @param sds Tibble from load_sds()
  #' @param min_identity Minimum sequence identity percentage (default: 90)
  #' @param min_length Minimum SD length in bp (default: 1000)
  #' @param same_chromosome If TRUE, keep only intra-chromosomal SDs (default: TRUE)
  #' @param orientation Filter by orientation: "inverted", "direct", or "both" (default: "inverted")
  #' @param remove_overlapping If TRUE, remove SD pairs where SD1 and SD2 overlap (default: FALSE)
  #' @param remove_identical If TRUE, remove SD pairs where SD1 and SD2 have same coordinates (default: TRUE)
  #' @param deduplicate If TRUE, remove reciprocal duplicate SD pairs (default: TRUE)
  #' @param normalize If TRUE, ensure SD1 is always upstream of SD2 (default: TRUE)
  #' @return Filtered and processed tibble of SDs
  
  # Validate orientation parameter
  if (!orientation %in% c("inverted", "direct", "both")) {
    stop("orientation must be 'inverted', 'direct', or 'both'")
  }
  
  # Track initial count
  n_initial <- nrow(sds)
  message(sprintf("Initial SD count: %d", n_initial))
  
  # Filter by identity
  n_before <- nrow(sds)
  filtered_sds <- sds %>%
    filter(identity >= min_identity)
  message(sprintf("After identity >= %d%%: %d (removed %d)", 
                  min_identity, nrow(filtered_sds), n_before - nrow(filtered_sds)))
  
  # Filter by length
  n_before <- nrow(filtered_sds)
  filtered_sds <- filtered_sds %>%
    filter(length >= min_length)
  message(sprintf("After length >= %d bp: %d (removed %d)", 
                  min_length, nrow(filtered_sds), n_before - nrow(filtered_sds)))
  
  # Filter by orientation
  if (orientation != "both") {
    n_before <- nrow(filtered_sds)
    filtered_sds <- filtered_sds %>%
      filter(orientation == !!orientation)
    message(sprintf("After orientation == '%s': %d (removed %d)", 
                    orientation, nrow(filtered_sds), n_before - nrow(filtered_sds)))
  }
  
  # Filter for same chromosome
  if (same_chromosome) {
    n_before <- nrow(filtered_sds)
    filtered_sds <- filtered_sds %>%
      filter(chr1 == chr2)
    message(sprintf("After same chromosome: %d (removed %d)", 
                    nrow(filtered_sds), n_before - nrow(filtered_sds)))
  }
  
  # Remove identical SD pairs
  if (remove_identical) {
    n_before <- nrow(filtered_sds)
    filtered_sds <- filtered_sds %>%
      filter(!(start1 == start2 & end1 == end2))
    message(sprintf("After removing identical pairs: %d (removed %d)", 
                    nrow(filtered_sds), n_before - nrow(filtered_sds)))
  }
  
  # Remove overlapping SD pairs
  if (remove_overlapping) {
    n_before <- nrow(filtered_sds)
    filtered_sds <- filtered_sds %>%
      filter(!(start1 <= end2 & start2 <= end1))
    message(sprintf("After removing overlapping pairs: %d (removed %d)", 
                    nrow(filtered_sds), n_before - nrow(filtered_sds)))
  }
  
  # Deduplicate reciprocal SD pairs
  if (deduplicate) {
    n_before <- nrow(filtered_sds)
    filtered_sds <- filtered_sds %>%
      mutate(
        pair_id = case_when(
          start1 <= start2 ~ paste(chr1, start1, end1, chr2, start2, end2, sep = "_"),
          TRUE ~ paste(chr2, start2, end2, chr1, start1, end1, sep = "_")
        )
      ) %>%
      distinct(pair_id, .keep_all = TRUE) %>%
      select(-pair_id)
    message(sprintf("After deduplication: %d (removed %d)", 
                    nrow(filtered_sds), n_before - nrow(filtered_sds)))
  }
  
  # Normalize so SD1 is always upstream of SD2
  # Normalize so SD1 is always upstream of SD2 (only for same chromosome)
  if (normalize) {
    n_flipped <- sum(filtered_sds$start1 > filtered_sds$start2 & filtered_sds$chr1 == filtered_sds$chr2)
    filtered_sds <- filtered_sds %>%
      mutate(
        needs_flip = start1 > start2 & chr1 == chr2,
        
        orig_start1 = start1,
        orig_end1 = end1,
        orig_start2 = start2,
        orig_end2 = end2,
        orig_strand1 = strand1,
        orig_strand2 = strand2,
        
        start1 = if_else(needs_flip, orig_start2, orig_start1),
        end1 = if_else(needs_flip, orig_end2, orig_end1),
        start2 = if_else(needs_flip, orig_start1, orig_start2),
        end2 = if_else(needs_flip, orig_end1, orig_end2),
        strand1 = if_else(needs_flip, orig_strand2, orig_strand1),
        strand2 = if_else(needs_flip, orig_strand1, orig_strand2)
      ) %>%
      select(-needs_flip, -starts_with("orig_"))
    message(sprintf("Normalized SD pairs (flipped %d intra-chromosomal)", n_flipped))
  }
  
  message(sprintf("Final SD count: %d (%.1f%% of initial)", 
                  nrow(filtered_sds), 100 * nrow(filtered_sds) / n_initial))
  
  return(filtered_sds)
} 
# run 
sds_filtered <- filter_sds(sds, 
                           min_identity = 0, 
                           min_length = 0, 
                           same_chromosome = F,
                           orientation = "both", 
                           remove_overlapping = F, 
                           remove_identical = T,                       
                           deduplicate = T, 
                           normalize = T)

#Initial SD count: 40758
#After identity >= 0%: 40758 (removed 0)
#After length >= 0 bp: 40758 (removed 0)
#After removing identical pairs: 40711 (removed 47)
#After deduplication: 39150 (removed 1561)
#Normalized SD pairs (flipped 1807 intra-chromosomal)
#Final SD count: 39150 (96.1% of initial)

# because the chromosome name in sds file has haplotype info, need to harmonize the chromosome name
sds_filtered_harmonized <- sds_filtered %>%
  mutate(
    haplotype = str_extract(chr1, "mat$|pat$"),
    sd_pair_id = row_number()
  )



# -----------------------------------------------------------------------------
# Plot A. INV: inversion size across each chromosome
# -----------------------------------------------------------------------------
# Note: FigS1A in Porubsky et al. (2022)

# note the inversion size is based on maternal haplotype
plot_inversion_sizes <- function(inv_data, categories) {
  
  # Check if chrom_id and category are consistent within each inv_id
  inconsistent <- inv_data %>%
    group_by(inv_id) %>%
    summarize(
      n_chrom = n_distinct(chrom_id),
      n_cat = n_distinct(category),
      .groups = "drop"
    ) %>%
    filter(n_chrom > 1 | n_cat > 1)
  
  if (nrow(inconsistent) > 0) {
    details <- inv_data %>%
      filter(inv_id %in% inconsistent$inv_id) %>%
      group_by(inv_id) %>%
      summarize(
        chrom_ids = paste(unique(chrom_id), collapse = ", "),
        categories = paste(unique(category), collapse = ", "),
        .groups = "drop"
      )
    warning("Inconsistent values found:\n", 
            paste(capture.output(print(details)), collapse = "\n"))
  }
  
  # Summarize by inv_id - average length, keep other attributes
  inversions <- inv_data %>%
    group_by(inv_id) %>%
    summarize(
      avg_length = mean(inv_length),
      chrom_id = chrom_id[1],
      category = category[1],
      .groups = "drop"
    ) %>%
    mutate(chrom_id = factor(chrom_id, levels = str_sort(unique(chrom_id), numeric = TRUE)))
  
  # Get count for label
  n_inv <- nrow(inversions)
  
  # Plot with color by category
  ggplot(inversions, aes(x = chrom_id, y = avg_length, color = category)) +
    geom_jitter(width = 0.07, alpha = 0.8, size = 2) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.3) +
    scale_y_log10(
      breaks = c(1e3, 1e4, 1e5, 1e6, 1e7),
      labels = c("1 kbp", "10 kbp", "100 kbp", "1 Mbp", "10 Mbp")
    ) +
    scale_color_manual(values = c(
      "dot" = "#D55E00",
      "macro" = "#0072B2",
      "micro" = "#F0E442",
      "chrZ" = "#CC79A7"
    )) +
    labs(x = "Chromosome", y = paste0("Inversion size (n=", n_inv, ")"), color = NULL) +
    theme_bw()
}

# run
pA <- plot_inversion_sizes(inversions_for_analysis, categories)
print(pA)





# -----------------------------------------------------------------------------
# Calculate flanking distance
# -----------------------------------------------------------------------------


# all inversions coordinates from small to large
any(inversions_for_analysis$start >= inversions_for_analysis$end)
# SD1 coordinates is always smaller than SD2 if they are on the same chromosome
any(sds_filtered_harmonized$start1 > sds_filtered_harmonized$start2 & 
      sds_filtered_harmonized$chr1 == sds_filtered_harmonized$chr2)



#' For each inversion, reports:
#' 1. Closest upstream SD (full SD pair info) - with tiebreaker for intra-chr and other SD proximity
#' 2. Closest downstream SD (full SD pair info) - with tiebreaker for intra-chr and other SD proximity
#' 3. Closest flanking SD pair(s) (best one selected by: closest distance, then inverted orientation)

calculate_inversion_flanking_sd <- function(inversion, sds) {
  #' Calculate distance from inversion breakpoints to flanking SDs
  #'
  #' For each inversion, reports:
  
  #' 1. Closest upstream SD (full SD pair info) - with tiebreaker for intra-chr and other SD proximity
  #' 2. Closest downstream SD (full SD pair info) - with tiebreaker for intra-chr and other SD proximity
  #' 3. Closest flanking SD pair (best one selected by: closest distance, then inverted orientation)
  #'
  #' @param inversion Data frame with columns: inv_id, chr, start, end, category
  #' @param sds Data frame of SD pairs with columns: chr1, start1, end1, chr2, start2, end2, orientation
  #'
  #' @return List with two data frames:
  #'   - wide: one row per inversion with all SD info in columns
  #'   - long: three rows per inversion (upstream, downstream, flanking) for visualization
  
  
  # Validate required columns
  required_inv_cols <- c("inv_id", "chr", "start", "end", "category")
  required_sd_cols <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "orientation")
  
  missing_inv <- setdiff(required_inv_cols, names(inversion))
  if (length(missing_inv) > 0) {
    stop("Missing columns in inversion data: ", paste(missing_inv, collapse = ", "))
  }
  
  missing_sd <- setdiff(required_sd_cols, names(sds))
  if (length(missing_sd) > 0) {
    stop("Missing columns in SD data: ", paste(missing_sd, collapse = ", "))
  }
  
  # ================================================================
  # Calculate flanking SDs (wide format)
  # ================================================================
  results_wide <- inversion %>%
    rowwise() %>%
    mutate(
      flanking_info = list({
        inv_5prime <- start
        inv_3prime <- end
        inv_chr <- chr
        
        # Get SDs where at least one member is on the inversion chromosome
        sd_with_inv_chr <- sds %>%
          filter(chr1 == inv_chr | chr2 == inv_chr)
        
        # Get SDs where both members are on the same chromosome (for flanking pairs)
        sd_same_chr <- sds %>%
          filter(chr1 == inv_chr & chr2 == inv_chr)
        
        # ============================================================
        # 1. Find closest UPSTREAM SD (at least partially upstream of 5' breakpoint)
        #    Tiebreaker: (1) intra-chromosomal pair, (2) other SD closer to inversion
        # ============================================================
        upstream_candidates <- bind_rows(
          # SD1 is on inversion chromosome and upstream
          sd_with_inv_chr %>%
            filter(chr1 == inv_chr & start1 < inv_5prime) %>%
            mutate(
              upstream_which = "SD1",
              upstream_dist = pmax(0, inv_5prime - end1)
            ),
          # SD2 is on inversion chromosome and upstream
          sd_with_inv_chr %>%
            filter(chr2 == inv_chr & start2 < inv_5prime) %>%
            mutate(
              upstream_which = "SD2",
              upstream_dist = pmax(0, inv_5prime - end2)
            )
        )
        
        if (nrow(upstream_candidates) > 0) {
          # Apply tiebreaker
          best_upstream <- upstream_candidates %>%
            mutate(
              # Is pair intra-chromosomal?
              intra_chr = chr1 == chr2,
              
              # Distance from "other" SD to nearest inversion breakpoint
              other_sd_dist_to_inv = case_when(
                upstream_which == "SD1" ~ pmin(
                  abs(start2 - inv_5prime),
                  abs(end2 - inv_5prime),
                  abs(start2 - inv_3prime),
                  abs(end2 - inv_3prime)
                ),
                upstream_which == "SD2" ~ pmin(
                  abs(start1 - inv_5prime),
                  abs(end1 - inv_5prime),
                  abs(start1 - inv_3prime),
                  abs(end1 - inv_3prime)
                )
              )
            ) %>%
            # Sort by: distance (primary), intra-chr (prefer TRUE), other SD distance (prefer closer)
            arrange(upstream_dist, desc(intra_chr), other_sd_dist_to_inv) %>%
            slice_head(n = 1)
          
          upstream_result <- list(
            upstream_which = best_upstream$upstream_which,
            upstream_dist = best_upstream$upstream_dist,
            upstream_pair_chr1 = best_upstream$chr1,
            upstream_pair_start1 = best_upstream$start1,
            upstream_pair_end1 = best_upstream$end1,
            upstream_pair_chr2 = best_upstream$chr2,
            upstream_pair_start2 = best_upstream$start2,
            upstream_pair_end2 = best_upstream$end2,
            upstream_pair_orientation = best_upstream$orientation
          )
        } else {
          upstream_result <- list(
            upstream_which = NA_character_,
            upstream_dist = NA_real_,
            upstream_pair_chr1 = NA_character_,
            upstream_pair_start1 = NA_integer_,
            upstream_pair_end1 = NA_integer_,
            upstream_pair_chr2 = NA_character_,
            upstream_pair_start2 = NA_integer_,
            upstream_pair_end2 = NA_integer_,
            upstream_pair_orientation = NA_character_
          )
        }
        
        # ============================================================
        # 2. Find closest DOWNSTREAM SD (at least partially downstream of 3' breakpoint)
        #    Tiebreaker: (1) intra-chromosomal pair, (2) other SD closer to inversion
        # ============================================================
        downstream_candidates <- bind_rows(
          # SD1 is on inversion chromosome and downstream
          sd_with_inv_chr %>%
            filter(chr1 == inv_chr & end1 > inv_3prime) %>%
            mutate(
              downstream_which = "SD1",
              downstream_dist = pmax(0, start1 - inv_3prime)
            ),
          # SD2 is on inversion chromosome and downstream
          sd_with_inv_chr %>%
            filter(chr2 == inv_chr & end2 > inv_3prime) %>%
            mutate(
              downstream_which = "SD2",
              downstream_dist = pmax(0, start2 - inv_3prime)
            )
        )
        
        if (nrow(downstream_candidates) > 0) {
          # Apply tiebreaker
          best_downstream <- downstream_candidates %>%
            mutate(
              # Is pair intra-chromosomal?
              intra_chr = chr1 == chr2,
              
              # Distance from "other" SD to nearest inversion breakpoint
              other_sd_dist_to_inv = case_when(
                downstream_which == "SD1" ~ pmin(
                  abs(start2 - inv_5prime),
                  abs(end2 - inv_5prime),
                  abs(start2 - inv_3prime),
                  abs(end2 - inv_3prime)
                ),
                downstream_which == "SD2" ~ pmin(
                  abs(start1 - inv_5prime),
                  abs(end1 - inv_5prime),
                  abs(start1 - inv_3prime),
                  abs(end1 - inv_3prime)
                )
              )
            ) %>%
            # Sort by: distance (primary), intra-chr (prefer TRUE), other SD distance (prefer closer)
            arrange(downstream_dist, desc(intra_chr), other_sd_dist_to_inv) %>%
            slice_head(n = 1)
          
          downstream_result <- list(
            downstream_which = best_downstream$downstream_which,
            downstream_dist = best_downstream$downstream_dist,
            downstream_pair_chr1 = best_downstream$chr1,
            downstream_pair_start1 = best_downstream$start1,
            downstream_pair_end1 = best_downstream$end1,
            downstream_pair_chr2 = best_downstream$chr2,
            downstream_pair_start2 = best_downstream$start2,
            downstream_pair_end2 = best_downstream$end2,
            downstream_pair_orientation = best_downstream$orientation
          )
        } else {
          downstream_result <- list(
            downstream_which = NA_character_,
            downstream_dist = NA_real_,
            downstream_pair_chr1 = NA_character_,
            downstream_pair_start1 = NA_integer_,
            downstream_pair_end1 = NA_integer_,
            downstream_pair_chr2 = NA_character_,
            downstream_pair_start2 = NA_integer_,
            downstream_pair_end2 = NA_integer_,
            downstream_pair_orientation = NA_character_
          )
        }
        
        # ============================================================
        # 3. Find closest FLANKING SD PAIR (same pair flanks both sides)
        #    Tiebreaker: (1) closest max_dist, (2) inverted orientation
        # ============================================================
        flanking_candidates <- sd_same_chr %>%
          mutate(
            config1_valid = (start1 < inv_5prime) & (end2 > inv_3prime),
            dist_5prime_config1 = if_else(config1_valid, pmax(0, inv_5prime - end1), NA_real_),
            dist_3prime_config1 = if_else(config1_valid, pmax(0, start2 - inv_3prime), NA_real_),
            max_dist_config1 = if_else(config1_valid, pmax(dist_5prime_config1, dist_3prime_config1), NA_real_),
            
            config2_valid = (start2 < inv_5prime) & (end1 > inv_3prime),
            dist_5prime_config2 = if_else(config2_valid, pmax(0, inv_5prime - end2), NA_real_),
            dist_3prime_config2 = if_else(config2_valid, pmax(0, start1 - inv_3prime), NA_real_),
            max_dist_config2 = if_else(config2_valid, pmax(dist_5prime_config2, dist_3prime_config2), NA_real_),
            
            best_config = case_when(
              config1_valid & config2_valid ~ if_else(max_dist_config1 <= max_dist_config2, 1L, 2L),
              config1_valid ~ 1L,
              config2_valid ~ 2L,
              TRUE ~ NA_integer_
            ),
            
            best_max_dist = case_when(
              best_config == 1L ~ max_dist_config1,
              best_config == 2L ~ max_dist_config2,
              TRUE ~ NA_real_
            ),
            best_dist_5prime = case_when(
              best_config == 1L ~ dist_5prime_config1,
              best_config == 2L ~ dist_5prime_config2,
              TRUE ~ NA_real_
            ),
            best_dist_3prime = case_when(
              best_config == 1L ~ dist_3prime_config1,
              best_config == 2L ~ dist_3prime_config2,
              TRUE ~ NA_real_
            )
          ) %>%
          filter(!is.na(best_config))
        
        if (nrow(flanking_candidates) > 0) {
          # Select best flanking pair: closest distance, then prefer inverted orientation
          best_flanking <- flanking_candidates %>%
            mutate(
              is_inverted = orientation %in% c("-", "inverted", "inv")
            ) %>%
            arrange(best_max_dist, desc(is_inverted)) %>%
            slice_head(n = 1)
          
          flanking_result <- list(
            flanking_config = best_flanking$best_config,
            flanking_dist_5prime = best_flanking$best_dist_5prime,
            flanking_dist_3prime = best_flanking$best_dist_3prime,
            flanking_max_dist = best_flanking$best_max_dist,
            flanking_pair_chr1 = best_flanking$chr1,
            flanking_pair_start1 = best_flanking$start1,
            flanking_pair_end1 = best_flanking$end1,
            flanking_pair_chr2 = best_flanking$chr2,
            flanking_pair_start2 = best_flanking$start2,
            flanking_pair_end2 = best_flanking$end2,
            flanking_pair_orientation = best_flanking$orientation
          )
        } else {
          flanking_result <- list(
            flanking_config = NA_integer_,
            flanking_dist_5prime = NA_real_,
            flanking_dist_3prime = NA_real_,
            flanking_max_dist = NA_real_,
            flanking_pair_chr1 = NA_character_,
            flanking_pair_start1 = NA_integer_,
            flanking_pair_end1 = NA_integer_,
            flanking_pair_chr2 = NA_character_,
            flanking_pair_start2 = NA_integer_,
            flanking_pair_end2 = NA_integer_,
            flanking_pair_orientation = NA_character_
          )
        }
        
        c(upstream_result, downstream_result, flanking_result)
      })
    ) %>%
    ungroup() %>%
    mutate(
      upstream_which = map_chr(flanking_info, "upstream_which"),
      upstream_dist = map_dbl(flanking_info, "upstream_dist"),
      upstream_pair_chr1 = map_chr(flanking_info, "upstream_pair_chr1"),
      upstream_pair_start1 = map_int(flanking_info, "upstream_pair_start1"),
      upstream_pair_end1 = map_int(flanking_info, "upstream_pair_end1"),
      upstream_pair_chr2 = map_chr(flanking_info, "upstream_pair_chr2"),
      upstream_pair_start2 = map_int(flanking_info, "upstream_pair_start2"),
      upstream_pair_end2 = map_int(flanking_info, "upstream_pair_end2"),
      upstream_pair_orientation = map_chr(flanking_info, "upstream_pair_orientation"),
      
      downstream_which = map_chr(flanking_info, "downstream_which"),
      downstream_dist = map_dbl(flanking_info, "downstream_dist"),
      downstream_pair_chr1 = map_chr(flanking_info, "downstream_pair_chr1"),
      downstream_pair_start1 = map_int(flanking_info, "downstream_pair_start1"),
      downstream_pair_end1 = map_int(flanking_info, "downstream_pair_end1"),
      downstream_pair_chr2 = map_chr(flanking_info, "downstream_pair_chr2"),
      downstream_pair_start2 = map_int(flanking_info, "downstream_pair_start2"),
      downstream_pair_end2 = map_int(flanking_info, "downstream_pair_end2"),
      downstream_pair_orientation = map_chr(flanking_info, "downstream_pair_orientation"),
      
      flanking_config = map_int(flanking_info, "flanking_config"),
      flanking_dist_5prime = map_dbl(flanking_info, "flanking_dist_5prime"),
      flanking_dist_3prime = map_dbl(flanking_info, "flanking_dist_3prime"),
      flanking_max_dist = map_dbl(flanking_info, "flanking_max_dist"),
      flanking_pair_chr1 = map_chr(flanking_info, "flanking_pair_chr1"),
      flanking_pair_start1 = map_int(flanking_info, "flanking_pair_start1"),
      flanking_pair_end1 = map_int(flanking_info, "flanking_pair_end1"),
      flanking_pair_chr2 = map_chr(flanking_info, "flanking_pair_chr2"),
      flanking_pair_start2 = map_int(flanking_info, "flanking_pair_start2"),
      flanking_pair_end2 = map_int(flanking_info, "flanking_pair_end2"),
      flanking_pair_orientation = map_chr(flanking_info, "flanking_pair_orientation")
    ) %>%
    select(-flanking_info)
  
  # ================================================================
  # Convert to long format (3 rows per inversion)
  # ================================================================
  results_long <- bind_rows(
    # Upstream (5' end)
    results_wide %>%
      select(inv_id, chr, start, end, category,
             sd_which = upstream_which,
             sd_dist = upstream_dist,
             sd_pair_chr1 = upstream_pair_chr1,
             sd_pair_start1 = upstream_pair_start1,
             sd_pair_end1 = upstream_pair_end1,
             sd_pair_chr2 = upstream_pair_chr2,
             sd_pair_start2 = upstream_pair_start2,
             sd_pair_end2 = upstream_pair_end2,
             sd_pair_orientation = upstream_pair_orientation) %>%
      mutate(sd_type = "upstream_5prime",
             flanking_dist_5prime = NA_real_,
             flanking_dist_3prime = NA_real_,
             flanking_config = NA_integer_),
    
    # Downstream (3' end)
    results_wide %>%
      select(inv_id, chr, start, end, category,
             sd_which = downstream_which,
             sd_dist = downstream_dist,
             sd_pair_chr1 = downstream_pair_chr1,
             sd_pair_start1 = downstream_pair_start1,
             sd_pair_end1 = downstream_pair_end1,
             sd_pair_chr2 = downstream_pair_chr2,
             sd_pair_start2 = downstream_pair_start2,
             sd_pair_end2 = downstream_pair_end2,
             sd_pair_orientation = downstream_pair_orientation) %>%
      mutate(sd_type = "downstream_3prime",
             flanking_dist_5prime = NA_real_,
             flanking_dist_3prime = NA_real_,
             flanking_config = NA_integer_),
    
    # Flanking pair (best one with tiebreaker)
    results_wide %>%
      select(inv_id, chr, start, end, category,
             sd_dist = flanking_max_dist,
             sd_pair_chr1 = flanking_pair_chr1,
             sd_pair_start1 = flanking_pair_start1,
             sd_pair_end1 = flanking_pair_end1,
             sd_pair_chr2 = flanking_pair_chr2,
             sd_pair_start2 = flanking_pair_start2,
             sd_pair_end2 = flanking_pair_end2,
             sd_pair_orientation = flanking_pair_orientation,
             flanking_config,
             flanking_dist_5prime,
             flanking_dist_3prime) %>%
      mutate(sd_type = "flanking",
             sd_which = case_when(
               flanking_config == 1L ~ "SD1_up_SD2_down",
               flanking_config == 2L ~ "SD2_up_SD1_down",
               TRUE ~ NA_character_
             ))
  ) %>%
    select(inv_id, chr, start, end, category, sd_type, sd_which, sd_dist,
           sd_pair_chr1, sd_pair_start1, sd_pair_end1,
           sd_pair_chr2, sd_pair_start2, sd_pair_end2,
           sd_pair_orientation,
           flanking_dist_5prime, flanking_dist_3prime, flanking_config) %>%
    arrange(inv_id, factor(sd_type, levels = c("upstream_5prime", "downstream_3prime", "flanking")))
  
  # ================================================================
  # Return both formats
  # ================================================================
  return(list(
    wide = results_wide,
    long = results_long
  ))
}

results <- calculate_inversion_flanking_sd(inversions_for_analysis, sds_filtered_harmonized)

# Access wide format (one row per inversion)
results$wide
# Access long format (three rows per inversion for visualization)
results$long






# -----------------------------------------------------------------------------
# Plot B. overview of different SD and inversion
# -----------------------------------------------------------------------------

# Categorize inversions based on SD proximity
categorize_inversions <- function(results_wide, distance_threshold = 10000) {
  
  results_wide %>%
    mutate(
      # Calculate minimum flanking distance (closer of the two breakpoints)
      flanking_min_dist = pmin(flanking_dist_5prime, flanking_dist_3prime, na.rm = TRUE),
      
      # Get minimum distance from upstream/downstream (whichever is closer)
      closest_sd_dist = pmin(upstream_dist, downstream_dist, na.rm = TRUE),
      
      # Determine if closest SD (upstream or downstream) is intra-chromosomal
      closest_is_upstream = case_when(
        is.na(upstream_dist) & is.na(downstream_dist) ~ NA,
        is.na(upstream_dist) ~ FALSE,
        is.na(downstream_dist) ~ TRUE,
        upstream_dist <= downstream_dist ~ TRUE,
        TRUE ~ FALSE
      ),
      
      # Check if closest SD pair is intra-chromosomal
      closest_sd_intra_chr = case_when(
        closest_is_upstream ~ upstream_pair_chr1 == upstream_pair_chr2,
        !closest_is_upstream ~ downstream_pair_chr1 == downstream_pair_chr2,
        TRUE ~ NA
      ),
      
      # Hierarchical categorization (mutually exclusive)
      sd_category = case_when(
        # 1. Flanking SD pair with BOTH ends <10kb (max dist < threshold)
        !is.na(flanking_max_dist) & flanking_max_dist < distance_threshold ~ 
          "Flanking SD pair (max dist <10kb)",
        
        # 2. Flanking SD pair with at least ONE end <10kb OR closest SD <10kb on same chr
        (!is.na(flanking_min_dist) & flanking_min_dist < distance_threshold) |
          (!is.na(closest_sd_dist) & closest_sd_dist < distance_threshold & closest_sd_intra_chr) ~ 
          "Nearest SD <10kb (pair on same chr)",
        
        # 3. Closest SD <10kb, inter-chromosomal
        !is.na(closest_sd_dist) & closest_sd_dist < distance_threshold & !closest_sd_intra_chr ~ 
          "Nearest SD <10kb (pair on diff chr)",
        
        # 4. No SD within threshold
        TRUE ~ "No SD <10kb"
      )
    ) %>%
    mutate(
      sd_category = factor(sd_category, levels = c(
        "Flanking SD pair (max dist <10kb)",
        "Nearest SD <10kb (pair on same chr)",
        "Nearest SD <10kb (pair on diff chr)",
        "No SD <10kb"
      ))
    )
}

# Apply categorization
inversions_categorized <- categorize_inversions(results$wide)

# get unique categorization
collapse_inversions_by_id <- function(inversions_categorized) {
  
  inversions_categorized %>%
    group_by(inv_id) %>%
    mutate(
      n_rows = n(),
      n_distinct_categories = n_distinct(sd_category),
      # Get all categories present for this inv_id
      all_categories = paste(sort(unique(as.character(sd_category))), collapse = "; ")
    ) %>%
    ungroup() %>%
    # For each inv_id, keep the row with the best (lowest factor level) category
    group_by(inv_id) %>%
    arrange(inv_id, as.numeric(sd_category)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    mutate(
      # Add note column
      haplotype_note = case_when(
        n_rows == 1 ~ NA_character_,
        n_distinct_categories == 1 ~ "Both haplotypes: same category",
        TRUE ~ paste0("Haplotypes differ: ", all_categories)
      )
    ) %>%
    select(-n_rows, -n_distinct_categories, -all_categories)
}

# Apply the function
inversions_collapsed <- collapse_inversions_by_id(inversions_categorized)



# Summarize counts
category_summary <- inversions_collapsed %>%
  count(sd_category, category) %>%
  mutate(
    category = factor(category, levels = c("macro", "micro", "dot", "chrZ"))
  )


# Generate stacked bar plot

# Define color palette for chromosome categories
chr_colors <- c(
  
  "macro" = "#0072B2",
  "micro" = "#F0E442",
  "dot" = "#D55E00",
  "chrZ" = "#CC79A7"
)

# Plot 1: Stacked bar chart
pB <- ggplot(category_summary, aes(x = sd_category, y = n, fill = category)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
  geom_text(
    aes(label = after_stat(y), group = sd_category),
    stat = "summary", fun = sum,
    vjust = -0.5, size = 3.5, fontface = "bold"
  ) +
  scale_fill_manual(
    values = chr_colors,
    #name = "Chromosome\nCategory"
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 12)) +
  labs(
    #title = "Inversion categorization by SD proximity (<10 kb threshold)",
    #subtitle = "Hierarchical assignment: flanking pair → nearest SD (same chr) → nearest SD (diff chr) → none",
    x = NULL,
    y = "Number of inversions"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.text.x = element_text(size = 10),
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  expand_limits(y = max(category_summary %>% group_by(sd_category) %>% summarise(total = sum(n)) %>% pull(total)) * 1.1) 

print(pB)





# -----------------------------------------------------------------------------
# Plot C. number of INV at different distance to nearest SD pair
# -----------------------------------------------------------------------------

plot_inversion_flanking_sd_distance <- function(results,
                                                distance_type = c("max", "min"),
                                                log_scale = FALSE,
                                                keep_all_haplotypes = FALSE) {
  #' Plot histogram of distance to flanking SD pairs
  
  #' Tests NAHR hypothesis: same SD pair flanks both inversion breakpoints
  #'
  #' @param results Output from calculate_inversion_flanking_sd() - a list with $wide and $long
  #' @param distance_type "max" (farther SD) or "min" (closer SD)
  #' @param log_scale If TRUE, use log10 scale for x-axis
  #' @param keep_all_haplotypes If TRUE, keep both mat/pat rows; if FALSE, keep only best per inv_id
  #'
  #' @return List with: plot, plot_data
  
  distance_type <- match.arg(distance_type)
  
  # Extract wide format results
  flanking_results <- results$wide
  
  # First calculate both distance metrics
  plot_data <- flanking_results %>%
    mutate(
      min_distance = pmin(flanking_dist_5prime, flanking_dist_3prime),
      max_distance = flanking_max_dist
    )
  
  # Filter and select based on the chosen distance type
  if (distance_type == "max") {
    plot_data <- plot_data %>%
      filter(!is.na(max_distance))
    
    if (!keep_all_haplotypes) {
      plot_data <- plot_data %>%
        group_by(inv_id) %>%
        slice_min(max_distance, n = 1, with_ties = FALSE) %>%
        ungroup()
    }
    
    plot_data <- plot_data %>%
      mutate(plot_distance = max_distance)
    x_label <- "Maximum distance to the closest flanking SD pair"
  } else {
    plot_data <- plot_data %>%
      filter(!is.na(min_distance))
    
    if (!keep_all_haplotypes) {
      plot_data <- plot_data %>%
        group_by(inv_id) %>%
        slice_min(min_distance, n = 1, with_ties = FALSE) %>%
        ungroup()
    }
    
    plot_data <- plot_data %>%
      mutate(plot_distance = min_distance)
    x_label <- "Minimum distance to the closest flanking SD pair"
  }
  
  # Set category factor levels
  plot_data <- plot_data %>%
    mutate(category = factor(category, levels = c("macro", "micro", "dot", "chrZ")))
  
  # Print summary
  n_total <- n_distinct(flanking_results$inv_id)
  n_with_flanking <- n_distinct(plot_data$inv_id)
  n_rows <- nrow(plot_data)
  
  message(sprintf("Unique inversions: %d", n_total))
  message(sprintf("Inversions with valid flanking SD pair: %d (%.1f%%)", 
                  n_with_flanking, 100 * n_with_flanking / n_total))
  
  if (keep_all_haplotypes) {
    message(sprintf("Total rows (including both haplotypes): %d", n_rows))
  }
  
  if (n_with_flanking > 0) {
    message(sprintf("\n%s distance - Median: %.1f bp, Mean: %.1f bp", 
                    str_to_title(distance_type),
                    median(plot_data$plot_distance), 
                    mean(plot_data$plot_distance)))
  }
  
  # Define colors
  category_colors <- c(
    "macro" = "#0072B2",
    "micro" = "#F0E442",
    "dot" = "#D55E00",
    "chrZ" = "#CC79A7"
  )
  
  # Create histogram
  if (log_scale) {
    plot_data <- plot_data %>%
      mutate(plot_distance = plot_distance + 1)  # Avoid log(0)
    
    p <- ggplot(plot_data, aes(x = plot_distance, fill = category)) +
      geom_histogram(color = "white") +
      geom_vline(xintercept = 10000, linetype = "dashed", color = "grey", linewidth = 0.8) +
      scale_x_log10(
        labels = scales::label_number(scale_cut = scales::cut_short_scale())
      ) +
      scale_y_continuous(breaks = scales::pretty_breaks()) +
      scale_fill_manual(values = category_colors) +
      labs(x = paste0(x_label, " (bp, log scale)"), 
           y = sprintf("Number of Inversions (n=%d)", n_with_flanking)) +
      theme_minimal() +
      theme(panel.grid.minor = element_blank(), legend.title = element_blank())
  } else {
    p <- ggplot(plot_data, aes(x = plot_distance / 1000, fill = category)) +
      geom_histogram(color = "white") +
      geom_vline(xintercept = 10000, linetype = "dashed", color = "grey", linewidth = 0.8) +
      scale_x_continuous(labels = scales::comma) +
      scale_y_continuous(breaks = scales::pretty_breaks()) +
      scale_fill_manual(values = category_colors) +
      labs(x = paste0(x_label, " (kbp)"), 
           y = sprintf("Number of Inversions (n=%d)", n_with_flanking)) +
      theme_minimal() +
      theme(panel.grid.minor = element_blank(), legend.title = element_blank())
  }
  
  return(list(
    plot = p,
    plot_data = plot_data
  ))
}


# keeps only one row per inv_id - the haplotype with the smaller distance. 
pC <- plot_inversion_flanking_sd_distance(results, 
                                         distance_type = "max", 
                                         log_scale = TRUE,
                                         keep_all_haplotypes = FALSE)
pC$plot





# -----------------------------------------------------------------------------
# Plot D. visualize inversion with flanking SDs 
# -----------------------------------------------------------------------------

# in one plot
plot_inversion_sd_locus <- function(plot_data, 
                                    n_show = 20,
                                    normalize = FALSE,
                                    dist_threshold = 10000,
                                    filter_by = c("max", "min"),  # NEW PARAMETER
                                    show_orientation = TRUE) {
  #' Plot inversions with flanking SD pairs as arrows
  #'
  #' @param plot_data Data frame from plot_inversion_flanking_sd_distance()$plot_data
  #' @param n_show Number of inversions to show (default 20, set to NULL for all)
  #' @param normalize If TRUE, normalize all inversions to same width for comparison
  
  #' @param dist_threshold Maximum distance to include (default 10000, set to NULL for no filtering)
  #' @param filter_by "max" to filter by max_distance, "min" to filter by min_distance
  #' @param show_orientation If TRUE, arrow direction indicates SD orientation
  #'
  #' @return ggplot object
  
  filter_by <- match.arg(filter_by)
  
  # Optionally filter by distance threshold
  if (!is.null(dist_threshold)) {
    if (filter_by == "max") {
      plot_data <- plot_data %>%
        filter(max_distance <= dist_threshold)
    } else {
      plot_data <- plot_data %>%
        filter(min_distance <= dist_threshold)
    }
  }
  
  if (nrow(plot_data) == 0) {
    stop("No inversions found after filtering")
  }
  
  # Select inversions to show (sorted by max_distance)
  plot_data <- plot_data %>%
    arrange(max_distance)
  
  if (!is.null(n_show)) {
    plot_data <- plot_data %>% head(n_show)
  }
  
  plot_data <- plot_data %>%
    mutate(
      # Calculate sizes
      inv_length = end - start,
      inv_half = inv_length / 2,
      inv_center = (start + end) / 2,
      sd1_length = flanking_pair_end1 - flanking_pair_start1,
      sd2_length = flanking_pair_end2 - flanking_pair_start2,
      
      # Positions relative to inversion center
      sd1_rel_start = flanking_pair_start1 - inv_center,
      sd1_rel_end = flanking_pair_end1 - inv_center,
      sd2_rel_start = flanking_pair_start2 - inv_center,
      sd2_rel_end = flanking_pair_end2 - inv_center,
      
      # Extract chromosome without haplotype suffix
      chr_clean = str_remove(chr, "_(pat|mat|hap1|hap2)$"),
      
      # Extract haplotype
      haplotype = str_extract(chr, "(pat|mat|hap1|hap2)$"),
      haplotype = if_else(is.na(haplotype), "unknown", haplotype),
      
      # Create label with inv_id, chr, AND haplotype
      inv_label = paste0(inv_id, " (", chr_clean, ", ", haplotype, ")")
    )
  
  if (normalize) {
    # Normalize to inversion half-length = 1
    plot_data <- plot_data %>%
      mutate(
        sd1_rel_start = sd1_rel_start / inv_half,
        sd1_rel_end = sd1_rel_end / inv_half,
        sd2_rel_start = sd2_rel_start / inv_half,
        sd2_rel_end = sd2_rel_end / inv_half,
        inv_half = 1
      )
  }
  
  # Determine arrow direction based on orientation
  if (show_orientation) {
    plot_data <- plot_data %>%
      mutate(
        is_inverted = flanking_pair_orientation %in% c("-", "inverted", "inv"),
        sd1_arrow_y = sd1_rel_start,
        sd1_arrow_yend = sd1_rel_end,
        sd2_arrow_y = if_else(is_inverted, sd2_rel_end, sd2_rel_start),
        sd2_arrow_yend = if_else(is_inverted, sd2_rel_start, sd2_rel_end)
      )
  } else {
    plot_data <- plot_data %>%
      mutate(
        sd1_arrow_y = sd1_rel_start,
        sd1_arrow_yend = sd1_rel_end,
        sd2_arrow_y = sd2_rel_start,
        sd2_arrow_yend = sd2_rel_end
      )
  }
  
  # Order by inversion size (small on top)
  plot_data <- plot_data %>%
    mutate(category = factor(category, levels = c("macro", "micro", "dot", "chrZ"))) %>%
    arrange(desc(inv_length)) %>%
    mutate(inv_label = factor(inv_label, levels = inv_label))
  
  # Define colors
  category_colors <- c(
    "macro" = "#0072B2",
    "micro" = "#F0E442",
    "dot" = "#D55E00",
    "chrZ" = "#CC79A7"
  )
  
  # Build plot
  p <- ggplot(plot_data, aes(x = inv_label)) +
    # Center line at 0
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
    # Inversion bar (horizontal line from -inv_half to +inv_half, centered at 0)
    geom_linerange(aes(ymin = -inv_half, ymax = inv_half, color = category), 
                   linewidth = 3, alpha = 0.8) +
    # SD1 arrow (upstream/5' side)
    geom_segment(aes(xend = inv_label, y = sd1_arrow_y, yend = sd1_arrow_yend),
                 color = "gray30", linewidth = 0.8, alpha = 0.6,
                 arrow = arrow(length = unit(0.08, "inches"), type = "closed"),
                 na.rm = TRUE) +
    # SD2 arrow (downstream/3' side)
    geom_segment(aes(xend = inv_label, y = sd2_arrow_y, yend = sd2_arrow_yend),
                 color = "gray30", linewidth = 0.8, alpha = 0.6,
                 arrow = arrow(length = unit(0.08, "inches"), type = "closed"),
                 na.rm = TRUE) +
    scale_color_manual(values = category_colors) +
    coord_flip() +
    labs(
      x = NULL,
      y = if (normalize) "Position relative to inversion (normalized)" else "Position relative to inversion center (bp)"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  if (!normalize) {
    p <- p + scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))
  }
  
  return(p)
}

# Filter by max_distance <= 10kb
pD <- plot_inversion_sd_locus(pC$plot_data, n_show = NULL, 
                        dist_threshold = 10000, 
                        filter_by = "max",  # NEW: filter by max_distance
                        normalize = FALSE)
print(pD)














# -----------------------------------------------------------------------------
# combine into one plot
# -----------------------------------------------------------------------------

library(cowplot)

# Remove legends from all plots
pA_noleg <- pA + theme(legend.position = "none",axis.title.x = element_blank()) + theme(plot.margin = margin(t = 5, r = 10, b = 10, l = 5))
pB_noleg <- pB + theme(legend.position = "none") + theme(plot.margin = margin(t = 5, r = 5, b = 10, l = 10))
pC_noleg <- pC$plot + theme(legend.position = "none") + theme(plot.margin = margin(t = 10, r = 10, b = 5, l = 5))
pD_noleg <- pD + theme(legend.position = "none") + theme(plot.margin = margin(t = 10, r = 5, b = 5, l = 10))


# Extract legend from one plot - use "right" position
legend <- get_legend(pA + theme(legend.position = "right"))


top_row <- plot_grid(pA_noleg, pB_noleg, labels = c("A", "B"), ncol = 2, align = "h")
bottom_row <- plot_grid(pD_noleg, pC_noleg, labels = c("C", "D"), ncol = 2, align = "h")

# Don't align between rows
plots_combined <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  align = "none"  # no vertical alignment between rows
)

combined_plot <- plot_grid(
  plots_combined,
  legend,
  nrow = 1,
  rel_widths = c(1, 0.1)
)

combined_plot
ggsave("Suuple_Figure.svg", combined_plot, width = 11, height = 8)
ggsave("Suuple_Figure.png", combined_plot, width = 11, height = 8)
# done





























