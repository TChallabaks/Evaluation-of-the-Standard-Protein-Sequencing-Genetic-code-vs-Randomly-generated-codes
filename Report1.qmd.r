# ===============================================================

# A Comprehensive Genetic Evaluation:

# A Comparison Between Protein Sequencing DNA Codes (SGC) & Randomly Generated Protein sequencing DNA codes

# Author: Olatunde Kudus Bakare

# Date: 2025-10-21

# ===============================================================

cat("Starting SGC Benchmark Analysis...\n")

suppressPackageStartupMessages({
  library(jsonlite)
  library(readr)
  library(tidyverse)
  library(readxl)
  library(cluster)
  library(factoextra)
  library(ggrepel)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(here)
})
# Define data paths (update this root path as needed)

data_path <- "~/Documents/CSB195 forked"

aa_prop_path      <- file.path(data_path, "dat", "aa_prop.json")
aa_feature_path   <- file.path(data_path, "dat", "aaFeatureSpace.4.1.Rds")
aa_sim_path       <- file.path(data_path, "dat", "aaSim.4.1.Rds")
aa_dist_path      <- file.path(data_path, "dat", "aaDist.4.1.Rds")
sgc_path          <- file.path(data_path, "dat", "SGC.csv")
breimann1_path    <- file.path(data_path, "dat", "Breimann_2024_Supplementary_Table_1.xlsx")
breimann2_path    <- file.path(data_path, "dat", "Breimann_2024_Supplementary_Table_3.xlsx")

# Auto-detect JSON or CSV

read_aa_table <- function(path) {
  txt <- read_file(path)
  txt_trim <- sub("^\s+", "", txt)
  if (nzchar(txt_trim) && (startsWith(txt_trim, "[") || startsWith(txt_trim, "{"))) {
    as_tibble(fromJSON(txt, flatten = TRUE))
  } else {
    suppressMessages(read_csv(path, show_col_types = FALSE))
  }
}

cat("Loading files...\n")

aa_prop           <- read_aa_table(aa_prop_path)
aaFeatureSpace.4.1 <- readRDS(aa_feature_path)
aa_dist           <- readRDS(aa_dist_path)
SGC               <- read_csv(sgc_path, show_col_types = FALSE)
breimann1         <- read_excel(breimann1_path)
breimann2         <- read_excel(breimann2_path)

cat("Files successfully loaded.\n")

#Normalize Amino Acid Table

normalize_aa_table <- function(df) {
  cols <- names(df)
  A_col   <- dplyr::case_when("A" %in% cols ~ "A", "one" %in% cols ~ "one", TRUE ~ NA_character_)
  Aaa_col <- dplyr::case_when("Aaa" %in% cols ~ "Aaa", "three" %in% cols ~ "three", TRUE ~ NA_character_)
  Amino_col <- dplyr::case_when("AminoAcid" %in% cols ~ "AminoAcid", "name" %in% cols ~ "name", TRUE ~ NA_character_)

  if (all(is.na(c(A_col, Aaa_col, Amino_col)))) stop("aa_prop missing expected columns.")

  out <- df
  if (!is.na(A_col))   out <- out |> mutate(A = toupper(.data[[A_col]]))
  if (!is.na(Aaa_col)) out <- out |> mutate(Aaa = toupper(.data[[Aaa_col]]))
  if (!is.na(Amino_col)) {
    out <- out |> mutate(AminoAcid = ifelse(toupper(.data[[Amino_col]]) == "STOP",
                                            "STOP",
                                            str_to_title(.data[[Amino_col]])))
  }
  out
}

aa_prop <- normalize_aa_table(aa_prop)
glimpse(aa_prop)

# Merge SGC and Amino Acid Properties
SGC <- SGC %>%
  mutate(
    Aaa = toupper(Aaa),
    A   = toupper(A),
    AminoAcid = ifelse(toupper(AminoAcid) == "STOP", "STOP", str_to_title(AminoAcid))
  )

AA_table <- SGC %>% left_join(aa_prop, by = "Aaa")

# Benchmark Computation (SGC vs Random Codes)

SGC_sense <- SGC %>%
  filter(toupper(AminoAcid) != "STOP") %>%
  mutate(Codon = toupper(str_replace_all(Codon, "\s+", "")),
         Codon = chartr("U", "T", Codon))

bases <- c("A","C","G","T")
mutate_pos <- function(codons, pos, to_base) {
  paste0(
    if (pos == 1) to_base else substr(codons, 1, 1),
    if (pos == 2) to_base else substr(codons, 2, 2),
    if (pos == 3) to_base else substr(codons, 3, 3)
  )
}

lookup <- SGC_sense %>% mutate(i = row_number()) %>% select(Codon, i)

edges <- list()
for (pos in 1:3) {
  from_base <- substr(SGC_sense$Codon, pos, pos)
  for (b in bases) {
    mask <- from_base != b
    if (any(mask)) {
      edges[[length(edges) + 1]] <- tibble(
        from = SGC_sense$Codon[mask],
        to   = mutate_pos(SGC_sense$Codon[mask], pos, b),
        pos  = pos
      )
    }
  }
}
edges_df <- bind_rows(edges)
edges_idx <- edges_df %>%
  inner_join(lookup, by = c("from" = "Codon")) %>%
  inner_join(lookup %>% rename(j = i), by = c("to" = "Codon"))

edges_labeled <- edges_idx %>%
  mutate(
    b_from = substr(from, pos, pos),
    b_to   = substr(to, pos, pos),
    ts     = ((b_from %in% c("A","G")) & (b_to %in% c("A","G")) & (b_from != b_to)) |
      ((b_from %in% c("C","T")) & (b_to %in% c("C","T")) & (b_from != b_to))
  )

A_ids <- SGC_sense$A
ai <- A_ids[edges_labeled$i]
aj <- A_ids[edges_labeled$j]
ri <- match(ai, rownames(aa_dist))
cj <- match(aj, colnames(aa_dist))
edges_labeled$dist <- aa_dist[cbind(ri, cj)]

D_ts <- sum(edges_labeled$dist[edges_labeled$ts],  na.rm = TRUE)
D_tv <- sum(edges_labeled$dist[!edges_labeled$ts], na.rm = TRUE)

SGC_target <- 9856.116
alpha <- (SGC_target - D_ts) / D_tv
SGC_cost <- D_ts + alpha * D_tv
cat(sprintf("✅ Benchmark SGC cost = %.3f (target = %.3f)\n", SGC_cost, SGC_target))

# Benchmark Computation (SGC vs Random Codes)

SGC_sense <- SGC %>%
  filter(toupper(AminoAcid) != "STOP") %>%
  mutate(Codon = toupper(str_replace_all(Codon, "\s+", "")),
         Codon = chartr("U", "T", Codon))

bases <- c("A","C","G","T")
mutate_pos <- function(codons, pos, to_base) {
  paste0(
    if (pos == 1) to_base else substr(codons, 1, 1),
    if (pos == 2) to_base else substr(codons, 2, 2),
    if (pos == 3) to_base else substr(codons, 3, 3)
  )
}

lookup <- SGC_sense %>% mutate(i = row_number()) %>% select(Codon, i)

edges <- list()
for (pos in 1:3) {
  from_base <- substr(SGC_sense$Codon, pos, pos)
  for (b in bases) {
    mask <- from_base != b
    if (any(mask)) {
      edges[[length(edges) + 1]] <- tibble(
        from = SGC_sense$Codon[mask],
        to   = mutate_pos(SGC_sense$Codon[mask], pos, b),
        pos  = pos
      )
    }
  }
}
edges_df <- bind_rows(edges)
edges_idx <- edges_df %>%
  inner_join(lookup, by = c("from" = "Codon")) %>%
  inner_join(lookup %>% rename(j = i), by = c("to" = "Codon"))

edges_labeled <- edges_idx %>%
  mutate(
    b_from = substr(from, pos, pos),
    b_to   = substr(to, pos, pos),
    ts     = ((b_from %in% c("A","G")) & (b_to %in% c("A","G")) & (b_from != b_to)) |
      ((b_from %in% c("C","T")) & (b_to %in% c("C","T")) & (b_from != b_to))
  )

A_ids <- SGC_sense$A
ai <- A_ids[edges_labeled$i]
aj <- A_ids[edges_labeled$j]
ri <- match(ai, rownames(aa_dist))
cj <- match(aj, colnames(aa_dist))
edges_labeled$dist <- aa_dist[cbind(ri, cj)]

D_ts <- sum(edges_labeled$dist[edges_labeled$ts],  na.rm = TRUE)
D_tv <- sum(edges_labeled$dist[!edges_labeled$ts], na.rm = TRUE)

SGC_target <- 9856.116
alpha <- (SGC_target - D_ts) / D_tv
SGC_cost <- D_ts + alpha * D_tv
cat(sprintf("✅ Benchmark SGC cost = %.3f (target = %.3f)\n", SGC_cost, SGC_target))

# Compare with Random Protein Sequencing Codes

set.seed(42)
n_random <- 1000
random_costs <- numeric(n_random)

edge_i <- edges_labeled$i
edge_j <- edges_labeled$j
edge_is_ts <- edges_labeled$ts

for (k in seq_len(n_random)) {
  perm <- sample(seq_along(A_ids))
  A_perm <- A_ids[perm]
  ai_p <- A_perm[edge_i]
  aj_p <- A_perm[edge_j]
  ri_p <- match(ai_p, rownames(aa_dist))
  cj_p <- match(aj_p, colnames(aa_dist))
  d_perm <- aa_dist[cbind(ri_p, cj_p)]
  Dts_p <- sum(d_perm[edge_is_ts], na.rm = TRUE)
  Dtv_p <- sum(d_perm[!edge_is_ts], na.rm = TRUE)
  random_costs[k] <- Dts_p + alpha * Dtv_p
}

hist(random_costs, breaks = 20, col = "lightblue", main = "Genetic Code Cost Distribution")
abline(v = SGC_cost, col = "red", lwd = 2)

#PCA & Clustering

aa_matrix <- as.data.frame(aaFeatureSpace.4.1)
aa_matrix <- aa_matrix[sapply(aa_matrix, is.numeric)]
pca_res <- prcomp(aa_matrix, scale. = TRUE)

fviz_pca_ind(pca_res, geom.ind = "point", fill.ind = "steelblue",
             pointsize = 3, repel = TRUE) +
  ggtitle("PCA of Standard Genetic Code Amino Acid Features") +
  theme_minimal()
hc_res <- hclust(dist(aa_matrix), method = "ward.D2")
fviz_dend(hc_res, k = 5, rect = TRUE) +
  ggtitle("Hierarchical Clustering of Amino Acids") +
  theme_minimal()

#Random Protein PCA Comparison

set.seed(42)
random_codons <- sample(SGC$Codon, size = 100000, replace = TRUE)
random_proteins <- tibble(RandomCodon = random_codons,
                          RandomAA = SGC$AminoAcid[match(random_codons, SGC$Codon)])

random_features <- matrix(rnorm(nrow(random_proteins) * ncol(aa_matrix)),
                          ncol = ncol(aa_matrix))
colnames(random_features) <- colnames(aa_matrix)
random_pca <- prcomp(random_features, scale. = TRUE)

std_scores <- as_tibble(pca_res$x[, 1:2]) %>% mutate(Type = "Standard Genetic Code")
rnd_scores <- as_tibble(random_pca$x[, 1:2]) %>% mutate(Type = "Random Code")

combined <- bind_rows(std_scores, rnd_scores)
ggplot(combined, aes(PC1, PC2, color = Type)) +
  geom_point(alpha = 0.6) +
  ggtitle("PCA Comparison: Standard vs Random Codes") +
  theme_minimal()


