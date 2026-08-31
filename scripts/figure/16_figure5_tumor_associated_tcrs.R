# ==============================================================================
# Script: 16_figure5_tumor_associated_tcrs.R
# Purpose: Generate Figure 5: Prioritization of Tumor-Enriched TCR Candidates
#
# Statistical design:
#   1) Tumor-enrichment tests are performed WITHIN each patient.
#   2) Primary comparison is matched Tumor vs NAT; Blood is not pooled into the
#      reference compartment.
#   3) Fisher exact P values are BH-adjusted across all tested clonotypes.
#   4) Statistical enrichment requires FDR < 0.05 and log2 OR > 1.
#   5) A separate abundance filter defines the experimental-priority shortlist.
#
# Biological claim discipline:
#   - TCRex output is treated as a database prediction, not functional proof.
#   - Database-unmapped tumor-enriched clones are "candidates", not proven
#     tumor-reactive TCRs.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(scales)
  library(stringr)
})

# Locate shared helper script under either common project layout.
.common_candidates <- c(
  "scripts/figure/11_figure_common.R",
  "scripts/11_figure_common.R",
  "11_figure_common.R"
)
.common_path <- .common_candidates[file.exists(.common_candidates)][1]
if (is.na(.common_path)) {
  stop(
    "Cannot find 11_figure_common.R. Expected one of: ",
    paste(.common_candidates, collapse = ", ")
  )
}
source(.common_path)
rm(.common_candidates, .common_path)

cat("\n============================================================\n")
cat(" Figure 5: Tumor-Enriched TCR Candidate Prioritization\n")
cat("============================================================\n\n")

# 1. Input / Output Paths & Analysis Parameters --------------------------------
master_path <- "data/processed/final_paired_tcell_metadata.rds"
tcrex_path  <- "data/processed/tcrex_results_score_gt_0.5.csv"
figure_dir  <- "results/figures"
table_dir   <- "results/tables"

TCREX_SCORE_THRESHOLD <- 0.5

# Statistical evaluability threshold: clones with fewer than this many cells
# across matched Tumor + NAT are not tested.
MIN_TN_CLONE_SIZE_FOR_TEST <- 3

# Statistical tumor-enrichment thresholds.
LOG2_OR_THRESHOLD <- 1
FDR_THRESHOLD <- 0.05

# Additional abundance thresholds used only for the experimental-priority
# shortlist. These do NOT change whether a clone is statistically enriched.
MIN_TUMOR_CELLS_FOR_PRIORITY <- 5
MIN_GLOBAL_CLONE_SIZE_FOR_PRIORITY <- 10

# Extremely small FDR values can compress the informative part of a volcano
# plot. Values above this -log10(FDR) are visually capped only for plotting.
VOLCANO_FDR_CAP <- 20

# Optional explicit override. Preferred for a final reproducible project.
# Example before running:
#   Sys.setenv(TCR_BETA_COLUMN = "CDR3b")
TCR_BETA_COLUMN <- Sys.getenv("TCR_BETA_COLUMN", unset = "")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(master_path)) {
  stop("Master dataset not found: ", master_path)
}
if (!file.exists(tcrex_path)) {
  stop("TCRex result file not found: ", tcrex_path)
}

master_df <- readRDS(master_path)
cat("Master dataset loaded. Total cells:", nrow(master_df), "\n")

# 2. Deterministic TCR-beta Column Resolution ----------------------------------
# Prefer beta-chain-specific names. A generic CDR3 field is accepted only with
# an explicit warning, rather than silently taking the first CDR3-like column.
resolve_tcr_beta_column <- function(df, explicit_column = "") {
  if (nzchar(explicit_column)) {
    if (!explicit_column %in% colnames(df)) {
      stop("Configured TCR_BETA_COLUMN not found in master metadata: ", explicit_column)
    }
    return(explicit_column)
  }
  
  beta_specific <- c(
    "CDR3b", "cdr3_b", "cdr3_beta", "cdr3b", "TRB_CDR3",
    "TRB_cdr3", "TRB_cdr3_aa", "junction_aa_beta", "cdr3_beta_aa"
  )
  beta_hits <- beta_specific[beta_specific %in% colnames(df)]
  
  if (length(beta_hits) == 1) {
    return(beta_hits[[1]])
  }
  if (length(beta_hits) > 1) {
    stop(
      "Multiple beta-chain CDR3 columns detected: ", paste(beta_hits, collapse = ", "),
      ". Set environment variable TCR_BETA_COLUMN explicitly."
    )
  }
  
  generic_candidates <- c("cdr3", "cdr3_aa", "junction_aa", "CDR3.aa")
  generic_hits <- generic_candidates[generic_candidates %in% colnames(df)]
  
  if (length(generic_hits) == 1) {
    warning(
      "No beta-specific CDR3 column name was found. Using generic column '",
      generic_hits[[1]],
      "'. Verify that this field is the TCR-beta CDR3 and set TCR_BETA_COLUMN explicitly for the final analysis."
    )
    return(generic_hits[[1]])
  }
  
  stop(
    "Unable to resolve a unique TCR-beta CDR3 column. Set TCR_BETA_COLUMN explicitly."
  )
}

tcr_beta_col <- resolve_tcr_beta_column(master_df, TCR_BETA_COLUMN)
cat("TCRex annotation will use beta-chain CDR3 column:", tcr_beta_col, "\n")

# 3. Robust but Auditable TCRex Loader -----------------------------------------
# TCRex exports can differ by version / download route (TSV, CSV, semicolon,
# whitespace-delimited, and occasionally weak/non-standard headers). We first
# resolve fields from explicit column names. Only if that fails do we use a
# content-based fallback with strict validation, and the detected mapping is
# printed to the console so the analysis never silently guesses columns.

normalize_colname <- function(x) {
  x %>%
    str_replace_all("\\ufeff", "") %>%
    str_trim() %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
}

resolve_named_column <- function(tbl, aliases, label, required = TRUE) {
  normalized_names <- normalize_colname(colnames(tbl))
  normalized_aliases <- normalize_colname(aliases)
  idx <- match(normalized_aliases, normalized_names, nomatch = 0)
  idx <- idx[idx > 0]
  
  if (length(idx) == 0) {
    if (required) {
      stop(
        "Could not find TCRex ", label, " column. Expected one of: ",
        paste(aliases, collapse = ", "),
        ". Available columns: ", paste(colnames(tbl), collapse = ", ")
      )
    }
    return(NA_character_)
  }
  
  colnames(tbl)[idx[[1]]]
}

standardize_tcrex_from_mapping <- function(tbl, cdr3_col, pathology_col, score_col, epitope_col = NA_character_) {
  out <- tibble(
    cdr3_aa = str_trim(as.character(tbl[[cdr3_col]])),
    pathology = if (!is.na(pathology_col)) str_trim(as.character(tbl[[pathology_col]])) else "Unspecified TCRex prediction",
    epitope = if (!is.na(epitope_col)) str_trim(as.character(tbl[[epitope_col]])) else NA_character_,
    score = suppressWarnings(as.numeric(as.character(tbl[[score_col]])))
  ) %>%
    filter(
      !is.na(cdr3_aa), cdr3_aa != "",
      !is.na(score), score > TCREX_SCORE_THRESHOLD
    ) %>%
    distinct(cdr3_aa, pathology, epitope, score, .keep_all = TRUE) %>%
    mutate(
      pathology_group = case_when(
        str_detect(pathology, regex("Melanoma|Tumou?r|Cancer|Myeloma|WT1|Leukemia|Lymphoma", ignore_case = TRUE)) ~
          "Tumor-associated TCRex prediction",
        str_detect(pathology, regex("Influenza|SARS|COVID|EBV|CMV|HIV|HCV|HBV|DENV|Yellow.?Fever|HSV|VZV", ignore_case = TRUE)) ~
          "Viral-associated TCRex prediction",
        TRUE ~ "Other TCRex prediction"
      )
    )
  
  if (nrow(out) == 0) {
    stop("TCRex table parsed successfully, but no valid predictions remained above score threshold > ", TCREX_SCORE_THRESHOLD)
  }
  out
}

standardize_tcrex_named <- function(tbl) {
  cdr3_col <- resolve_named_column(
    tbl,
    c(
      "cdr3", "cdr3_aa", "cdr3 sequence", "cdr3_sequence", "cdr3.beta", "cdr3b",
      "cdr3_beta", "cdr3b_aa", "trb_cdr3", "trb_cdr3_aa", "junction_aa"
    ),
    "CDR3"
  )
  pathology_col <- resolve_named_column(
    tbl,
    c(
      "pathology", "disease", "condition", "pathology_name", "disease_name",
      "specificity", "antigen", "antigen_source", "source_pathology"
    ),
    "pathology"
  )
  score_col <- resolve_named_column(
    tbl,
    c(
      "score", "prediction_score", "probability", "tcrex_score", "prediction probability",
      "prediction_probability", "tcr_probability", "tcr_probability_score"
    ),
    "score"
  )
  epitope_col <- resolve_named_column(
    tbl,
    c("epitope", "peptide", "antigen_epitope", "peptide_sequence", "epitope_sequence"),
    "epitope",
    required = FALSE
  )
  
  cat(
    "TCRex named-column mapping: CDR3=", cdr3_col,
    "; pathology=", pathology_col,
    "; score=", score_col,
    if (!is.na(epitope_col)) paste0("; epitope=", epitope_col) else "; epitope=<not available>",
    "\n",
    sep = ""
  )
  
  standardize_tcrex_from_mapping(tbl, cdr3_col, pathology_col, score_col, epitope_col)
}

infer_tcrex_columns <- function(tbl) {
  if (ncol(tbl) < 3) {
    stop("Content-based fallback requires at least 3 parsed columns.")
  }
  
  char_cols <- lapply(tbl, function(x) str_trim(as.character(x)))
  nonempty <- lapply(char_cols, function(x) x[!is.na(x) & x != ""])
  
  # Beta-chain CDR3 amino-acid strings: strongly enriched for CASS/CS motifs and
  # amino-acid-only sequences of plausible length. We choose a column only when
  # the evidence is unique enough to be auditable.
  cdr3_strength <- vapply(nonempty, function(v) {
    if (length(v) == 0) return(0)
    mean(str_detect(v, "^C[A-Z]{6,30}$") & str_detect(v, "^(CAS|CS|CA)"))
  }, numeric(1))
  
  cdr3_candidates <- which(cdr3_strength >= 0.20)
  if (length(cdr3_candidates) == 0) {
    # Small files can have only a handful of predictions; permit >=2 strong hits.
    cdr3_hits <- vapply(nonempty, function(v) sum(str_detect(v, "^C[A-Z]{6,30}$") & str_detect(v, "^(CAS|CS|CA)")), numeric(1))
    cdr3_candidates <- which(cdr3_hits >= 2)
  }
  if (length(cdr3_candidates) != 1) {
    stop(
      "Unable to uniquely infer the TCRex CDR3 column from content. Candidate columns: ",
      if (length(cdr3_candidates)) paste(colnames(tbl)[cdr3_candidates], collapse = ", ") else "none"
    )
  }
  cdr3_idx <- cdr3_candidates[[1]]
  
  # Score: predominantly numeric and bounded in [0, 1].
  score_strength <- vapply(nonempty, function(v) {
    if (length(v) == 0) return(0)
    z <- suppressWarnings(as.numeric(v))
    valid <- !is.na(z)
    if (mean(valid) < 0.7) return(0)
    mean(z[valid] >= 0 & z[valid] <= 1)
  }, numeric(1))
  score_candidates <- setdiff(which(score_strength >= 0.90), cdr3_idx)
  if (length(score_candidates) != 1) {
    stop(
      "Unable to uniquely infer the TCRex score column from content. Candidate columns: ",
      if (length(score_candidates)) paste(colnames(tbl)[score_candidates], collapse = ", ") else "none"
    )
  }
  score_idx <- score_candidates[[1]]
  
  # Pathology/disease: prefer a column containing recognizable disease/source
  # terms. This is used only for annotation categories, not the enrichment test.
  pathology_regex <- regex(
    "Influenza|SARS|COVID|EBV|CMV|HIV|HCV|HBV|DENV|Yellow.?Fever|HSV|VZV|Melanoma|Tumou?r|Cancer|Myeloma|WT1|Leukemia|Lymphoma",
    ignore_case = TRUE
  )
  pathology_hits <- vapply(nonempty, function(v) sum(str_detect(v, pathology_regex)), numeric(1))
  pathology_hits[c(cdr3_idx, score_idx)] <- -1
  pathology_idx <- which.max(pathology_hits)
  if (pathology_hits[pathology_idx] <= 0) {
    stop(
      "Unable to infer a pathology/disease column from TCRex content. ",
      "Please keep the original TCRex header or rename the relevant column to 'pathology'."
    )
  }
  
  # Epitope is optional. Choose a peptide-like amino-acid column, excluding CDR3.
  epitope_strength <- vapply(nonempty, function(v) {
    if (length(v) == 0) return(0)
    mean(str_detect(v, "^[ACDEFGHIKLMNPQRSTVWY]{8,15}$"))
  }, numeric(1))
  epitope_strength[c(cdr3_idx, score_idx, pathology_idx)] <- 0
  epitope_idx <- if (max(epitope_strength, na.rm = TRUE) >= 0.20) which.max(epitope_strength) else NA_integer_
  
  cdr3_col <- colnames(tbl)[cdr3_idx]
  score_col <- colnames(tbl)[score_idx]
  pathology_col <- colnames(tbl)[pathology_idx]
  epitope_col <- if (!is.na(epitope_idx)) colnames(tbl)[epitope_idx] else NA_character_
  
  cat(
    "TCRex validated content-based mapping: CDR3=", cdr3_col,
    "; pathology=", pathology_col,
    "; score=", score_col,
    if (!is.na(epitope_col)) paste0("; epitope=", epitope_col) else "; epitope=<not available>",
    "\n",
    sep = ""
  )
  
  standardize_tcrex_from_mapping(tbl, cdr3_col, pathology_col, score_col, epitope_col)
}

# Native TCRex text export parser ------------------------------------------------
# Some TCRex exports use a whitespace-separated header but tab-separated data.
# readr::read_delim/read_table can interpret this inconsistently. We therefore
# recognize the canonical 7-column TCRex layout directly from raw lines and
# split on ANY run of whitespace (spaces or tabs). Pathology is allowed to be
# multi-word by anchoring the first four and last two fields.
read_tcrex_native_export <- function(path) {
  lines <- read_lines(path)
  lines <- lines[str_trim(lines) != ""]
  
  if (length(lines) < 2) {
    stop("TCRex file contains fewer than two non-empty lines.")
  }
  
  header_tokens <- str_split(str_trim(lines[[1]]), "\\s+")[[1]]
  header_norm <- normalize_colname(header_tokens)
  canonical_norm <- normalize_colname(c("TRBV", "CDR3", "TRBJ", "Epitope", "Pathology", "Score", "BPR"))
  
  if (!all(canonical_norm %in% header_norm)) {
    stop(
      "Raw TCRex header does not match the canonical TRBV/CDR3/TRBJ/Epitope/Pathology/Score/BPR layout. Header: ",
      paste(header_tokens, collapse = " | ")
    )
  }
  
  parsed_rows <- lapply(seq_along(lines[-1]), function(i) {
    line_no <- i + 1
    tokens <- str_split(str_trim(lines[[line_no]]), "\\s+")[[1]]
    
    # Canonical rows contain at least seven fields:
    # TRBV | CDR3 | TRBJ | Epitope | Pathology [possibly multi-word] | Score | BPR
    if (length(tokens) < 7) {
      stop(
        "Malformed TCRex row ", line_no, ": expected >=7 whitespace-delimited fields, found ",
        length(tokens), ". Row: ", lines[[line_no]]
      )
    }
    
    pathology_start <- 5
    pathology_end <- length(tokens) - 2
    
    tibble(
      TRBV = tokens[[1]],
      CDR3 = tokens[[2]],
      TRBJ = tokens[[3]],
      Epitope = tokens[[4]],
      Pathology = paste(tokens[pathology_start:pathology_end], collapse = " "),
      Score = tokens[[length(tokens) - 1]],
      BPR = tokens[[length(tokens)]]
    )
  })
  
  out <- bind_rows(parsed_rows)
  
  if (nrow(out) == 0) {
    stop("Canonical TCRex raw parser produced zero rows.")
  }
  
  cat(
    "Detected native TCRex mixed-whitespace export: ", nrow(out),
    " raw prediction rows.\n",
    sep = ""
  )
  
  out
}

read_tcrex_table <- function(path) {
  # First handle the native TCRex export observed in this project: the header
  # is space-separated while data rows are tab-separated. Splitting raw lines
  # on arbitrary whitespace is deterministic for this canonical layout.
  native <- tryCatch(
    read_tcrex_native_export(path),
    error = function(e) NULL
  )
  if (!is.null(native)) {
    standardized <- tryCatch(
      standardize_tcrex_named(native),
      error = function(e) {
        stop("Native TCRex export was parsed, but semantic mapping failed: ", conditionMessage(e))
      }
    )
    cat(
      "TCRex table parsed using native mixed-whitespace parser with ",
      nrow(standardized), " predictions above score threshold.\n",
      sep = ""
    )
    return(standardized)
  }
  
  # Fallback: try common TCRex/download formats, with and without headers. Parsing success
  # requires >=3 columns; semantic mapping then has to pass named or validated
  # content-based resolution.
  parser_specs <- list(
    TSV = list(type = "delim", delim = "\t"),
    CSV = list(type = "delim", delim = ","),
    SEMICOLON = list(type = "delim", delim = ";"),
    PIPE = list(type = "delim", delim = "|"),
    WHITESPACE = list(type = "table", delim = NA_character_)
  )
  
  errors <- character()
  
  for (has_header in c(TRUE, FALSE)) {
    for (parser_name in names(parser_specs)) {
      spec <- parser_specs[[parser_name]]
      parsed <- tryCatch({
        if (spec$type == "delim") {
          read_delim(
            path,
            delim = spec$delim,
            col_names = has_header,
            trim_ws = TRUE,
            show_col_types = FALSE,
            progress = FALSE
          )
        } else {
          read_table(
            path,
            col_names = has_header,
            trim_ws = TRUE,
            show_col_types = FALSE,
            progress = FALSE
          )
        }
      }, error = function(e) NULL)
      
      attempt_label <- paste0(parser_name, if (has_header) "[header]" else "[no-header]")
      
      if (is.null(parsed) || ncol(parsed) < 3 || nrow(parsed) == 0) {
        errors <- c(errors, paste0(attempt_label, ": parse failed / <3 columns / 0 rows"))
        next
      }
      
      # Prefer explicit column names whenever a header was read.
      if (has_header) {
        named <- tryCatch(
          standardize_tcrex_named(parsed),
          error = function(e) {
            errors <<- c(errors, paste0(attempt_label, " named mapping: ", conditionMessage(e)))
            NULL
          }
        )
        if (!is.null(named)) {
          cat("TCRex table parsed as ", attempt_label, " with ", nrow(named), " predictions above score threshold.\n", sep = "")
          return(named)
        }
      }
      
      inferred <- tryCatch(
        infer_tcrex_columns(parsed),
        error = function(e) {
          errors <<- c(errors, paste0(attempt_label, " content mapping: ", conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(inferred)) {
        cat("TCRex table parsed as ", attempt_label, " with ", nrow(inferred), " predictions above score threshold.\n", sep = "")
        return(inferred)
      }
    }
  }
  
  # Give an immediately useful diagnostic instead of a generic parsing error.
  preview <- tryCatch(read_lines(path, n_max = 3), error = function(e) character())
  stop(
    "Unable to parse/map the TCRex table after trying TSV, CSV, semicolon, pipe, and whitespace formats.\n",
    "First non-empty lines of the file:\n",
    paste(preview, collapse = "\n"),
    "\n\nDetailed attempts:\n- ", paste(errors, collapse = "\n- ")
  )
}

tcrex_clean <- read_tcrex_table(tcrex_path)

# 4. Clean Single-Cell Metadata & Build Canonical Clone Metadata ---------------
required_columns <- c("patient", "source", "patient_clonotype_clean", tcr_beta_col)
missing_cols <- setdiff(required_columns, colnames(master_df))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

df <- master_df %>%
  # Remove any upstream derived columns that are recomputed canonically below.
  # This prevents dplyr joins from creating .x/.y suffixes and silently drifting definitions.
  select(-any_of(c("expansion_status", "global_expansion_status", "global_clone_size", "is_clonally_expanded"))) %>%
  mutate(
    patient = as.character(patient),
    source = as.character(source),
    patient_clonotype_clean = as.character(patient_clonotype_clean),
    cdr3_beta = as.character(.data[[tcr_beta_col]])
  ) %>%
  filter(
    !is.na(patient),
    !is.na(source),
    source %in% TISSUE_LEVELS,
    !is.na(patient_clonotype_clean),
    patient_clonotype_clean != "",
    patient_clonotype_clean != "NA"
  ) %>%
  mutate(source = factor(source, levels = TISSUE_LEVELS))

clone_cdr3_check <- df %>%
  group_by(patient, patient_clonotype_clean) %>%
  summarise(
    n_beta_cdr3 = n_distinct(cdr3_beta[!is.na(cdr3_beta) & cdr3_beta != ""]),
    .groups = "drop"
  )

if (any(clone_cdr3_check$n_beta_cdr3 > 1)) {
  stop(
    "Detected patient-clonotypes mapped to multiple beta-chain CDR3 sequences. ",
    "Resolve upstream clonotype construction before TCRex annotation."
  )
}

clone_meta <- df %>%
  group_by(patient, patient_clonotype_clean) %>%
  summarise(
    global_clone_size = n(),
    cdr3_beta = {
      x <- unique(cdr3_beta[!is.na(cdr3_beta) & cdr3_beta != ""])
      if (length(x) == 0) NA_character_ else x[[1]]
    },
    .groups = "drop"
  ) %>%
  mutate(global_expansion_status = assign_expansion_tier(global_clone_size))

# 5. TCRex Annotation: Highest-Score Row Kept Consistently ---------------------
tcrex_annotation <- clone_meta %>%
  select(patient, patient_clonotype_clean, cdr3_beta) %>%
  left_join(tcrex_clean, by = c("cdr3_beta" = "cdr3_aa")) %>%
  group_by(patient, patient_clonotype_clean) %>%
  arrange(desc(score), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    patient,
    patient_clonotype_clean,
    is_tcrex_matched = !is.na(score),
    top_pathology_group = ifelse(is_tcrex_matched, pathology_group, "Unmapped"),
    top_pathology = ifelse(is_tcrex_matched, pathology, "Unmapped"),
    top_epitope = ifelse(is_tcrex_matched, epitope, NA_character_),
    top_tcrex_score = score
  )

# 6. Within-Patient Tumor vs NAT Fisher Exact Tests ----------------------------
# The tested clone is compared with all OTHER T cells from the SAME patient,
# using matched Tumor and NAT compartments only.
patient_tn_totals <- df %>%
  filter(source %in% c("Tumor", "NAT")) %>%
  count(patient, source, name = "total_cells") %>%
  pivot_wider(names_from = source, values_from = total_cells, values_fill = 0) %>%
  filter(Tumor > 0, NAT > 0) %>%
  rename(patient_tumor_total = Tumor, patient_nat_total = NAT)

clone_tn_counts <- df %>%
  filter(patient %in% patient_tn_totals$patient, source %in% c("Tumor", "NAT")) %>%
  count(patient, patient_clonotype_clean, source, name = "cell_count") %>%
  pivot_wider(names_from = source, values_from = cell_count, values_fill = 0) %>%
  mutate(tn_clone_size = Tumor + NAT) %>%
  left_join(patient_tn_totals, by = "patient") %>%
  left_join(
    clone_meta %>% select(patient, patient_clonotype_clean, global_clone_size, global_expansion_status),
    by = c("patient", "patient_clonotype_clean")
  )

run_within_patient_fisher <- function(a, b, c, d) {
  #               Clone   Other cells
  # Tumor           a         c
  # NAT             b         d
  mat <- matrix(c(a, c, b, d), nrow = 2, byrow = TRUE)
  res <- fisher.test(mat)
  tibble(
    p_value = res$p.value,
    fisher_odds_ratio = unname(res$estimate[[1]])
  )
}

tested_clones <- clone_tn_counts %>%
  filter(tn_clone_size >= MIN_TN_CLONE_SIZE_FOR_TEST) %>%
  mutate(
    a = Tumor,
    b = NAT,
    c = patient_tumor_total - Tumor,
    d = patient_nat_total - NAT
  )

if (nrow(tested_clones) > 0) {
  fisher_results <- bind_rows(lapply(seq_len(nrow(tested_clones)), function(i) {
    run_within_patient_fisher(
      tested_clones$a[i], tested_clones$b[i],
      tested_clones$c[i], tested_clones$d[i]
    )
  }))
  
  tested_clones <- tested_clones %>%
    bind_cols(fisher_results) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      log2_odds_ratio = log2_or_ha(a, b, c, d),
      neg_log10_fdr = -log10(pmax(p_adj, 1e-300)),
      is_tumor_enriched =
        p_adj < FDR_THRESHOLD &
        log2_odds_ratio > LOG2_OR_THRESHOLD
    )
} else {
  tested_clones <- tested_clones %>%
    mutate(
      p_value = numeric(), fisher_odds_ratio = numeric(), p_adj = numeric(),
      log2_odds_ratio = numeric(), neg_log10_fdr = numeric(),
      is_tumor_enriched = logical()
    )
}

# 7. Unified Clone-Level Classification ----------------------------------------
clonotype_results <- clone_meta %>%
  left_join(
    tcrex_annotation,
    by = c("patient", "patient_clonotype_clean")
  ) %>%
  left_join(
    tested_clones %>%
      select(
        patient, patient_clonotype_clean,
        Tumor, NAT, tn_clone_size,
        p_value, p_adj, log2_odds_ratio, neg_log10_fdr,
        is_tumor_enriched
      ),
    by = c("patient", "patient_clonotype_clean")
  ) %>%
  mutate(
    was_tested = !is.na(p_adj),
    is_tumor_enriched = replace_na(is_tumor_enriched, FALSE),
    
    # Statistical enrichment and experimental prioritization are deliberately
    # separated. A clone can be statistically tumor-enriched yet remain below
    # the abundance threshold chosen for a practical validation shortlist.
    passes_priority_abundance =
      was_tested &
      replace_na(Tumor, 0) >= MIN_TUMOR_CELLS_FOR_PRIORITY &
      global_clone_size >= MIN_GLOBAL_CLONE_SIZE_FOR_PRIORITY,
    is_prioritized_candidate =
      is_tumor_enriched & passes_priority_abundance,
    
    candidate_category = case_when(
      !was_tested ~ "Not tested for tumor enrichment",
      is_tumor_enriched & !is_tcrex_matched ~
        "Tumor-enriched, TCRex-unmapped",
      is_tumor_enriched & top_pathology_group == "Tumor-associated TCRex prediction" ~
        "Tumor-enriched + TCRex tumor-associated",
      is_tumor_enriched & top_pathology_group == "Viral-associated TCRex prediction" ~
        "Tumor-enriched + TCRex viral-associated",
      is_tumor_enriched & is_tcrex_matched ~
        "Tumor-enriched + other TCRex prediction",
      TRUE ~ "Tested, not tumor-enriched"
    ),
    candidate_category = factor(
      candidate_category,
      levels = c(
        "Tumor-enriched, TCRex-unmapped",
        "Tumor-enriched + TCRex tumor-associated",
        "Tumor-enriched + TCRex viral-associated",
        "Tumor-enriched + other TCRex prediction",
        "Tested, not tumor-enriched",
        "Not tested for tumor enrichment"
      )
    )
  )

novel_candidates <- clonotype_results %>%
  filter(is_prioritized_candidate, !is_tcrex_matched) %>%
  arrange(p_adj, desc(global_clone_size))

write_csv(
  novel_candidates,
  file.path(table_dir, "Figure5_tumor_enriched_tcrex_unmapped_candidates.csv")
)
write_csv(
  clonotype_results,
  file.path(table_dir, "Figure5_tcrex_annotated_clonotypes.csv")
)

cat("\nSummary of clone classification:\n")
print(table(clonotype_results$candidate_category, useNA = "ifany"))
cat("\n")

# 8. Panel A: Within-Patient Tumor-Enrichment Volcano Plot ---------------------
volcano_df <- clonotype_results %>%
  filter(was_tested) %>%
  mutate(
    volcano_category = case_when(
      is_tumor_enriched & !is_tcrex_matched ~ "TCRex-unmapped",
      is_tumor_enriched & top_pathology_group == "Tumor-associated TCRex prediction" ~ "TCRex tumor-associated",
      is_tumor_enriched & top_pathology_group == "Viral-associated TCRex prediction" ~ "TCRex viral-associated",
      is_tumor_enriched & is_tcrex_matched ~ "Other TCRex prediction",
      TRUE ~ "Not significant"
    ),
    plot_neg_log10_fdr = pmin(neg_log10_fdr, VOLCANO_FDR_CAP)
  )

volcano_colors <- c(
  "TCRex-unmapped" = TISSUE_COLORS[["Tumor"]],
  "TCRex tumor-associated" = "#e7298a",
  "TCRex viral-associated" = "#7570b3",
  "Other TCRex prediction" = "#1b9e77",
  "Not significant" = "grey82"
)

p5a <- ggplot(volcano_df, aes(x = log2_odds_ratio, y = plot_neg_log10_fdr, color = volcano_category)) +
  geom_vline(
    xintercept = c(-LOG2_OR_THRESHOLD, LOG2_OR_THRESHOLD),
    linetype = "dashed", color = "grey50", linewidth = 0.5
  ) +
  geom_hline(
    yintercept = -log10(FDR_THRESHOLD),
    linetype = "dashed", color = "grey50", linewidth = 0.5
  ) +
  geom_point(aes(size = global_clone_size), alpha = 0.82) +
  scale_color_manual(
    values = volcano_colors,
    name = "TCR/TCRex category",
    drop = TRUE
  ) +
  scale_size_continuous(range = c(1.5, 6.5), name = "Global Clone Size") +
  scale_y_continuous(
    limits = c(0, VOLCANO_FDR_CAP),
    breaks = seq(0, VOLCANO_FDR_CAP, by = 5),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  labs(
    title = "A. Within-patient tumor enrichment and TCRex annotation",
    subtitle = paste0(
      "Tumor vs matched NAT Fisher tests; BH-FDR; y capped at ", VOLCANO_FDR_CAP
    ),
    x = "Log2 odds ratio (Tumor vs NAT)",
    y = paste0("-Log10(BH-FDR), capped at ", VOLCANO_FDR_CAP)
  ) +
  project_theme() +
  theme(legend.position = "right")

# 9. Panel B: Top Tumor-Enriched Clones ----------------------------------------
top20_tumor_clones <- clonotype_results %>%
  filter(is_prioritized_candidate) %>%
  arrange(p_adj, desc(global_clone_size)) %>%
  slice_head(n = 20) %>%
  mutate(
    clonotype_label = paste0(patient, ": ", substring(patient_clonotype_clean, 1, 14)),
    # Keep plot labels short; full epitope/score information remains in the CSV.
    annotation_tag = case_when(
      !is_tcrex_matched ~ "Unmapped",
      TRUE ~ top_pathology
    ),
    bar_category = case_when(
      !is_tcrex_matched ~ "TCRex-unmapped",
      top_pathology_group == "Tumor-associated TCRex prediction" ~ "TCRex tumor-associated",
      top_pathology_group == "Viral-associated TCRex prediction" ~ "TCRex viral-associated",
      TRUE ~ "Other TCRex prediction"
    )
  )

p5b <- ggplot(
  top20_tumor_clones,
  aes(x = reorder(clonotype_label, global_clone_size), y = global_clone_size, fill = bar_category)
) +
  geom_col(color = "black", linewidth = 0.25, width = 0.72) +
  geom_text(
    aes(label = annotation_tag),
    hjust = -0.06, size = 2.5, fontface = "bold", color = "grey20"
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = volcano_colors[names(volcano_colors) != "Not significant"],
    name = "TCR/TCRex category",
    drop = TRUE
  ) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.28))) +
  labs(
    title = "B. Top prioritized tumor-enriched clonotypes",
    subtitle = paste0(
      "Priority: BH-FDR < ", FDR_THRESHOLD,
      ", log2 OR > ", LOG2_OR_THRESHOLD,
      ", Tumor cells >= ", MIN_TUMOR_CELLS_FOR_PRIORITY,
      ", global clone size >= ", MIN_GLOBAL_CLONE_SIZE_FOR_PRIORITY
    ),
    x = "Tumor-enriched patient-clonotypes",
    y = "Global clone size (cell count)"
  ) +
  project_theme() +
  theme(
    legend.position = "right",
    plot.margin = margin(5.5, 32, 5.5, 5.5)
  )

# 10. Panel C: Tumor-Enrichment Rate Among Evaluable Clonotypes ---------------
# Panel C is intentionally simplified relative to Panels A/B. For the summary
# composition, TCRex-unmapped enriched clones and enriched clones with an
# "other" (non-tumor/non-viral) TCRex prediction are combined into one orange
# category: neither carries a tumor-associated or viral-associated TCRex label.
# This avoids a visually negligible green sliver / legend entry while preserving
# the biologically important distinction between tumor-associated, viral-associated,
# and all remaining tumor-enriched candidates.
panel_c_levels <- c(
  "Tumor-enriched, no tumor/viral TCRex annotation",
  "Tumor-enriched + TCRex tumor-associated",
  "Tumor-enriched + TCRex viral-associated",
  "Tested, not tumor-enriched"
)

panel_c_colors <- c(
  "Tumor-enriched, no tumor/viral TCRex annotation" = TISSUE_COLORS[["Tumor"]],
  "Tumor-enriched + TCRex tumor-associated" = "#e7298a",
  "Tumor-enriched + TCRex viral-associated" = "#7570b3",
  "Tested, not tumor-enriched" = "grey75"
)

panel_c_df <- clonotype_results %>%
  filter(was_tested) %>%
  mutate(
    panel_c_category = case_when(
      candidate_category %in% c(
        "Tumor-enriched, TCRex-unmapped",
        "Tumor-enriched + other TCRex prediction"
      ) ~ "Tumor-enriched, no tumor/viral TCRex annotation",
      as.character(candidate_category) == "Tumor-enriched + TCRex tumor-associated" ~
        "Tumor-enriched + TCRex tumor-associated",
      as.character(candidate_category) == "Tumor-enriched + TCRex viral-associated" ~
        "Tumor-enriched + TCRex viral-associated",
      TRUE ~ "Tested, not tumor-enriched"
    ),
    panel_c_category = factor(panel_c_category, levels = panel_c_levels)
  ) %>%
  count(global_expansion_status, panel_c_category, name = "n", .drop = TRUE) %>%
  filter(n > 0, !is.na(panel_c_category)) %>%
  group_by(global_expansion_status) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()

# Only show legend categories that are genuinely present after the biologically
# motivated collapse above.
panel_c_observed_categories <- panel_c_levels[
  panel_c_levels %in% unique(as.character(panel_c_df$panel_c_category))
]
panel_c_colors_used <- panel_c_colors[panel_c_observed_categories]

panel_c_summary <- clonotype_results %>%
  group_by(global_expansion_status) %>%
  summarise(
    n_tested = sum(was_tested),
    n_enriched = sum(was_tested & is_tumor_enriched),
    enrichment_rate = ifelse(n_tested > 0, n_enriched / n_tested * 100, NA_real_),
    .groups = "drop"
  ) %>%
  complete(
    global_expansion_status = factor(EXPANSION_LEVELS, levels = EXPANSION_LEVELS),
    fill = list(n_tested = 0L, n_enriched = 0L, enrichment_rate = NA_real_)
  ) %>%
  mutate(
    global_expansion_status = factor(global_expansion_status, levels = EXPANSION_LEVELS),
    summary_label = ifelse(
      n_tested > 0,
      sprintf("tested n=%d\n%.1f%% enriched", n_tested, enrichment_rate),
      "Not evaluable"
    )
  )

write_csv(
  panel_c_summary,
  file.path(table_dir, "Figure5_tested_tier_enrichment_summary.csv")
)

p5c <- ggplot(
  panel_c_df,
  aes(x = global_expansion_status, y = percentage, fill = panel_c_category)
) +
  geom_col(color = "black", linewidth = 0.18, width = 0.72) +
  geom_text(
    data = panel_c_summary,
    aes(x = global_expansion_status, y = 105, label = summary_label),
    inherit.aes = FALSE,
    size = 2.5, fontface = "bold", lineheight = 0.95, color = "grey20"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(
    values = panel_c_colors_used,
    breaks = panel_c_observed_categories,
    limits = panel_c_observed_categories,
    name = "Classification among tested clonotypes",
    drop = TRUE
  ) +
  scale_y_continuous(
    breaks = c(0, 25, 50, 75, 100),
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(0, 112), clip = "off") +
  labs(
    title = "C. Tumor-enrichment rate among statistically evaluable clonotypes",
    subtitle = paste0(
      "Denominator = tested clonotypes within each tier (Tumor + NAT clone size >= ",
      MIN_TN_CLONE_SIZE_FOR_TEST,
      "); orange combines TCRex-unmapped and other non-tumor/non-viral predictions"
    ),
    x = "Global clone expansion tier",
    y = "Percentage of tested clonotypes (%)"
  ) +
  project_theme() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right",
    plot.margin = margin(8, 12, 5.5, 5.5)
  )

# 11. Combine & Export Figure 5 -------------------------------------------------
figure5_3panel <- (p5a | p5b) / p5c +
  plot_layout(heights = c(1.08, 0.92)) +
  plot_annotation(
    title = "Figure 5. Prioritization of tumor-enriched TCR clonotypes using tissue enrichment and TCRex annotation",
    subtitle = "Within-patient Tumor-vs-NAT statistics identify candidates; TCRex matches are treated as database predictions rather than functional proof",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10, color = "grey30")
    )
  )

paths <- save_figure_pair(
  figure5_3panel,
  file.path(figure_dir, "Figure5_tumor_associated_tcrs"),
  width = 13.5,
  height = 10.8
)

cat("============================================================\n")
cat(" Figure 5 completed successfully.\n")
cat(" PNG:", paths[["png"]], "\n")
cat(" SVG:", paths[["svg"]], "\n")
cat(" Candidate table:", file.path(table_dir, "Figure5_tumor_enriched_tcrex_unmapped_candidates.csv"), "\n")
cat("============================================================\n\n")
