#!/usr/bin/env Rscript

# Harmonize BCAC East Asian and BioBank Japan breast cancer GWAS summary
# statistics, then combine them using fixed-effect inverse-variance
# meta-analysis.
#
# Usage:
#   Rscript EAS_BCAC_BBJ_MetaAnalysis.R \
#     all.meta.icogs.oncoarray.asian.overall.Apr2016.1.txt \
#     GWASsummary_BrC_Japanese_SakaueKanai2020.auto.txt \
#     EAS_BCAC_BBJ_meta_sumdata.txt

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3L) {
  stop(
    paste(
      "Usage: Rscript EAS_BCAC_BBJ_MetaAnalysis.R",
      "<BCAC_EAS_summary_statistics>",
      "<BBJ_summary_statistics>",
      "<output_file>"
    ),
    call. = FALSE
  )
}

bcac_path <- args[[1L]]
bbj_path <- args[[2L]]
output_path <- args[[3L]]

suppressPackageStartupMessages(library(data.table))

check_columns <- function(data, required_columns, label) {
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      label,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

make_sorted_variant_id <- function(chr_pos, allele1, allele2) {
  allele1 <- toupper(allele1)
  allele2 <- toupper(allele2)
  first_allele <- ifelse(allele1 < allele2, allele1, allele2)
  second_allele <- ifelse(allele1 < allele2, allele2, allele1)
  paste(chr_pos, first_allele, second_allele, sep = "_")
}

bcac_raw <- fread(bcac_path)
bbj_raw <- fread(bbj_path)

check_columns(
  bcac_raw,
  c(
    "MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr",
    "P-value"
  ),
  "BCAC input"
)
check_columns(
  bbj_raw,
  c(
    "v", "SNPID", "Allele1", "Allele2", "AF.Controls", "BETA", "SE",
    "p.value"
  ),
  "BBJ input"
)

if (any(!grepl("^[^_]+_[^_]+_", bcac_raw$MarkerName))) {
  stop("Some BCAC MarkerName values do not contain chromosome and position.",
       call. = FALSE)
}
if (any(!grepl("^[^:]+:[^:]+:", bbj_raw$v))) {
  stop("Some BBJ v values do not contain chromosome and position.",
       call. = FALSE)
}

# BCAC Allele1 is the effect allele.
bcac_chr_pos <- sub("^([^_]+_[^_]+)_.*", "\\1", bcac_raw$MarkerName)
bcac_effect_allele <- toupper(bcac_raw$Allele1)
bcac_non_effect_allele <- toupper(bcac_raw$Allele2)

bcac <- data.table(
  unique_SNP_id = make_sorted_variant_id(
    bcac_chr_pos,
    bcac_effect_allele,
    bcac_non_effect_allele
  ),
  effect_allele_BCAC = bcac_effect_allele,
  non_effect_allele_BCAC = bcac_non_effect_allele,
  Freq_effect_BCAC = bcac_raw$Freq1,
  BETA_BCAC = bcac_raw$Effect,
  SE_BCAC = bcac_raw$StdErr,
  P_BCAC = bcac_raw[["P-value"]],
  N_eff_BCAC = 1 / (
    2 * bcac_raw$Freq1 * (1 - bcac_raw$Freq1) * bcac_raw$StdErr^2
  )
)

# BBJ Allele2 is the effect allele in the source file.
bbj_chr_pos <- sub("^([^:]+):([^:]+).*", "\\1_\\2", bbj_raw$v)
bbj_effect_allele <- toupper(bbj_raw$Allele2)
bbj_non_effect_allele <- toupper(bbj_raw$Allele1)

bbj <- data.table(
  unique_SNP_id = make_sorted_variant_id(
    bbj_chr_pos,
    bbj_effect_allele,
    bbj_non_effect_allele
  ),
  effect_allele_BBJ = bbj_effect_allele,
  non_effect_allele_BBJ = bbj_non_effect_allele,
  Freq_effect_BBJ = bbj_raw$AF.Controls,
  BETA_BBJ = bbj_raw$BETA,
  SE_BBJ = bbj_raw$SE,
  P_BBJ = bbj_raw$p.value,
  N_eff_BBJ = as.integer(
    1 / (
      2 * bbj_raw$AF.Controls * (1 - bbj_raw$AF.Controls) * bbj_raw$SE^2
    )
  ),
  SNPID = bbj_raw$SNPID
)

if (anyDuplicated(bcac$unique_SNP_id)) {
  warning("BCAC input contains duplicate harmonized variant identifiers.",
          call. = FALSE)
}
if (anyDuplicated(bbj$unique_SNP_id)) {
  warning("BBJ input contains duplicate harmonized variant identifiers.",
          call. = FALSE)
}

combined <- merge(
  bbj,
  bcac,
  by = "unique_SNP_id",
  all = TRUE,
  sort = FALSE
)

# BBJ defines the meta-analysis allele orientation when it is available. For
# BCAC-only variants, the original BCAC orientation is retained.
combined[, effect_allele_meta := ifelse(
  !is.na(effect_allele_BBJ),
  effect_allele_BBJ,
  effect_allele_BCAC
)]
combined[, non_effect_allele_meta := ifelse(
  !is.na(non_effect_allele_BBJ),
  non_effect_allele_BBJ,
  non_effect_allele_BCAC
)]

# flip_sign records whether the original BCAC effect allele differed from the
# BBJ effect allele. Reorient every BCAC effect-allele field together. The BCAC
# beta is flipped once, and the effect-allele frequency is complemented once.
combined[, flip_sign :=
  !is.na(effect_allele_BBJ) &
  !is.na(effect_allele_BCAC) &
  effect_allele_BBJ != effect_allele_BCAC]

combined[flip_sign == TRUE, `:=`(
  BETA_BCAC = -BETA_BCAC,
  Freq_effect_BCAC = 1 - Freq_effect_BCAC
)]
combined[!is.na(effect_allele_BCAC), `:=`(
  effect_allele_BCAC = effect_allele_meta,
  non_effect_allele_BCAC = non_effect_allele_meta
)]

if (any(!is.na(combined$BETA_BBJ) &
        (!is.finite(combined$SE_BBJ) | combined$SE_BBJ <= 0))) {
  stop("BBJ contains a nonpositive or nonfinite SE for an analyzed variant.",
       call. = FALSE)
}
if (any(!is.na(combined$BETA_BCAC) &
        (!is.finite(combined$SE_BCAC) | combined$SE_BCAC <= 0))) {
  stop("BCAC contains a nonpositive or nonfinite SE for an analyzed variant.",
       call. = FALSE)
}

combined[, w_BBJ := ifelse(!is.na(BETA_BBJ), 1 / SE_BBJ^2, NA_real_)]
combined[, w_BCAC := ifelse(!is.na(BETA_BCAC), 1 / SE_BCAC^2, NA_real_)]

both_studies <- !is.na(combined$w_BBJ) & !is.na(combined$w_BCAC)

combined[, BETA_meta := ifelse(
  both_studies,
  (w_BBJ * BETA_BBJ + w_BCAC * BETA_BCAC) / (w_BBJ + w_BCAC),
  ifelse(is.na(w_BBJ), BETA_BCAC, BETA_BBJ)
)]
combined[, SE_meta := ifelse(
  both_studies,
  sqrt(1 / (w_BBJ + w_BCAC)),
  ifelse(is.na(w_BBJ), SE_BCAC, SE_BBJ)
)]
combined[, P_meta := 2 * pnorm(-abs(BETA_meta / SE_meta))]
combined[, N_eff_meta := ifelse(
  !is.na(N_eff_BBJ) & !is.na(N_eff_BCAC),
  N_eff_BBJ + N_eff_BCAC,
  ifelse(is.na(N_eff_BBJ), N_eff_BCAC, N_eff_BBJ)
)]

bcac_present <- !is.na(combined$BETA_BCAC)
bbj_present <- !is.na(combined$BETA_BBJ)

stopifnot(
  !any(
    bcac_present &
      (combined$effect_allele_BCAC != combined$effect_allele_meta |
       combined$non_effect_allele_BCAC != combined$non_effect_allele_meta),
    na.rm = TRUE
  ),
  !any(
    bbj_present &
      (combined$effect_allele_BBJ != combined$effect_allele_meta |
       combined$non_effect_allele_BBJ != combined$non_effect_allele_meta),
    na.rm = TRUE
  ),
  !any(
    !is.na(combined$Freq_effect_BCAC) &
      (combined$Freq_effect_BCAC < 0 | combined$Freq_effect_BCAC > 1)
  )
)

output_columns <- c(
  "unique_SNP_id",
  "effect_allele_BBJ", "non_effect_allele_BBJ", "Freq_effect_BBJ",
  "BETA_BBJ", "SE_BBJ", "P_BBJ", "N_eff_BBJ", "SNPID",
  "effect_allele_BCAC", "non_effect_allele_BCAC", "Freq_effect_BCAC",
  "BETA_BCAC", "SE_BCAC", "P_BCAC", "N_eff_BCAC", "flip_sign",
  "w_BBJ", "w_BCAC", "BETA_meta", "SE_meta", "P_meta",
  "effect_allele_meta", "non_effect_allele_meta", "N_eff_meta"
)

setcolorder(combined, output_columns)

output_directory <- dirname(output_path)
if (output_directory != ".") {
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
}

# The historical release format is whitespace-delimited.
write.table(
  combined,
  file = output_path,
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

message("Wrote ", format(nrow(combined), big.mark = ","), " variants to ",
        output_path)
message("BCAC records reoriented to BBJ: ",
        format(sum(combined$flip_sign), big.mark = ","))
