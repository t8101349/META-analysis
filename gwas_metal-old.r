args <- commandArgs(trailingOnly = TRUE)
file_mvp <- args[1]
file_gwas <- args[2]

library(data.table)

# Read in both datasets
mvp <- fread(file_mvp, header = TRUE)
adj <- fread(file_gwas, header = TRUE)

# Standardize column names for merging
setnames(mvp, c("SNP_ID", "chrom", "pos", "ref", "alt", "af", "num_samples", "or", "pval"), 
              c("ID", "CHROM", "POS", "REF", "ALT", "MAF", "OBS_CT", "OR", "P"))
setnames(adj, c("ID", "#CHROM", "POS", "REF", "ALT", "A1_FREQ", "OBS_CT", "OR", "LOG(OR)_SE", "P"), 
              c("ID", "CHROM", "POS", "REF", "ALT", "MAF", "OBS_CT", "OR", "LOG_OR_SE", "P"))

# Check MAF and flip if necessary

mvp[!is.na(MAF) & MAF > 0.5,
    c("MAF", "REF", "ALT") := .(1 - MAF, ALT, REF)
]

adj[!is.na(MAF) & MAF > 0.5,
    c("MAF", "REF", "ALT") := .(1 - MAF, ALT, REF)
]


# 假設有 CI 欄位，例如 "ci" 格式為 "lower,upper"
# 例如: "0.775901,1.2075"
extract_se_from_ci <- function(or, ci) {
  if (is.na(or) || is.na(ci) || ci == "") return(NA_real_)
  ci <- gsub("[()\\[\\] ]", "", ci)
  ci <- as.numeric(strsplit(ci, ",")[[1]])
  if (length(ci) != 2 || any(is.na(ci)) || any(ci <= 0) || or <= 0) {
    return(NA_real_)
  }
  log_l <- log(ci[1])
  log_u <- log(ci[2])
  # 95% CI: log(OR) ± 1.96*SE
  se <- (log_u - log_l) / (2 * 1.96)
  return(se)
}

# 對 mvp 資料計算 LOG_OR_SE
mvp[, LOG_OR_SE := mapply(extract_se_from_ci, OR, ci)]

# Merge datasets on SNP ID
meta_data <- merge(mvp, adj, by = "ID", suffixes = c("_mvp", "_adj"))

# Prepare METAL-style input for both studies
metal_input <- rbind(
  mvp[, .(
    MARKER = ID,
    CHROM = CHROM,
    POS = POS,
    REF = REF,
    ALT = ALT,
    OBS_CT = OBS_CT,
    ALT_FREQS = MAF,
    EFFECT = log(OR),
    LOG_OR_SE = LOG_OR_SE, # If SE available, use it
    P = P
  )],
  adj[, .(
    MARKER = ID,
    CHROM = CHROM,
    POS = POS,
    REF = REF,
    ALT = ALT,
    OBS_CT = OBS_CT,
    ALT_FREQS = MAF,
    EFFECT = log(OR),
    LOG_OR_SE = LOG_OR_SE,
    P = P
  )]
)

# Write out for METAL or further analysis
fwrite(metal_input, "metal_input.txt", sep = "\t")

# Meta-analysis (inverse-variance weighted, fixed effect)
meta_results <- meta_data[!is.na(LOG_OR_SE_adj) & !is.na(OR_mvp), {
  # Use log(OR) and SE from both studies if available
  log_or1 <- log(OR_mvp)
  log_or2 <- log(OR_adj)
  se1 <- LOG_OR_SE_mvp
  se2 <- LOG_OR_SE_adj
  # If SE missing, skip meta for that SNP
  w1 <- 1 / se1^2
  w2 <- 1 / se2^2
  meta_log_or <- (log_or1 * w1 + log_or2 * w2) / (w1 + w2)
  meta_se <- sqrt(1 / (w1 + w2))
  meta_z <- meta_log_or / meta_se # z-test
  meta_p <- 2 * pnorm(-abs(meta_z))
  .(MARKER = ID, META_LOG_OR = meta_log_or, META_SE = meta_se, META_P = meta_p)
}]


fwrite(meta_results, "meta_results.txt", sep = "\t")

# Rscript /home/Weber/Metal/gwas_metal.r MVP_R4.1000G_AGR.Phe_185.META.GIA.dbGaP.txt.gz GWAS/ISPROSTATECANCER_TPMI_imputed_adjGWAS.glm.logistic