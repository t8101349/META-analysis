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

# Meta-analysis (inverse-variance weighted, fixed effect) + Heterogeneity statistics
meta_results <- meta_data[!is.na(LOG_OR_SE_adj) & !is.na(OR_mvp) & !is.na(LOG_OR_SE_mvp), {
  # log(OR) and SE from both studies
  log_or <- c(log(OR_mvp), log(OR_adj))
  se <- c(LOG_OR_SE_mvp, LOG_OR_SE_adj)
  w <- 1 / se^2
  meta_log_or <- sum(w * log_or) / sum(w)
  meta_se <- sqrt(1 / sum(w))
  #z-test
  meta_z <- meta_log_or / meta_se
  meta_p <- 2 * pnorm(-abs(meta_z))
  
  # Q-test (Cochran's Q)
  Q <- sum(w * (log_or - meta_log_or)^2)
  k <- length(log_or)
  # I² statistic
  I2 <- if (Q > (k - 1)) ((Q - (k - 1)) / Q) * 100 else 0
  
  .(META_LOG_OR = meta_log_or,
    META_SE = meta_se,
    META_P = meta_p,
    Q = Q,
    I2 = I2)
}, by = ID]


fwrite(meta_results, "meta_fixed_results.txt", sep = "\t")

# 隨機效應模型 (DerSimonian-Laird)
meta_random_results <- meta_data[!is.na(LOG_OR_SE_adj) & !is.na(OR_mvp) & !is.na(LOG_OR_SE_mvp), {
  log_or <- c(log(OR_mvp), log(OR_adj))
  se <- c(LOG_OR_SE_mvp, LOG_OR_SE_adj)
  w <- 1 / se^2
  meta_log_or_fixed <- sum(w * log_or) / sum(w)
  k <- length(log_or)
  Q <- sum(w * (log_or - meta_log_or_fixed)^2)
  # DerSimonian-Laird tau^2
  tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
  w_star <- 1 / (se^2 + tau2)
  meta_log_or_random <- sum(w_star * log_or) / sum(w_star)
  meta_se_random <- sqrt(1 / sum(w_star))
  meta_z_random <- meta_log_or_random / meta_se_random
  meta_p_random <- 2 * pnorm(-abs(meta_z_random))
  I2 <- if (Q > (k - 1)) ((Q - (k - 1)) / Q) * 100 else 0
  .(META_LOG_OR = meta_log_or_random,
    META_SE = meta_se_random,
    META_P = meta_p_random,
    Q = Q,
    I2 = I2)
}, by = ID]

fwrite(meta_random_results, "meta_random_results.txt", sep = "\t")

# Rscript /home/Weber/Metal/gwas_metal.r MVP_R4.1000G_AGR.Phe_185.META.GIA.dbGaP.txt.gz GWAS/ISPROSTATECANCER_TPMI_imputed_adjGWAS.glm.logistic