
library(data.table)
library(metafor)

# 讀檔
dt <- fread("metal_input.txt")

# 建立空的結果表
meta_res <- dt[, .(MARKER, META_BETA = NA_real_, META_SE = NA_real_, META_P = NA_real_)]

# 對每個 SNP 計算 fixed effect
for (i in seq_len(nrow(dt))) {
  res <- rma(yi = dt$EFFECT[i], sei = dt$LOG_OR_SE[i], method = "FE")
  meta_res$META_BETA[i] <- res$beta
  meta_res$META_SE[i]   <- res$se
  meta_res$META_P[i]    <- res$pval
}

# 輸出結果
fwrite(meta_res, "meta_results_metafor.txt", sep = "\t")
