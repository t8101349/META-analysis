##
# Rscript RunTopR2.R --input Pheno_TPMI_imputed_adjGWAS.glm.logistic --thresh 1e-5 --annotate 5e-8 --ymax 20 --output Aneurysm
##
library("topr")
library("data.table")
library("dplyr")
library("ggplot2")
library("argparser")
# 解析參數
parser <- arg_parser("Process manhattan_plot by topr")
parser <- add_argument(parser, "--input", type = 'character', help = "read input meta_results file")
parser <- add_argument(parser, "--meta", type = 'character', help = "read meta_input file")
parser <- add_argument(parser, "--thresh", type = 'numeric', help = "thresh for gene list", default = 1e-5)
parser <- add_argument(parser, "--annotate", type = 'numeric', help = "annotate gene thresh in plot", default = 5e-8)
parser <- add_argument(parser, "--ymax", type='character', help="ymax", default = "NA")
parser <- add_argument(parser, "--output", type = 'character', help = "output name and title")
parser <- add_argument(parser, "--subtitle", type = 'character', help = "subtitle for plot", default = "")
args <- parse_args(parser, commandArgs(trailingOnly = TRUE))


# -- 引數賦值 --
i <- args$input #meta_results.txt的路徑
meta <- args$meta #metal_input.txt的路徑
out <- args$output  #檔名
thresh <- args$thresh #找基因時SNP的最大值
annotate = args$annotate #畫圖時標基因的最大值
ymax <- args$ymax  #畫圖的Y最大值， NULL表示自動
subtitle_text <- args$subtitle

if (args$ymax == "NA" || args$ymax == "NULL") {
  ymax <- NULL
} else {
  ymax <- as.numeric(args$ymax)
}


# Step 1: 合併 metal_input.txt 與 metal_results.txt
metal_input <- fread(i)      # MARKER, CHROM, POS, REF, ALT, OBS_CT, ALT_FREQS, EFFECT, LOG_OR_SE, P
metal_results <- fread(meta)  # ID, META_LOG_OR, META_SE, META_P

# 內連結合
merge_dataframe <- merge(metal_input, metal_results, by.x = "MARKER", by.y = "ID", all = FALSE)


# Step 2: 調整欄位名稱與順序，符合 output 格式
output_df <- merge_dataframe %>%
  transmute(
    # output格式: #CHROM, POS, ID, REF, ALT, PROVISIONAL_REF?, A1, OMITTED, A1_FREQ, FIRTH?, TEST, OBS_CT, OR, LOG(OR)_SE, L95, U95, Z_STAT, P, ERRCODE
    # 這裡只填入可取得的欄位，其餘補NA或預設值
    `#CHROM` = CHROM,
    POS = POS,
    ID = MARKER,
    REF = REF,
    ALT = ALT,
    `PROVISIONAL_REF?` = NA,
    A1 = ALT,
    OMITTED = NA,
    A1_FREQ = ALT_FREQS,
    `FIRTH?` = NA,
    TEST = "META",
    OBS_CT = OBS_CT,
    OR = exp(META_LOG_OR),
    `LOG(OR)_SE` = META_SE,
    L95 = exp(META_LOG_OR - 1.96 * META_SE),
    U95 = exp(META_LOG_OR + 1.96 * META_SE),
    Z_STAT = META_LOG_OR / META_SE,
    P = META_P,
    ERRCODE = "."
  )

# Step 3: 輸出合併後檔案
last <- sub(".*_", "", out)
output_name <- paste0("finaloutput_", last, "_input.txt")
fwrite(output_df, output_name, sep = "\t", na = "NA")

# Step 4: 使用 topr 套件進行後續分析
# 讀取輸入檔案
cmuh_data <- fread(output_name, sep="\t", header=TRUE)
names(cmuh_data)[names(cmuh_data) == "#CHROM"] <- "CHROM"

# 確保欄位型別正確
cmuh_data[, CHROM := as.integer(CHROM)]
cmuh_data[, POS := as.integer(POS)]

# 找 lead SNP 並加上基因註解
CMUH_SNP <- get_lead_snps(cmuh_data,thresh=thresh) %>% annotate_with_nearest_gene()

# 繪製 Manhattan 圖
print(paste("Generating TIFF file:", paste0("manhattan_plot_", out, ".tiff")))
tiff_filename <- paste0("manhattan_plot_", out, ".tiff")
# subtitle_text <- paste0("Case: ", case_count, ", Control: ", control_count)

manhattan_plot <- manhattan(
  cmuh_data,
  annotate = annotate,
  color = c("#2d42ba"),
  title = out,
  ymax = ymax
) +
  ggtitle(out, subtitle = subtitle_text) +
  theme(
    axis.text.x = element_text(size = 7)
  )

ggsave(
  filename = paste0("manhattan_plot_", out, ".tiff"),
  plot = manhattan_plot,
  device = "tiff",
  width = 9,
  height = 4,
  units = "in",
  dpi = 600
)


tiff_filename1 <- paste0("QQ_plot_", out, ".tiff")
qq_plot <- qqtopr(cmuh_data)

ggsave(filename = tiff_filename1, plot = qq_plot, device = "tiff", width = 4, height = 4, units = "in", dpi = 600)

fwrite(CMUH_SNP, paste0(out, "_SNP_Gene.txt"), sep="\t", na="NA")






