# This script is used to draw metagene plots of DNC intervals (Fig. 3A, Fig. 5B, SFig. 3A-C, SFig. 5A-B);

library(tidyverse)
library(rtracklayer)

scripts <- c("batch_read_track_data", "metagene_matrix", "draw_metagene_plot", "merge_and_normalize_GRanges")
for (script in scripts) {
  paste0("https://github.com/Maxim-Ivanov/Utility_functions/blob/main/", script, ".R?raw=TRUE") %>% devtools::source_url()
}

# Load the confirmed DNC loci:
dnc_dir <- "." # change to the directory with DNC_loci.RDS file returned by 04-Find_DNC_loci.R
dnc <- file.path(dnc_dir, "DNC_loci.RDS") %>% readRDS()

# Make windows [-100, +500 bp] around TSS of protein-coding genes in DNC loci:
host_win <- mcols(dnc)$host_coord %>% resize(500, "start") %>% resize(600, "end")
# Make windows [-100, +500 bp] around TSS of DNC transcripts in DNC loci:
dnc_win <- dnc %>% resize(500, "start") %>% resize(600, "end")

# Load our NET-seq data for wt, hda1 and hda3 yeast strains:
netseq_dir <- "." # change to the directory with NET-seq Bedgraph files returned by the 02-Alignment_of_novel_yeast_NET-seq_data.sh script
bg_files <- list.files(netseq_dir, pattern = "^NETseq_s\\d{2}.*rep[12]_norm1M\\.bedgraph\\.gz$")
netseq_data <- batch_read_track_data(bg_files, dir = netseq_dir, format = "bedGraph", seqinfo = seqinfo(dnc))
names(netseq_data) <- names(netseq_data) %>% str_replace("_norm1M.bedgraph.gz", "") %>% str_replace("^NETseq_s\\d{2}_", "")

# Merge replicates and normalize to 1M tags:
idx <- seq(1, length(netseq_data), by = 2)
netseq_data_merged <- lapply(idx, function(x) { merge_and_normalize_GRanges(netseq_data[c(x, x + 1)]) })
names(netseq_data_merged) <- names(netseq_data)[idx] %>% str_replace("_rep1", "")

##### Draw NET-seq metagene plots for Fig. 3A --------------------------------------------------------------------

draw_dnc_metagene_plot <- function(data, dnc_win, host_win, ttl_suffix) {
  ml1 <- lapply(data, metagene_matrix, intervals = dnc_win, scaling = FALSE, matrix.length = 600, anchor = "start", skip.zeros = FALSE, skip.outliers = FALSE, skip.top.obs = TRUE, shrink = TRUE)
  draw_metagene_plot(ml1, x.axis = seq(-19, 100), vline = 0, title = paste0("DNC windows (n=", length(dnc_win), ") - ", ttl_suffix), 
                   xlabel = "Windows [-100, +500] relative to DNC TSS (5 bp bins)", width = 10, height = 8, units = "in")
  ml2 <- lapply(data, metagene_matrix, intervals = host_win, scaling = FALSE, matrix.length = 600, anchor = "start", skip.zeros = FALSE, skip.outliers = FALSE, skip.top.obs = TRUE, shrink = TRUE)
  ml2 <- lapply(ml2, function(x) { return(x[, seq(ncol(x), 1)]) }) # mirror the matrices horizontally
  draw_metagene_plot(ml2, x.axis = seq(-99, 20), vline = 0, title = paste0("Coding windows (n=", length(host_win), ") - ", ttl_suffix), 
                   xlabel = "Windows [-100, +500] relative to coding TSS (5 bp bins)", width = 10, height = 8, units = "in")
}

draw_dnc_metagene_plot(netseq_data_merged, dnc_win, host_win, "hda1 hda2 vs WT")


##### Draw NET-seq metagene plots for SFig. 3A-C ----------------------------------------------------------------

# Stratify DNC loci by NET-seq FPKM of the DNC transcript:
grp <- ifelse(mcols(dnc)$fpkm < 10, "Low", ifelse(mcols(dnc)$fpkm < 30, "Medium", "High"))

# Draw metagene plots stratified by expression level of the DNC transcript:
draw_dnc_metagene_plot(netseq_data_merged, dnc_win[grp == "Low"], host_win[grp == "Low"], paste0("hda1 hda2 vs WT (low expr n=", sum(grp == "Low"), ")"))
draw_dnc_metagene_plot(netseq_data_merged, dnc_win[grp == "Medium"], host_win[grp == "Medium"], paste0("hda1 hda2 vs WT (medium expr n=", sum(grp == "Medium"), ")"))
draw_dnc_metagene_plot(netseq_data_merged, dnc_win[grp == "High"], host_win[grp == "High"], paste0("hda1 hda2 vs WT (high expr n=", sum(grp == "High"), ")"))


##### Draw ChIP-seq metagene plots (Fig. 5B, SFig. 5A-B) ---------------------------------------------------------

# Generate trios of windows ([DNC TSS - 1000 bp, DNC TSS], [DNC TSS, host TSS], [host TSS + 500 bp] on the host strand):
win1 <- dnc %>% resize(ifelse(width(.) > 1000, 1000, width(.)), "start")
strand(win1) <- ifelse(strand(win1) == "+", "-", "+")
win2 <- pgap(dnc, mcols(dnc)$host_coord, ignore.strand = TRUE)
strand(win2) <- ifelse(strand(win2) == "+", "-", "+")
win3 <- mcols(dnc)$host_coord %>% resize(ifelse(width(.) > 500, 500, width(.)), "start")

# Download ChIP-seq data from Ha et al., (PMID 31537788):
url <- list("H3ac_wt_rep1" = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3445776&format=file&file=GSM3445776%5FYSB787%2D1%2DH3ac%5FnoDup%2Enorm%2Ebw", 
            "H3ac_wt_rep2" = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3445777&format=file&file=GSM3445777%5FYSB787%2D2%2DH3ac%5FnoDup%2Enorm%2Ebw", 
            "H3_wt_rep1" = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3445780&format=file&file=GSM3445780%5FYSB787%2D1%2DH3%5FnoDup%2Enorm%2Ebw", 
            "H3_wt_rep2" = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3445781&format=file&file=GSM3445781%5FYSB787%2D2%2DH3%5FnoDup%2Enorm%2Ebw", 
            "H3ac_hda1_rep1" = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3445782&format=file&file=GSM3445782%5FYTK113%2D1%2DH3ac%5FnoDup%2Enorm%2Ebw", 
            "H3ac_hda1_rep2" = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3445783&format=file&file=GSM3445783%5FYTK113%2D2%2DH3ac%5FnoDup%2Enorm%2Ebw", 
            "H3_hda1_rep1" = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3445786&format=file&file=GSM3445786%5FYTK113%2D1%2DH3%5FnoDup%2Enorm%2Ebw", 
            "H3_hda1_rep2" = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3445787&format=file&file=GSM3445787%5FYTK113%2D2%2DH3%5FnoDup%2Enorm%2Ebw")

for (i in seq_along(url)) {
  download.file(url[[i]], paste0("Ha2019_", names(url)[[i]], ".bw"), method = "curl")
}

# Load ChIP-seq tracks into R:
bw_dir <- "." # change to the directory with downloaded ChIP-seq files
bw_files <- list.files(bw_dir, pattern = "^Ha2019.*bw$")
bw_data <- batch_read_track_data(bw_files, dir = bw_dir, format = "BigWig")
names(bw_data) <- names(bw_data) %>% str_replace("^Ha2019_", "") %>% str_replace(".bw", "")

# Fix the chromosome names:
bw_data <- lapply(bw_data, function(x) { seqlevels(x) <- seqlevels(x) %>% str_replace("Mito", "M") %>% paste0("chr", .); x <- sortSeqlevels(x); return(x) })

# Merge replicates:
idx <- seq(1, length(bw_data), by = 2)
bw_data_merged <- lapply(idx, function(x) { merge_and_normalize_GRanges(bw_data[c(x, x + 1)]) })
names(bw_data_merged) <- names(bw_data)[idx] %>% str_replace("_rep1", "")

draw_chipseq_metagenes <- function(data, win1, win2, win3, title) {
  ml1 <- lapply(data, metagene_matrix, intervals = win1, scaling = FALSE, matrix.length = 1000, anchor = "end", skip.zeros = FALSE, skip.outliers = FALSE, skip.top.obs = TRUE, shrink = TRUE)
  ml2 <- lapply(data, metagene_matrix, intervals = win2, scaling = TRUE, matrix.length = 50, skip.zeros = FALSE, skip.outliers = FALSE, skip.top.obs = TRUE)
  ml3 <- lapply(data, metagene_matrix, intervals = win3, scaling = FALSE, matrix.length = 500, anchor = "start", skip.zeros = FALSE, skip.outliers = FALSE, skip.top.obs = TRUE, shrink = TRUE)
  ml <- mapply(cbind, ml1, ml2, ml3, SIMPLIFY = FALSE)
  ml <- lapply(ml, function(x) { return(x[, seq(ncol(x), 1)]) }) # flip the plot horizontally
  draw_metagene_plot(ml, x.axis = seq(-149, 200), vline = c(-50, 0), title = title, xlabel = "[-1Kb, DNC TSS] (5 bp bins) + [DNC TSS, host TSS] (50 bins) + [host TSS, +500bp] (5 bp bins)", width = 10, height = 8, units = "in")
}

# Draw SFig. 5A:
draw_chipseq_metagenes(bw_data_merged[c("H3ac_wt", "H3ac_hda1")], win1, win2, win3, paste0("H3ac in hda1 vs WT (n=", length(win1), ")"))

# Draw SFig. 5B:
draw_chipseq_metagenes(bw_data_merged[c("H3_wt", "H3_hda1")], win1, win2, win3, paste0("H3 in hda1 vs WT (n=", length(win1), ")"))

# "Normalize" H3ac signal by H3:
calculate_log2_ratio_of_two_gr <- function(gr1, gr2, pseudocount = 0.1) {
  stopifnot(class(gr1) == "GRanges" && !is.null(score(gr1)))
  stopifnot(class(gr2) == "GRanges" && !is.null(score(gr2)))
  # (if there was any strand info in the input GRanges, it is lost)
  cov1 <- coverage(gr1, weight = "score") + pseudocount
  cov2 <- coverage(gr2, weight = "score") + pseudocount
  ratio <- log2(cov1 / cov2)
  out <- bindAsGRanges(score = ratio)
}

bw_data_ratio <- list("wt" = calculate_log2_ratio_of_two_gr(bw_data_merged[["H3ac_wt"]], bw_data_merged[["H3_wt"]]), 
                      "hda1" = calculate_log2_ratio_of_two_gr(bw_data_merged[["H3ac_hda1"]], bw_data_merged[["H3_hda1"]]))

# Draw Fig. 5B:
ml1 <- lapply(bw_data_ratio, metagene_matrix, intervals = win1, scaling = FALSE, matrix.length = 1000, anchor = "end", skip.zeros = FALSE, skip.outliers = FALSE, skip.top.obs = TRUE, shrink = TRUE)
ml2 <- lapply(bw_data_ratio, metagene_matrix, intervals = win2, scaling = TRUE, matrix.length = 50, skip.zeros = FALSE, skip.outliers = FALSE, skip.top.obs = TRUE)
ml3 <- lapply(bw_data_ratio, metagene_matrix, intervals = win3, scaling = FALSE, matrix.length = 500, anchor = "start", skip.zeros = FALSE, skip.outliers = FALSE, skip.top.obs = TRUE, shrink = TRUE)
ml <- mapply(cbind, ml1, ml2, ml3, SIMPLIFY = FALSE)
ml <- lapply(ml, function(x) { return(x[, seq(ncol(x), 1)]) }) # flip the plot horizontally
draw_metagene_plot(ml, x.axis = seq(-199, 150), vline = c(0, 50), title = paste0("Ratio of H3ac to H3 in hda1 vs WT (n=", length(win1), ")"), 
                   xlabel = "[-1Kb, DNC TSS] (5 bp bins) + [DNC TSS, host TSS] (50 bins) + [host TSS, +500bp] (5 bp bins)", 
                   width = 10, height = 8, units = "in", ylabel = "log2(H3ac / H3)")
