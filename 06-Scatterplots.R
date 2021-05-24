# This script uses flow cytometry data (provided by 05-Load_flow_cytometry_data.R) to draw the following scatterplots:
# 1) Directionality scores in ORC2 and GCG1 screens (Fig. 1D-E);
# 2) Directionality scores in ORC2 screen vs GCG1 screen (Fig. 2B);
# 3) NET-seq coding/non-coding ratios in cac2 vs wt (SFig. 1C);

library(tidyverse)
library(rtracklayer)

scripts <- c("batch_read_track_data", "get_overlapping_scores", "trim_by_down_or_upstream_features")
for (script in scripts) {
  paste0("https://github.com/Maxim-Ivanov/Utility_functions/blob/main/", script, ".R?raw=TRUE") %>% devtools::source_url()
}

# Load all SacCer3 genes:
genes <- readRDS("Improved_annotation_for_Scerevisiae.RDS") # file produced by the 03-Correct_and_expand_SacCer3_annotation.R script
names(genes) <- mcols(genes)$name

##### Scatterplots of directionality scores --------------------------------------------------------------

# Load dataframes returned by the Fortessa script:
df_dir <- "." # change to the directory with dataframes produced by the 05-Load_flow_cytometry_data.R script
orc <- file.path(df_dir, "df_orc2.RDS") %>% readRDS()
gcg <- file.path(df_dir, "df_gcg1.RDS") %>% readRDS()

# Exclude duplicated plates from the ORC2 screen:
orc <- filter(orc, !str_detect(Plate, "_uthra"))

# Skip one unsuccessful plate from the ORC2 screen:
orc <- filter(orc, Plate != "170223_3-1x49-71")

# Deduplicate on gene names:
orc <- orc %>% mutate(Name = str_to_upper(Name)) %>% filter(!duplicated(Name))
gcg <- gcg %>% mutate(Name = str_to_upper(Name)) %>% filter(!duplicated(Name))

# Skip control samples:
orc <- orc %>% filter(Type == "S" & !Name %in% c("-", "BY4741") & !is.na(YFP_well))
gcg <- gcg %>% filter(Type == "S" & !Name %in% c("-", "BY4741") & !is.na(YFP_well))

# Export directionality scores as Tab-delimited files:
#orc_tbl <- orc %>% as_tibble() %>% select(Name, Geom_dist, Pval_Geom) %>% dplyr::rename(Geom_dist_ORC2 = Geom_dist, Pval_Geom_ORC2 = Pval_Geom)
#gcg_tbl <- gcg %>% as_tibble() %>% select(Name, Geom_dist, Pval_Geom) %>% dplyr::rename(Geom_dist_GCG1 = Geom_dist, Pval_Geom_GCG1 = Pval_Geom)
#orc_gcg_tbl <- full_join(orc_tbl, gcg_tbl, by = "Name")
#temp_tbl <- mcols(genes) %>% as_tibble() %>% select(name, name_2) %>% dplyr::rename(Name = name, Name_2 = name_2)
#orc_gcg_tbl <- left_join(orc_gcg_tbl, temp_tbl, by = "Name") %>% select(Name, Name_2, everything())
#write_tsv(orc_gcg_tbl, "ORC2_GCG1_directionality_scores.txt")

# Add the favorite genes:
fav <- read_tsv("Genes_to_highlight_on_scatterplots.txt")
orc <- left_join(orc, fav, by = "Name")
gcg <- left_join(gcg, fav, by = "Name")

# Sort tibbles for better visibility of chosen genes:
orc <- orc %>% mutate(alpha = ifelse(is.na(Complex), 0.01, 1), Complex_2 = ifelse(is.na(Label), NA, Complex), alpha_2 = ifelse(is.na(Label), 0.01, alpha)) %>% arrange(alpha, alpha_2)
gcg <- gcg %>% mutate(alpha = ifelse(is.na(Complex), 0.01, 1), Complex_2 = ifelse(is.na(Label), NA, Complex), alpha_2 = ifelse(is.na(Label), 0.01, alpha)) %>% arrange(alpha, alpha_2)

# Define colors:
palette = list("Rpd3 HDAC" = "#000000", "Hda1C" = "#E69F00", "Histones" = "#56B4E9", "Protein Urmylation" = "#009E73", 
               "CAF-I" = "#F0E442", "SAGA HAT" = "#0072B2", "SAGA deubiquitinase" = "#D55E00", "ISW2 chromatin remodeller" = "#673191", "SWI/SNF" = "#D15E81")

# ORC2 plot:
ttl <- "ORC2 scatterplot (Fig. 1E)"
coef <- lm(orc$YFP_well_norm ~ orc$mCh_well_norm)$coefficients
ggplot(orc, aes(x = mCh_well_norm, y = YFP_well_norm, fill = Complex_2, colour = Complex_2, alpha = alpha_2)) + 
  geom_point(shape = 21, size = 2) + 
  ggtitle(ttl) + theme_bw() + 
  geom_abline(intercept = coef[[1]], slope = coef[[2]]) + 
  coord_cartesian(xlim = c(0, 3), ylim = c(0.7, 1.5)) + 
  geom_text(aes(label = Label), nudge_y = 0.02) + 
  scale_fill_manual(values = palette, na.value = "#999999") + 
  scale_colour_manual(values = palette, na.value = "#999999") + 
  xlab("mCherry") + ylab("YFP") + 
  theme(legend.position = "none")
for (ext in c(".png", ".pdf")) {
  ggsave(paste0(ttl, ext), width = 6, height = 6, units = "in")
}

# GCG1 plot:
ttl <- "GCG1 scatterplot (Fig. 1D)"
coef <- lm(gcg$YFP_well_norm ~ gcg$mCh_well_norm)$coefficients
ggplot(gcg, aes(x = mCh_well_norm, y = YFP_well_norm, fill = Complex_2, colour = Complex_2, alpha = alpha_2)) + 
  geom_point(shape = 21, size = 2) + 
  ggtitle(ttl) + theme_bw() + 
  geom_abline(intercept = coef[[1]], slope = coef[[2]]) + 
  coord_cartesian(xlim = c(-0.2, 4.5), ylim = c(0.8, 1.5)) + 
  geom_text(aes(label = Label), nudge_y = 0.02) + 
  scale_fill_manual(values = palette, na.value = "#999999") + 
  scale_colour_manual(values = palette, na.value = "#999999") + 
  xlab("mCherry") + ylab("YFP") + 
  theme(legend.position = "none")
for (ext in c(".png", ".pdf")) {
  ggsave(paste0(ttl, ext), width = 6, height = 6, units = "in")
}

# Join dataframes by gene names:
orc_cut <- orc %>% select(Name, Geom_dist) %>% dplyr::rename(ORC2 = Geom_dist)
gcg_cut <- gcg %>% select(Name, Geom_dist, Complex, Label, alpha) %>% dplyr::rename(GCG1 = Geom_dist)
joined <- inner_join(orc_cut, gcg_cut, by = "Name")

# ORC2 vs GCG1 plot:
ttl <- "ORC2 vs GCG1 scatterplot (Fig. 2B)"
ggplot(joined, aes(x = GCG1, y = ORC2, fill = Complex, colour = Complex, alpha = alpha)) + 
  geom_point(shape = 21, size = 2) +
  ggtitle(ttl) + theme_bw() + 
  coord_cartesian(xlim = c(-0.47, 0.2), ylim = c(-0.45, 0.2)) + 
  scale_fill_manual(values = palette, na.value = "#999999") + 
  scale_colour_manual(values = palette, na.value = "#999999") + 
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_text(aes(label = Label), nudge_y = 0.02) + 
  guides(alpha = FALSE) + 
  theme(legend.title = element_blank())
for (ext in c(".png", ".pdf")) {
  ggsave(paste0(ttl, ext), width = 8, height = 6, units = "in")
}


##### Scatterplot of NET-seq coding/non-coding ratios ------------------------------------------------------

# Load NET-seq data for wt and cac2 samples from Marquardt et al., 2014 (PMID 24949978):
bg_dir <- "." # change to the directory with NET-seq Bedgraph files produced by the 01-Remapping_published_NET-seq_datasets.sh script
bg_files <- list.files(bg_dir, pattern = "NETseq_Marquardt2014.*rep[12]_first_base.bedgraph.gz$")
bg_data <- batch_read_track_data(bg_files, bg_dir, format = "bedGraph", seqinfo = seqinfo(genes))
names(bg_data) <- names(bg_data) %>% str_replace("_first_base.bedgraph.gz", "") %>% str_replace("NETseq_Marquardt2014_", "")

# Add GCG1 and ORC2 names to genes:
mapping <- tibble(id = c("YER163C", "YBR060C"), name2 = c("GCG1", "ORC2"))
mcols(genes)$name2 <- tibble("id" = mcols(genes)$name) %>% left_join(mapping, by = "id") %>% .$name2

# Find all coding windows ([-100, +500 bp] from TSS of protein-coding gene) and DNC windows ([-500, + 100 bp], flipped to the opposite strand):
npcd <- genes[mcols(genes)$type == "gene" & seqnames(genes) != "chrM"] # nuclear protein-coding genes
genes_bad <- genes[!mcols(genes)$type %in% c("CUTs", "SUTs", "ncRNA_gene", "HC", "MC", "LC", "Nascent")] # intervals to be excluded from DNC windows
win_coding <- suppressWarnings(flank(npcd, 100) %>% trim() %>% resize(width(.) + 500, "start"))
win_dnc <- suppressWarnings(flank(npcd, 500) %>% trim() %>% resize(width(.) + 100, "start"))
strand(win_dnc) <- ifelse(strand(win_dnc) == "+", "-", "+")
win_coding <- trim_by_down_or_upstream_features(win_coding, genes_bad)
win_dnc <- trim_by_down_or_upstream_features(win_dnc, genes_bad)

# Retain only pairs of windows longer than 200 bp:
good <- width(win_coding) >= 200 & width(win_dnc) >= 200
win_coding <- win_coding[good]
win_dnc <- win_dnc[good]

# Quantify NET-seq tags on the windows:
cm_coding <- get_overlapping_scores(win_coding, bg_data, value = "count_matrix")
cm_dnc <- get_overlapping_scores(win_dnc, bg_data, value = "count_matrix")
good2 <- rowSums(cm_coding) > 0 & rowSums(cm_dnc) > 0
cm_coding <- cm_coding[good2, ]
cm_dnc <- cm_dnc[good2, ]
win_coding <- win_coding[good2]
win_dnc <- win_dnc[good2]
ratio <- (cm_coding + 1) / (cm_dnc + 1)
# Calculate values for the X and Y axis (assuming that the columns 1:2 are cac2 and columns 3:4 are wt replicates):
x <- ratio[, 3:4] %>% rowMeans() %>% log2() %>% unname() 
y <- log2(rowMeans(ratio[, 1:2]) / rowMeans(ratio[, 3:4])) %>% unname()

tbl <- tibble(x = x, y = y, name = mcols(win_dnc)$name,  label = mcols(win_dnc)$name2, alpha = ifelse(is.na(label), 0.2, 1)) %>% arrange(alpha)

ggplot(tbl, aes(x= x, y = y, colour = label, fill = label, alpha = alpha)) + geom_point(shape = 21, size = 2) + 
  geom_hline(yintercept = 0, linetype = "dotted") + geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = -log2(1.5), colour = "red", linetype = "dotted") + xlim(-6, 12) + ylim(-4, 2) + theme_bw() + 
  geom_text(aes(label = label), nudge_y = 0.1) + 
  theme(legend.position = "none")
for (ext in c(".png", ".pdf")) {
  ggsave(paste0("SFig. 1C", ext), width = 7, height = 7, units = "in")
}

