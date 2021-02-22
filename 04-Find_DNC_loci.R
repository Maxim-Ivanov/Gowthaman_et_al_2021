# This script uses the improved annotation of yeast transcriptome for unbiased detection of DNC loci;
# DNC (Divergent Non-Coding transcription) loci were defined as promoters of known protein-coding genes 
# with experimental evidence for transcription in the divergent orientation;
# Transcription in coding and non-coding direction must start from the same NFR (Nucleosome-Free Region);

library(GenomicRanges)
library(tidyverse)
library(readxl)
library(devtools)

scripts <- c("trim_by_down_or_upstream_features")
for (script in scripts) {
  paste0("https://github.com/Maxim-Ivanov/Utility_functions/blob/main/", script, ".R?raw=TRUE") %>% devtools::source_url()
}

all_genes <- readRDS("Improved_annotation_for_Scerevisiae.RDS") # file produced by the 03-Correct_and_expand_SacCer3_annotation.R script
mcols(all_genes)$orig_coord <- NULL

##### Find CUTs, SUTs, novel HC/MC/LC genes and nascent-only intervals 
##### starting within 500 bp from TSS of a known protein-coding gene on the opposite strand ------------------------------

pcd_genes <- all_genes %>% `[`(mcols(.)$type == "gene" & seqnames(.) != "chrM")
good <- mcols(all_genes)$type %in% c("CUTs", "SUTs", "HC", "MC", "LC", "Nascent")
good_genes <- all_genes[good]
bad_genes <- all_genes[!good]

win <- suppressWarnings(flank(pcd_genes, 500) %>% trim())
mcols(win)$host <- granges(pcd_genes)
strand(win) <- ifelse(strand(win) == "+", "-", "+")
win_trim <- TranscriptomeReconstructoR:::trim_by_down_or_upstream_features(win, bad_genes, mode = "down") %>% `[`(width(.) > 0)

good_start <- resize(good_genes, 1, "start")
hits <- findOverlaps(win_trim, good_start)
win_par <- win_trim[queryHits(hits)]
good_par <- good_start[subjectHits(hits)]
start_diff <- pgap(resize(win_par, 0, "start"), resize(good_par, 0, "start")) %>% width()
best <- tapply(start_diff, queryHits(hits), TranscriptomeReconstructoR:::find_uniq_min, simplify = FALSE) %>% unlist() %>% unname()
hits <- hits[best]
win_par <- win_trim[queryHits(hits)]
dnc <- good_genes[subjectHits(hits)]

mcols(dnc)$host_name <- mcols(win_par)$name
mcols(dnc)$host_coord <- mcols(win_par)$host
mcols(dnc)$host_fpkm <- mcols(win_par)$fpkm

# Load NFR from Chereji et al., 2018 (PMID 29426353):
download.file("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5807854/bin/13059_2018_1398_MOESM2_ESM.xlsx", "Chereji2018_NFR.xlsx", method = "curl")
nfr_tbl <- read_excel("Chereji2018_NFR.xlsx")
nfr <- GRanges(seqnames = nfr_tbl$Chr, IRanges(start = nfr_tbl$`NDR Center`, width = 0), strand = "*", seqinfo = seqinfo(dnc)) %>% 
  resize(nfr_tbl$`NDR Width`, "center") %>% unique() %>% GenomicRanges::reduce()

# Load nucleosome positions from Jiang et al., 2009 (PMID 19814794):
# (download file https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2009-10-10-r109/MediaObjects/13059_2009_2259_MOESM1_ESM.ZIP and unzip it)
nps <- read_excel("Additional data file 1.xls", sheet = "Measured nucleosomes", range = "B31:D61141")
from <- nps$chrom %>% unique() %>% sort()
to <- seqlevels(dnc)[1:16]
for (i in seq_along(from)) {
  nps$chrom <- str_replace(nps$chrom, from[[i]], to[[i]])
}
nps <- suppressWarnings(GRanges(seqnames = nps$chrom, ranges=IRanges(nps$start, end=nps$end), seqinfo=seqinfo(dnc)) %>% trim())

# Skip weak nucleosomes which are found completely within NFRs:
nps <- nps[!overlapsAny(nps, nfr, type = "within")]

# Count strong nucleosomes in the interval between DNC TSS and host TSS:
between <- pgap(dnc, mcols(dnc)$host_coord, ignore.strand = TRUE)
hits <- findOverlaps(nps, between, type = "within")
tbl_1 <- hits %>% as_tibble() %>% group_by(subjectHits) %>% summarize(nucl_count = n()) %>% dplyr::rename(idx = subjectHits)
tbl_2 <- tibble(idx = 1:length(dnc))
tbl <- left_join(tbl_2, tbl_1, by = "idx") %>% mutate(nucl_count = ifelse(is.na(nucl_count), 0, nucl_count))

# Retain only DNC loci without intervening nucleosomes:
dnc <- dnc[tbl$nucl_count == 0]

# Save the results for future use:
saveRDS(dnc, "DNC_loci.RDS")
