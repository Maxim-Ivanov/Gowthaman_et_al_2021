# This script uses published ONT Direct RNA-seq, CAGE-seq, 3'READS and NET-seq data to:
# 1) Correct borders of SacCer3 genes by the experimental TSS and PAS data;
# 2) Call novel (previously unannotated) genes and lncRNAs using the TranscriptomeReconstructoR package;

library(rtracklayer)
library(GenomicAlignments)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
library(tidyverse)
library(readxl)
library(devtools)

if (!"TranscriptomeReconstructoR" %in% rownames(installed.packages())) {
  devtools::install_github("Maxim-Ivanov/TranscriptomeReconstructoR", build_vignettes = TRUE, ref = "main")
}
library(TranscriptomeReconstructoR)

scripts <- c("batch_read_track_data", "get_overlapping_scores", "find_all_introns")
for (script in scripts) {
  paste0("https://github.com/Maxim-Ivanov/Utility_functions/blob/main/", script, ".R?raw=TRUE") %>% devtools::source_url()
}

##### Define custom functions ---------------------------------------------------------------------------------

adjust_tx <- function(tx, genes, bad, mode) {
  stopifnot(class(tx) == "GRanges" && class(genes) == "GRanges")
  stopifnot(class(bad) == "logical")
  stopifnot(mode %in% c("tss", "pas"))
  stopifnot(length(tx) == length(bad) && length(tx) == length(genes))
  out0 <- tx[!bad]
  tx2 <- tx[bad]
  genes2 <- genes[bad]
  if (mode == "tss") {
    new <- punion(resize(genes2, 1, "start"), resize(tx2, 1, "end"), fill.gap = TRUE)
  } else {
    new <- punion(resize(genes2, 1, "end"), resize(tx2, 1, "start"), fill.gap = TRUE)
  }
  names(new) <- names(tx2)
  ranges(tx2) <- ranges(new)
  out <- c(out0, tx2) %>% `[`(order(as.integer(names(.))))
  return(out)
}

trim_windows <- function(win, genes, idx, mode) {
  stopifnot(class(win) == "GRanges" && class(genes) == "GRanges")
  stopifnot(class(idx) == "integer")
  stopifnot(length(win) == length(idx) && length(win) == length(genes))
  stopifnot(mode %in% c("up", "down"))
  na <- is.na(idx)
  out0 <- win[na]
  #message(length(out0), " host genes have no neighbor;")
  win2 <- win[!na]
  idx2 <- idx[!na]
  neighbor <- genes[idx2]
  over <- poverlaps(win2, neighbor) %>% as.logical()
  out1 <- win2[!over]
  #message(length(out1), " windows did not overlap with the neighbor;")
  win2 <- win2[over]
  neighbor <- neighbor[over]
  win2_names <- names(win2) 
  if (mode == "up") {
    new_ranges <- pgap(neighbor, resize(win2, 0, "end")) %>% ranges()
  } else {
    new_ranges <- pgap(neighbor, resize(win2, 0, "start")) %>% ranges()
  }
  names(new_ranges) <- win2_names # otherwise the names of new_ranges are inherited from neighbor (which are shifted by 1)
  ranges(win2) <- new_ranges
  #message(length(win2), " windows were trimmed by the neighbor;")
  out <- c(out0, out1, win2) %>% `[`(order(as.integer(names(.))))
  return(out)
}

find_tc <- function(genes, win, tc, mode) {
  stopifnot(length(win) == length(genes))
  stopifnot(mode %in% c("tss", "pas"))
  old_names <- names(genes)
  names(genes) <- 1:length(genes)
  over <- win %over% tc
  # Add dummy TC to genes which do not overlap with any TC:
  g1 <- genes[!over]
  if (mode == "tss") {
    dummy <- resize(g1, 1, "start") %>% granges() %>% unname()
  } else {
    dummy <- resize(g1, 1, "end") %>% granges() %>% unname()
  }
  mcols(g1)[, mode] <- dummy
  mcols(g1)[, paste0(mode, "_summit")] <- dummy
  mcols(g1)[, paste0(mode, "_real")] <- FALSE
  # Find the strongest TC for the other genes:
  g2 <- genes[over]
  w2 <- win[over]
  hits <- findOverlaps(w2, tc)
  scores <- score(tc)[subjectHits(hits)]
  best <- tapply(scores, queryHits(hits), TranscriptomeReconstructoR:::find_uniq_max, simplify = FALSE) %>% unlist() %>% unname()
  hits <- hits[best]
  tc_best <- tc[subjectHits(hits)]
  mcols(g2)[, mode] <- granges(tc_best) %>% unname()
  mcols(g2)[, paste0(mode, "_summit")] <- mcols(tc_best)$thick %>% granges() %>% unname()
  mcols(g2)[, paste0(mode, "_real")] <- TRUE
  out <- c(g1, g2)
  out <- out[order(as.integer(names(out)))]
  names(out) <- old_names
  return(out)
}

adjust_known_genes <- function(genes, long_reads, tss, pas, adjust = "summit", offset = 100) {
  stopifnot(grepl("GRangesList", class(long_reads)))
  stopifnot(class(genes) == "GRanges")
  stopifnot(class(tss) == "GRanges")
  stopifnot(!is.null(score(tss)))
  stopifnot(class(pas) == "GRanges")
  stopifnot(!is.null(score(pas)))
  stopifnot(adjust %in% c("summit", "outer_border"))
  stopifnot(is.numeric(offset) && length(offset) == 1 && offset >= 0)
  old_names <- names(genes)
  names(genes) <- 1:length(genes)
  nano <- long_reads %>% range() %>% unlist()
  # For each Nanopore read, find its best mate among the called genes:
  hits <- findOverlaps(nano, genes)
  par1 <- nano[queryHits(hits)]
  par2 <- genes[subjectHits(hits)]
  overlap <- width(pintersect(par1, par2)) / width(par1)
  best_mate <- tapply(overlap, queryHits(hits), TranscriptomeReconstructoR:::find_uniq_max, simplify = FALSE) %>% unlist() %>% unname()
  hits <- hits[best_mate]
  # Group best mate reads by gene and merge them into transcribed intervals:
  tbl <- tibble(a = subjectHits(hits), b = queryHits(hits)) %>% arrange(a)
  h2 <- Hits(from = tbl$a, to = tbl$b, nLnode = length(genes), nRnode = length(nano))
  tx <- nano[subjectHits(h2)] %>% split(queryHits(h2)) %>% range() %>% unlist() %>% unname()
  names(tx) <- genes[unique(queryHits(h2))] %>% names()
  # Add dummy intervals for the genes w/o Nanopore coverage:
  tx <- genes[-unique(queryHits(h2))] %>% granges() %>% c(tx) %>% `[`(order(as.integer(names(.)))) # tx is parallel to genes
  # If annotated TSS/PAS is located within another gene, force tx to the annotated gene border:
  bad_tss <- countOverlaps(resize(genes, 1, "start"), genes) > 1
  bad_pas <- countOverlaps(resize(genes, 1, "end"), genes) > 1
  tx <- adjust_tx(tx, genes, bad_tss, mode = "tss")
  tx <- adjust_tx(tx, genes, bad_pas, mode = "pas")
  # Generate windows between annotated and empirical gene starts (ends):
  win_up <- punion(resize(genes, 0, "start"), resize(tx, 0, "start"), fill.gap = TRUE)
  win_down <- punion(resize(genes, 0, "end"), resize(tx, 0, "end"), fill.gap = TRUE)
  # Extend the windows by <offset> bp each side (to increase the chance of hitting the relevant TSS/PAS):
  win_up_ext <- suppressWarnings(trim(win_up %>% resize(width(.) + offset*2, "center")))
  win_down_ext <- suppressWarnings(trim(win_down %>% resize(width(.) + offset*2, "center")))
  # If the TSS and PAS windows of the same gene do overlap, resize them to the middle of their overlapping area:
  over <- poverlaps(win_up_ext, win_down_ext)
  w1_up <- win_up_ext[!over]
  w1_down <- win_down_ext[!over]
  w2_up <- win_up_ext[over]
  w2_down <- win_down_ext[over]
  area <- pintersect(w2_up, w2_down)
  middle <- resize(area, 1, "center")
  w2_up <- w2_up %>% resize(0, "start") %>% pgap(middle)
  w2_down <- w2_down %>% resize(0, "end") %>% pgap(middle)
  win_up_ext <- c(w1_up, w2_up) %>% `[`(order(as.integer(names(.))))
  win_down_ext <- c(w1_down, w2_down) %>% `[`(order(as.integer(names(.))))
  # Trim windows by the nearest up- and downstream genes:
  idx_up <- follow(genes, genes)
  idx_down <- precede(genes, genes)
  win_up_ext <- trim_windows(win_up_ext, genes, idx_up, mode = "up")
  win_down_ext <- trim_windows(win_down_ext, genes, idx_down, mode = "down")
  # Search for TSS (PAS) in the extended windows, choose the strongest one (if any):
  genes <- find_tc(genes, win_up_ext, tss, mode = "tss")
  genes <- find_tc(genes, win_down_ext, pas, mode = "pas")
  # Adjust annotated genes to the TSS/PAS:
  mcols(genes)$orig_coord <- granges(genes) %>% unname()
  if (adjust == "outer_border") {
    new_rng <- punion(mcols(genes)$tss, mcols(genes)$pas, fill.gap = TRUE) %>% ranges()
  } else {
    new_rng <- punion(mcols(genes)$tss_summit, mcols(genes)$pas_summit, fill.gap = TRUE) %>% ranges()
  }
  ranges(genes) <- new_rng
  names(genes) <- old_names
  return(genes)
}

add_tss_or_pas_info <- function(genes, tc, mode) {
  stopifnot(mode %in% c("tss", "pas"))
  old_names <- names(genes)
  names(genes) <- 1:length(genes)
  term_base <- resize(genes, 1, ifelse(mode == "tss", "start", "end")) %>% granges()
  over <- term_base %over% tc
  out1 <- genes[!over]
  if (length(out1) > 0) {
    mcols(out1)[, mode] <- out1 %>% granges() %>% unname() %>% resize(1, ifelse(mode == "tss", "start", "end"))
    mcols(out1)[, paste0(mode, "_summit")] <- mcols(out1)[, mode]
  }
  out2 <- genes[over]
  if (length(out2) > 0) {
    tc_par <- tc[findOverlaps(out2, tc, select = "first")]
    mcols(out2)[, mode] <- tc_par %>% granges() %>% unname()
    mcols(out2)[, paste0(mode, "_summit")] <- mcols(tc_par)$thick
  }
  out <- c(out1, out2) %>% `[`(order(as.integer(names(.))))
  names(out) <- old_names
  return(out)
}


##### Load SacCer3 genes from SGD -----------------------------------------------------------------------

# (the annotation available from Bioconductor library TxDb.Scerevisiae.UCSC.sacCer3.sgdGene is not good, because it lacks ncRNA genes)
download.file("http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/archive/saccharomyces_cerevisiae.20170114.gff.gz", "SacCer3_SGD.gff.gz")
gff <- import("SacCer3_SGD.gff.gz", format = "GFF3")
mcols(gff) <- mcols(gff)[, c("type", "Name")]
names(mcols(gff)) <- names(mcols(gff)) %>% str_replace("Name", "name")
gff <- gff %>% sortSeqlevels() %>% sort()
seqinfo(gff, new2old = as.integer(1:17)) <- seqinfo(txdb)
genes_gff <- gff[grepl("gene", mcols(gff)$type)]
mcols(genes_gff)$type <- mcols(genes_gff)$type %>% droplevels()

# Load CUTs and SUTs from Xu et al. 2009 (PMID 19169243):
# (download archive from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2766638/bin/NIHMS137055-supplement-supZip.zip and unzip it)
xu <- read_excel("SI_Tab3.xls", sheet = 1) %>% dplyr::select(chr, start, end, strand, type, name) %>% filter(type %in% c("CUTs", "SUTs")) # Supplementary Table 3
xu <- GRanges(seqnames = xu$chr, IRanges(xu$start, end = xu$end), strand = xu$strand, type = xu$type, name = xu$name)
xu <- xu %>% sortSeqlevels() %>% sort()
seqinfo(xu, new2old = c(1:16, NA)) <- seqinfo(txdb)

# Update coordinates from sacCer2 to sacCer3:
# (download chain file from http://hgdownload.soe.ucsc.edu/goldenPath/sacCer2/liftOver/sacCer2ToSacCer3.over.chain.gz and gunzip it)
ch <- import.chain("sacCer2ToSacCer3.over.chain")
xu_upd <- liftOver(xu, ch) %>% range() %>% unlist()
mcols(xu_upd) <- mcols(xu)

# Combine CUT/SUTs with SacCer3 genes:
genes <- c(genes_gff, xu_upd) %>% sort()
names(genes) <- mcols(genes)$name


##### Load and pre-process the experimental data --------------------------------------------------------------

# Call TSS from yeast CAGE-seq data:
# (to produce CAGE-seq Begraph files, run https://github.com/Maxim-Ivanov/Reanalysis_of_published_datasets/blob/main/CAGE-seq/Lu_2019_(PMID_31076411).sh)
cage_dir <- "." # change to the directory with CAGE-seq Bedgraph files returned by the Lu_2019_(PMID_31076411).sh script
cage_files <- list.files(cage_dir, pattern = "^Lu2019.*rep.\\.bedgraph\\.gz$")
cage_data <- batch_read_track_data(cage_files, dir = cage_dir, format = "bedGraph", seqinfo = seqinfo(txdb))
names(cage_data) <- names(cage_data) %>% str_replace(".bedgraph.gz", "")
tss <- TranscriptomeReconstructoR::call_TCs(cage_data, min_support = 2)

# Call PAS from yeast 3'READS data:
# (to produce 3'READS Begraph files, run https://github.com/Maxim-Ivanov/Reanalysis_of_published_datasets/blob/main/3'READS/Liu_2018_(PMID_28916539).sh)
pas_dir <- "." # change to the directory with 3'READS Bedgraph files returned by the Liu_2018_(PMID_28916539).sh script
pas_files <- list.files(pas_dir, pattern = "^Liu2018.*rep.\\.bedgraph\\.gz$")
pas_data <- batch_read_track_data(pas_files, dir = pas_dir, format = "bedGraph", seqinfo = seqinfo(txdb))
names(pas_data) <- names(pas_data) %>% str_replace(".bedgraph.gz", "")
pas <- TranscriptomeReconstructoR::call_TCs(pas_data, min_support = 2)

# Call transcribed intervals from the merged NET-seq data (coverage of whole reads):
netseq_dir_1 <- "." # change to the directory with NET-seq Bedgraph files returned by the 01-Remapping_published_NET-seq_datasets.sh script
wr_file <- "NETseq_merged_wt_whole_read.bedgraph.gz"
netseq_wr_data <- batch_read_track_data(list(wr_file), dir = netseq_dir_1, format = "bedGraph", seqinfo = seqinfo(txdb))
names(netseq_wr_data) <- names(netseq_wr_data) %>% str_replace(".bedgraph.gz", "")
netseq <- call_transcribed_intervals(netseq_wr_data[[1]], min_signal = 10, max_gapwidth = 100, min_width = 200)
transcribed <- netseq[[1]]
transcribed_gaps <- netseq[[2]]

# Also load NET-seq Bedgraph file with reads truncated to the first bases (for FPKM calculations):
# i) Merged NET-seq WT tracks from the published studies:
merged_fb_file <- "NETseq_merged_wt_first_base_norm1M.bedgraph.gz"
netseq_merged_fb_data <- batch_read_track_data(list(merged_fb_file), dir = netseq_dir_1, format = "bedGraph", seqinfo = seqinfo(txdb))
names(netseq_merged_fb_data) <- names(netseq_merged_fb_data) %>% str_replace("_norm1M.bedgraph.gz", "")
# ii) Individual replicates of our NET-seq tracks (wt, hda1-3, gcn5, ngg1, ada2, ubp8, sgf73):
netseq_dir_2 <- "." # change to the directory with NET-seq Bedgraph files returned by the 02-Alignment_of_novel_yeast_NET-seq_data.sh script
bg_files <- list.files(netseq_dir_2, pattern = "^NETseq_s\\d{2}.*rep[12].*norm1M.bedgraph.gz$")
our_netseq_data <- batch_read_track_data(bg_files, dir = netseq_dir_2, format = "bedGraph", seqinfo = seqinfo(txdb))
names(our_netseq_data) <- names(our_netseq_data) %>% str_replace("_norm1M.bedgraph.gz", "") %>% str_replace("NETseq_s\\d{2}_", "")

# Load ONT Direct RNA-seq reads:
# (to produce RDS files, run Garalde_2018_(PMID_29334379).sh followed by Filter_Direct_RNA-seq_BAM_files.R)
# (both scripts are available from https://github.com/Maxim-Ivanov/Reanalysis_of_published_datasets/blob/main/ONT_Direct_RNA-seq/)
rds_dir <- "." # change to the directory with RDS files returned by the Filter_Direct_RNA-seq_BAM_files.R script
rds_files <- list.files(rds_dir, pattern = "^Garalde2018.*rep.*RDS$")
ont_ga <- file.path(rds_dir, rds_files) %>% lapply(readRDS) %>% Reduce(c, .) %>% sort()

# Find all splice sites:
ebg <- exonsBy(txdb, by = "gene")
introns <- find_all_introns(ebg)
introns <- introns[elementNROWS(introns) > 0]
introns_unl <- unlist(introns)
donor <- resize(introns_unl, 1, "start") %>% resize(6, "center")
acc <- resize(introns_unl, 1, "end") %>% resize(6, "center")

# Skip split ONT reads that do not match known introns:
unspl <- njunc(ont_ga) == 0
ont_spl <- ont_ga[!unspl] %>% grglist()
gaps <- psetdiff(unlist(range(ont_spl)), ont_spl)
group <- rep(1:length(gaps), times = elementNROWS(gaps))
gap_start <- resize(gaps, 1, "start") %>% unlist()
gap_end <- resize(gaps, 1, "end") %>% unlist()
over_donor <- gap_start %over% donor
over_acc <- gap_end %over% acc
good_start <- tapply(over_donor, group, all, simplify = FALSE) %>% as.logical()
good_end <- tapply(over_acc, group, all, simplify = FALSE) %>% as.logical()
unspl[!unspl] <- good_start & good_end
ont_ga <- ont_ga[unspl]

# Convert GAlignments to GRangesList:
ont <- ont_ga %>% unname() %>% grglist()

# Also skip fusion ONT reads (covering at least 2 disjoint known genes by more than 50% each):
fusion_idx <- ont %>% range() %>% unlist() %>% TranscriptomeReconstructoR:::find_fusion_tx(genes, skip = FALSE)
if (length(fusion_idx) > 0) {
  message(length(fusion_idx), " fusion transcripts skipped;")
  ont <- ont[-fusion_idx]
}


##### Correct and extend the SacCer3/CUT/SUT annotation by the experimental data -----------------------------

# Adjust known genes to summits of nearby TSS and PAS (connected by ONT reads):
genes_adj <- adjust_known_genes(genes, ont, tss, pas)

# Skip Nanopore reads strongly overlapping with adjusted genes:
ont_free <- ont %>% TranscriptomeReconstructoR:::add_mcols_to_grl() %>% unlist(use.names = FALSE) %>% 
  TranscriptomeReconstructoR:::find_free_reads(genes_adj, max_overlap_called = 0.1, min_read_width = 1000) %>% 
  split(mcols(.)$read_id)

# Call novel genes from the remaining reads and the unused TSS/PAS:
tss_free <- tss[tss %outside% resize(genes_adj, 1, "start")]
pas_free <- pas[pas %outside% resize(genes_adj, 1, "end")]
out <- ont_free %>% TranscriptomeReconstructoR::extend_long_reads_to_TSS_and_PAS(tss_free, pas_free) %>% 
  TranscriptomeReconstructoR::call_transcripts_and_genes()
hml_genes <- out[[1]]
reads_free <- out[[5]]

# Combine adjusted and novel genes:
names(genes_adj) <- NULL
mcols(genes_adj) <- mcols(genes_adj) %>% subset(select = c("type", "name", "orig_coord", "tss", "tss_summit", "pas", "pas_summit"))
score(genes_adj) <- 0
mcols(hml_genes)$orig_coord <- granges(hml_genes)
hml_genes <- hml_genes %>% 
  add_tss_or_pas_info(tss, mode = "tss") %>% 
  add_tss_or_pas_info(pas, mode = "pas")
all_genes <- c(genes_adj, hml_genes) %>% sort()

# Augment genes with NET-seq data:
out2 <- TranscriptomeReconstructoR::process_nascent_intervals(all_genes, transcribed, tss, pas, reads_free, transcribed_gaps)
all_genes <- out2[[1]]
nascent_tx <- out2[[3]]

# Skip nascent-only transcripts called within used TSS or PAS:
tss_used <- tss[tss %over% resize(all_genes, 1, "start")]
pas_used <- pas[pas %over% resize(all_genes, 1, "end")]
nascent_tx <- nascent_tx %>% 
  TranscriptomeReconstructoR:::trim_by_down_or_upstream_features(tss_used, mode = "down", offset = 20) %>% 
  TranscriptomeReconstructoR:::trim_by_down_or_upstream_features(pas_used, mode = "up", offset = 20) %>% 
  `[`(width(.) >= 500)

# Combine all_genes and nascent_tx:
mcols(all_genes) <- mcols(all_genes) %>% subset(select = c("type", "name", "orig_coord"))
mcols(nascent_tx) <- mcols(nascent_tx) %>% subset(select = c("type", "name"))
mcols(nascent_tx)$orig_coord <- nascent_tx %>% granges() %>% unname() # dummy column for compatibility with all_genes
final_genes <- c(all_genes, nascent_tx) %>% sort()

# Count NET-seq FPKM on the merged samples from all published studies:
mcols(final_genes)$fpkm <- get_overlapping_scores(final_genes, netseq_merged_fb_data, value = "count_matrix") %>% 
  as.numeric() %>% `/`(width(final_genes)) %>% `*`(1e03) %>% round(1)

# Save the new yeast annotation for future use by 04-Find-DNC_loci.R:
saveRDS(final_genes, "Improved_annotation_for_Scerevisiae.RDS")

# Also count FPKM values on our NET-seq data (for Fig. 2D, Fig. 3B-C, SFig. 2B):
cm <- get_overlapping_scores(final_genes, our_netseq_data, value = "count_matrix") %>% apply(2, `/`, width(final_genes)) %>% `*`(1e03) %>% round(1)
mcols(final_genes)$orig_coord <- NULL
final_genes %>% as_tibble() %>% bind_cols(as_tibble(cm)) %>% write_tsv("NETseq_FPKM_all_genes.txt")
