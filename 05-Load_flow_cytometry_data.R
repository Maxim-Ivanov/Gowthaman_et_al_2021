# This script imports raw flow cytometry data, calculate some statistics (e.g. the directionality scores) and returns a data frame;
# It used the flowCore library to read FCS files which correspond to individual wells on 96- or 384-well plates;
# The FCS files are expected to be placed in subdirectories (the name of subdirectory is regarded as the plate name);
# Each subdirectory is expected to contain a sample information file;
# Sample info file is a Tab-delimited text file with three columns: i) well number (e.g. A01); ii) human-readable sample name (e.g. BY4741, or WT); iii) sample type.
# The allowed sample types are: S  for regular samples, N for negative controls, F for FPR controls, C for other controls;

if (!requireNamespace("flowCore", quietly = TRUE)) {
  BiocManager::install("flowCore")
  library(flowCore)
}

fc_dir <- "."
orc2_dir <- file.path(fc_dir, "ORC2_screen")
gcg1_dir <- file.path(fc_dir, "GCG1_screen")

##### Define custom functions ----------------------------------------------------------------------------------------------

# This function detects rows of matrix or data frame which contain outlier values:
quantile_filtering <- function(df, columns, quant) {
  # columns should be a character vector of column names, quant should be a list of lower and upper quantile value
  q1 <- min(quant[[1]], quant[[2]])
  q2 <- max(quant[[1]], quant[[2]])
  results <- matrix(nrow = nrow(df), ncol = length(columns))
  for (i in seq_along(columns)) {
    column <- columns[[i]]
    data <- df[, column]
    low <- quantile(data, q1, na.rm = TRUE)
    high <- quantile(data, q2, na.rm = TRUE)
    good <- data >= low & data <= high
    results[, i] <- good
  }
  out <- apply(results, MARGIN = 1, all)
  return(out)
}

# This function calculates well-level statistics from a single flowFrame object:
process_flowframe <- function(ff, cell_shape_filter, min_ncells, ch1, ch2) {
  mat <- exprs(ff) # convert flowFrame object to a matrix
  fcs_filename <- basename(keyword(ff)$FIL)
  well <- strsplit(fcs_filename, "_")[[1]][[4]] # extract the well number from the original filename
  # Filter cells by their shape (as determined by the front scatter and side scatter values):
  # (we assume that most of the cells are normal, thus we skip cells which have outlier FSC or SSC values)
  if (is.numeric(cell_shape_filter)) {
    good_shape <- quantile_filtering(mat, columns = c("FSC.A", "SSC.A"), quant = c(cell_shape_filter, 1 - cell_shape_filter))
    mat <- mat[good_shape, ] # matrix can become empty at this step...
  }
  # Skip cells which have non-meaningful YFP or mCherry fluorescence values:
  good_values <- mat[, ch1] > 0 & mat[, ch2] > 0 & !is.na(mat[, ch1]) & !is.na(mat[, ch2])
  mat <- mat[good_values, ]
  # If the number of cells after filtering is below the provided threshold, skip the well:
  if (is.numeric(min_ncells) && nrow(mat) < min_ncells) {
    return(data.frame("Well" = well, "N_cells" = nrow(mat), "YFP_well" = NA, "mCh_well" = NA, "r_well" = NA)) 
  }
  mcherry <- log10(mat[, ch1]) # use log10-transformed fluorescence values
  yfp <- log10(mat[, ch2])
  df <- data.frame("Well" = well, "N_cells" = nrow(mat), "YFP_well" = median(yfp), "mCh_well" = median(mcherry),  "r_well" = cor(yfp, mcherry, method = "spearman"))
  return(df)
}

# This function background normalizes the fluorescence values by negative control samples
# (then the negative control samples are removed from further consideration)
process_negative_controls <- function(df) {
  df_neg <- df[!df$Non_neg & !is.na(df$r_well), ]
  if (nrow(df_neg) > 0) {
    df$YFP_well <- df$YFP_well - mean(df_neg$YFP_well) # the subtraction can produce values below zero in negative samples
    df$mCh_well <- df$mCh_well - mean(df_neg$mCh_well)
    # Check that all fluorescence values in non-negative samples remain positive after subtraction:
    damaged <- !is.na(df$r_well) & df$Non_neg & (df$YFP_well < 0 | df$mCh_well < 0)
    if (any(damaged)) {
      warning("\t", sum(damaged), " non-negative samples were damaged by the background correction!")
    }
  } else {
    message("\tCannot subtract negative controls...")
  }
  return(df)
}

# This function calculates normal-based p-values:
add_pval_column <- function(df, rows, data_column, pval_column) {
  data <- df[rows, data_column]
  col_mean <- mean(data, na.rm = TRUE)
  col_sd <- sd(data, na.rm = TRUE)
  pvals <- pnorm(df[, data_column], col_mean, col_sd) # P[X <= x]
  df[, pval_column] <- ifelse(df[, data_column] <= col_mean, pvals, 1 - pvals)
  return(df)
}

# This function calculate plate-level statistics:
calc_plate_stats <- function(df, mode, regression_filter, min_fpr) {
  df$Ratio_well <- df$YFP_well / df$mCh_well
  chosen_rows <- !is.na(df$Ratio_well) & df$Non_neg & df$Ann # by default choose all wells which are not negative controls, contain sufficient number of cells and are annotated
  other_ctrl <- !is.na(df$Type) & df$Type == "C"
  chosen_rows <- chosen_rows & !other_ctrl # also exclude "Other controls", if any
  using_all <- TRUE
  if (mode == "all_wells") {
    message("\t\tUsing all samples in the group;")
  } else if (mode == "FPR_only") {
    fpr <- chosen_rows & !is.na(df$Type) & df$Type == "F" # limit the choice of rows to FPR samples only
    # Check that the number of FPR samples is enough to calculate the plate-level stats: 
    if (sum(fpr) >= min_fpr) {
      message("\t\tUsing FPR samples only;")
      chosen_rows <- fpr
      using_all <- FALSE
    } else {
      warning("\t\tOnly ", sum(fpr), " valid FPR samples found! Using all samples in the group instead;")
    }
  }
  # Calculate Group 3 stats:
  df$YFP_plate <- median(df$YFP_well[chosen_rows], na.rm = TRUE)
  df$mCh_plate <- median(df$mCh_well[chosen_rows], na.rm = TRUE)
  df$r_plate <- cor(df$YFP_well[chosen_rows], df$mCh_well[chosen_rows], method = "spearman", use = "na.or.complete")
  df$Ratio_plate <- df$YFP_plate / df$mCh_plate
  # Calculate the "geometric" distance (the "directionality score"):
  # (if using all samples, then skip the potential outliers which may disturb the regression line)
  clean_rows <- chosen_rows
  if (isTRUE(using_all) && is.numeric(regression_filter) && length(regression_filter) == 1) {
    decent_rows <- quantile_filtering(df[chosen_rows, ], columns = c("YFP_well", "mCh_well"), quant = c(regression_filter, 1 - regression_filter))
    clean_rows[chosen_rows] <- decent_rows
  }
  regression_line <- coefficients(lm(formula = mCh_well ~ YFP_well, data = df[clean_rows, ]))
  intercept <- regression_line[[1]]
  slope <- regression_line[[2]]
  df$Clean <- clean_rows
  # Dot product distance between a point and a line (signed):
  # (a negative Geom_dist value means that data point is located below the regression line, i.e. in the YFP half-space) 
  df$Geom_dist <- (slope * df$YFP_well - df$mCh_well + intercept) / sqrt(slope ^ 2 + 1)
  # Calculate other distance measures:
  df$Ratio_well_norm <- df$Ratio_well / df$Ratio_plate
  df$YFP_well_norm <- df$YFP_well / df$YFP_plate
  df$mCh_well_norm <- df$mCh_well / df$mCh_plate
  ## Calculate p-values (Group 5 stats):
  # P-value for the "geometric" distance (assuming Normal distribution of Geom_dist values around the regression line):
  df <- add_pval_column(df, rows = clean_rows, data_column = "Geom_dist", pval_column = "Pval_Geom")
  # P-value for the difference between Ratio_well and the plate mean of Ratio_well
  df <- add_pval_column(df, rows = clean_rows, data_column = "Ratio_well", pval_column = "Pval_Ratio")
  # P-value for the difference between YFP_well and its plate mean:
  df <- add_pval_column(df, rows = clean_rows, data_column = "YFP_well", pval_column = "Pval_YFP")
  # The same for mCherry:
  df <- add_pval_column(df, rows = clean_rows, data_column = "mCh_well", pval_column = "Pval_mCh")
  return(df)
}

read_sample_info <- function(sample_file) {
  out <- read.table(sample_file, sep = "\t", header = FALSE)
  out <- out[, 1:3] # skip extra columns if any
  names(out) <- c("Well", "Name", "Type") # the first three columns in the Tab-delimited sample file are expected to look like: "B08 snf5 S"
  return(out)
}

#' Load FCS data and calculate directionality scores
#'
#' @param dir path to the directory which contains FCS files (or subdirectories with FCS files).
#' @param subdirs If set to TRUE, the script looks for FCS files in subdirectories of \code{dir}. The names of subdirectories are considered as names of individual plates.
#' @param sample_info_file Name of the annotation file (should be placed into the directory with relevant FCS files).
#' @param cell_shape_filter The stringency of removal of odd-shaped cells based on the distribution of front scatter (FSC) and side scatter (SSC) values.
#' @param min_ncells A well is considered valid if the number of cells (after filtering by cell shape) does not drop below this number.
#' @param mode if \code{mode == "all_wells"} (default), the plate-level stats will be calculated on all wells (except the negative controls). 
#' Otherwise, if \code{mode == "FPR_only"}, the plate-level stats will be calculated only on samples annotated as FPR controls. 
#' This can be useful if the analyzed yeast samples were heavily selected (e.g. they are the best hits from a big screening).
#' @param subtract_negative If TRUE (default), then all fluorescense values are background corrected by the samples annotated as negative controls (N).
#' @param regression_filter Skip potential outlier observations when calculating the regression line. 
#' By default, \code{regression_filter = 0.05} which means that a reliable data point should have both YFP_well and mCh_well values between their 5% and 95% quantiles. 
#' When the regression line is calculated on FPR samples only (\code{mode == "FPR_only"}), the value of regression_filter is ignored.
#' To suppress the regression filter in the default "all_wells" mode, set \code{regression_filter ==  NULL}.
#' @param min_fpr The minimal number of valid FPR control samples to calculate the plate stats in the "FPR_only" mode. Ignored when \code{mode != "FPR_only"}.
#' @return Data frame with the following columns:
#' \enumerate {
#'     \item Identifiers and general information:
#'     \itemize {
#'         \item \code{Plate} Plate name (i.e. name of the directory containing FCS files);
#'         \item \code{Well} Well name (e.g. A01) extracted from the FCS file name;
#'         \item \code{Name} Human-readable name of the yeast sample (if available from the annotation file, otherwise NA);
#'         \item \code{Type} Sample type (S for regular samples, N for negative controls, F for FPR controls) if available from the annotation file (otherwise NA);
#'         \item \code{N_cells} Number of cells in given well (after filtering for proper cell shape and positive fluorescence values);
#'     }
#'     \item Statistics on individual wells:
#'     \itemize {
#'         \item \code{YFP_well} Median YFP value in given well;
#'         \item \code{mCh_well} Median mCherry value in given well;
#'         \item \code{r_well} Spearman correlation coefficient between YFP and mCherry signals within the well;
#'         \item \code{Ratio_well} Ratio between YFP_well and mCh_well;
#'     }
#'     \item Plate-level statistics:
#'     (observe that if \code{mode == "FPR_only"}, then plate-level stats will be calculated only on the FPR control samples, not on all wells)
#'     \itemize {
#'         \item \code{YFP_plate} Median YFP value on given plate (calculated as the median value of YFP_well);
#'         \item \code{mCh_plate} Median mCherry value on given plate (calculated as the median value of mCh_well);
#'         \item \code{r_plate} Spearman correlation coefficient between YFP_well and mCh_well values within the plate;
#'         \item \code{Ratio_plate} Ratio between YFP_plate and mCh_plate;
#'     }
#'     \item Directionality measures:
#'     \itemize {
#'         \item \code{Geom_dist} "Geometric" distance from given well to the regression line (the "Directionality score").
#'         (observe that if \code{mode == "FPR_only"}, then the regression line is calculated only on the FPR control samples, not on all wells);
#'         \item \code{Ratio_well_norm} = Ratio_well / Ratio_plate (i.e. the YFP/mCherry ratio for each well is expressed in units of the plate median YFP/mCherry ratio);
#'         \item \code{YFP_well_norm} = YFP_well / YFP_plate;
#'         \item \code{mCh_well_norm} = mCh_well / mCh_plate;
#'     }
#'     \item P-values (assuming Normal distribution):
#'     \itemize {
#'         \item \code{Pval_Geom} Statistical significance of \code{Geom_dist} (the Directionality score) assuming its Normal distributuion around the regression line;
#'         \item \code{Pval_Ratio} Significance of the difference between given \code{Ratio_well} and the plate mean of \code{Ratio_well} values;
#'         \item \code{Pval_YFP} Significance of the difference between given \code{YFP_well} and the plate mean of \code{YFP_well};
#'         \item \code{Pval_mCh} Significance of the difference between given \code{mCh_well} and the plate mean of \code{mCh_well};
#'     }
#'     \item Other columns:
#'     \itemize {
#'         \item \code{Non_neg} TRUE if the well was NOT annotated as negative control;
#'         \item \code{Clean} TRUE if the well passes the "regression filter" and thus can be safely used to calculate the regression line;
#'     }
#' }
load_plate_data <- function(dir = ".", subdirs = FALSE, sample_info_file = NULL, cell_shape_filter = 0.1, 
                            min_ncells = 100, mode = "all_wells", subtract_negative = TRUE, min_fpr = 5, 
                            regression_filter = 0.05, ch1 = "mCherry.H", ch2 = "B535.345.H") {
  # Count FCS files either in the provided directory ("dir"), or in its subdirectories (depending on value of "subdirs"):
  if (isTRUE(subdirs)) {
    subf <- list.dirs(dir, full.names = FALSE)
    subf <- subf[nchar(subf) > 0] # skip the first record which is always empty
    nfiles <- lapply(subf, function(x) { return(length(list.files(file.path(dir, x), pattern = "fcs$"))) }) # read names of FCS files in each subdirectory
    names(nfiles) <- subf
  } else {
    nfiles <- list(length(list.files(dir, pattern = "fcs$"))) # if subdirs == FALSE, nfiles is a list of length 1
  }
  nfiles <- nfiles[nfiles > 0] # skip directories which do not contain FCS files
  if (length(nfiles) == 0) {
    stop("There are no FCS files to load!")
  } else {
    message("There are ", sum(as.integer(nfiles)), " FCS files in ", length(nfiles), " directories;")
  }
  final_df_list <- vector("list", length(nfiles)) # if subdirs == FALSE, final_df_list is a list of length 1
  # Process all FCS files one by one:
  for (i in seq_along(nfiles)) {
    if (!is.null(names(nfiles))) {
      plate_name <- names(nfiles)[[i]] # define the plate name
      folder <- file.path(dir, plate_name) # define path to the FCS directory
    } else {
      if (dir == ".") {
        plate_name <- basename(getwd())
      } else {
        plate_name <- basename(dir)
      }
      folder <- dir
    }
    message("Reading ", nfiles[[i]], " FCS files from plate ", plate_name, ";")
    # Batch read FCS files in given directory:
    fs <- as(read.flowSet(path = folder, pattern = "fcs$", alter.names = TRUE), "list")
    # Extract well-level statistics from all FCS files in the directory (using custom function process_flowframe() defined above):
    message("\tCalculating well-level stats;")
    df_list <- unname(lapply(fs, process_flowframe, cell_shape_filter = cell_shape_filter, min_ncells = min_ncells, ch1 = ch1, ch2 = ch2))
    df <- do.call(rbind, df_list)
    # Add sample info, if available:
    if (is.character(sample_info_file) && length(sample_info_file) == 1) {
      sample_file <- file.path(folder, sample_info_file)
      if (file.exists(sample_file)) {
        message("\tAdding sample info;")
        mapping <- read_sample_info(sample_file)
        df <- merge(df, mapping, by = "Well", all.x = TRUE)
      } else {
        warning("\tCannot find sample info file ", sample_info_file, " in directory ", folder, "!")
        df$Name <- NA; df$Type <- NA
      }
    } else {
      df$Name <- NA; df$Type <- NA
    }
    df$Plate <- plate_name # add the plate name
    # Mark negative control samples:
    neg_wells <- !is.na(df$Type) & df$Type == "N"
    df$Non_neg <- !neg_wells
    # Subtract and mark negative control samples, if any:
    if (sum(neg_wells) > 0 && isTRUE(subtract_negative)) {
      message("\tBackground correction by ", sum(neg_wells), " negative samples;")
      df <- process_negative_controls(df)
    }
    # Mark annotated samples (found in the sample info file):
    df$Ann <- !is.na(df$Name)
    # Reorder columns in the data frame (just for better readability):
    before <- c("Plate", "Well", "Name", "Type")
    after <- c("N_cells", "YFP_well", "mCh_well", "r_well", "Non_neg", "Ann")
    if (is.null(df$Screen)) {
      df <- df[, c(before, after)]
    } else {
      df <- df[, c(before, "Screen", after)]
    }
    # Calculate plate-levels stats:
    if (is.null(df$Screen)) {
      message("\tCalculating plate-level stats on all wells;")
      df <- calc_plate_stats(df, mode = mode, min_fpr = min_fpr, regression_filter = regression_filter)
    } else {
      scr <- df$Ann & df$Screen != ""
      if (sum(scr) == 0) {
        warning("\t\tNo wells with known reporter name!")
      }
      uniq_scr <- unique(df$Screen[scr])
      if (length(uniq_scr) > 1) {
        last <- uniq_scr[[length(uniq_scr)]]
        other <- uniq_scr[1:(length(uniq_scr) - 1)]
        message("\tCalculating plate-level stats separately on ", paste(other, sep = ","), " and ", last, " reporters;")
        if (length(!scr) > 0) {
          df_noscr <- df[!scr, ]
          fake <- matrix(nrow = nrow(df_noscr), ncol = 14)
          df_noscr <- cbind(df_noscr, as.data.frame(fake))
        }
        df <- df[scr, ]
        by_obj <- by(df, list(df$Screen), calc_plate_stats, mode = mode, min_fpr = min_fpr, regression_filter = regression_filter, simplify = FALSE)
        df <- do.call(rbind, by_obj)
        rownames(df) <- NULL
        if (length(!scr) > 0) {
          colnames(df_noscr) <- colnames(df)
          df <- rbind(df, df_noscr)
        }
      } else {
        message("\tCalculating plate-level stats on all wells;")
        df <- calc_plate_stats(df, mode = mode, min_fpr = min_fpr, regression_filter = regression_filter)
      }
    }
    # Sort data frame by well names:
    well_letter <- substr(df$Well, 1, 1)
    well_number <- substr(df$Well, 2, 3)
    df <- df[order(well_letter, well_number), ]
    # Append the plate-specific data frame to final_df_list:
    final_df_list[[i]] <- df
  }
  # Combine results from different plates (n >= 1):
  out <- do.call(rbind, final_df_list)
  return(out)
}


##### Use the functions defines above to import the raw FCS files ---------------------------------------------------------------

df_orc2 <- load_plate_data(dir = orc2_dir, subdirs = TRUE, sample_info_file = "Info.txt")
df_gcg1 <- load_plate_data(dir = gcg1_dir, subdirs = TRUE, sample_info_file = "Info.txt")

# Save the data frames for future use:
saveRDS(df_orc2, "df_orc2.RDS")
saveRDS(df_gcg1, "df_gcg1.RDS")
