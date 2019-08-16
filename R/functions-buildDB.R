# buildDB ---------------------------------------------------------------------------------------------------------
#' @title Build a chemical reference database table using rda files
#'
#' @description Function generates a chemical reference database table from rda files in the selected directory (rda_dir).
#'
#' @param rda_dir \code{character} with absolute path to the directory with rda files.
#' @param out_dir \code{character} with absolute path to the directory where database file should be written.
#' @param peakgr_thr \code{numeric} specifying intensity threshold (100\% BPI) for peaks in a chemical spectrum. If set to more than 0, peaks below will be removed.
#' @param compound_thr \code{numeric} specifying intensity threshold (100\% BPI) for peaks in all compound's chemical spectra. If set to more than 0, peaks below will be removed.
#'
#' @details Each rda file must be written for each chemical compond separately and contain a \code{chem.file} list, comprising the following four entries:
#' \itemize{
#' \item \code{header} listing experimental details.
#' \item \code{id} specifying chemical compound identity.
#' \item \code{reference}
#' \item \code{analytical} listing all detected peaks, which were built into chemical spectra (i.e. peak groups).
#' }
#'
#' The database table contains columns: "peakid", "mz", "rt", "into", "chemid", "dbid", "dbname".
#'
#' @return Function returns generated database table and writes it to a csv file in the selected directory (out_dir).
#'
#' @export
#'
#' @seealso \code{\link{annotateDS}}
#'
buildDB <- function(
                    rda_dir = NULL,
                    out_dir = NULL,
                    peakgr_thr = 0,
                    compound_thr = 0) {
  if (is.null(rda_dir)) {
    stop("'rda_dir' is required!")
  }
  if (!dir.exists(rda_dir)) {
    stop("incorrect filepath for 'rda_dir' provided")
  }
  if (is.null(out_dir)) {
    stop("'out_dir' is required!")
  }
  if (!dir.exists(out_dir)) {
    stop("incorrect filepath for 'out_dir' provided")
  }
  if (peakgr_thr != 0) {
    warning("peakgr_thr value above 0 is selected. Recommended value is 0")
  }
  if (compound_thr != 0) {
    warning("compound_thr above 0 is selected. Recommended value is 0")
  }
  db_files <- list.files(path = rda_dir, pattern = "*.rda", full.names = TRUE)
  pb <- utils::txtProgressBar(min = 0, max = length(db_files), style = 3)
  db <- NULL ## initiate empty database

  #### ---- for each DB compound
  for (chem in 1:length(db_files)) {

    ## load and check rda file
    load(file = db_files[chem])
    if (!exists("chem.file")) {
      stop("rda file is empty: ", db_files[chem])
    }
    if (length(chem.file) != 4 |
      any(names(chem.file) != c("header", "id", "reference", "analytical"))) {
      stop("rda file is corrupted: ", db_files[chem])
    }

    #### ---- for each peak group in the compound's file
    compound <- NULL
    peakgrs <- length(chem.file$analytical)
    for (gr in 1:peakgrs) {
      ## normalise all features IT within the single peak group to %BPI and remove features with IT < peakgr_thr
      pkg <- chem.file$analytical[[gr]]
      pkg_spec <- pkg$spec
      into_norm <- ((pkg_spec[, "IT"]) / max(pkg_spec[, "IT"])) * 100
      peakgr <- setNames(pkg_spec[which(into_norm > peakgr_thr), c("MZ", "IT")], nm = c("mz", "into"))
      if (nrow(peakgr) > 1) {
        # add RT value for the peakgroup
        peakgr$rt <- pkg$rt
        # add peakgroup number
        peakgr$peakgr <- gr
        # append to final compound matrix
        compound <- rbind(compound, peakgr)
      } else {
        next()
      }
    }

    if (is.null(compound)) {
      ## omit a compound with no peak-groups with more than 1 peak
      next()
    }

    ## normalise all features IT within compound (all peak groups) to %BPI and remove spectrum features with IT < compound_thr
    into_norm <- ((compound[, "into"]) / max(compound[, "into"])) * 100
    compound <- compound[which(into_norm > compound_thr), ]

    ## retain the most intense peakgroup
    peakgrs_ints <- sapply(unique(compound$peakgr), function(x) {
      sum(compound[which(compound$peakgr == x), "into"])
    })
    peakgr_top <- unique(compound$peakgr)[which(peakgrs_ints == max(peakgrs_ints))]
    compound <- compound[compound$peakgr == peakgr_top, ]
    compound <- compound[, -c(match("peakgr", names(compound)))]

    ## add chemical DB  details
    compound$chemid <- chem
    compound$dbid <- chem.file$id$NPCID
    compound$dbname <- chem.file$id$name
    db <- rbind(db, compound)

    # remove table from the environment before loading next compound in next loop
    rm(chem.file)
    # update progress bar
    utils::setTxtProgressBar(pb, chem)
  }
  ## arrange peaks by chemid and into, and get unique peak numbers
  db <- db[order(db$chemid, -db$into), ]
  db$rt <- db$rt * 60
  db$peakid <- 1:nrow(db)
  db <- db[, c("peakid", "mz", "rt", "into", "chemid", "dbid", "dbname")]
  write.csv(db, file.path(out_dir, "database.csv"), quote = T, row.names = F)
  message("\n Database was created.")
}
