# ============================================================
# HotSpot Null Simulation - SLURM Array Version
# Each SLURM array task runs ONE null iteration
# ============================================================

cat("Job started at:", format(Sys.time()), "\n")

# ---------------------------
# RNG setup (array-specific)
# ---------------------------
RNGkind("L'Ecuyer-CMRG")

array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(array_id)) {
  stop("SLURM_ARRAY_TASK_ID not found. Are you running as a SLURM array job?")
}

set.seed(100000 + array_id)

cat("Running array task:", array_id, "\n")

# ---------------------------
# Load data
# ---------------------------
load("Locs_1-22-26.rda")   # must contain object: Locs

# ---------------------------
# Prepare base data
# ---------------------------
Locs_base <- Locs[, 1:2]

# ---------------------------
# Helper functions
# ---------------------------
assign2bins <- function(dat, bins) {
  binborders <- c(bins$start, bins$end)
  binborders <- binborders[order(binborders)]

  assigned_bins <- cut(
    dat,
    c(-Inf, binborders, Inf),
    labels = c(
      rbind("none", rownames(bins))[1:(2 * nrow(bins))],
      "none"
    )
  )

  return(assigned_bins)
}

FindHotSpots <- function(locations,
                         WindowSize = 500,
                         EffGenomeSize = 2182075994) {

  out <- NULL

  Lambda <- nrow(locations) * WindowSize / EffGenomeSize
  Threshold <- qpois(10^-5, Lambda, lower.tail = FALSE)

  for (chr in unique(locations$chr)) {

    xx <- locations$position[locations$chr == chr]

    BinnedLocs <- NULL

    for (offset in 0:(WindowSize - 1)) {
      BinnedLocs <- c(
        BinnedLocs,
        floor((xx - offset) / WindowSize) * WindowSize + offset
      )
    }

    BinnedLocs <- table(BinnedLocs)
    SigBins <- BinnedLocs[BinnedLocs >= Threshold]

    if (length(SigBins) > 0) {

      # Merge overlapping significant bins
      SigLocs <- as.numeric(names(SigBins))
      ends <- which(diff(SigLocs) >= WindowSize)
SigLocs <- data.frame(
        chr = chr,
        start = SigLocs[c(1, ends + 1)],
        end = SigLocs[c(ends, length(SigLocs))] + WindowSize - 1
      )

      # Trim hot spots to insertion range
      assigned_bins <- assign2bins(xx, SigLocs)

      SigLocs[, 2:3] <- t(
        matrix(
          unlist(by(xx, assigned_bins, range)[rownames(SigLocs)]),
          nrow = 2
        )
      )

      SigLocs$insertions <- table(assigned_bins)[rownames(SigLocs)]
      rownames(SigLocs) <- paste(chr, rownames(SigLocs), sep = "_hs")

      out <- rbind(out, SigLocs)
    }
  }

  return(out)
}

# ---------------------------
# Run ONE null simulation
# ---------------------------

LocsNULL <- Locs_base

LocsNULL$position <- Locs_base$position +
  round((runif(nrow(Locs_base)) - 0.5) * 1e6)

HSnull <- FindHotSpots(LocsNULL)

# ---------------------------
# Save output (unique file)
# ---------------------------

outfile <- paste0("HotSpotNULL_", array_id, ".rda")
save(HSnull, file = outfile)

cat("Saved:", outfile, "\n")
cat("Job finished at:", format(Sys.time()), "\n")
