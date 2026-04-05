run_pls_pipeline <- function(
    base,
    ahba_file,
    pheno_file,
    pheno_name,
    spins_file,
    pathways_file,
    ncomp = 5,
    nperm = 1000,
    nboot = 1000,
    MinSizeGO = 15,
    MaxSizeGO = 500
) {
  
  library(tidyverse)
  library(data.table)
  library(xlsx)
  library(tidyr)
  library(svMisc)
  library(plsdepot)
  library(fgsea)
  
  # -----------------------------------------------
  message("Step 1: Loading AHBA, phenotype, and spins")
  
  AHBA <- read.csv(paste0(base, ahba_file))
  cat("Number of rows in AHBA:", nrow(AHBA), "\n")
  
  # Identify rows with all NA except ROI
  removed_row_numbers <- which(rowSums(is.na(AHBA[, -1])) >= (ncol(AHBA) - 1))
  AHBA <- AHBA[rowSums(is.na(AHBA[, -1])) < (ncol(AHBA) - 1), ]
  cat("Number of removed rows in AHBA:", length(removed_row_numbers), "\n")
  
  data <- read.csv(paste0(base, pheno_file))
  cat("Number of rows in data:", nrow(data), "\n")
  
  spins <- read.csv(paste0(base, spins_file))
  max_index <- nrow(spins) + 1
  message("Max index in spins: ", max_index)
  
  if (colnames(data)[1] != "ROI") {
    stop("Error: First column of 'data' must be 'ROI'")
  }
  if (colnames(AHBA)[1] != "ROI") {
    stop("Error: First column of 'AHBA' must be 'ROI'")
  }

  missing_in_data <- setdiff(AHBA$ROI, data$ROI)
  missing_in_ahba <- setdiff(data$ROI, AHBA$ROI)

  if (length(missing_in_data) == 0 && length(missing_in_ahba) == 0) {
    cat("All ROI values match between data and AHBA.\n")
  } else {
    cat("ROIs in AHBA but not in data:\n")
    print(missing_in_data)
  
    cat("ROIs in data but not in AHBA:\n")
    print(missing_in_ahba)
  
    stop("Error: ROI values do not match between data and AHBA. Execution stopped.")
  }
  
  PLS_data <- inner_join(data, AHBA, by = "ROI")
  
  PLS_data <- inner_join(data, AHBA, by = "ROI")
  Y <- PLS_data[, pheno_name, drop = FALSE]
  X <- PLS_data[, !(colnames(PLS_data) %in% c("ROI", pheno_name))]
  
  # -----------------------------------------------
  message("Step 2: Running PLS model")
  
  pls.model <- plsreg1(
    predictors = X,
    response = Y,
    comps = ncomp,
    crosval = TRUE
  )
  
  variance <- as.data.frame(pls.model$R2)
  variance$Component <- 1:ncomp
  names(variance)[1] <- "R2"
  
  # Realign components
  R1 <- cor(pls.model$x.scores, X)
  for (c in 1:ncomp) {
    if (abs(max(R1[c,])) < abs(min(R1[c,]))) {
      pls.model$x.scores[, c] <- -pls.model$x.scores[, c]
      pls.model$mod.wgs[, c] <- -pls.model$mod.wgs[, c]
    }
  }
  
  # -----------------------------------------------
  message("Step 3: Spin permutations")
  
  perm.var.expl <- NULL
  
  for (col in colnames(spins)[1:nperm]) {
    
    permuted_indices <- spins[[col]][spins[[col]] < max_index]
    
    permuted_indices <- setdiff(permuted_indices, removed_row_numbers)
    
    perm.model <- plsreg1(
      X,
      data[permuted_indices, pheno_name],
      comps = ncomp,
      crosval = TRUE
    )
    
    perm.var.expl <- rbind(perm.model$R2, perm.var.expl)
  }
  
  Spin.pvalue <- sapply(1:ncomp, function(c)
    sum(pls.model$R2[c] < perm.var.expl[, c]) / nperm
  )
  
  Results_PLS_Spin <- cbind(variance, Spin.pvalue)
  write.csv(Results_PLS_Spin,
            paste0(base, "Results_PLS_Spin.csv"),
            row.names = FALSE)
  
  # -----------------------------------------------
  message("Step 4: Bootstrapping PLS weights")
  
  npred <- ncol(X)
  boot.weights <- array(NA, dim = c(npred, ncomp, nboot))
  
  for (b in 1:nboot) {
    set.seed(b)
    bootorder <- sample(1:nrow(Y), nrow(Y), TRUE)
    Xrand <- X[bootorder, ]
    Yrand <- Y[bootorder, ]
    
    boot.model <- plsreg1(Xrand, Yrand, comps = ncomp, crosval = FALSE)
    
    for (c in 1:ncomp) {
      if (cor(boot.model$mod.wgs[, c], pls.model$mod.wgs[, c]) < 0) {
        boot.model$mod.wgs[, c] <- -boot.model$mod.wgs[, c]
      }
    }
    
    boot.weights[, , b] <- boot.model$mod.wgs
  }
  
  corr.weights <- matrix(NA, nrow = npred, ncol = ncomp)
  for (c in 1:ncomp) {
    boot.sd <- apply(boot.weights[, c, ], 1, sd)
    corr.weights[, c] <- pls.model$mod.wgs[, c] / boot.sd
  }
  
  Results_PLS_cWeights <- data.frame(
    name = colnames(X),
    corr.weights
  )
  colnames(Results_PLS_cWeights)[-1] <- as.character(1:ncomp)
  
  write.csv(Results_PLS_cWeights,
            paste0(base, "Results_PLS_cWeights.csv"),
            row.names = FALSE)
  
  # -----------------------------------------------
  message("Step 5: FGSEA")

  # GO should be in .gmx.txt format
  Pathways <- read.csv(
    paste0(base, pathways_file),
    header = TRUE,
    sep = "\t"
  )
  
  # first PLS weights are used for  FGSEA
  Stats <- setNames(Results_PLS_cWeights[, 2], Results_PLS_cWeights[, 1])
  
  Results_FGSEA <- fgsea(
    pathways = Pathways,
    stats = Stats,
    minSize = MinSizeGO,
    maxSize = MaxSizeGO
  )
  
  Results_FGSEA <- as.data.frame(Results_FGSEA)
  Results_FGSEA$leadingEdge <- sapply(
    Results_FGSEA$leadingEdge,
    paste,
    collapse = ","
  )
  
  write.csv(Results_FGSEA,
            paste0(base, "Results_FGSEA.csv"),
            row.names = FALSE)
  
  message("Pipeline complete")
}
