# Functions for preprocessing Seurat objects

# - Fix duplicated cell names -------------------------------------------------
dedup_barcodes <- function(mtx) {
  if (sum(duplicated(colnames(mtx))) > 0) {
    barcodes <- colnames(mtx)
    dup_codes <- which(duplicated(barcodes))
    barcodes[dup_codes] <- paste0(barcodes[dup_codes], "x")
    colnames(mtx) <- barcodes
    cat(paste(length(dup_codes), "barcodes replaced"))
  }
  return(mtx)
}

# - Remove gene names that have lower expression than duplicate name ---------
dedup_genes <- function(mtx) {
  # calculate rowSums
  genes <- data.frame(
    "name" = rownames(mtx),
    "sum" = rowSums(mtx),
    "position" = seq(nrow(mtx))
  )
  
  # arrange by decreasing expression
  genes_to_keep <- genes %>% 
    arrange(desc(sum)) %>%
    filter(!duplicated(name)) %>%
    .$position
  
  # keep only positions left in genes df
  cat(paste(length(genes) - length(genes_to_keep), "genes removed"))
  mtx2 <- mtx[genes, ]
  return(mtx2)
}

# - Turn counts into TPM ------------------------------------------------------
make_tpm <- function(mtx) {
  mtx <- apply(mtx, 2, function(x) log2(x / sum(x) * 10^6 + 1))
  mtx <- Matrix(mtx, sparse = TRUE)
  return(mtx)
}

# - Wipe clean metadata slot --------------------------------------------------
wipe_metadata <- function(object) {
  object@meta.data <- data.frame("nGene" = rep(NA, ncol(object@raw.data)))
  rownames(object@meta.data) <- colnames(object@raw.data)
  return(object)
}


# - Calculate percent mitochondrial ------------------------------------------
percent_mito <- function(object) {
  mito_genes <- str_detect(rownames(object@raw.data), "^mt")
  mito <- 
    Matrix::colSums(object@raw.data[mito_genes, ]) /
    Matrix::colSums(object@raw.data)
  object@meta.data$percent_mito <- mito
  return(object)
}
