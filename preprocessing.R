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
  genes <- genes %>% 
    arrange(desc(sum)) %>%
    filter(!duplicated(name)) %>%
    .$position
  
  # keep only positions left in genes df
  cat(paste(length(genes), "genes removed"))
  mtx2 <- mtx[genes, ]
  return(mtx2)
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