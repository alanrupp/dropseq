library(tidyverse)
library(Seurat)

# find highly expressed genes -----------------------------------------------
find_expressed <- function(object, min_cells = 10, min_counts = 1) {
  calc <- rowSums(object@raw.data >= min_counts) >= min_cells
  calc <- calc[calc == TRUE]
  return(names(calc))
}


# Test multiple clusters and genes ------------------------------------------
de_tx_test <- function(object, clusters = NULL, genes = NULL, 
                       tx1, tx2, test) {
  # ensure test is available
  if (!(test %in% c("MAST", "Wilcoxon", "DESeq2"))) {
    cat(paste("test not available:", test, "\n"))
    stop()
  }
  
  # function to grab cells by cluster and treatment
  grab_cells <- function(cluster, treatment) {
    intersect(
      names(object@ident)[object@ident == cluster], 
      rownames(object@meta.data)[object@meta.data$treatment == treatment]
    )
  }
  
  # filter clusters
  if (is.null(clusters)) {
    clusters <- levels(object@ident)
  } else {
    clusters <- clusters[clusters %in% object@ident]
  }
  
  # filter genes
  if (is.null(genes)) {
    genes <- rownames(object@data)
  } else {
    genes <- genes[genes %in% rownames(object@data)]
  }
  
  # check that all clusters have cells with that treatment
  no_samples <- function() {
    tx_table <- table(object@ident, object@meta.data$treatment) %>%
      as.data.frame()
    tx_table <- filter(tx_table, Var1 %in% clusters & Var2 %in% c(tx1, tx2))
    tx_table %>% filter(Freq == 0) %>% .$Var1
  }
  
  # remove clusters with no treated samples
  clusters <- clusters[!(clusters %in% no_samples())]
  
  # - TESTS --------------------
  if (test == "MAST") {
    results <- 
      map(clusters,
          ~ MASTDETest(object,
                       cells.1 = grab_cells(.x, tx1),
                       cells.2 = grab_cells(.x, tx2),
                       genes.use = genes)
      )
  } else if (test == "DESeq2") {
    results <- 
      map(clusters,
          ~ DESeq2DETest(object,
                         cells.1 = grab_cells(.x, tx1),
                         cells.2 = grab_cells(.x, tx2),
                         genes.use = genes)
    )
  } else if (test == "Wilcoxon") {
    results <- 
      map(clusters,
          ~ WilcoxDETest(object,
                         cells.1 = grab_cells(.x, tx1),
                         cells.2 = grab_cells(.x, tx2),
                         genes.use = genes)
      )
  }
  
  # genes to column
  results <- map(results, ~ rownames_to_column(.x, "gene"))
  
  # add cluster info
  results <- map(1:length(results), 
                 ~mutate(results[[.x]], "cluster" = clusters[.x]))
  
  # correct p values
  results <- bind_rows(results) %>%
    mutate(q_val = p.adjust(p_val, method = "BH"))
  
  return(results)
}


# plot differences by treatment ---------------------------------------------
treatment_violin <- function(object, cluster, genes, tx1, tx2,
                             jitter = TRUE, tpm = FALSE) {
  # check that genes are in object
  genes <- genes[genes %in% rownames(object@data)]
  
  # grab cells by treatment
  grab_cells <- function(cluster, tx) {
    intersect(
      names(arc@ident)[arc@ident == cluster], 
      rownames(arc@meta.data)[arc@meta.data$treatment == tx]
    )
  }
  cells1 <- grab_cells(cluster, tx1)
  cells2 <- grab_cells(cluster, tx2)
  
  # grab data using input cells
  if (length(genes) == 1) {
    df <- object@data[genes, union(cells1, cells2)] %>%
      as.matrix() %>% t() %>% as.data.frame()
    rownames(df) <- genes
  } else {
    df <- object@data[genes, union(cells1, cells2)] %>%
      as.matrix() %>% as.data.frame()
  }
  
  # add treatment info
  df <- df %>%
    rownames_to_column("gene") %>%
    gather(-gene, key = "cell", value = "Counts") %>%
    mutate("Treatment" = case_when(
      cell %in% cells1 ~ tx1,
      cell %in% cells2 ~ tx2
    )) %>%
    mutate(Treatment = factor(Treatment, levels = c(tx1, tx2)))
  
  # plot
  plt <- ggplot(df, aes(x = Treatment, y = Counts, fill = Treatment)) +
    geom_violin(show.legend = FALSE) +
    theme_bw() +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.grid = element_blank())
  if (jitter == TRUE) {
    plt <- plt + geom_jitter(show.legend = FALSE, height = 0)
  }
  if (length(genes) > 1) {
    plt <- plt + facet_wrap(~gene)
  }
  
  plt
}

# - Finding coexpression of 2 genes -------------------------------------------
two_keys <- function(object, genes = NULL, clusters = NULL) {
  if (is.null(genes)) {
    genes <- rownames(object@data)
  } else {
    genes <- genes[genes %in% rownames(object@data)]
  }
  
  if (is.null(clusters)) {
    clusters <- unique(sort(object@ident))
  } else {
    clusters <- clusters[clusters %in% object@ident]
  }
  
  # remove duplicates
  cat("Generating search grid")
  search_grid <- expand.grid(genes, genes) %>% set_names("gene1", "gene2")
  dedup <- function(grid_df) {
    grid2 <- apply(grid_df, 1, sort) %>% t() %>% 
      as.data.frame() %>% mutate('combined' = paste(V1, V2))
    grid_df <- grid_df[!duplicated(grid2$combined), ]
    return(grid_df)
  }
  search_grid <- dedup(search_grid)
  search_grid <- search_grid[search_grid$gene1 != search_grid$gene2, ]
  
  # find coexpression rate by clusters vs. other clusters
  coexpression <- function(cluster, genes) {
    tbl <- table(object@ident, 
                 object@data[genes[1], ] > 0 & object@data[genes[2], ] > 0)
    tbl <- as.data.frame(tbl)
    tbl <- tbl %>% 
      group_by(Var1) %>%
      summarize("Prop" = sum(Freq[Var2 == TRUE]) / sum(Freq)) %>%
      mutate("Group" = ifelse(Var1 == cluster, "target", "other")) %>%
      select(Group, Prop) %>%
      group_by(Group) %>%
      summarize(Prop = max(Prop)) %>%
      spread(key = "Group", value = "Prop") %>%
      mutate("Gene1" = genes[1], "Gene2" = genes[2]) %>%
      mutate("Cluster" = cluster)
    return(tbl)
  }
  
  # apply by cluster then put together later
  by_cluster <- function(cluster) {
    cat(paste("Calculating cluster", cluster, "\n"))
    result <- apply(search_grid, 1, function(x) coexpression(cluster, c(x[1], x[2])))
    return(bind_rows(result))
  }
  
  results <- map(clusters, ~by_cluster(.x)) %>% bind_rows()
  results <- mutate(results, "Diff" = target - other) %>%
    filter(Diff > 0) %>%
    group_by(Cluster) %>%
    arrange(desc(Diff))
  return(results)
}
