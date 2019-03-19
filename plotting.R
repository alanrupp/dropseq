
# - Violin with treatment -----------------------------------------------------
violin_plot <- function(object, clusters = NULL, genes, tx = NULL,
                        jitter = TRUE) {
  if (is.null(clusters)) {
    clusters <- sort(unique(object@ident))
  }
  
  # check that genes & clusters are in object
  genes <- genes[genes %in% rownames(object@data)]
  clusters <- clusters[clusters %in% object@ident]
  
  # function to grab cells by treatment
  grab_cells <- function(clusters, tx = NULL) {
    by_cluster <- names(arc@ident)[arc@ident %in% clusters]
    if (!is.null(tx)) {
      by_tx <- rownames(object@meta.data)[object@meta.data$treatment %in% tx]
      by_cluster <- intersect(by_cluster, by_tx)
    }
    return(by_cluster)
  }
  
  cells1 <- grab_cells(clusters, tx[1])
  cells2 <- grab_cells(clusters, tx[1])
  
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

# - Column plot of mean gene expression ---------------------------------------
column_plot <- function(object, genes, clusters = NULL) {
  if (is.null(clusters)) {
    clusters <- levels(object@ident)
  }
  
  genes <- genes[genes %in% rownames(object@data)]
  
  # grab data
  if (length(genes) == 0) {
    return(NULL)
  } else if (length(genes) == 1) {
    df <- object@data[genes, ] %>% as.data.frame() %>% set_names(genes) %>% t()
  } else {
    df <- object@data[genes, ] %>% as.matrix() %>% as.data.frame()
  }
  
  # make tidy
  df <- df %>%
    rownames_to_column("gene") %>%
    gather(-gene, key = "barcode", value = "counts")
  
  # add cluster info
  df <- left_join(df, data.frame("barcode" = names(object@ident),
                                 "cluster" = object@ident), by = "barcode")
  # plot
  plt <- 
    df %>%
    group_by(gene, cluster) %>%
    summarize(mean = mean(counts)) %>%
    ggplot(aes(x = cluster, y = mean, fill = cluster)) +
    geom_col(show.legend = FALSE) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(y = "Mean expression", x = element_blank()) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  if (length(genes) > 1) {
    plt <- plt + facet_wrap(~gene, scales = "free_y")
  }
  return(plt)
}

# - Proportion plot -----------------------------------------------------------
proportion_plot <- function(object, genes, clusters = NULL) {
  if (is.null(clusters)) {
    clusters <- sort(unique(object@ident))
  }
  
  get_table <- function(gene) {
    table(object@ident, object@data[gene, ] > 0) %>%
      as.data.frame() %>%
      mutate("Gene" = gene)
  }
  
  df <- map(genes, ~get_table(.x)) %>% bind_rows()
  
  plt <-
    ggplot(df, aes(x = Var1, y = Freq, fill = Var2)) + 
    geom_col(show.legend = FALSE) + 
    scale_fill_manual(values = c('gray', 'red')) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Cluster", y = "Cells") +
    theme_bw() +
    theme(panel.grid = element_blank())
  if (length(genes) > 1) {
    plt <- plt + facet_wrap(~Gene)
  }
  return(plt)
}


# - Correlation plot ---------------------------------------------------------