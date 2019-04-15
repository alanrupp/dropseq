# - Summarize data ------------------------------------------------------------
summarize_data <- function(object, genes, clusters = NULL) {
  if (is.null(clusters)) {
    clusters <- unique(object@ident)
  } else {
    clusters <- clusters[clusters %in% object@ident]
  }
  
  # keep valid genes
  genes <- genes[genes %in% rownames(object@data)]
  
  # grab data
  if (length(genes) == 0) {
    return(NULL)
  } else if (length(genes) == 1) {
    df <- object@data[genes, ] %>% as.data.frame() %>% t() %>% as.data.frame()
    rownames(df) <- genes
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
  
  # filter by selected clusters
  df <- filter(df, cluster %in% clusters)
  
  # summarize
  df <- df %>%
    group_by(gene, cluster) %>%
    summarize(avg = mean(counts), prop = sum(counts > 0)/n()) %>%
    ungroup()
  
  return(df)
}


# - Heatmap plot --------------------------------------------------------------
order_clusters <- function(object, cells = NULL, genes = NULL, scale = TRUE,
                          col_names = FALSE, row_names = FALSE, 
                          heatmap_legend = FALSE) {
  if (is.null(cells)) {
    cells <- colnames(object@data)
  }
  if (is.null(genes)) {
    genes <- rownames(object@data)
  } 
  
  cells <- cells[cells %in% colnames(object@data)]
  genes <- genes[genes %in% rownames(object@data)]
  
  if (scale == TRUE) {
    matrix <- object@scale.data[genes, ]
  } else {
    matrix <- object@data[genes, ]
  }
  
  # functions to grab cells and calculate mean expression
  mean_expression <- function(matrix, cells) {
    rowMeans(matrix[, cells])
  }
  cells_use <- function(object, ident) {
    names(object@ident)[object@ident == ident]
  }
  
  # calculate mean expression by cluster
  mean_values <-
    map(unique(object@ident),
        ~mean_expression(matrix, cells = cells_use(object, .x))) %>%
    as.data.frame() %>%
    set_names(unique(object@ident)) %>%
    t()
  
  # auto clip high and low to center at 0
  auto_clip <- function(mtx) {
    values <- as.matrix(mtx)
    if (abs(min(values)) < max(values)) {
      clip <- abs(min(values))
      mtx <- apply(mtx, c(1,2), function(x) ifelse(x > clip, clip, x))
    } else if (abs(min(values)) < max(values)) {
      clip <- max(values)
      mtx <- apply(mtx, c(1,2), function(x) ifelse(x < -clip, -clip, x))
    }
    return(mtx)
  }
  
  if (scale == TRUE) {
    mean_values <- auto_clip(mean_values)
  }
  
  dists <- dist(mean_values)
  clusts <- hclust(dists)
  
  return(clusts$labels[clusts$order])
}


heatmap_plot <- function(object, genes = NULL, cells = NULL, scale = TRUE,
                         cluster_order = NULL,
                         col_names = FALSE, row_names = FALSE, 
                         heatmap_legend = FALSE) {
  
  if (is.null(cells)) {
    cells <- colnames(object@data)
  }
  if (is.null(genes)) {
    genes <- rownames(object@data)
  } 
  
  cells <- cells[cells %in% colnames(object@data)]
  genes <- genes[genes %in% rownames(object@data)]
  
  if (scale == TRUE) {
    matrix <- object@scale.data[genes, ]
    # make red --> white --> blue color palette
    endcolors <- c("firebrick3", "white", "dodgerblue3")
    color_pal <- c(colorRampPalette(c(endcolors[1], endcolors[2]))(50), 
                   colorRampPalette(c(endcolors[2], endcolors[3]))(51)[-1])
  } else {
    matrix <- object@data[genes, ]
    # make white --> blue color palette
    endcolors <- c("white", "dodgerblue3")
    color_pal <- colorRampPalette(c(endcolors[1], endcolors[2]))(100)
  }
  
  # functions to grab cells and calculate mean expression
  mean_expression <- function(matrix, cells) {
    rowMeans(matrix[, cells])
  }
  cells_use <- function(object, ident) {
    names(object@ident)[object@ident == ident]
  }
  
  # calculate mean expression by cluster
  mean_values <-
    map(unique(object@ident),
        ~mean_expression(matrix, cells = cells_use(object, .x))) %>%
    as.data.frame() %>%
    set_names(unique(object@ident)) %>%
    t()
  
  # auto clip high and low to center at 0
  auto_clip <- function(mtx) {
    values <- as.matrix(mtx)
    if (abs(min(values)) < max(values)) {
      clip <- abs(min(values))
      mtx <- apply(mtx, c(1,2), function(x) ifelse(x > clip, clip, x))
    } else if (abs(min(values)) < max(values)) {
      clip <- max(values)
      mtx <- apply(mtx, c(1,2), function(x) ifelse(x < -clip, -clip, x))
    }
    return(mtx)
  }
  
  if (scale == TRUE) {
    mean_values <- auto_clip(mean_values)
  }
  
  # plot
  plt <-
    pheatmap::pheatmap(mean_values, 
                       show_colnames = col_names,
                       treeheight_col = 0,
                       legend = heatmap_legend,
                       color = color_pal)
  return(plt)
}


# - Violin plot ---------------------------------------------------------------
violin_plot <- function(object, genes, tx = NULL, clusters = NULL, 
                        jitter = TRUE, stacked = FALSE, order_genes = FALSE,
                        ncol = NULL) {
  if (is.null(clusters)) {
    clusters <- sort(unique(object@ident))
  }
  
  # check that genes & clusters are in object
  genes <- genes[genes %in% rownames(object@data)]
  clusters <- clusters[clusters %in% object@ident]
  
  # function to grab cells by treatment
  grab_cells <- function(clusters) {
    names(object@ident)[object@ident %in% clusters]
  }
  
  # grab data using input cells
  if (length(genes) == 1) {
    df <- object@data[genes, grab_cells(clusters)] %>%
      as.matrix() %>% t() %>% as.data.frame()
    rownames(df) <- genes
  } else {
    df <- object@data[genes, grab_cells(clusters)] %>%
      as.matrix() %>% as.data.frame()
  }
  
  # make tidy
  df <- df %>%
    rownames_to_column("gene") %>%
    gather(-gene, key = "cell", value = "Counts")
  
  if (order_genes == TRUE) {
    df <- mutate(df, gene = factor(gene, levels = genes))
  }
  
  # add cluster info
  cluster_df <- data.frame("cell" = names(object@ident),
                           "Cluster" = object@ident)
  df <- left_join(df, cluster_df, by = "cell")
  
  # plot ------------
  # add treatment info
  if (!is.null(tx)) {
    tx_df <- data.frame("Treatment" = object@meta.data$treatment) %>%
      rownames_to_column("cell") %>%
      filter(Treatment %in% tx)
    df <- inner_join(df, tx_df, by = "cell") %>%
      mutate(Treatment = factor(Treatment, levels = tx))
    plt <- ggplot(df, aes(x = Treatment, y = Counts, fill = Treatment)) +
      facet_wrap(~Cluster)
  } else {
    plt <- ggplot(df, aes(x = Cluster, y = Counts, fill = Cluster))
  }
  
  # generic plot attributes
  plt <- plt + 
    geom_violin(show.legend = FALSE, scale = "width") +
    theme_bw() +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.grid = element_blank())
  
  # add jitter
  if (jitter == TRUE) {
    plt <- plt + geom_jitter(show.legend = FALSE, height = 0, alpha = 0.2,
                             stroke = 0)
  }
  
  # facet wrap for multiple genes or multiple genes + treatments
  if (length(genes) > 1) {
    if (!is.null(tx)) {
      plt <- plt + facet_wrap(~gene + Cluster, scales = "free_y")
    } else {
      plt <- plt + facet_wrap(~gene, scales = "free_y")
    }
  }
  return(plt)
}


# - Column plot of mean gene expression ---------------------------------------
column_plot <- function(object, genes, clusters = NULL, 
                        reorder_clusters = FALSE) {
  df <- summarize_data(object, genes, clusters)
  
  # plot
  plt <- 
    ggplot(df, aes(x = cluster, y = counts, fill = cluster)) +
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
  } else {
    clusters <- clusters[clusters %in% object@ident]
  }
  
  get_table <- function(gene) {
    table(object@ident, object@data[gene, ] > 0) %>%
      as.data.frame() %>%
      mutate("Gene" = gene)
  }
  
  df <- map(genes, ~get_table(.x)) %>% bind_rows()
  
  # keep only specified clusters
  df <- filter(df, Var1 %in% clusters)
  
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


# - Flamemap plot ------------------------------------------------------------
flamemap <- function(object, genes, cells = NULL, n_bars = 100,
                     order_genes = FALSE, cluster_labels = TRUE,
                     icicle = FALSE) {
  genes <- genes[genes %in% rownames(object@data)]
  if (is.null(cells)) {
    cells <- colnames(object@data)
  } else {
    cells <- cells[cells %in% colnames(object@data)]
  }
  
  if (length(genes) == 0) {
    print("Invalid gene input")
    exit()
  }
  
  # grab cluster information and arrange by cluster name
  clusters <- data.frame(cell = names(object@ident), cluster = object@ident)
  clusters <- arrange(clusters, cluster, cell)
  clusters <- filter(clusters, cell %in% cells)
  
  # grab cluster info for plotting ticks on x axis
  scale_factor <- nrow(clusters)/n_bars
  cluster_stats <- clusters %>% mutate("pos" = seq(n())) %>% 
    group_by(cluster) %>% 
    summarize(max = max(pos), mid = mean(pos)) %>%
    mutate_if(is.numeric, ~ .x / scale_factor)
    
  # grab relevant data
  mtx <- object@data[genes, filter(clusters, cell %in% cells)$cell]
  
  # reshape
  df <- 
    mtx %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    gather(-gene, key = "cell", value = "counts") %>%
    mutate(cell = factor(cell, levels = filter(clusters, cell %in% cells)$cell)) %>%
    group_by(gene) %>%
    mutate("bar" = ntile(cell, n_bars)) %>%
    group_by(gene, bar) %>%
    arrange(desc(counts)) %>%
    ungroup() %>%
    group_by(gene, bar) %>%
    mutate("height" = as.numeric(counts > 1) / n())
  
  # order genes by user-defined order
  if (order_genes == TRUE) {
    df <- df %>%
      mutate(gene = factor(gene, levels = genes))
  }
  
  if (icicle == TRUE) {
    df <- mutate(df, counts = -counts, height = -height)
  }
  
  # plot
  plt <-
    ggplot(df, aes(x = bar, y = height, fill = counts)) +
    geom_col(show.legend = FALSE) +
    geom_hline(aes(yintercept = 0)) +
    scale_x_continuous(expand = c(0, 0), breaks = c(1, cluster_stats$max+0.5)) +
    labs(y = element_blank(), x = element_blank()) +
    facet_wrap(~gene) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.text.x = element_blank(),
          axis.line.x = element_blank())
  
  if (icicle == TRUE) {
    plt <- plt + scale_fill_gradient(low = "royalblue", high = "aquamarine") +
      geom_text(data = cluster_stats,
                aes(x = mid, y = 0.015, label = cluster), inherit.aes = FALSE) +
      scale_y_continuous(expand = c(0, 0), limits = c(-1, 0.03),
                         labels = seq(1, 0, by = -0.25))
  } else {
    plt <- plt + scale_fill_gradient(low = "yellow", high = "red") +
      geom_text(data = cluster_stats,
                aes(x = mid, y = -0.015, label = cluster), inherit.aes = FALSE)  +
      scale_y_continuous(expand = c(0, 0), limits = c(-0.03, 1))
  }
  
  return(plt)
}


# - UMAP plot ----------------------------------------------------------------
umap_plot <- function(object, genes = NULL, cells = NULL, clusters = NULL, 
                      legend = FALSE, cluster_label = FALSE,
                      ncol = NULL, xlim = NULL, ylim = NULL) {
  # pull UMAP data
  umap <- data.frame(
    UMAP1 = object@dr$umap@cell.embeddings[, 1],
    UMAP2 = object@dr$umap@cell.embeddings[, 2]
  ) %>% rownames_to_column("barcode")
  
  if (is.null(clusters)) {
    clusters <- sort(unique(object@ident))
  } else {
    cluster_bool <- object@ident %in% clusters
    clusters <- clusters[cluster_bool]
    umap <- umap[cluster_bool, ]
  }
  
  pull_data <- function(genes) {
    if (length(genes) == 1) {
      df <- object@data[genes, ] %>% as.data.frame() %>% set_names(genes)
    } else {
      df <- object@data[genes, ] %>% as.matrix() %>% t() %>% as.data.frame()
    }
    return(df)
  }
  
  results <- pull_data(genes) %>%
    rownames_to_column("barcode") %>%
    left_join(., umap, by = "barcode") %>%
    select(-barcode) %>%
    gather(-starts_with("UMAP"), key = "gene", value = "value")
  
  # plot
  plt <-
    ggplot(results, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = value), show.legend = legend, stroke = 0) +
    scale_color_gradient(low = "gray90", high = "navyblue") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    facet_wrap(~gene, ncol = ncol)
  if (cluster_label == TRUE) {
    cluster_center <- 
      left_join(umap, data.frame("barcode" = names(object@ident),
                                 "cluster" = object@ident), by = "barcode") %>%
      select(-barcode)
    cluster_center <- cluster_center %>%
      group_by(cluster) %>%
      summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))
    plt <- plt + geom_text(data = cluster_center, aes(label = cluster))
  }
  if (!is.null(xlim)) {
    plt <- plt + xlim(xlim)
  }
  if (!is.null(ylim)) {
    plt <- plt + ylim(ylim)
  }
  return(plt)
}


# - Dot plot ------------------------------------------------------------------
dot_plot <- function(object, genes, clusters = NULL, 
                     gene_order = FALSE,
                     cluster_order = FALSE,
                     cluster_labels = FALSE) {
  
  # get summary data
  df <- summarize_data(object, genes, clusters)
  
  if (cluster_order == TRUE) {
    df <- df %>%
      mutate(cluster = factor(cluster, levels = clusters))
  }
  
  if (gene_order == TRUE) {
    df <- df %>%
      mutate(gene = factor(gene, levels = genes))
  }
  
  # plot
  plt <- 
    ggplot(df, aes(x = gene, y = fct_rev(cluster), size = avg,
                   color = prop)) +
    geom_point(show.legend = FALSE) +
    scale_x_discrete(position = "top") +
    scale_radius(limits = c(0.01, NA)) +
    scale_color_gradient(low = "blue", high = "red") +
    theme_void()
  
  if (cluster_labels == TRUE) {
    plt <- plt + theme(axis.text.x = element_text(color = "black", angle = 89, 
                                                  vjust = 0.5, hjust = 0.5),
                       axis.text.y = element_text())
  } else {
    plt <- plt + theme(axis.text.x = element_text(color = "black", angle = 89, vjust = 1,
                                     hjust = ))
  }
  
  return(plt)
}


# - Correlation plot ---------------------------------------------------------
correlation_plot <- function(object, genes, clusters = NULL, 
                             view_by_cluster = FALSE,
                             fit_line = FALSE) {
  if (length(genes) != 2) {
    cat("Please specify exactly 2 genes")
    return(NULL)
  }
  
  if (is.null(clusters)) {
    clusters <- sort(unique(object@ident))
  } else {
    clusters <- clusters[clusters %in% object@ident]
  }
  
  cluster_df <- data.frame("barcode" = names(object@ident),
                           "Cluster" = object@ident) %>%
    filter(Cluster %in% clusters)
  
  df <- object@data[genes, ] %>% as.matrix() %>% t() %>% as.data.frame() %>%
    rownames_to_column("barcode")
  
  # add cluster ids
  df <- inner_join(df, cluster_df, by = "barcode") %>% select(-barcode)
  
  #return(select(df, !!genes[1]))
  
  # plot
  plt <- ggplot(df, aes(x = !!sym(genes[1]), y = !!sym(genes[2]))) +
    geom_point(aes(color = Cluster)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    theme(panel.grid = element_blank())
  if (fit_line == TRUE) {
    plt <- plt +  geom_smooth(method = "lm", se = FALSE, 
                              linetype = "dashed", color = "gray")
  }
  if (view_by_cluster == TRUE) {
    plt <- plt + facet_wrap(~Cluster)
  }
  return(plt)
}


# - Coexpression plot --------------------------------------------------------
coexpression_plot <- function(object, genes, clusters = NULL, title = TRUE) {
  if (is.null(clusters)) {
    clusters <- unique(object@ident)
  } else {
    clusters <- clusters[clusters %in% object@ident]
  }
  
  get_table <- function(object, genes) {
    table(object@ident, 
          object@data[genes[1],] > 0 & object@data[genes[2],] > 0) %>%
      as.data.frame()
  }
  
  df <- get_table(object, genes) %>%
    group_by(Var1) %>%
    summarize("Proportion" = Freq[Var2 == TRUE] / sum(Freq))
  
  plt <- ggplot(df, aes(x = Var1, y = "Coexpression", fill = Proportion)) +
    geom_tile() +
    scale_fill_viridis_c() + 
    labs(x = "Cluster", y = paste(genes[1], "+", genes[2])) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  if (title == TRUE) {
    plt <- plt + ggtitle(paste("Coexpression of", genes[1], "and", genes[2]))
  }
  return(plt)
}
