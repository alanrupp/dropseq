# convert gene_id to gene_name
id_to_name <- function(gene_ids) {
  features %>%
    filter(X1 %in% gene_ids) %>%
    .$X2
}

# convert gene_id to gene_name
name_to_id <- function(gene_names) {
  features %>%
    filter(X2 %in% gene_names) %>%
    .$X1
}
  
# filter markers dataframe for top genes
sig_markers <- function(df, n = NULL, p = 0.05) {
  genes <-
    df %>%
    filter(p_val_adj < p) %>%
    arrange(p_val_adj) %>%
    filter(!duplicated(gene)) %>%
    arrange(cluster)
  if (!is.null(n)) {
    genes <- 
      genes %>%
      arrange(p_val_adj) %>%
      slice(1:n)
  }
  return(genes$gene)
}



cells_by_cluster <- function(object, cluster) {
  names(object@ident)[object@ident == cluster]
}

# - Plotting ------------------------------------------------------------------
heatmap <- function(object, cells, genes, palette = c("seq", "YlGnBu")) {
  groups <- unique(object@ident)
  map(groups, ~rowMeans(object@data[genes, cells])) %>%
  as.data.frame() %>%
  set_names(groups) %>%
  t() %>%
  pheatmap::pheatmap(., 
                     show_colnames = FALSE,
                     color = RColorBrewer::brewer.pal(palette),
                     legend = FALSE)
}



# heatmap function
heatmap <- function(object, cells, genes, scale = TRUE) {
  if (scale == TRUE) {
    matrix <- object@scale.data
  } else {
    matrix <- object@data
  }
  
  mean_expression <- function(matrix, cells, genes) {
    rowMeans(matrix[genes, cells])
  }
  
  cells_use <- function(object, ident) {
    names(object@ident)[object@ident == ident]
  }
  
  map(unique(object@ident),
      ~mean_expression(matrix[, cells_use(object, .x)])) %>%
    as.data.frame() %>%
    set_names(unique(object@ident)) %>%
    t() %>%
    pheatmap::pheatmap(., show_colnames = FALSE,
                       legend = FALSE,
                       color = RColorBrewer::brewer.pal("div", "RdYlBu"))
    
  
}

genes_on_umap <- function(object, genes, legend = FALSE) {
  ids <- name_to_id(genes)
  # use gene_name, not gene_id
  object@data[ids, ] %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    #set_names(genes) %>%
    bind_cols(., data.frame("UMAP1" = object@dr$umap@cell.embeddings[, 1],
                            "UMAP2" = object@dr$umap@cell.embeddings[, 2])) %>%
    gather(-starts_with("UMAP"), key = "gene", value = "counts") %>%
    ggplot(aes(x = UMAP1, y = UMAP2, color = counts)) +
    geom_point(show.legend = legend, stroke = 0) +
    theme_classic() +
    scale_color_continuous(low = "gray90", high = "navyblue") +
    facet_wrap(~gene)
}

