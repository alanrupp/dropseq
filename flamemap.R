flamemap <- function(object, genes, cells, n_bars = 100) {
  genes <- genes[genes %in% rownames(object@data)]
  cells <- cells[cells %in% colnames(object@data)]

  if (length(genes) == 0) {
    print("Invalid gene input")
    exit()
  }

  # grab cluster information and arrange by cluster name
  clusters <- data.frame(cell = names(object@ident), cluster = object@ident)
  clusters <- arrange(clusters, cluster, cell)

  # grab relevant data
  mtx <- object@data[genes, filter(clusters, cell %in% cells)$cell]

  # plot
  mtx %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    gather(-gene, key = "cell", value = "counts") %>%
    group_by(gene) %>%
    mutate("bar" = ntile(cell, n_bars)) %>%
    group_by(gene, bar) %>%
    arrange(desc(counts)) %>%
    ggplot(aes(x = bar, y = counts, fill = counts)) +
    geom_col() +
    scale_fill_gradient(low = "yellow", high = "red") +
    facet_wrap(~gene) +
    theme_bw() +
    theme(panel.grid = element_blank())
}
