treatment_violin <- function(object, cluster, genes, tx1, tx2,
                             jitter = TRUE) {
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