library(tidyverse)
library(Seurat)

grab_areas <- function(areas = NULL) {
  # bring in all available areas
  all_areas <- list.files(path = "~/Desktop/Allen_Institute", 
                          pattern = "*.csv",
                          full.names = TRUE)
  area_names <- str_extract(all_areas, "[A-Za-z]+(?=_v_)")
  
  if (is.null(areas)) {
    areas <- area_names
  }
  
  # catch if entered area is not in names of available areas
  map(areas,
      ~ ifelse(.x %in% area_names, "", 
               print(paste("Warning:", .x, "not available"))))
  areas <- areas[areas %in% area_names]
  return(areas)
}

# from selected areas, grab enrichment scores
grab_area_enrichment <- function(areas) {
  map(areas, ~read_csv(paste0("~/Desktop/Allen_Institute/", .x, "_v_HY.csv")))
}


# get average scaled difference for all genes within cluster
enrichment_score <- function(object, cluster, areas = NULL) {

  # grab area enrichment scores
  available_areas <- grab_areas(areas)
  enrichment <- grab_area_enrichment(available_areas)[[1]]
  
  genes <- unique(enrichment$"gene-symbol")
  genes <- genes[genes %in% rownames(object@scale.data)]
  cells <- names(object@ident)[object@ident == cluster]
  cells <- cells[cells %in% colnames(object@scale.data)]
  
  # calculate mean expression
  rowMeans(object@scale.data[genes, cells])

}

# score for multiple areas
score_multiple <- function(object, cluster, areas = NULL) {
  if (is.null(areas)) {
    areas <- grab_areas()
  }
  score <- map(areas, ~enrichment_score(object, cluster, .x))
  score <- map(score, ~as.data.frame(.x))
  score <- map(score, ~rename(.x, "score" = ".x"))
  map(1:length(score),
      ~mutate(score[[.x]], "area" = areas[.x])) %>%
    bind_rows()
}

# get scores across all clusters and areas
full_scores <- function(object, clusters = NULL, areas) {
  if (is.null(clusters)) {
    clusters = unique(object@ident)
    clusters = sort(clusters)
  }
  result <- map(clusters, ~ score_multiple(neurons, .x))
  result <- map(1:length(result), ~ mutate(result[[.x]], "cluster" = clusters[.x]))
  return(result)
}

# summarize the scorse by mean z-score
summarize_scores <- function(scores) {
  scores %>%
    bind_rows() %>%
    group_by(area, cluster) %>%
    summarize(score = mean(score)) %>%
    arrange(desc(score)) %>%
    ungroup()
}

# get list of brain areas above 0 z-score
predict_brain_area <- function(scores_summary) {
  scores_summary %>%
    filter(score > 0) %>%
    group_by(cluster) %>%
    summarize(areas = list(area)) %>%
    mutate(areas = unlist(areas))
}

# plot heatmap
heatmap_plot <- function(object, genes = NULL, cells = NULL, scores,
                         scale_clip_low = -3, scale_clip_high = 3) {
  if (is.null(genes)) {
    genes <- object@var.genes
  }
  if (is.null(cells)) {
    cells <- sample(colnames(object@data), 1000)
  }
  
  # calculate mean z-score by gene and cluster
  cells_use <- function(object, cluster) {
    cell_names <- names(object@ident)[object@ident == cluster]
    cell_names <- cell_names[cell_names %in% cells]
  }
  
  mean_scale_data <-
    map(unique(object@ident),
        ~ rowMeans(object@scale.data[genes, cells_use(object, .x)])
        ) %>%
    as.data.frame() %>%
    set_names(unique(object@ident))
  
  # clip data with big outliers
  mean_scale_data <- apply(mean_scale_data, c(1,2),
                           function(x) ifelse(x > scale_clip_high, scale_clip_high, x))
  mean_scale_data <- apply(mean_scale_data, c(1,2),
                           function(x) ifelse(x < scale_clip_low, scale_clip_low, x))
  
  # cluster using euclidean distance
  distances <- dist(t(mean_scale_data))
  clusters <- hclust(distances)
  
  cluster_order <- clusters$label[clusters$order]

  dendro <- 
    ggdendro::ggdendrogram(clusters, rotate = TRUE) +
    theme_void()
  
  # make heatmap
  tiles <- 
    mean_scale_data %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    gather(-gene, key = "cluster", value = "expr") %>%
    mutate(cluster = factor(cluster, levels = cluster_order)) %>%
    ggplot(aes(x = gene, y = cluster, fill = expr)) +
    geom_tile(show.legend = FALSE) +
    scale_fill_distiller(palette = "RdYlBu", direction = 1) +
    theme_void() +
    theme(axis.text.y = element_text())
  
  # add score plot
  score_plot <-
    scores %>%
    filter(score > 0) %>%
    group_by(cluster) %>%
    mutate(position = seq(1:n())) %>%
    ungroup() %>%
    complete(cluster) %>%
    mutate(cluster = factor(cluster, levels = cluster_order)) %>%
    ggplot(aes(x = position, y = cluster, size = score, label = area)) +
    geom_text(show.legend = FALSE) +
    scale_x_reverse() +
    theme_void()

  cowplot::plot_grid(score_plot, tiles, dendro, 
                     ncol = 3,
                     rel_widths = c(0.4, 0.4, 0.2))
}

# - Run program -------------------------------------------------------------
results <- full_scores(neurons)
scores_summary <- summarize_scores(results)

ggsave("~/Desktop/VMH.png", 
       heatmap_plot(neurons, 
             genes = filter(markers, p_val_adj < 0.05 & !duplicated(gene))$gene, 
             scores = scores_summary, 
             scale_clip_low = -2, 
             scale_clip_high = 2), 
       width = 9, height = 5, units = "in", dpi = 300)