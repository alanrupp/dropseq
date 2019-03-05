umap_plot <- function(object, cells = NULL, clusters = NULL, 
                      features = NULL, legend = FALSE) {
  
  umap <- data.frame(
    UMAP1 = object@dr$umap@cell.embeddings[, 1],
    UMAP2 = object@dr$umap@cell.embeddings[, 2]
  )
  
  pull_feature <- function(feature) {
    if (feature %in% colnames(object@meta.data)) {
      return(object@meta.data[, feature])
    } else if (feature %in% rownames(object@data)) {
      return(object@data[feature, ])
    } else if (feature == "ident") {
      return(object@ident)
    }
  }
  
  results <- map(features, ~pull_feature(.x)) %>% 
    bind_cols() %>%
    set_names(features) %>%
    bind_cols(umap)
  
  results %>%
    gather(-starts_with("UMAP"), key = "feature", value = "value") %>%
    mutate(value = as.numeric(value)) %>%
    ggplot(aes(x = UMAP1, y = UMAP2, color = value)) +
    geom_point(show.legend = legend, stroke = 0) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    facet_wrap(~feature)
}