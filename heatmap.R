# heatmap function
heatmap <- function(object, 
                    cells = NULL, 
                    genes = NULL, 
                    scale = TRUE,
                    col_names = FALSE, 
                    row_names = FALSE, 
                    heatmap_legend = TRUE) {
  
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
    color_pal <- RColorBrewer::brewer.pal("div", "RdYlBu")
  } else {
    matrix <- object@data[genes, ]
    color_pal <- RColorBrewer::brewer.pal("seq", "YlGnBu")
  }
  
  mean_expression <- function(matrix, cells) {
    rowMeans(matrix[, cells])
  }
  
  cells_use <- function(object, ident) {
    names(object@ident)[object@ident == ident]
  }
  
  map(unique(object@ident),
      ~mean_expression(matrix, cells = cells_use(object, .x))) %>%
    as.data.frame() %>%
    set_names(unique(object@ident)) %>%
    t() %>%
    pheatmap::pheatmap(., 
                       show_colnames = col_names,
                       treeheight_col = 0,
                       legend = heatmap_legend,
                       color = color_pal)
}
