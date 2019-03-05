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


