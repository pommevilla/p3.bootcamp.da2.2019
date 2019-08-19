write.cor_to_edgelist <- function(cor_matrix, fout){
  cor_matrix %>% 
    melt() %>% 
    filter(value != 0) %>% 
    write.csv(fout, quote = FALSE, row.names = FALSE)
}

get_expr_data <- function(prot, treats, genes){
  e.data <- prot %>%
    column_to_rownames(var = genes) %>%
    dplyr::select(treats - 1) %>%
    t()

  return(e.data)
}

plot_cor_heatmap <- function(cor_matrix, triangle = TRUE){

  sig_count <- cor_matrix %>%
    melt() %>%
    filter(value != 0) %>%
    count()

  if (triangle == TRUE) {
    cor_matrix[upper.tri(cor_matrix)] <- NA
  }


  p <- ggplot(na.omit(melt(cor_matrix)), aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    coord_equal() +
    labs(x = "", y = "", fill = "Correlation") +
    scale_fill_gradient2(low = "red", mid = "black", high = "steelblue", space = "Lab",
                         breaks = c(-1, 0, 1), labels = c(-1, 0, 1),
                         limits = c(-1, 1)) +
    labs(subtitle = paste0("# Interactions: ", sig_count)) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(), panel.background = element_blank())
  return(p)
}

bootstrap_cor <- function(prot, treats, genes = "UNIQID", B = 1000, method = "pearson") {

  library(data.table)

  dat.observed <- prot %>%
    column_to_rownames(var = genes) %>%
    dplyr::select(treats - 1) %>%
    t()

  n <- ncol(dat.observed)
  boots = data.table(c_score = numeric(), counts = numeric())

  for (i in 1:B){
    dat.perm <- data.table(dat.observed[, sample(n, n, replace = TRUE)])

    cor.perm <- data.table(c(round(cor(dat.perm, method = method), 3)))
    cor.perm[, counts := .N, by = .(V1)]
    cor.perm <- unique(cor.perm)

    boots <- rbindlist(list(boots, cor.perm))[, lapply(.SD, sum, na.rm = TRUE),
                                              by = .(c_score)]

  }

  boots <- boots %>% arrange(-c_score) %>% data.table()

  return(boots)
}

bootstrap_quantiles <- function(boots, p) {
  boots[, prop := counts/sum(counts)]

  quantiles <- boots[, list(
    upper = c_score[sum(cumsum(prop) <= (p / 2))],
    lower = c_score[sum(cumsum(prop) <= (1 - (p / 2)))]
  )]

  return(quantiles)
}

threshold_cor_matrix <- function(cor.matrix, quantiles){

  cor.matrix[cor.matrix >= quantiles$upper] <- 1
  cor.matrix[cor.matrix <= quantiles$lower] <- -1
  cor.matrix[cor.matrix < quantiles$upper & cor.matrix > quantiles$lower] <- 0

  # Additionally set the diagonal to 0.
  diag(cor.matrix) <- 0

  return(cor.matrix)

}

prepare_network <- function(network) {
  network <- delete_vertices(network,
                             V(network)[degree(network) == 0])
  V(network)$size = degree(network)
  V(network)$btwn <- betweenness(network, weights = rep(1, length = length(E(network))))

  return(network)
}

get_treatment_cor <- function(prot, treats, genes = "UNIQID", method = "pearson"){
  e.data <- get_expr_data(prot, treats, genes)

  cor.data <- cor(e.data, method = method)

  return(cor.data)
}

plot_bootstrap <- function(boots, quantiles) {
  ggplot(boots, aes(x = c_score, y = counts)) +
    geom_col() +
    scale_y_continuous(limits = c(0, max(boots$counts))) +
    labs(x = "Pearson correlation", y = "Counts",
         title = "Simulated correlation values (1000 iterations)") +
    geom_vline(xintercept = c(quantiles$lower, quantiles$upper), color = "#BB0000", linetype = 'dashed') +
    geom_text(aes(x = quantiles$lower , label = paste(quantiles$lower), y = max(boots$counts) * 0.75),
              colour = "#BB0000", size = 3, nudge_x = 0.05) +
    geom_text(aes(x = quantiles$upper , label = paste(quantiles$upper), y = max(boots$counts) * 0.75),
              colour = "#BB0000", size = 3, nudge_x = -0.05)
}

scale_and_merge_prot_transc <- function(p, t, p.treats = 8:16, t.treats = 6:14) {
  prot[p.treats] <- scale(prot[p.treats])
  trans[t.treats] <- scale(trans[t.treats])
  pt_merged <- merge(prot, trans, by.x = "UNIQID", by.y = "Unique_ID") 
  genotype_info <- pt_merged %>% 
    dplyr::select(p.treats, t.treats + (dim(proteomics)[2] - 1)) %>% 
    t()
  colnames(genotype_info) <-pt_merged$UNIQID
  return(genotype_info)
}

pca_nth_component_histogram <- function(pca, n, thresh = c(0.05, 0.95)){
  loadings <- pca$loadings$X
  component_of_interest <- loadings[, n]
  qs <- quantile(component_of_interest, thresh)
  data.frame(val = component_of_interest) %>%
    ggplot(aes(val)) +
    geom_histogram() +
    geom_vline(xintercept = qs, color = "#BB0000", linetype = 'dashed') +
    labs(x = "Coefficient", y = "Counts",
         title = paste0("Distribution of coefficients for principal component: ", as.character(substitute(n))))
}

get_top_loadings <- function(pca, component, thresh = c(0.05, 0.95)) {
  component_of_interest <- pca$loadings$X[, component]
  qs <- quantile(component_of_interest, thresh)
  top_loadings <- component_of_interest[component_of_interest > qs[2] | 
                                          component_of_interest < qs[1]] %>% 
    sort(decreasing = TRUE) 
  
  return(top_loadings)
}

plot_dd <- function(g){
  g.dd <- data.frame(0:(length(degree_distribution(g)) - 1), 
                     degree_distribution(g, cumulative = FALSE))
  colnames(g.dd) <-c("k", "p_k")
  ggplot(g.dd, aes(x = k, y = p_k)) + geom_bar(stat = "identity") +
    labs(x = "Degree (k)", y = "Proportion of nodes of degree k (p(k))", title = "Degree distribution: k vs. p(k)", 
         subtitle = paste("Graph:", deparse(substitute(g)))) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
}