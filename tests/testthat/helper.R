
make_dummy_neighborhoods_data.frame <- function(nrow = 10, n_cells = 100, cell_names = TRUE){
  cells <- if(cell_names){
    paste0("cell_", seq_len(n_cells))
  }else{
    seq_len(n_cells)
  }

  neighborhoods <- lapply(seq_len(nrow), \(idx){
    if(idx == 1){
      integer(0L)
    }else if(idx == nrow){
      seq_len(n_cells)
    }else{
      size <- sample.int(n_cells, size = 1)
      sample.int(n_cells, size, replace = FALSE)
    }
  })
  if(cell_names){
    names <- paste0("cell_", seq_len(n_cells))
    neighborhoods <- lapply(neighborhoods, \(nei) names[nei])
  }

  data.frame(name = paste0("gene_", seq_len(nrow)),
             neighborhood = I(neighborhoods),
             n_cells = lengths(neighborhoods),
             sel_statistic = rnorm(nrow, mean = 10)^2,
             pval = runif(nrow),
             adj_pval = runif(nrow),
             f_statistic = rnorm(nrow, mean = 10)^2,
             df1 = 2,
             df2 = rnorm(1, mean = 3)^2,
             lfc = rnorm(nrow),
             did_pval = runif(nrow),
             did_adj_pval = runif(nrow),
             did_lfc = rnorm(nrow))
}
