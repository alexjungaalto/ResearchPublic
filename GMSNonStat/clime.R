library("fastclime")
library("matrixcalc")

# adjacency matrix of chain graph
chain_adj <- function(dimens) {
  first_row <- rep(0, dimens)
  first_row[2] = 1
  return (toeplitz(first_row))
}

# neighborhood as a binary vector
find_neighbor <- function(node_idx, adj_matrix) {
  return (adj_matrix[node_idx,])
}

generate_precison_matrix <- function(rho_min, adj_matrix, dimens) {
  C_0 <- (1/3)*rho_min*(matrix(1, dimens,dimens)-diag(dimens)) + diag(dimens)
  return (C_0 * (adj_matrix + diag(dimens))) # add diagonal for the covariance matrix
}

neighbor_by_clime <-function(query_node, blk_len, dimens, cov_matrix) {
  samples <- mvrnorm(n = blk_len, rep(0, dimens), cov_matrix)
  #apply CLIME
  out1 = fastclime(samples, 0.1)
  out2 = fastclime.selector(out1$lambdamtx, out1$icovlist, 0.2)
  return (out2$adaj[query_node,])
}

write_error_to_csv <- function(x_vals, y_vals, filename, dimens) {
  filename = paste(sprintf(filename, dimens), 'csv', sep='.')
  df = data.frame(x_vals, y_vals)
  write.table(df,file=filename, sep = "," , row.names = FALSE, col.names = c("x","y"))
}

exp_clime <-function() {
  # some parameter settings
  blklen_vals <- seq(10,100, by=10)
  dim_vals <- c(64, 128, 256, 512)
  rho_min_fac_vals <- c(1, 0.9, 0.8, 1.2)
  error_rate <- matrix(0, length(blklen_vals), length(dim_vals)) 
  
  num_blks <- 4
  num_tests <- 100
  sparsity <- 2
  
  # Repeat the estimation for each value of dimension and blk length
  
  #for each dimension value
  for (idx_dim in 1:length(dim_vals)) {
    dimens = dim_vals[idx_dim]
    # initialize empty precision and covariance matrices
    C = array(0, dim=c(dimens,dimens, num_blks))
    
    # the underlying matrices of the whole process without block structure
    adj_matrix <- chain_adj(dimens)
    K_0 <- generate_precison_matrix(rho_min_fac_vals[idx_dim], adj_matrix, dimens) 
    
    # The graph varies from blk to blk by remove an edge from the chain
    for (idx_blk in 1:num_blks) {
      K <- K_0
      K[idx_blk,idx_blk+1] = 0
      K[idx_blk+1,idx_blk] = 0
      C[,,idx_blk]=chol2inv(chol(K))
      
    }
    # for each blk length value
    query_node = 2
    true_neighbors = find_neighbor(query_node, adj_matrix)
    for (idx_blklen in 1:length(blklen_vals)) {
      blk_len = blklen_vals[idx_blklen]
      num_samples = num_blks * blk_len
      acc = 0
      #repeat the tes num_tests times
      for (idx in 1:num_tests) {
        
        estimated_neighbors = rep(0, dimens)
        #run clime for each blk and take the union of the neighbor
        for (idx_blk in 1:num_blks) {
          estimated_neighbors = estimated_neighbors + neighbor_by_clime(query_node, blk_len, dimens, C[,,idx_blk])
        }      
        #neighbors?
        estimated_neighbors = as.numeric(estimated_neighbors>0)
        if (norm(true_neighbors-estimated_neighbors, type="2")<0.1) {
          acc = acc + 1
        }
      }
      #the error rate for each pair of blk length and dimension
      error_rate[idx_blklen, idx_dim] = 1 - acc/num_tests
    }
  }
  
  # write the error rate into csv file
  for (idx_dim in 1:length(dim_vals)) {
    x_vals = num_blks * blklen_vals
    y_vals = error_rate[,idx_dim]
    write_error_to_csv(x_vals, y_vals, 'CLIMERawSample%02d', dim_vals[idx_dim])
    
    scale_x_vals = blklen_vals*num_blks*rho_min_fac_vals[idx_dim]/log(dim_vals[idx_dim]) 
    write_error_to_csv(scale_x_vals, y_vals, 'CLIMEScaledSample%02d', dim_vals[idx_dim])
  }
}


