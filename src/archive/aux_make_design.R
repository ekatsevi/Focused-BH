library(MASS)

make_cov = function(design_options){
  stopifnot(!is.null(design_options$cov_type))
  cov_type = design_options$cov_type
  stopifnot(cov_type %in% c("AR", "block_AR"))
  if(cov_type == "AR"){
    stopifnot(!is.null(design_options$rho))
    rho = design_options$rho
    Sigma = make_AR_cov(rho, m)
  }
  if(cov_type == "block_AR"){
    block_sizes = design_options$block_sizes
    rho_values = design_options$rho_values
    Sigma = make_block_AR_cov(block_sizes, rho_values)
  }
  return(Sigma)
}

make_AR_cov = function(rho, n){
  Sigma = rho^(abs(outer(1:n, 1:n, "-")))
  return(Sigma)
}

make_block_AR_cov = function(block_sizes, rho_values){
  n = sum(block_sizes)  
  num_blocks = length(block_sizes)
  Sigma = matrix(0, n, n)
  idx = 1
  for(block in 1:num_blocks){
    block_size = block_sizes[block]
    rho = rho_values[block]
    block_idx = idx:(idx +block_size-1)
    Sigma[block_idx, block_idx] = make_AR_cov(rho, block_size)
    idx = idx + block_size
  }
  return(Sigma)
}

