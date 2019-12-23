AR_cov = function(p, rho){
  Sigma = outer(1:p, 1:p, function(i, j)(rho^(abs(i-j))))
  return(Sigma)
}

conditional_AR = function(Sigma, group){
  Sigma_inv = solve(Sigma[-group, -group, drop = FALSE])
  # Sigma_inv[abs(Sigma_inv) < 1e-10] = 0
  mean_vec = Sigma[group, -group, drop = FALSE]%*%Sigma_inv
  cov_mat = Sigma[group, group,drop = FALSE] - 
    Sigma[group, -group, drop = FALSE]%*%Sigma_inv%*%Sigma[-group, group, drop = FALSE]
  out = c()
  out$mean_vec = mean_vec
  out$cov_mat = cov_mat
  return(out)
}

# Sigma_inv is precision matix
# k is variable to condition on the rest
conditional_gaussian_fast = function(Sigma_inv, k){
  out = c()
  out$mean_vec = -1/Sigma_inv[k,k]*Sigma_inv[k, -k,drop=FALSE]
  out$cov_mat = 1/Sigma_inv[k,k]
  return(out)
}

soft_threshold = function(t, lambda){
  return(sign(t)*max(0, abs(t)-lambda))
}

select_lambda = function(err_estimates, num_increasing){
  num_lambda = length(err_estimates)
  err_estimates_next_few = numeric(num_lambda)
  err_estimates_next_few[1:(num_lambda-1)] = sapply(1:(num_lambda-1), 
                                                    function(idx)(min(err_estimates[(idx+1):min(idx+num_increasing, num_lambda)])))
  err_estimates_next_few[num_lambda] = Inf
  return(min(which(err_estimates <= err_estimates_next_few)))
}
