make_design_two_groups = function(m, N_cases, N, A, design_options){
  stopifnot(!is.null(design_options$design_type))
  design_type = design_options$design_type
  
  stopifnot(design_type %in% c("Gaussian"))

  if(design_type == "Gaussian"){
    Sigma = make_cov(design_options)
    identity_cov = FALSE
    if(design_options$cov_type == "AR"){
      if(design_options$rho == 0){
        identity_cov = TRUE
      }
    }
    if(identity_cov){
      X = matrix(rnorm(N*m), N, m)
    } else{
      X = mvrnorm(N, numeric(m), Sigma)      
    }
  }
  return(X)
}

make_design_regression = function(m, N, design_options){
  
  stopifnot(!is.null(design_options$design_type))
  design_type = design_options$design_type
  
  stopifnot(design_type %in% c("orthogonal", "Gaussian"))

  if(design_type == "orthogonal"){
    X = matrix(0, N, m)
    X[1:m, 1:m] = diag(m)
  }
  if(design_type == "Gaussian"){
    Sigma = make_cov(design_options)
    X = withSeed(mvrnorm(N, numeric(m), Sigma), seed = 12345)
    X = scale(X)/sqrt(N-1)
  }
  return(X)
}

make_beta = function(m, G, beta_options){
  stopifnot(!is.null(beta_options$beta_type))
  beta_type = beta_options$beta_type
  
  stopifnot(beta_type %in% c("random_locations", "line", "first_k", 
                             "forest_of_lines", "custom", "node_based")) 
  
  beta = matrix(0, m, 1)
  
  if(beta_type == "random_locations"){
    stopifnot(!is.null(beta_options$k))

    k = beta_options$k

    beta = make_beta_random_locations(m, k)
  }
  
  if(beta_type == "node_based"){
    stopifnot(!is.null(beta_options$non_null_nodes))
    non_null_nodes = beta_options$non_null_nodes
    beta = matrix(0, m, 1)
    for(non_null_node in non_null_nodes){
      beta[G$groups[[non_null_node]]] = 1
    }
  }
  
  if(beta_type == "line"){
    D = beta_options$D
    root_node_size = beta_options$root_node_size
    leaf_node_size = beta_options$leaf_node_size
    non_nulls_per_depth = beta_options$non_nulls_per_depth
    beta = make_beta_line(D, root_node_size, leaf_node_size, non_nulls_per_depth)
  }
  
  if(beta_type == "first_k"){
    k = beta_options$k
    beta[1:k] = 1
  }
  
  if(beta_type == "custom"){
    non_nulls = beta_options$non_nulls
    beta[non_nulls] = 1
  }
  
  
  if(beta_type == "forest_of_lines"){
    D = beta_options$D
    W = beta_options$W
    root_node_size = beta_options$root_node_size
    leaf_node_size = beta_options$leaf_node_size
    first_non_null = beta_options$first_non_null
    num_non_nulls = beta_options$num_non_nulls
    
    beta = make_beta_forest_of_lines(D, W, root_node_size, leaf_node_size, first_non_null, num_non_nulls)
  }
  
  return(beta)
}

make_DAG = function(groups, children){
  n = length(groups)
  m = length(unique(unlist(groups)))

  group_sizes = sapply(groups, length)
  ### get children, if necessary ###
  if(is.null(children)){
    # reorder groups by size for convenience

    ord = order(group_sizes)
    groups_ordered = list()
    for(k in 1:n){
      groups_ordered[[k]] = groups[[ord[k]]]
    }
    groups = groups_ordered
    group_sizes = group_sizes[ord]
    
    
    children = vector("list", n)
    for(k in n:2){
      # j is candidate child
      for(j in (k-1):1){
        if(group_sizes[j] < group_sizes[k]){ # don't both to check if they are the same size; assumes no identical groups
          # check if group j is a subset of group k
          if(length(union(groups[[j]], groups[[k]])) == length(groups[[k]])){
            # if so, check if group j is a subset of any other child of group k
            new_child = TRUE
            for(existing_child in children[[k]]){
              if(length(union(groups[[j]], groups[[existing_child]])) == length(groups[[existing_child]])){
                new_child = FALSE
                break
              }
            }
            if(new_child){
              children[[k]] = union(children[[k]], j) 
            }
          }
        }
      }
    }
  }
  
  # get parents
  parents = vector("list", n)
  for(i in 1:n){
    for(j in children[[i]]){
      parents[[j]] = c(parents[[j]], i)
    }
  }
  
  # get which groups contain which genes
  groups_containing = vector("list", m)
  for(k in 1:n){
    for(j in groups[[k]]){
      groups_containing[[j]] = c(groups_containing[[j]], k)
    }
  }

  depths = numeric(n)
  for(i in order(group_sizes, decreasing = TRUE)){
    if(length(parents[[i]]) == 0){
      depths[i] = 1
    } else{
      depths[i] = 1 + max(depths[parents[[i]]])
    }
  }

  
  G = c()
  G$n = n
  G$m = m
  G$groups = groups
  G$C = children
  G$Pa = parents
  G$groups_containing = groups_containing
  G$depths = depths
  return(G)
}

make_base_weights = function(G, base_weights_type){
  stopifnot(base_weights_type %in% c("constant", "IC", "soft"))
  if(base_weights_type == "constant"){
    base_weights = rep(1, G$n)
  }
  if(base_weights_type == "IC"){
      group_sizes = sapply(G$groups, length)
      base_weights = -log(group_sizes/G$m)
      base_weights = base_weights/max(base_weights)
  }
  if(base_weights_type == "soft"){
    base_weights = numeric(G$n)
    group_sizes = sapply(G$groups, length)
    for(i in order(group_sizes, decreasing = FALSE)){
      if(length(G$C[[i]]) == 0){
        base_weights[i] = 1
      } else{
        base_weights[i] = min(sapply(G$C[[i]], function(child)(base_weights[child]/length(G$Pa[[child]]))))
      }
    }    
  }
  return(base_weights)  
}

make_p_values = function(X, y, G, p_vals_options){
  stopifnot(!is.null(p_vals_options$data_type))
  stopifnot(!is.null(p_vals_options$p_vals_type))
  
  data_type = p_vals_options$data_type
  p_vals_type = p_vals_options$p_vals_type
  
  switch(data_type,
         "regression" = {
           stopifnot(p_vals_type %in% c("Goeman_global", "hypergeometric")) 
           if(p_vals_type == "Goeman_global"){
             library(globaltest)
             colnames(X) = 1:ncol(X)
             gt_output=gt(y,X,subsets=G$groups) 
             P=p.value(gt_output)
           }
           if(p_vals_type == "hypergeometric"){
             marginal_associations = abs(cor(X, y))
             
             thresholding_type = p_vals_options$thresholding_type
             switch(thresholding_type,
                    "best_K" = {
                      K = p_vals_options$K
                      ord = order(marginal_associations, decreasing = TRUE)
                      top_genes = logical(G$m)
                      top_genes[ord[1:K]] = TRUE
                    },
                    
                    "p_value_threshold" = {
                      print("Implement later...")
                    }
             )
             process_term = function(j){
               conting_table = matrix(0, 2, 2)
               conting_table[1,1] = sum(top_genes & G$groups_containing[j,])
               conting_table[1,2] = sum(!top_genes & G$groups_containing[j,])
               conting_table[2,1] = sum(top_genes & !G$groups_containing[j,])
               conting_table[2,2] = sum(!top_genes & !G$groups_containing[j,])
               output = fisher.test(conting_table, alternative = "greater")
               return(output$p.value)
             }
             P = sapply(1:n, process_term)
           }
         },
         
         "two_groups" = {
           stopifnot(p_vals_type %in% c("Simes", "Fisher")) # could also implement hypergeometric, Fisher
           if(p_vals_type == "Simes"){
             P_gene = pnorm(t(X)%*%y, lower.tail = FALSE)
             P = sapply(G$groups, function(group)(min(sort(P_gene[group])*length(group)/(1:length(group)))))
           }
           if(p_vals_type == "Fisher"){
             P_gene = pnorm(t(X)%*%y, lower.tail = FALSE)
             P = sapply(G$groups, function(group)(pchisq(-2*sum(log(P_gene[group])), df = 2*length(group), lower.tail = FALSE)))
           }
         }
  )
  return(P)
}
