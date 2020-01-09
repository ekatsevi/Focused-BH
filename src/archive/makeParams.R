if(params_version == 1){
  # parameter being varied
  parameter_name = "SNR"
  data_type = "two_groups"
  N_cases = 100
  N_controls = 100
  N = N_cases + N_controls
  DAG_version = 1
  q = 0.1
  p_vals_options = list()
  p_vals_options$p_vals_type = "Simes"
  p_vals_options$data_type = data_type
  params = list()
  for(param_index in 1:5){
    p = c()
    p$SNR = seq(0, 0.5, length.out = 5)[param_index] # signal-to-noise ratio
    
    # design options
    p$design_options = list()
    p$design_options$design_type = "Gaussian" 
    p$design_options$cov_type = "AR"
    p$design_options$rho = 0
    
    # beta options
    p$beta_options = list()
    p$beta_options$beta_type = "custom"
#    p$beta_options$non_nulls = c(17, 109, 145, 146, 203, 268, 366, 452, 470, 596, 599, 604, 657, 685) # 118
    p$beta_options$non_nulls = c(17, 109, 145, 146) # 118
    
    params[[param_index]] = p
  }
}

if(params_version == 2){
  # parameter being varied
  parameter_name = "SNR"
  data_type = "two_groups"
  N_cases = 100
  N_controls = 100
  N = N_cases + N_controls
  DAG_version = 2
  q = 0.1
  p_vals_options = list()
  p_vals_options$p_vals_type = "Fisher"
  p_vals_options$data_type = data_type
  params = list()
  for(param_index in 1:11){
    p = c()
    p$SNR = seq(0, 5, length.out = 11)[param_index] # signal-to-noise ratio
    
    # design options
    p$design_options = list()
    p$design_options$design_type = "Gaussian" 
    p$design_options$cov_type = "AR"
    p$design_options$rho = 0
    
    # beta options
    p$beta_options = list()
    p$beta_options$beta_type = "random_locations"
    p$beta_options$k = 4
    
    params[[param_index]] = p
  }
}

if(params_version == 3){
  # parameter being varied
  parameter_name = "SNR"
  data_type = "two_groups"
  N_cases = 100
  N_controls = 100
  N = N_cases + N_controls
  DAG_version = 1
  q = 0.2
  p_vals_options = list()
  p_vals_options$p_vals_type = "Fisher"
  p_vals_options$data_type = data_type
  params = list()
  for(param_index in 1:5){
    p = c()
    p$SNR = seq(0, 0.5, length.out = 5)[param_index] # signal-to-noise ratio
    
    # design options
    p$design_options = list()
    p$design_options$design_type = "Gaussian" 
    p$design_options$cov_type = "AR"
    p$design_options$rho = 0
    
    # beta options
    p$beta_options = list()
    p$beta_options$beta_type = "custom"
    p$beta_options$non_nulls = c(17, 109, 145, 146, 203, 268, 366, 452, 470, 596, 599, 604, 657, 685) # 118
#    p$beta_options$non_nulls = c(17, 109, 145, 146) # 118
    
    params[[param_index]] = p
  }
}

if(params_version == 4){ # DAG, Simes
  # parameter being varied
  parameter_name = "A"
  data_type = "two_groups"
  N_cases = 100
  N_controls = 100
  N = N_cases + N_controls
  DAG_version = 1
  q = 0.1
  p_vals_options = list()
  p_vals_options$p_vals_type = "Simes"
  p_vals_options$data_type = data_type
  params = list()
  for(param_index in 1:11){
    p = c()
    p$A = seq(0, 1, length.out = 11)[param_index] # signal-to-noise ratio
    
    # design options
    p$design_options = list()
    p$design_options$design_type = "Gaussian" 
    p$design_options$cov_type = "AR"
    p$design_options$rho = 0
    
    # beta options
    p$beta_options = list()
    p$beta_options$beta_type = "node_based"
#    p$beta_options$non_nulls = c(17, 109, 145, 146, 203, 268, 366, 452, 470, 596, 599, 604, 657, 685) # 118
#    p$beta_options$non_nulls = c(17, 109, 145, 146, 203, 268, 366, 452, 470, 596, 599, 604, 657, 685,
                               #  40, 89, 213, 244) # 116 & 100
    p$beta_options$non_null_nodes = c(71, 138, 167)
    #    p$beta_options$non_nulls = c(17, 109, 145, 146) # 118
    
    params[[param_index]] = p
  }
}


if(params_version == 5){ # DAG, Fisher
  # parameter being varied
  parameter_name = "A"
  data_type = "two_groups"
  N_cases = 100
  N_controls = 100
  N = N_cases + N_controls
  DAG_version = 1
  q = 0.1
  p_vals_options = list()
  p_vals_options$p_vals_type = "Fisher"
  p_vals_options$data_type = data_type
  params = list()
  for(param_index in 1:11){
    p = c()
    p$A = seq(0, 1, length.out = 11)[param_index] # signal-to-noise ratio
    
    # design options
    p$design_options = list()
    p$design_options$design_type = "Gaussian" 
    p$design_options$cov_type = "AR"
    p$design_options$rho = 0
    
    # beta options
    p$beta_options = list()
    p$beta_options$beta_type = "node_based"
    #    p$beta_options$non_nulls = c(17, 109, 145, 146, 203, 268, 366, 452, 470, 596, 599, 604, 657, 685) # 118
    #    p$beta_options$non_nulls = c(17, 109, 145, 146, 203, 268, 366, 452, 470, 596, 599, 604, 657, 685,
    #  40, 89, 213, 244) # 116 & 100
    p$beta_options$non_null_nodes = c(71, 138, 167)
    params[[param_index]] = p
  }
}

if(params_version == 6){ # for tree, signals scattered, Simes
  # parameter being varied
  parameter_name = "A"
  data_type = "two_groups"
  N_cases = 100
  N_controls = 100
  N = N_cases + N_controls
  DAG_version = 2
  q = 0.1
  p_vals_options = list()
  p_vals_options$p_vals_type = "Simes"
  p_vals_options$data_type = data_type
  params = list()
  for(param_index in 1:11){
    p = c()
    p$A = seq(0, 1, length.out = 11)[param_index] # signal-to-noise ratio
    
    # design options
    p$design_options = list()
    p$design_options$design_type = "Gaussian" 
    p$design_options$cov_type = "AR"
    p$design_options$rho = 0
    
    # beta options
    p$beta_options = list()
    p$beta_options$beta_type = "random_locations"
    p$beta_options$k = 4
    
    params[[param_index]] = p
  }
}

if(params_version == 7){ # for tree, signals clustered, Simes
  # parameter being varied
  parameter_name = "A"
  data_type = "two_groups"
  N_cases = 100
  N_controls = 100
  N = N_cases + N_controls
  DAG_version = 2
  q = 0.1
  p_vals_options = list()
  p_vals_options$p_vals_type = "Simes"
  p_vals_options$data_type = data_type
  params = list()
  for(param_index in 1:5){ # change from 5 to 11
    p = c()
    p$A = seq(0, 10, length.out = 5)[param_index] # signal-to-noise ratio
    
    # design options
    p$design_options = list()
    p$design_options$design_type = "Gaussian" 
    p$design_options$cov_type = "AR"
    p$design_options$rho = 0
    
    # beta options
    p$beta_options = list()
    p$beta_options$beta_type = "first_k"
    p$beta_options$k = 4
    
    params[[param_index]] = p
  }
}

if(params_version == 8){ # for tree, signals clustered, Fisher
  # parameter being varied
  parameter_name = "A"
  data_type = "two_groups"
  N_cases = 100
  N_controls = 100
  N = N_cases + N_controls
  DAG_version = 2
  q = 0.1
  p_vals_options = list()
  p_vals_options$p_vals_type = "Fisher"
  p_vals_options$data_type = data_type
  params = list()
  for(param_index in 1:5){ # change from 5 to 11
    p = c()
    p$A = seq(0, 10, length.out = 5)[param_index] # signal-to-noise ratio
    
    # design options
    p$design_options = list()
    p$design_options$design_type = "Gaussian" 
    p$design_options$cov_type = "AR"
    p$design_options$rho = 0
    
    # beta options
    p$beta_options = list()
    p$beta_options$beta_type = "first_k"
    p$beta_options$k = 4
    
    params[[param_index]] = p
  }
}

if(params_version == 9){ # for tree, signals scattered, Fisher
  # parameter being varied
  parameter_name = "A"
  data_type = "two_groups"
  N_cases = 100
  N_controls = 100
  N = N_cases + N_controls
  DAG_version = 2
  q = 0.1
  p_vals_options = list()
  p_vals_options$p_vals_type = "Fisher"
  p_vals_options$data_type = data_type
  params = list()
  for(param_index in 1:11){
    p = c()
    p$A = seq(0, 1, length.out = 11)[param_index] # signal-to-noise ratio
    
    # design options
    p$design_options = list()
    p$design_options$design_type = "Gaussian" 
    p$design_options$cov_type = "AR"
    p$design_options$rho = 0
    
    # beta options
    p$beta_options = list()
    p$beta_options$beta_type = "random_locations"
    p$beta_options$k = 4
    
    params[[param_index]] = p
  }
}

if(params_version == 10){ # for tree, signals scattered, Simes
  # parameter being varied
  parameter_name = "A"
  data_type = "two_groups"
  N_cases = 100
  N_controls = 100
  N = N_cases + N_controls
  DAG_version = 2
  q = 0.1
  p_vals_options = list()
  p_vals_options$p_vals_type = "Simes"
  p_vals_options$data_type = data_type
  params = list()
  for(param_index in 1:11){
    p = c()
    p$A = seq(0, 10, length.out = 11)[param_index] # signal-to-noise ratio
    
    # design options
    p$design_options = list()
    p$design_options$design_type = "Gaussian" 
    p$design_options$cov_type = "AR"
    p$design_options$rho = 0
    
    # beta options
    p$beta_options = list()
    p$beta_options$beta_type = "random_locations"
    p$beta_options$k = 4
    
    params[[param_index]] = p
  }
}

if(params_version == 11){ # for tree, signals scattered, Simes
  # parameter being varied
  parameter_name = "A"
  data_type = "two_groups"
  N_cases = 100
  N_controls = 100
  N = N_cases + N_controls
  DAG_version = 5
  q = 0.1
  p_vals_options = list()
  p_vals_options$p_vals_type = "Simes"
  p_vals_options$data_type = data_type
  params = list()
  for(param_index in 1:11){
    p = c()
    p$A = seq(0, 1, length.out = 11)[param_index] # signal-to-noise ratio
    
    # design options
    p$design_options = list()
    p$design_options$design_type = "Gaussian" 
    p$design_options$cov_type = "AR"
    p$design_options$rho = 0
    
    # beta options
    p$beta_options = list()
    p$beta_options$beta_type = "custom"
#    p$beta_options$non_nulls = c(1,2,3,4,13,14,17,18,21,22,34,42)
#    p$beta_options$non_nulls = c()1 + 16*(0:9)
    p$beta_options$non_nulls = c(1,2,3,4,17, 18, 21,22, 1 + 16*(2:14))
    
    params[[param_index]] = p
  }
}

if(params_version == 12){ # for tree, signals scattered, Simes
  # parameter being varied
  parameter_name = "A"
  data_type = "two_groups"
  N_cases = 100
  N_controls = 100
  N = N_cases + N_controls
  DAG_version = 5
  q = 0.1
  p_vals_options = list()
  p_vals_options$p_vals_type = "Simes"
  p_vals_options$data_type = data_type
  params = list()
  for(param_index in 1:11){
    p = c()
    p$A = seq(0, 1, length.out = 11)[param_index] # signal-to-noise ratio
    
    # design options
    p$design_options = list()
    p$design_options$design_type = "Gaussian" 
    p$design_options$cov_type = "AR"
    p$design_options$rho = 0
    
    # beta options
    p$beta_options = list()
    p$beta_options$beta_type = "custom"
    #    p$beta_options$non_nulls = c(1,2,3,4,13,14,17,18,21,22,34,42)
    #    p$beta_options$non_nulls = c()1 + 16*(0:9)
    p$beta_options$non_nulls = c(1,2,3,4,17, 18, 21,22, 1 + 32*(2:14))
    
    params[[param_index]] = p
  }
}


construct_SH = function (G) 
{
    graph = construct(G$groups)
    
    if(length(graph@parents) != G$n){
        graph = new("DAGstructure", parents = G$Pa, children = G$C, sets = G$groups, twoway = FALSE)
    }

    return(graph)
}

source("makeDAGparams.R", local = TRUE)
G = make_DAG(groups, children)
G_Structured_Holm = construct_SH(G)

roots = which(sapply(G$Pa, length) == 0)
G_Y = G
G_Y$n = G_Y$n + 1
G_Y$C[[G_Y$n]] = roots
G_Y$Pa[[G_Y$n]] = integer(0)
for(root in roots){
  G_Y$Pa[[root]] = G_Y$n
}
G_Y$depths = c(G_Y$depths, 0)

G_Yekutieli = t(do.call(cbind, 
                        lapply(1:G_Y$n, 
                               function(i)(if(length(G_Y$C[[i]]) > 0) sapply(G_Y$C[[i]], 
                                                                                   function(j)(c(as.character(i), as.character(j)))) else matrix(0,2,0)))))
