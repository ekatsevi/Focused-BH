switch(data_type,
  "regression" = {
    # construct design matrix
    X = make_design_regression(m, N, design_options)

    if(is.null(A)){
        ## get signal amplitude
        A = as.numeric(sqrt(N*SNR/(t(X%*%beta)%*%(X%*%beta))))
    }
  },
    
  "two_groups" = {
    # construct response vector
    y = numeric(N)
    y[1:N_cases] = 1/N_cases
    y[(N_cases+1):N] = -1/N_controls
    y = y/sqrt(1/N_cases + 1/N_controls)
    
    ## get signal amplitude
    if(is.null(A)){
        A = sqrt(SNR*m*(1/N_cases + 1/N_controls)/sum(beta*beta))
    }
  }
)

