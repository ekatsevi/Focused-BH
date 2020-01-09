# set seed
set.seed(100000*params_version + 1000*param_idx + rep)

switch(data_type,
       "regression" = {
         # generate response
         y = X%*%(A*beta) + rnorm(N)
       },
       
       "two_groups" = {
         # generate X
         X = make_design_two_groups(m, N_cases, N, A, design_options)
         cases_indicator = matrix(0, N, 1)
         cases_indicator[1:N_cases] = 1
         X = X + A*(cases_indicator %*% t(beta))
       }   
)
P = make_p_values(X, y, G, p_vals_options)