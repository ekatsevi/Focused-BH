##### run the different methods #####

num_methods = nrow(methods)
rejections = matrix(FALSE, num_methods, n)
t_stars = numeric(num_methods)

for(method in 1:num_methods){
  method_type = methods$method_type[method]
  method_options = methods$method_options[method]
  
  switch(method_type,
         BH = {
           R = p.adjust(P, method = "fdr") <= q
         },
         
         Storey_BH = {
           R = qvalue(P, lambda = q, fdr.level = q)$significant
         },
         
         filtered_BH = {

             ## can hopefully delete this later, as t should be packaged in as part of the
             ## precomputations
             t = unique(c(seq(0, 0.01, length.out = 1000),
                          seq(0.01, 0.1, length.out = 1000),
                          seq(0.1, 1, length.out = 1000)))
             
             switch(method_options,
                    "conservative" = {
                        gamma = function(s)(n*s)
                    },
                    
                    "permutation" = {
                        gamma = approxfun(x = t, y = gamma_permutation[,param_idx], 
                                          method = "constant", yleft = 0, rule = 2, f = 1, ties = mean)
                    },
                    
                    "oracle" = {
                        gamma = approxfun(x = t, y = gamma_oracle[,param_idx], 
                                          method = "constant", yleft = 0, rule = 2, f = 1, ties = mean)
                    },
                    
                    "Storey" = {
                        lambda = q
                        n_0 = (1+sum(P > lambda))/(1-lambda)
                        gamma = function(s)(n_0*s)
                    }
                    )
             R = filtered_BH(P, G, q, gamma, filter_type, base_weights)
         },
         
         DAGGER = {
           R = DAGGER(P, G, q)
         },
         
         Yekutieli = {
           adjust = method_options
           if(adjust){
               ## q_Yekutieli = q/(2*1.44*max(G$depths))
               q_Yekutieli = q/(2*max(G$depths))
               
           } else{
             q_Yekutieli = q
           }
           unadj.p.values = c(P,0)
           names(unadj.p.values) = sapply(1:(n+1), function(i)(as.character(i)))
           hyp.tree = hFDR.adjust(unadj.p.values, G_Yekutieli, q_Yekutieli)
           R = hyp.tree@p.vals[[2]] <= q_Yekutieli
           R[is.na(R)] = FALSE
           R = R[1:n]
         },
         
         Structured_Holm = {
           if(min(P) >= q/n){
             R = rep(FALSE, n)
           } else{
             R = structuredHolm(G_Structured_Holm, alpha_max = q, pvalues=P)@rejected
           }
         }
  )
  rejections[method,] = R
    if(any(R)){
    t_stars[method] = max(P[R])
  } else{
    t_stars[method] = 0
  }
}