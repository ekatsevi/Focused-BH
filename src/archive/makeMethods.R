if(methods_version == 1){ # change for new method sets
  method_type = numeric(0) # method type; BH, filtered_BH, Structured_Holm, filtered_BH_oracle, Yekutieli
  method_name = numeric(0) # method type; BH, filtered_BH, Structured_Holm, filtered_BH_oracle, Yekutieli
  method_options = numeric(0) # any options for these methods
  color = numeric(0)
  linetype = numeric(0)
  methods = data.frame(method_type, method_name, method_options, color, linetype)
  
  methods[1,] = c("Storey_BH", "Storey-BH", NA, "firebrick1", "solid")
  methods[2,] = c("Structured_Holm", "Structured Holm", NA, "blue", "solid")
  methods[3,] = c("filtered_BH", "Focused BH", "conservative", "dodgerblue", "solid")
  methods[4,] = c("filtered_BH", "Focused BH (permutation)", "permutation", "dodgerblue", "longdash")
  methods[5,] = c("filtered_BH", "Focused BH (oracle)", "oracle", "cyan", "solid")
  methods[6,] = c("Yekutieli", "Yekutieli", TRUE, "purple", "solid")
}

num_methods = nrow(methods)