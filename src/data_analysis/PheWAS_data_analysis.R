rm(list = ls())
machine = "local"
source("setup.R")

######### DATA ANALYSIS PARAMETERS ###################
allele = "B_2705"              # HLA allele to test
q = 0.05                       # FDR target level
min_affected = 50              # consider only diseases with 
                               #   at least this many cases
filter_name = "outer_nodes"    # filter to use
methods = c("Focused_BH",      # methods to apply
            "Structured_Holm", 
            "BH", 
            "Leaf_BH", 
            "Yekutieli")

############# PROCESSED FILES ###################
# processed ICD tree file
graph_file = sprintf("%s/data/processed/biobank/ICD_graph.Rda", base_dir)
# processed diseases file
HES_affected_file = sprintf("%s/data/processed/biobank/HES_affected.Rda", base_dir)
# file storing number of cases per disease
num_affected_file = sprintf("%s/data/processed/biobank/HES_number_affected.tsv", base_dir)
# processed HLA file
HLA_split_file = sprintf("%s/data/processed/biobank/HLA_split.tsv", base_dir)
# PheWAS p-values 
pvals_filename = sprintf("%s/data/processed/biobank/pvals_HES_%s_m%d.tsv", base_dir, allele, min_affected)
# reduced graph 
reduced_graph_filename = sprintf("%s/data/processed/biobank/G_at_least_%d_cases.Rda", base_dir, min_affected)
# final sets of rejections for each method
rejections_filename = sprintf("%s/data/processed/biobank/PheWAS_rejections.tsv", base_dir)

######### CREATE P-VALUES ###############
if(!file.exists(pvals_filename)){
  # read in processed data
  HLA_split_df = read_tsv(HLA_split_file, col_types = cols(.default = col_double()))
  df_num_affected = read_tsv(num_affected_file, col_types = "ci")
  load(HES_affected_file)
  load(graph_file)
  
  # restrict attention to the allele of interest
  genotypes = round(HLA_split_df[[allele]])
  
  # Carry out chi-square test for association with HLA allele
  # for each disease with at least 50 cases
  num_affected = df_num_affected$num_affected
  m = length(num_affected)
  P = numeric(m)
  P[num_affected < min_affected] = 1
  for(k in which(num_affected >= min_affected)){
    affecteds = affected[[as.character(G$info$node_ids[k])]]
    aff = c(sum(genotypes[affecteds] == 0, na.rm = TRUE), 
            sum(genotypes[affecteds] == 1, na.rm = TRUE), 
            sum(genotypes[affecteds] == 2, na.rm = TRUE))
    unaff = c(sum(genotypes == 0, na.rm = TRUE), 
              sum(genotypes == 1, na.rm = TRUE), 
              sum(genotypes == 2, na.rm = TRUE)) - aff
    cont_table = data.frame(aff, unaff)
    P[k] = chisq.test(cont_table)$p.value
  }
  
  # save to file
  write_tsv(data.frame(P), path = pvals_filename)
}

####### EXCLUDE DISEASES WITH FEWER THAN 50 CASES FROM TREE ###############
if(!file.exists(reduced_graph_filename)){
  # load full graph and number affected per disease
  load(graph_file)
  df_num_affected = read_tsv(num_affected_file, col_types = "ci")
  
  # remove diseases with too few cases
  num_affected = df_num_affected$num_affected
  to_remove = num_affected < min_affected
  G_reduced = subset_graph(to_remove, G)
  
  # convert to Structured Holm and Yekutieli formats
  G_SH_reduced = new("DAGstructure", parents = G_reduced$Pa, children = G_reduced$C, 
                     sets = G_reduced$sets, twoway = FALSE)
  G_Y_reduced = convert_to_Yekutieli_format(G_reduced)
  
  # save to file
  save(G_reduced, G_SH_reduced, G_Y_reduced, file = reduced_graph_filename)
}

####### APPLY EACH HIERARCHICAL TESTING METHOD ###############
if(!file.exists(rejections_filename)){
  # read in the p-values, reduced graph, and cases per disease
  P = read_tsv(pvals_filename, col_types = "d") %>% pull()
  load(reduced_graph_filename)
  df_num_affected = read_tsv(num_affected_file, col_types = "ci")
  num_affected = df_num_affected$num_affected
  
  # remove p-values from diseases with too few cases
  to_remove = num_affected < min_affected
  P = P[!to_remove]
  
  # set up array to store rejections
  m_reduced = G_reduced$m
  num_methods = length(methods)
  rejections = array(FALSE, 
                     dim = c(2,num_methods, m_reduced), 
                     dimnames = list(c("before_filter", "after_filter"), 
                                     methods, 
                                     G_reduced$node_names))
  
  # apply BH
  rejections["before_filter", "BH",] = p.adjust(P, "fdr") <= q
  
  # apply Leaf BH 
  leaves = sapply(G_reduced$C, length) == 0
  rejections["before_filter", "Leaf_BH",leaves] = p.adjust(P[leaves], "fdr") <= q
  
  # apply Structured Holm
  if(min(P) >= q/m_reduced){
    rejections["before_filter", "Structured_Holm",] = rep(FALSE, m_reduced)
  } else{
    rejections["before_filter", "Structured_Holm",] = structuredHolm(G_SH_reduced, 
                                                                     alpha_max = q, 
                                                                     pvalues=P)@rejected
  }
  
  # apply Yekutieli
  q_Yekutieli = q/(2*max(get_depths(G_reduced$Pa)))
  rejections["before_filter", "Yekutieli",] = Yekutieli(P, G_Y_reduced, q_Yekutieli)
  
  # apply Focused BH
  rejections["before_filter", "Focused_BH",] = FocusedBH(filter_name, P, G_reduced, q)
  
  # Filter the results
  for(method in methods){
    R = rejections["before_filter", method,]
    rejections["after_filter", method, ] = run_filter(P, R, G_reduced, filter_name)
  }
  
  # massage into data frame
  tree_file = sprintf("%s/data/raw/biobank/coding19.tsv", base_dir)
  tree = read_tsv(tree_file, col_types = c("cciic"))
  df_rejections = melt(rejections)
  names(df_rejections) = c("filter_status", "method", "node", "rejected")
  df_rejections = as_tibble(df_rejections) %>% mutate(node = as.character(node))
  node_to_meaning = tree %>% select(coding, meaning) %>% rename(node = coding)
  df_rejections = df_rejections %>% left_join(node_to_meaning, by = "node")
  
  # save to file
  write_tsv(df_rejections, rejections_filename)
  
  # add provenance information to output file
  call_info = "#\n#\n### FUNCTION CALL: ###\n# source(\"PheWAS_data_analysis.R\")"
  write(call_info, rejections_filename, append = TRUE)
  add_git_info(rejections_filename)
  add_R_session_info(rejections_filename)
}
