rm(list = ls())
machine = "local"
source("setup.R")

######### DATA ANALYSIS PARAMETERS ###################
q = 0.1                        # FDR target level
filter_name = "REVIGO"         # filter to use
methods = c("Focused_BH",      # methods to apply
            "Structured_Holm", 
            "BH")

############# PROCESSED FILES ###################
# enrichment p-values from GOrilla
GOrilla_output = sprintf("%s/data/raw/GO/GO.xls", base_dir)
# preprocessed GO graph
reduced_graph_file = sprintf("%s/data/processed/GO/GO_reduced_graph.Rda", base_dir)
# final sets of rejections for each method
rejections_filename = sprintf("%s/data/processed/GO/REVIGO_rejections.tsv", base_dir)

############# READ IN P-VALUES AND GRAPH STRUCTURE, APPLY METHODS ###################
if(!file.exists(rejections_filename)){
  # read in GOrilla data, graph structure
  GO <- read_excel(GOrilla_output)
  load(reduced_graph_file)
  
  # define a set of p-values for each node in the graph
  ids = G$node_names
  m = G$m
  P = rep(1, m)
  significant_ids = match(GO$`GO Term`, ids)
  cat(sprintf("%d out of %d of the significant GO terms are not in the older ontology available to REVIGO.\n", 
              sum(is.na(significant_ids)), length(significant_ids)))
  P[significant_ids[!is.na(significant_ids)]] = GO$`P-value`[!is.na(significant_ids)]
  names(P) = ids

  # array to store rejections
  num_methods = length(methods)
  rejections = array(FALSE, 
                     dim = c(2,num_methods, m), 
                     dimnames = list(c("before_filter", "after_filter"), 
                                     methods, 
                                     G$node_names))
  
  # apply BH
  rejections["before_filter", "BH",] = p.adjust(P, "fdr") <= q
  
  # apply Structured Holm
  if(min(P) >= q/m){
    rejections["before_filter", "Structured_Holm",] = rep(FALSE, m)
  } else{
    rejections["before_filter", "Structured_Holm",] = structuredHolm(G_SH, 
                                                                     alpha_max = q, 
                                                                     pvalues=P)@rejected
  }
  
  # apply Focused BH
  rejections["before_filter", "Focused_BH",] = FocusedBH(filter_name, P, G, q)
  
  # Filter the results
  for(method in methods){
    R = rejections["before_filter", method,]
    rejections["after_filter", method, ] = run_filter(P, R, G, filter_name)
  }
  
  # massage results into data frame
  df_rejections = melt(rejections)
  names(df_rejections) = c("filter_status", "method", "node", "rejected")
  df_rejections = as_tibble(df_rejections) %>% mutate(node = as.character(node))
  
  # save results to file
  write_tsv(df_rejections, rejections_filename)
  
  # add provenance information to output file
  call_info = "#\n#\n### FUNCTION CALL: ###\n# source(\"REVIGO_data_analysis.R\")"
  write(call_info, rejections_filename, append = TRUE)
  add_git_info(rejections_filename)
  add_R_session_info(rejections_filename)
}