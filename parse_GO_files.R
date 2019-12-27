# raw data filenames
ontology_file = sprintf("%s/data/raw/go-basic.obo", base_dir)
annotations_file = sprintf("%s/data/raw/goa_human.gaf", base_dir)
gene_list_file = sprintf("%s/data/raw/gene_list.txt", base_dir)
older_terms_file = sprintf("%s/data/raw/go_201410-termdb-tables/term.txt", base_dir)

# processed data filenames
GO_graph_file = sprintf("%s/data/processed/GO_graph.Rda", base_dir)
GO_reduced_graph_file = sprintf("%s/data/processed/GO_reduced_graph.Rda", base_dir)
GO_reduced_graph_100_file = sprintf("%s/data/processed/GO_reduced_graph_100.Rda", base_dir)

# convert the GO graph to my format
if(!file.exists(GO_graph_file)){
  cat(sprintf("Parsing Gene Ontology (this will take several minutes)...\n"))
  
  # read the ontology structure
  ontology_raw = get_ontology(ontology_file, 
                              propagate_relationships = c("is_a", "part_of"), 
                              extract_tags = "everything")
  # read the gene annotations
  annotations_raw = read_tsv(annotations_file, comment = "!", 
                             col_types = "--c-c-----c------", 
                             col_names = FALSE)
  names(annotations_raw)  = c("gene", "id", "aliases")
  
  # extract terms in the biological process ontology
  biological_processes = unname(ontology_raw$id[ontology_raw$namespace == "biological_process"])
  
  # propagate annotations upward through the ontology
  annotations = vector("list", length(biological_processes))
  names(annotations) = biological_processes
  for(idx in 1:length(biological_processes)){
    GO_id = biological_processes[idx]
    genes = annotations_raw %>% filter(id == GO_id) %>% pull(gene)
    for(id in c(GO_id, ontologyIndex::get_ancestors(ontology_raw, GO_id))){
      annotations[[id]] = c(annotations[[id]], genes)
    }
  }
  annotations = sapply(annotations, unique)
  
  # index the GO terms
  C = ontology_raw$children[biological_processes]
  Pa = ontology_raw$parents[biological_processes]
  m = length(C)
  GO_ids = names(C)
  stopifnot(all(names(Pa) == GO_ids) & length(Pa) == m)
  term_index = c()
  for(i in 1:m){
    term_index[[GO_ids[i]]] = i
  }
  Pa_indexed = unname(sapply(Pa, function(parents_list)(unname(term_index[parents_list]))))
  C_indexed = unname(sapply(C, function(children_list)(unname(term_index[children_list]))))
  
  # index the genes
  gene_names = unique(unlist(annotations))
  num_genes = length(gene_names)
  gene_index = c()
  for(i in 1:num_genes){
    gene_index[[gene_names[i]]] = i
  }
  gene_sets = unname(sapply(GO_ids, function(id)(unname(gene_index[unlist(annotations[id])]))))
  
  # create a DAG structure in my format and save to file
  G = c()
  G$m = m
  G$num_items = num_genes
  G$Pa = Pa_indexed
  G$C = C_indexed
  G$sets = gene_sets
  G$node_names = GO_ids
  G$item_names = gene_names
  save(G, file = GO_graph_file)
}
 
# create GO graph for real data analysis, removing terms with no annotations and terms that
# were not present in the version of GO used for REVIGO
if(!file.exists(GO_reduced_graph_file)){
  cat(sprintf("Preparing reduced GO graph for data analysis...\n"))
  
  # load the full GO graph
  load(GO_graph_file) 
  
  # list of terms available to REVIGO
  older_terms = read_tsv(older_terms_file, col_names = FALSE, col_types = "iccciii")
  older_terms = older_terms %>% filter(grepl("GO", X4)) %>% dplyr::select(X4) %>% pull()

  # nodes that are not in REVIGO or that have no annotations  
  to_remove = !(G$node_names %in% older_terms) | sapply(G$sets, length) == 0

  # restrict the graph to the remaining nodes
  G = subset_graph(to_remove, G)
  
  # convert graph to Structured Holm format
  G_SH = new("DAGstructure", parents = G$Pa, children = G$C, sets = G$sets, twoway = FALSE)
  
  # save to file
  save(G, G_SH, file = GO_reduced_graph_file)
}

# create GO graph for numerical experiments, removing terms with over 100 annotations
if(!file.exists(GO_reduced_graph_100_file)){
  cat(sprintf("Preparing reduced GO graph for numerical experiments...\n"))
  
  # load graph from data analysis
  load(GO_reduced_graph_file)
  
  # remove terms with over 100 annotations
  to_remove = sapply(G$sets, length) > 100
  G = subset_graph(to_remove, G)
  
  # convert graph to Structured Holm format
  G_SH = new("DAGstructure", parents = G$Pa, children = G$C, sets = G$sets, twoway = FALSE)

  # save to file
  save(G, G_SH, file = GO_reduced_graph_100_file)
}