ontology_file = sprintf("%s/data/raw/go-basic.obo", base_dir)
annotations_file = sprintf("%s/data/raw/goa_human.gaf", base_dir)
gene_list_file = sprintf("%s/data/raw/gene_list.txt", base_dir)
older_terms_file = sprintf("%s/data/raw/go_201410-termdb-tables/term.txt", base_dir)

annotations_processed_file = sprintf("%s/data/processed/annotations.Rda", base_dir)
reduced_graph_file = sprintf("%s/data/processed/reduced_graph.Rda", base_dir)
reduced_graph_50_file = sprintf("%s/data/processed/reduced_graph_50.Rda", base_dir)

if(file.exists(annotations_file) & file.exists(reduced_graph_file)){
  cat(sprintf("GO files already parsed!"))
} else{
ontology_raw = get_ontology(ontology_file, 
                        propagate_relationships = c("is_a", "part_of"), 
                        extract_tags = "everything")
annotations_raw = read_tsv(annotations_file, comment = "!", col_types = "--c-c-----c------", col_names = FALSE)
names(annotations_raw)  = c("gene", "id", "aliases")

biological_processes = unname(ontology_raw$id[ontology_raw$namespace == "biological_process"])

if(!file.exists(annotations_processed_file)){
  annotations = vector("list", length(biological_processes))
  names(annotations) = biological_processes
  for(idx in 1:length(biological_processes)){
    print(idx)
    GO_id = biological_processes[idx]
    genes = annotations_raw %>% filter(id == GO_id) %>% pull(gene)
    for(id in c(GO_id, ontologyIndex::get_ancestors(ontology_raw, GO_id))){
      annotations[[id]] = c(annotations[[id]], genes)
    }
  }
  
  save(annotations, file = annotations_processed_file)
}

if(!file.exists(reduced_graph_file)){
  load(annotations_file)
  
  gene_list = read_tsv(gene_list_file, col_names = FALSE)
  gene_list = gene_list$X1
  older_terms = read_tsv(older_terms_file, col_names = FALSE, col_types = "iccciii")
  older_terms = older_terms %>% filter(grepl("GO", X4)) %>% select(X4) %>% pull()
  
  G_full = c()
  G_full$C = ontology_raw$children[biological_processes]
  G_full$Pa = ontology_raw$parents[biological_processes]
  to_remove = biological_processes[ontology_raw$obsolete[ontology_raw$namespace == "biological_process"] | 
                                     sapply(annotations, length) == 0 | 
                                     !(biological_processes %in% older_terms)]
  
  G = subset_graph(to_remove, G_full)
  remaining_nodes = names(G$C)
  G$gene_sets = sapply(annotations[remaining_nodes], unique)
  save(G, file = reduced_graph_file)
}

if(!file.exists(reduced_graph_50_file)){
  load(reduced_graph_file)
  num_annotations = sapply(G$gene_sets, length)
  to_remove = names(which(num_annotations < 50 | num_annotations > 10000))
  
  G = subset_graph(to_remove, G)
  remaining_nodes = names(G$C)
  G$gene_sets = annotations[remaining_nodes]
  save(G, file = reduced_graph_50_file)
}

}