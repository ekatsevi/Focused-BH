library(tidyverse)
library(ontologyIndex)

options(tibble.print_max = Inf)
data_dir="/home/ekatsevi/project-files/filtered_BH/DAG/data/raw"
ontology_file = sprintf("%s/go-basic.obo", data_dir)
annotations_file = sprintf("%s/goa_human.gaf", data_dir)
gene_list_file = sprintf("%s/gene_list.txt", data_dir)
older_terms_file = "/home/ekatsevi/Dropbox/Research/Projects/HierTest/code/REVIGO/go_201410-termdb-tables/term.txt"

ontology_raw = get_ontology(ontology_file, 
                        propagate_relationships = c("is_a", "part_of"), 
                        extract_tags = "everything")
annotations_raw = read_tsv(annotations_file, comment = "!", col_types = "--c-c-----c------", col_names = FALSE)
names(annotations_raw)  = c("gene", "id", "aliases")

biological_processes = unname(ontology_raw$id[ontology_raw$namespace == "biological_process"])

annotations = vector("list", length(biological_processes))
names(annotations) = biological_processes
for(idx in 1:length(biological_processes)){
  print(idx)
  GO_id = biological_processes[idx]
  genes = annotations_raw %>% filter(id == GO_id) %>% pull(gene)
  for(id in c(GO_id, get_ancestors(ontology_raw, GO_id))){
    annotations[[id]] = c(annotations[[id]], genes)
  }
}

save(annotations, file = "/home/ekatsevi/project-files/filtered_BH/DAG/data/processed/annotations.Rda")

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

G = subset_graph_2(to_remove, G_full)
remaining_nodes = names(G$C)
G$gene_sets = annotations[remaining_nodes]
save(G, file = "/home/ekatsevi/project-files/filtered_BH/DAG/data/processed/reduced_graph.Rda")