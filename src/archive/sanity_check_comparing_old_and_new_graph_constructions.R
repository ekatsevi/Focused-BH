load(sprintf("%s/data/processed/GO_reduced_graph.Rda", base_dir))
G_new = G
load(sprintf("%s/data/processed/reduced_graph.Rda", base_dir))
G_old = G

length(G_old$C) == G_new$m

all(sapply(1:G_new$m, function(i)(all(G_new$node_names[G_new$Pa[[i]]] == G_old$Pa[[i]]))))

all(sapply(1:G_new$m, function(i)(all(G_new$node_names[G_new$C[[i]]] == G_old$C[[i]]))))

identical(G_new$sets, unname(G_old$gene_sets))

load(sprintf("%s/data/processed/GO_reduced_graph_100.Rda", base_dir))
G_new = G
load(sprintf("%s/data/processed/reduced_graph_100.Rda", base_dir))
G_old = G

length(G_old$C) == G_new$m

all(sapply(1:G_new$m, function(i)(all(G_new$node_names[G_new$Pa[[i]]] == G_old$Pa[[i]]))))

all(sapply(1:G_new$m, function(i)(all(G_new$node_names[G_new$C[[i]]] == G_old$C[[i]]))))

identical(G_new$sets, unname(G_old$gene_sets))

######################################################

load(sprintf("%s/data/processed/GO_reduced_graph.Rda", base_dir))
G_new = G
load(sprintf("%s/data/processed/old_GO_reduced_graph.Rda", base_dir))
G_old = G

G_old$m == G_new$m

identical(G_old$Pa, G_new$Pa)

identical(G_old$C, G_new$C)

all(sapply(1:G_new$m, function(i)(all(G_old$sets[[i]] == G_new$item_names[G_new$sets[[i]]]))))

load(sprintf("%s/data/processed/GO_reduced_graph_100.Rda", base_dir))
G_new = G
load(sprintf("%s/data/processed/old_GO_reduced_graph_100.Rda", base_dir))
G_old = G

G_old$m == G_new$m

identical(G_old$Pa, G_new$Pa)

identical(G_old$C, G_new$C)

all(sapply(1:G_new$m, function(i)(setequal(G_old$sets[[i]], G_new$item_names[G_new$sets[[i]]]))))