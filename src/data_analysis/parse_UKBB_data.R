############# RAW FILES ###################
# raw ICD tree file (downloaded from https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=19)
tree_file = sprintf("%s/data/raw/biobank/coding19.tsv", base_dir)
# Names of HLA alleles (downloaded from biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/ukb_hla_v2.txt)
allele_file = sprintf("%s/data/raw/biobank/ukb_hla_v2.txt", base_dir)
# Raw phenotype file (requires UKBB application)
phen_file = sprintf("%s/data/raw/biobank/ukb25261.csv", base_dir)

############# PROCESSED FILES ###################
# processed ICD tree file
graph_file = sprintf("%s/data/processed/biobank/ICD_graph.Rda", base_dir)
# processed diseases file
disease_file = sprintf("%s/data/processed/biobank/diagnosed_disease.tab", base_dir)
# further processed diseases file
HES_affected_file = sprintf("%s/data/processed/biobank/HES_affected.Rda", base_dir)
# file storing number of cases per disease
num_affected_file = sprintf("%s/data/processed/biobank/HES_number_affected.tsv", base_dir)
# processed HLA file
HLA_file = sprintf("%s/data/processed/biobank/HLA.tab", base_dir)
# further processed HLA file
HLA_split_file = sprintf("%s/data/processed/biobank/HLA_split.tsv", base_dir)

############# READ IN AND PARSE ICD-10 GRAPH STRUCTURE ###################
if(!file.exists(graph_file)){
  cat(sprintf("Parsing ICD tree...\n"))
  
  # read the raw tree structure
  tree = read_tsv(tree_file, col_types = c("cciic"))
  m = nrow(tree)
  
  # create an index of the ICD-10 codes
  id2idx = list()
  for(j in 1:m){
    id2idx[as.character(tree$node_id[j])] = j
  }
  
  # populate children and parents of each node
  C = vector("list", m)
  Pa = vector("list", m)
  for(j in 1:m){
    node_id = tree$node_id[j]
    parent_id = tree$parent_id[j]
    if(parent_id != 0){
      parent_idx = id2idx[[as.character(parent_id)]]
      Pa[[j]] = c(Pa[[j]], parent_idx)
      C[[parent_idx]] = c(C[[parent_idx]], j)
    }
  }
  
  # extract coding for each ICD-10 term
  node_names = tree$coding
  
  # extract leaves and treat each leaf as an item
  leaves = sapply(C, length) == 0
  item_names = node_names[leaves]
  num_items = sum(leaves)
  
  # index the leaf nodes and assign leaf nodes
  # to their ICD-10 codes
  item_index = integer(G$num_items)
  item_index[leaves] = 1:num_items
  sets = vector("list", m)
  for(node in 1:m){
    R_delta = logical(m)
    R_delta[node] = TRUE
    sets[[node]] = item_index[get_descendants(R_delta, C) & leaves]
  }
  
  # create graph structure for Focused BH
  G = c()
  G$m = m
  G$num_items = num_items
  G$C = C
  G$Pa = Pa
  G$sets = sets
  G$node_names = node_names
  G$item_names = item_names
  
  # convert graph to Structured Holm and Yekutieli formats
  G_SH = new("DAGstructure", parents = G$Pa, children = G$C, sets = G$sets, twoway = FALSE)
  G_Yekutieli = convert_to_Yekutieli_format(G)
  
  # save to file
  save(G, G_SH, G_Yekutieli, file = graph_file)
}

############# EXTRACT DISEASES AND HLA GENOTYPES ###################
if(!file.exists(disease_file) | !file.exists(HLA_file)){
  # read in the raw phenotype file
  phenotypes_raw = read_csv(phen_file, guess_max=500000)
  
  # extract fields for diagnosed diseases and save to file
  fields = sapply(0:379, function(i)(paste("f.41202.0.", i, sep = "")))
  HES = phenotypes_raw %>% select(fields)
  write_tsv(HES, path = disease_file)
  
  # extract fields for HLA and save to file
  fields = "f.22182.0.0"
  HLA = phenotypes_raw %>% select(fields)
  write_tsv(HLA, path = HLA_file)
}

############# PARSE DISEASE INFORMATION ###################
if(!file.exists(HES_affected_file)){
  # read in disease information and tree
  HES = read_tsv(disease_file, col_types = cols(.default = col_character()))
  tree = read_tsv(tree_file, col_types = c("cciic"))
  
  # convert data to list of affected individuals per disease
  n = nrow(HES)
  p = ncol(HES)
  diseases = unique(unlist(sapply(1:p, function(j)(unique(HES[[j]])))))
  diseases = diseases[!is.na(diseases)]
  m0 = length(diseases)
  affected = vector("list", m0)
  names(affected) = diseases
  for(i in 1:n){
    if(i %% 10000 == 0){
      print(i)
    }
    for(j in 1:p){
      disease = HES[i,j]
      if(!is.na(disease)){
        affected[[as.character(disease)]] = c(affected[[as.character(disease)]], i)
      } else{
        break
      }
    }
  }
  
  # re-index diseases based on their node ids
  new_names = sapply(names(affected), function(name)(tree$node_id[tree$coding == name]))
  names(affected) = new_names
  
  # propagate disease status up the tree
  for(j in which(tree$selectable == "N")){
    affecteds = c()
    to_explore = G$C[[j]] # should be non-empty
    while(length(to_explore) > 0){
      k = to_explore[1]
      if(tree$selectable[k] == "Y"){
        affecteds = c(affecteds, affected[[as.character(tree$node_id[k])]]) 
      } else{
        to_explore = c(to_explore, G$C[[k]])
      }
      if(length(to_explore) == 1){
        to_explore = numeric(0)
      } else{
        to_explore = to_explore[2:length(to_explore)]
      }
    }
    affected[[as.character(tree$node_id[j])]] = affecteds
  }
  
  # calculate number of cases per disease
  num_affected = sapply(1:m, function(k)(length(affected[[as.character(tree$node_id[k])]])))
  df_num_affected = tibble(coding = tree$coding, num_affected = num_affected)
  
  # save to file
  save(affected, file = HES_affected_file)
  write_tsv(df_num_affected, path = num_affected_file)
}

############# PARSE HLA INFORMATION ###################
if(!file.exists(HLA_split_file)){
  # read HLA data
  HLA = read_tsv(HLA_file)
  HLA_names = read_tsv(allele_file, col_names = FALSE) %>% pull()
  
  # split the genotypes for each individual
  n = nrow(HLA)
  p = length(strsplit(HLA$f.22182.0.0[1], split = ",")[[1]])
  HLA_split = matrix(0, n, p)
  for(i in 1:n){
    if(i %% 10000 == 0){
      print(i)
    }
    HLA_split[i,] = as.numeric(strsplit(HLA$f.22182.0.0[i], split = ",")[[1]])
  }
  HLA_split_df = data.frame(HLA_split)
  
  # label with the names of the alleles
  names(HLA_split_df) = HLA_names
  
  # save to file
  write_tsv(HLA_split_df, path = HLA_split_file, col_names = TRUE)
}