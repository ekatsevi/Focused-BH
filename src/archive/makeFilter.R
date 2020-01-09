if(filter_version == 0){
  filter_type = "trivial"
  base_weights_type = "constant"
}

if(filter_version == 1){
  filter_type = "soft_outer_nodes"
  base_weights_type = "constant"
}

if(filter_version == 2){
  filter_type = "soft_outer_nodes"
  base_weights_type = "IC"
}

if(filter_version == 3){
  filter_type = "hard_outer_nodes"
  base_weights_type = "constant"
}

if(filter_version == 4){
  filter_type = "hard_outer_nodes"
  base_weights_type = "IC"
}

if(filter_version == 5){
  filter_type = "hard_outer_nodes"
  base_weights_type = "soft"
}