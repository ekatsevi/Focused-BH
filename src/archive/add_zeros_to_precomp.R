for(file in list.files("/home/ekatsevi/project-files/focused-bh/precomp/GO_Simes")){
  filename = sprintf("/home/ekatsevi/project-files/focused-bh/precomp/GO_Simes/%s", file)
  df = read_tsv(filename)
  df = df %>% add_row(pvalue = 0, V_oracle = 0, V_permutation = 0, .before = 1)
  write_tsv(df, filename)
}