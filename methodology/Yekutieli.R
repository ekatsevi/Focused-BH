# Wrapper around Yekutieli method
# P:           vector of m p-values
# G_Yekutieli: graph in Yekutieli format
# q:           target FDR level
Yekutieli = function(P, G_Yekutieli, q){
  unadj.p.values = c(P, 0)
  m = length(P)
  names(unadj.p.values) = sapply(1:(m+1), as.character)
  hyp.tree = my.hFDR.adjust(unadj.p.values, G_Yekutieli, q)
  R = hyp.tree@p.vals[[2]] <= q
  R[is.na(R)] = FALSE
  return(R[1:m])
}

# Apply Yekutieli method (same code as from StructSSI, but 
# with small modification due to bug in multtest package)
# unadjp:           vector of m p-values
# tree.el:          graph in Yekutieli format
# alpha:            target FDR level
my.hFDR.adjust = function(unadjp, tree.el, alpha = 0.05) 
{
  if (is.null(names(unadjp))) {
    names(unadjp) <- 1:length(unadjp)
  }
  if (!all(names(unadjp) %in% unique(as.vector(tree.el)))) {
    stop("Names of elements in unadjp do not match names of tree nodes")
  }
  p.vals <- data.frame(unadjp, adjp = NA)
  hyp.tree.unadjusted <- new("hypothesesTree", tree = tree.el, 
                             p.vals = p.vals, alpha = alpha)
  root <- FindRoot(hyp.tree.unadjusted@tree)
  if (hyp.tree.unadjusted@p.vals[root, "unadjp"] > alpha) {
    warning("Root hypothesis p-value equal to ", 
            hyp.tree.unadjusted@p.vals[root, "unadjp"], 
            ". Fail to reject any hypotheses, terminating procedure.")
    hyp.tree.unadjusted@p.vals[, "adj.significance"] <- "-"
    return(hyp.tree.unadjusted)
  }
  # made change here, because of bug in multtest package
  hyp.tree <- my.hFDR.internal(hyp.tree.unadjusted)
  hyp.tree@p.vals[root, "adjp"] <- hyp.tree@p.vals[root, "unadjp"]
  hyp.tree@p.vals[, "adj.significance"] <- SignificanceStars(alpha, 
                                                             hyp.tree@p.vals[, "adjp"])
  return(hyp.tree)
}

# Yekutieli method helper (same code as from StructSSI, but 
# with small modification due to bug in multtest package)
# hyp.tree: StructSSI object containing adjustted and unadjusted p-values
my.hFDR.internal = function (hyp.tree) 
{
  tree <- graph.edgelist(hyp.tree@tree)
  tree.el.tmp <- data.frame(parent = hyp.tree@tree[, 1], 
                            child = hyp.tree@tree[, 2], stringsAsFactors = F)
  root <- FindRoot(tree.el.tmp)
  children <- tree.el.tmp[which(tree.el.tmp$parent == root), "child"]
  children.p.vals <- hyp.tree@p.vals[children, ]
  # made change here, due to small bug in multtest package
  if(length(children) == 1){
    children.p.vals$adjp = children.p.vals$unadjp
  } else{
    adjust <- mt.rawp2adjp(children.p.vals$unadjp, "BH")
    children.p.vals$adjp <- adjust$adjp[order(adjust$index), 
                                        "BH"]
  }
  hyp.tree@p.vals[children, "adjp"] <- children.p.vals$adjp
  rejected <- rownames(children.p.vals)[which(children.p.vals$adjp < 
                                                hyp.tree@alpha)]
  for (child in rejected) {
    subcomp <- subcomponent(tree, child, "out")
    if (length(subcomp) > 1) {
      subtree.igraph <- induced.subgraph(graph = tree, 
                                         vids = subcomp)
      subtree.names <- get.vertex.attribute(subtree.igraph, 
                                            "name")
      subtree <- new("hypothesesTree", alpha = hyp.tree@alpha, 
                     tree = get.edgelist(subtree.igraph))
      subtree@p.vals <- hyp.tree@p.vals[subtree.names, 
                                        ]
      hyp.tree@p.vals[subtree.names, ] <- my.hFDR.internal(subtree)@p.vals
    }
  }
  return(hyp.tree)
}

# convert graph structure to Yekutieli format
convert_to_Yekutieli_format = function(G){
  roots = which(sapply(G$Pa, length) == 0)
  # add an imaginary root node, whose roots are current root nodes
  G_augmented = G
  G_augmented$m = G_augmented$m + 1
  G_augmented$C[[m+1]] = roots
  G_augmented$Pa[[m+1]] = integer(0)
  for(root in roots){
    G_augmented$Pa[[root]] = m + 1
  }
  
  G_Yekutieli = t(do.call(cbind, 
                          lapply(1:G_augmented$m, 
                                 function(i)(if(length(G_augmented$C[[i]]) > 0) sapply(G_augmented$C[[i]], 
                                                                                       function(j)(c(as.character(i), as.character(j)))) else matrix(0,2,0)))))
  return(G_Yekutieli)
}
