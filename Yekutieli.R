Yekutieli = function(P, G_Yekutieli, q){
  unadj.p.values = c(P, 0)
  m = length(P)
  names(unadj.p.values) = sapply(1:(m+1), as.character)
  hyp.tree = my.hFDR.adjust(unadj.p.values, G_Yekutieli, q)
  R = hyp.tree@p.vals[[2]] <= q
  R[is.na(R)] = FALSE
  rejections["Yekutieli",] = R[1:m]
}

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
  hyp.tree <- my.hFDR.internal(hyp.tree.unadjusted)
  hyp.tree@p.vals[root, "adjp"] <- hyp.tree@p.vals[root, "unadjp"]
  hyp.tree@p.vals[, "adj.significance"] <- SignificanceStars(alpha, 
                                                             hyp.tree@p.vals[, "adjp"])
  return(hyp.tree)
}

my.hFDR.internal = function (hyp.tree) 
{
  tree <- graph.edgelist(hyp.tree@tree)
  tree.el.tmp <- data.frame(parent = hyp.tree@tree[, 1], 
                            child = hyp.tree@tree[, 2], stringsAsFactors = F)
  root <- FindRoot(tree.el.tmp)
  children <- tree.el.tmp[which(tree.el.tmp$parent == root), "child"]
  children.p.vals <- hyp.tree@p.vals[children, ]
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
