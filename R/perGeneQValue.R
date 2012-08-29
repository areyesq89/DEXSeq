##--------------------------------------------------
## get q-value for each gene
##--------------------------------------------------
perGeneQValue = function(ecs, p = "pvalue", method = perGeneQValueExact) {
  wTest     = which(fData(ecs)$testable)
  ## use only those exons that were testable
  pvals     = fData(ecs)[[p]][wTest]
  ## 'factor' removes ununsed levels
  geneID    = factor(fData(ecs)$geneID[wTest])
  geneSplit = split(seq(along=geneID), geneID)

  ## summarise p-values of exons for one gene: take the minimum
  pGene = sapply(geneSplit, function(i) min(pvals[i]))
  stopifnot(all(is.finite(pGene)))

  ## Determine the thetas to be used
  theta = unique(sort(pGene))

  ## compute q-values associated with each theta
  q = method(pGene, theta, geneSplit)

  ## return a named vector of q-values per gene
  res        = rep(NA_real_, length(pGene))
  names(res) = names(geneSplit)
  res        = q[match(pGene, theta)]
  stopifnot(!any(is.na(res)))
  return(res)
}

##----------------------------------------------------------------------
## For each value of theta, determine how many minima of random per-exon
## p-values are smaller (using the number of exons per gene)
##----------------------------------------------------------------------
perGeneQValueBySimulation = function(pGene, theta, geneSplit, nperm = 24) {
  nr = sum(listLen(geneSplit))
  pRand = apply(matrix(runif(nr*nperm), nrow=nr, ncol=nperm), 2,
    function(p) sapply(geneSplit, function(i) min(p[i])))

  ## check that the apply/sapply stuff worked as intended
  stopifnot(nrow(pRand)==length(pGene), ncol(pRand)==nperm)

  hTest   = hist(pGene, breaks=c(theta,+Inf), plot=FALSE)
  hRand   = hist(pRand, breaks=c(theta,+Inf), plot=FALSE)
  stopifnot(sum(hTest$counts)==length(pGene),
            sum(hRand$counts)==length(pRand))

  numPos       = cumsum(hTest$counts)
  numFalsePos  = cumsum(hRand$counts)/nperm

  return(numFalsePos/numPos)
}

##--------------------------------------------------
## Exact computation - see methods part of the paper
##---------------------------------------------------
perGeneQValueExact = function(pGene, theta, geneSplit) {
  stopifnot(length(pGene)==length(geneSplit))

  ## Compute the numerator \sum_{i=1}^M 1-(1-theta)^{n_i}
  ## Below we first identify the summands which are the same
  ## (because they have the same n_i), then do the sum via the
  ## mapply
  numExons     = listLen(geneSplit)
  tab          = tabulate(numExons)
  notZero      = (tab>0)
  numerator    = mapply(function(m, n) m * (1 - (1-theta)^n),
                            m = tab[notZero],
                            n = which(notZero))
  numerator    = rowSums(numerator)

  ## Compute the denominator: for each value of theta, the number
  ## of genes with pGene <= theta[i].
  ## Note that in cut(..., right=TRUE), the intervals are
  ## right-closed (left open) intervals.
  bins   = cut(pGene, breaks=c(-Inf, as.vector(theta)), right = TRUE, include.lowest = TRUE)
  counts = tabulate(bins, nbins = nlevels(bins))
  denom  = cumsum(counts)
  stopifnot(denom[length(denom)]==length(pGene))

  return(numerator/denom)
}

