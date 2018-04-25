#import stats
generateLen = function(vec,filt)
{
  if(length(vec)!=sum(!filt))
  {
    stop("Wrong lengths detected")
  }
  out = rep(NA,length(filt))
  out[!filt] = vec
  out
}
gqtl = function(bulkA,bulkB = NULL,geneticMap,splitVar,getDistortion=T,kern = function(d){(1-d^3)^3 / sum((1-d^3)^3)},W=25,lim= 0.9)
{
  if(length(geneticMap)!=length(splitVar))
  {
    stop("geneticMap and splitVar must have same length")
  }
  if(!getDistortion %in% c(F,T))
    stop("getDistortion has to be either FALSE or TRUE")
  if(ncol(bulkA)!=2)
    stop("bulkA is a 2 column matrix, each column has counts of an allele")
  if(!is.null(bulkB))
  {
    if(nrow(bulkB)!=nrow(bulkA))
      stop("bulkB and bulkA must have the same number of rows (variants)")
    if(ncol(bulkB)!=2)
      stop("bulkB should be a 2 column matrix, each column has counts for an allele")
    n1 = bulkB[,1]
    n3 = bulkB[,2]
  }
  n2 = bulkA[,1]
  n4 = bulkA[,2]
  filt = getq(n2,n4)
  filteredOut = filt>=lim | filt<=1-lim
  print(sprintf("Removed %s variants based on allele frequency criterion",sum(filteredOut)))
  n2 = n2[!filteredOut]
  n4 = n4[!filteredOut]
  geneticMapOut = geneticMap
  geneticMap = geneticMap[!filteredOut]
  splitVarOut = splitVar
  splitVar = splitVar[!filteredOut]
  GstatA = Gstat(n2,n4)
  if(any(is.na(GstatA)))
  {
    stop(sprintf("Identified %d NA values for G statistic based on bulkA. Make sure to remove rows with all 0 values!",sum(is.na(GstatA))))
  }
  res = list()
  GsmoothA = NULL
  for(v in unique(splitVar))
  {
    curInd= splitVar==v
    GsmoothA = c(GsmoothA,smoothQTL(GstatA[curInd],geneticMap[curInd],Kern="kern",W = W))
  }
  paramA = getPars(GsmoothA)
  res = list()
  res["GsmoothA"] = list(generateLen(GsmoothA,filteredOut))
  res["paramA"] = list(paramA)
  res["GstatA"] = list(generateLen(GstatA,filteredOut))
  res["pvalA"] = list(generateLen(stats::plnorm(GsmoothA,meanlog = paramA[1],sdlog = paramA[2],lower.tail = F),filteredOut))
  if(!is.null(bulkB))
  {
    n1 = n1[!filteredOut]
    n3 = n3[!filteredOut]
    if(getDistortion)
    {
      qRawA = getq(n2,n4)
      qRawB = getq(n1,n3)
      qSmoothA = NULL
      qSmoothB = NULL
      for(v in unique(splitVar))
      {
        qSmoothA = c(qSmoothA,smoothQTL(qRawA[splitVar==v],geneticMap[splitVar==v],Kern="kern",W = W))
        qSmoothB = c(qSmoothB,smoothQTL(qRawB[splitVar==v],geneticMap[splitVar==v],Kern="kern",W = W)) 
      }
      res["qSmoothA"] = list(generateLen(qSmoothA,filteredOut))
      res["qSmoothB"] = list(generateLen(qSmoothB,filteredOut))
    }
    else
    {
      qSmoothA = 0.5
    }
    GstatB = Gstat(n2,n4,n1,n3,q = qRawA)
    GsmoothB = NULL
    for(v in unique(splitVar))
    {
      GsmoothB = c(GsmoothB,smoothQTL(G = GstatB[splitVar==v],Map = geneticMap[splitVar==v],Kern = "kern",W = W))
    }
    paramB = getPars(GsmoothB)
    pvalB = stats::plnorm(GsmoothB,meanlog = paramB[1],sdlog = paramB[2],lower.tail = F)
    res["removed"] = list(which(filteredOut))
    res["GstatB"] = list(generateLen(GstatB,filteredOut))
    res["GsmoothB"] = list(generateLen(GsmoothB,filteredOut))
    res["paramB"] = list(paramB)
    res["pvalB"] = list(generateLen(pvalB,filteredOut))
    res["geneticMap"] = list(geneticMapOut)
    res["splitVar"] = list(splitVarOut)
  }
  return(res)
}

gqtlDF = function(res)
{
  df = with(res,data.frame(GstatA,GsmoothA,pvalA,GstatB,GsmoothB,pvalB,qSmoothA = qSmoothA,qSmoothB = qSmoothB))
}