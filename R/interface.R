#import stats

#' Title Internal funcion to turn values into NA based on filter
#'
#' @param vec vector
#' @param filt filter vector (boolean)
#' @return
#' output with filtered values turned to NA. This function is basic, only convenience to avoid repeating code.
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

#' Title Run X-QTL analysis
#'
#' @param bulkA A two column matrix with the allele counts in the low/control group 
#' @param bulkB A two column matrix with the allele counts in the high/treatment group 
#' @param geneticMap A vector of the genetic distances. All chromosomes are combined into one vector, and within each chromosome the values should monotonically increase. The function is compatible with the genetic map vector produced by r/qtl. has to be the same length as the row number of bulkA and bulkB
#' @param splitVar A vector used to split the genetic map and the allele counts to chromosomes. Has to be the same length as the row number of bulkA and bulkB
#' @param getDistortion Should segregration distortion be calculated from the data? If FALSE, X-QTL will run assuming baseline allele frequencies are all 0.5. 
#' @param kern The kernel function for the smoother. This is currently a placeholder since the function uses a C++ implementation of the smoother with the kernel baked in.
#' @param W The window size for the smoother. When smoothing the value for an SNV, a window of size W will be taken around it. See the Magwene et al 2011 paper for discussion of values. Defaults to 25, the low bound of the values they recommend 
#' @param lim A filter for SNVs which show distortion. This is calculated based on the low/control group. The G statistic becomes highly unstable when considering SNVs that are nearly fixed, so they are filtered out.
#'
#' @return
#' A list with the following values:
#' GstatA: The G statistic for the control/low group (compared to a baseline of 0.5/0.5). Can be used to determine segregation distortion in the baseline
#' GstatB: The G statistic for the high/treatment group (compared to the control group). This is the statistic used for the X-QTL
#' GsmoothA: The smoothed G statistic for the control/low group. 
#' GsmoothB: The smoothed G statistic for the treatment/high group. Used directly to calculate P-values.
#' PvalA: P-values for segregation distortion in the baseline. Beware! If large parts of the genome show distortion (as in C. elegans) This will not work, since it derives the null distribution robustly from the values themselves.
#' PvalB: XQTL P-values.
#' qSmoothA: The smoothed allele frequencies in the control/low group. Can be useful to observe segregation distortion in the baseline
#' qSmoothB: The smoothed allele frequences in the treatment/high group. Together with qSmoothA, can be used to determine direction of observed QTL
#' geneticMap: The genetic map supplied to the function. This is returned for convenience/sanity check
#' splitVar: The chromosome vector supplied to the function. This is returned for convenience/sanity check
#' @export
gqtl = function(bulkA,bulkB = NULL,geneticMap,splitVar,getDistortion=T,kern = function(d){(1-d^3)^3 / sum((1-d^3)^3)},W=25,lim= 0.99)
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
    GsmoothA = c(GsmoothA,smoothQTL(GstatA[curInd],geneticMap[curInd],W = W))
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
        qSmoothA = c(qSmoothA,smoothQTL(qRawA[splitVar==v],geneticMap[splitVar==v],W = W))
        qSmoothB = c(qSmoothB,smoothQTL(qRawB[splitVar==v],geneticMap[splitVar==v],W = W)) 
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
      GsmoothB = c(GsmoothB,smoothQTL(G = GstatB[splitVar==v],Map = geneticMap[splitVar==v],W = W))
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


#' Title Transform XQTL list to data frame
#'
#' @param res Results from the gqtl function
#'
#' @return
#' Convenience function that simply converts the list returned by gqtl to data frame
#' @export
gqtlDF = function(res)
{
  df = with(res,data.frame(GstatA,GsmoothA,pvalA,GstatB,GsmoothB,pvalB,qSmoothA = qSmoothA,qSmoothB = qSmoothB))
}