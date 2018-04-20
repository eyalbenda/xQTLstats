#import robust
#import stats

getq = function(n1,n3)
{
  q = n1/(n1+n3)
  q[q==1] = 0.99
  q[q==0] = 0.01
  q
}
Gstat = function(n2,n4,n1=NULL,n3=NULL,q=0.5)
{
  if(is.null(n1))
  {
    if(!is.null(n3))
      stop("if n1 is specified, n3 must also be specified")
  }
  if(is.null(n3))
    if(!is.null(n1))
      stop("if n3 is specified, n1 must also be specified")
  if(is.null(n1))
  {
    print("unspecified n1 and n3, simulating based on q")
    cov1 = n2 + n4
    n1 = cov1*q
    n3 = cov1*(1-q)
  }
  C = apply(cbind(n1+n3,n2+n4),1,mean)
  ((1-q)*(n2-n1) + q*(n3-n4))^2/(2*C*q*(1-q))
}


 
smoothQTLwrap = function(G,Map,Kern="kern",W=20)
{
  GsmoothRaw = smoothQTL(G = G[!is.na(G)],Map = Map[!is.na(G)],Kern = Kern,W = W)
  if(any(is.na(G)))
  {
    Gsmoothed = rep(NA,length(G))
    Gsmoothed[!is.na(G)] = GsmoothRaw
  } else
  {
    Gsmoothed = GsmoothRaw
  }
  return(Gsmoothed)
}
#Your responsibility to run per chrom
Gsmooth = function(G,kern = function(d){(1-d^3)^3 / sum((1-d^3)^3)},map,W=20)
{
  if(length(G)!=length(map))
    stop("Length of G statistic vector must be the same as length of map vector")
  if((max(map)-min(map))< W*2)
    stop("For safety, the window W can't be less than 2x the span of the map")
  cutoff = W/2
  Gsmoothed = rep(NA,length(G))
  firstEdge = findInterval(cutoff,map-map[1])+1
  secondEdge = rev(which((map[length(map)]-map) > (cutoff)))[1]
  G.padded = c(G[firstEdge:2],G,G[(length(G)-1):secondEdge])
  map.padded = c(map[1]-map[firstEdge:2]+ map[1],map,map[length(map)] + map[length(map)]-map[(length(map)-1):secondEdge])
  k = 1
  for(i in (firstEdge):(firstEdge+length(G)-1))
  {
    startEnd = findInterval(c(cutoff,-cutoff),map.padded-map.padded[i])
    startInd = startEnd[1]+1
    endInd = startEnd[2]+1
    distS = map.padded[startInd:endInd] - map.padded[i]
    distS = abs(distS/max(distS))
    Gsmoothed[k] = kern(d = distS)%*%G.padded[startInd:endInd]
    k = k+1
  }
  Gsmoothed
}

getPars = function(Gsmoothed)
{
  robust::fitdstnRob(Gsmoothed,"lognorm")$estimate
}

