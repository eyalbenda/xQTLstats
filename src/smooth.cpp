#include "RcppArmadillo.h"
#include <algorithm>


using namespace Rcpp;
//' Find interval
//' 
//' This takes the fast findInterval2 implementation from Hadley Wickham (http://adv-r.had.co.nz/Rcpp.html), so all credit to him
//' @param x the cutoff points to be found
//' @param breaks the vector to look for.
//' @export
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
IntegerVector findInterval2(NumericVector x, NumericVector breaks) {
  IntegerVector out(x.size());
  NumericVector::iterator it, pos;
  IntegerVector::iterator out_it;
  NumericVector::iterator brbeg = breaks.begin();
  NumericVector::iterator brend = breaks.end();
  for(it = x.begin(), out_it = out.begin(); it != x.end(); 
  ++it, ++out_it) {
    pos = std::upper_bound(brbeg, brend ,*it);
    *out_it = std::distance(brbeg, pos);
  }
  
  return out;
}

//' Find interval implementation in armadillo
//' 
//' This takes the fast findInterval2 implementation from Hadley Wickham (http://adv-r.had.co.nz/Rcpp.html) and just replaces everything with RcppArmadillo objects. Since smoothQTL is implemented in armadillo this saves some conversions
//' @param x the cutoff points to be found
//' @param breaks the vector to look for.
//' @export
// [[Rcpp::export]]
arma::vec findIntervalArmadillo(arma::vec& x, arma::vec& breaks) {
  arma::vec out(x.size());
  arma::vec::iterator it, pos;
  arma::vec::iterator out_it;
  arma::vec::iterator brbeg = breaks.begin();
  arma::vec::iterator brend = breaks.end();
  for(it = x.begin(), out_it = out.begin(); it != x.end(); 
  ++it, ++out_it) {
    pos = std::upper_bound(brbeg, brend ,*it);
    *out_it = std::distance(brbeg, pos);
  }
  
  return out;
}


arma::vec tricube(arma::vec D)
{
  arma::vec denom = pow((1-pow(D,3)),3);
  arma::vec tricube = denom / sum(denom);
  return tricube;
}

//' Calculate smoothed G statistic
//' 
//' This function implements the smoothing procedure outlined in Magwene et al 2011.
//'
//' @param G the G statistic
//' @param Map the genetic map (see gqtl function)
//' @param W the window size for the smoothing (see gqtl function)
//' @export
// [[Rcpp::export]]
NumericVector smoothQTL(arma::vec G,arma::vec Map,NumericVector W) 
{
  if(G.size()!=Map.size())
  {
    throw std::invalid_argument("G and Map must be of sampe length");
  }
  if((max(Map)-min(Map))< (int)W[0]*2)
  {
    throw std::invalid_argument("For safety, the window W can't be less than 2x the span of the map"); 
  }
  arma::vec::fixed<1> cutoff;
  cutoff[0] = (double)W[0]/2;
  arma::vec MapDiff = Map-Map[0];
  arma::vec firstEdge = findIntervalArmadillo(cutoff,MapDiff);
  arma::vec MapDiffSecond = Map-Map[Map.size()-1];
  arma::vec::fixed<1> cutoffNeg;
  cutoffNeg = -1*cutoff[0];
  arma::vec secondEdge = findIntervalArmadillo(cutoffNeg,MapDiffSecond);
  arma::vec Gpadded(firstEdge[0]+G.size() + Map.size() - secondEdge[0]);
  arma::vec mapPadded(Gpadded.size());
  for(int i=0;i<firstEdge[0];i++)
  {
    Gpadded(i) = G(firstEdge[0]-i); 
    mapPadded(i) = Map[0]-Map(firstEdge[0]-i);
  }
  
  for(int i=0;i<G.size();i++)
  {
    Gpadded(i+firstEdge[0]) = G(i); 
    mapPadded(i+firstEdge[0]) = Map(i);
  }
  for(int i=0;i<G.size()-secondEdge[0];i++)
  {
    Gpadded(i+firstEdge[0]+G.size()) = G(G.size()-i-2);
    mapPadded(i+firstEdge[0]+G.size()) = Map(G.size()-1)+ Map(G.size()-1)-Map(G.size()-i-2);
  }
  arma::vec Gsmoothed(G.size());
  int k = 0;
  for(int i=firstEdge[0];i<firstEdge[0]+G.size();i++)
  {
    arma::vec curMapDiff = mapPadded-mapPadded(i);
    int startInd = findIntervalArmadillo(cutoffNeg,curMapDiff)[0];
    int endInd = findIntervalArmadillo(cutoff,curMapDiff)[0];
    arma::uvec indices = arma::regspace<arma::uvec>(startInd,endInd);
    arma::vec distS = mapPadded.elem(indices) - mapPadded(i);
    arma::vec distSnorm = arma::abs(distS/arma::max(distS));
    arma::rowvec GpaddedTrans = Gpadded.elem(indices).t();
    Gsmoothed(k) = dot(tricube(distSnorm),GpaddedTrans);
    k+=1;
  }
  return wrap(Gsmoothed);
}
