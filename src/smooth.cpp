#include "RcppArmadillo.h"
#include <algorithm>


using namespace Rcpp;

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

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
arma::colvec tricube(arma::colvec D)
{
  arma::vec denom = pow((1-pow(D,3)),3);
  arma::vec tricube = denom / sum(denom);
  return tricube;
}
// [[Rcpp::export]]
NumericVector smoothQTL(arma::colvec G,arma::colvec Map,CharacterVector Kern,NumericVector W) 
  {
    if(G.size()!=Map.size())
    {
      throw std::invalid_argument("G and Map must be of sampe length");
    }
    if((max(Map)-min(Map))< (int)W[0]*2)
    {
      throw std::invalid_argument("For safety, the window W can't be less than 2x the span of the map"); 
    }
    NumericVector cutoff = wrap((double)W[0]/2);
    Rcpp::IntegerVector firstEdge = findInterval2(cutoff,wrap((Map-Map[1])));
    Rcpp::IntegerVector secondEdge = findInterval2(wrap(-1 * (double)cutoff[0]),wrap((Map-Map[Map.size()-1])));
    
    arma::colvec Gpadded(firstEdge[0]+G.size() + Map.size() - secondEdge[0]);
    arma::colvec mapPadded(Gpadded.size());
    
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
    
    arma::colvec Gsmoothed(G.size());
    int k = 0;
    for(int i=firstEdge[0];i<firstEdge[0]+G.size();i++)
    {
      int startInd = findInterval2(wrap(-1 *(double)cutoff[0]),wrap(mapPadded-mapPadded(i)))[0];
      int endInd = findInterval2(wrap(cutoff[0]),wrap(mapPadded-mapPadded(i)))[0];
      arma::uvec indices = arma::regspace<arma::uvec>(startInd,endInd);
      arma::colvec distS = mapPadded.elem(indices) - mapPadded(i);
      arma::colvec distSnorm = arma::abs(distS/arma::max(distS));
      arma::rowvec GpaddedTrans = Gpadded.elem(indices).t();
      Gsmoothed(k) = dot(tricube(distSnorm),GpaddedTrans);
      k+=1;
    }
    return wrap(Gsmoothed);
  }
