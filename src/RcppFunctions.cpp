#include <RcppArmadillo.h>
#include <SpatialIndex.h>
#include "LAS.h"
using namespace Rcpp;

// [[Rcpp::export(rng = false)]]
IntegerVector C_count_in_disc(NumericVector X, NumericVector Y, NumericVector x, NumericVector y, double radius, int ncpu)
{
  unsigned int n = x.length();
  IntegerVector output(n);

  lidR::GridPartition tree(X,Y);

  #pragma omp parallel for num_threads(ncpu)
  for(unsigned int i = 0 ; i < n ; i++)
  {
    lidR::Circle disc(x[i], y[i], radius);
    std::vector<lidR::PointXYZ> pts;
    tree.lookup(disc, pts);

    #pragma omp critical
    {
      output[i] = pts.size();
    }
  }

  return output;
}

// [[Rcpp::export(rng = false)]]
IntegerVector C_li2012_auto(S4 las, double dt1, double dt2, double Zu, double th_tree, double radius)
{
  LAS pt(las);
  return pt.segment_trees_auto(dt1, dt2, Zu, th_tree, radius);
}
