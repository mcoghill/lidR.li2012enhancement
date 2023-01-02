#ifndef LAS_H
#define LAS_H

#include <RcppArmadillo.h>
#define NDEBUG 1
using namespace Rcpp;

class LAS
{
public:
  S4 las;
  void new_filter(LogicalVector b);
  void filter_local_maxima(NumericVector ws, double min_height, bool circular);
  NumericVector X;
  NumericVector Y;
  NumericVector Z;
  NumericVector T;
  IntegerVector I;
  unsigned int ncpu;
  unsigned int npoints;
  std::vector<bool> filter;
  std::vector<bool> skip;

public:
  LAS(S4 las, int npcu = 1);
  IntegerVector segment_trees_auto(double dt1, double dt2, NumericVector R, double Zu, double th_tree, double radius);

private:
  unsigned int sensor;
  enum TYPES {UKN = 0, ALS = 1, TLS = 2, UAV = 3, DAP = 4, MLS = 5};
};

#endif //LAS_H