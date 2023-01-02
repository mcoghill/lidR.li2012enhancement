#include "LAS.h"
#include "Progress.h"
#include "myomp.h"
#include <SpatialIndex.h>

using namespace lidR;

LAS::LAS(S4 las, int ncpu)
{
  Rcpp::List index = las.slot("index");
  this->sensor = index["sensor"];

  this->las = las;

  DataFrame data = as<DataFrame>(las.slot("data"));
  this->X = data["X"];
  this->Y = data["Y"];
  this->Z = data["Z"];

  if (data.containsElementNamed("Intensity"))
    this->I = data["Intensity"];

  if (data.containsElementNamed("gpstime"))
    this->T = data["gpstime"];

  this->npoints = X.size();
  this->ncpu = ncpu;
  this->filter.resize(npoints);
  std::fill(filter.begin(), filter.end(), false);
  this->skip.resize(npoints);
  std::fill(skip.begin(), skip.end(), false);
}

void LAS::filter_local_maxima(NumericVector ws, double min_height, bool circular)
{
  bool abort = false;
  bool vws = ws.length() > 1;

  SpatialIndex tree(las);
  Progress pb(npoints, "Local maximum filter: ");

#pragma omp parallel for num_threads(ncpu)
  for (unsigned int i = 0 ; i < npoints ; i++)
  {
    if (abort) continue;
    if (pb.check_interrupt()) abort = true;
    pb.increment();

    double hws = (vws) ? ws[i]/2 : ws[0]/2;

    if (Z[i] < min_height)
      continue;

    // Get the points within a windows centered on the current point
    std::vector<PointXYZ> pts;
    if(!circular)
    {
      Rectangle rect(X[i]-hws, X[i]+hws, Y[i]-hws, Y[i]+hws);
      tree.lookup(rect, pts);
    }
    else
    {
      Circle circ(X[i], Y[i], hws);
      tree.lookup(circ, pts);
    }

    // Initialize the highest point using the central point
    PointXYZ p(X[i], Y[i], Z[i], i);
    double zmax = Z[i];
    bool is_lm = true;

    // Search if one is higher
    for (auto pt : pts)
    {
      double z = pt.z;

      // Found one higher, it is not a LM
      if(z > zmax)
      {
        is_lm = false;
        break;
      }

      // Found one equal. If this one was already tagged LM we can't have two lm
      // The first tagged has the precedence
      if (z == zmax && filter[pt.id])
      {
        is_lm = false;
        break;
      }
    }

#pragma omp critical
{
  filter[i] = is_lm;
}
  }

  if (abort) throw Rcpp::internal::InterruptedException();

  return;
}

IntegerVector LAS::segment_trees_auto(double dt1, double dt2, NumericVector R, double Zu, double th_tree, double radius)
{
  double xmin = min(X);
  double ymin = min(Y);

  unsigned int ni = npoints;            // Number of points
  unsigned int n  = ni;                 // Number of remaining points
  unsigned int k  = 1;                  // Current tree ID

  // The ID of each point (returned object)
  IntegerVector idtree(ni);
  std::fill(idtree.begin(), idtree.end(), NA_INTEGER);

  // Square distance to speed up computation (dont need sqrt)
  radius = radius * radius;
  dt1 = dt1 * dt1;
  dt2 = dt2 * dt2;

  /* =====================
   * LI ET AL ALGORITHHM *
   ======================*/

  // Li, W., Guo, Q., Jakubowski, M. K., & Kelly, M. (2012). A New Method for Segmenting Individual
  // Trees from the Lidar Point Cloud. Photogrammetric Engineering & Remote Sensing, 78(1), 75–84.
  // https://doi.org/10.14358/PERS.78.1.75

  // Find if a point is a local maxima within an R windows
  LogicalVector is_lm;
  if (R.length() > 1 || R[0] > 0)
  {
    filter_local_maxima(R, 0, true);
    is_lm = Rcpp::wrap(filter);
  }
  else
  {
    is_lm = LogicalVector(ni);
    std::fill(is_lm.begin(), is_lm.end(), true);
  }

  // A progress bar and abort options
  Progress p(ni, "Tree segmentation: ");

  // U the points to be segmented (see Li et al. page 78)
  std::vector<PointXYZ*> U(ni);
  for (unsigned int i = 0 ; i < ni ; ++i)
    U[i] = new PointXYZ(X[i], Y[i], Z[i], i);

  // N and P groups (see Li et al. page 78)
  std::vector<PointXYZ*> P,N;
  P.reserve(100);
  N.reserve(100);

  // A dummy point out of the dataset (see Li et al. page 79)
  PointXYZ* dummy = new PointXYZ(xmin-100,ymin-100,0,-1);

  // Z-sort the point cloud U
  std::sort(U.begin(), U.end(), ZSort<PointXYZ>());

  while(n > 0)
  {
    PointXYZ* u = U[0];
    std::vector<bool> inN(n);

    // Stop the algo is the highest point u, which is the target tree top, is below a threshold
    // Addition from original algo to limit over segmentaton
    if (u->z < th_tree)
    {
      p.update(ni);
    }
    else
    {
      if (p.check_interrupt())
      {
        for (unsigned int i = 0 ; i < U.size() ; i++) delete U[i]; // # nocov
        delete dummy; // # nocov
        p.exit(); // # nocov
      }

      p.update(ni-n);

      // Initial step no point in P or N
      P.clear();
      N.clear();

      // element 0 is the current highest point and is in P (target tree)
      P.push_back(u);
      idtree[u->id] = k;

      // Add the dummy point in N
      N.push_back(dummy);

      // Compute the distance between the current point u and all the &other points of U
      // This is not in the original algo. This is an optimisation to reduce the computation
      // time (see line 136).
      std::vector<double> d = sqdistance(U, *u);

      // Loop over each point of U (but the global maximum that is alredy in P)
      for (unsigned int i = 1 ; i < n ; ++i)
      {
        u = U[i];

        // If d > radius this point u is far and thus it is not in the current segmented tree
        // We don't need to apply the li et al. rules. This speed up a lot the computations
        if(d[i] > radius)
        {
          inN[i] = true;
        }
        // If d <= radius classify the point u based on Li et al. rules
        else
        {
          std::vector<double> dP = sqdistance(P, *u);
          std::vector<double> dN = sqdistance(N, *u);

          double dmin1 = *std::min_element(dP.begin(), dP.end());
          double dmin2 = *std::min_element(dN.begin(), dN.end());
          double dt    = (u->z > Zu) ? dt2 : dt1;

          if(is_lm[u->id]) // if u is a local maximum
          {
            if (dmin1 > dt || (dmin1 < dt && dmin1 > dmin2))
            {
              inN[i] = true;
              N.push_back(u);
            }
            else
            {
              P.push_back(u);
              idtree[u->id] = k;
            }
          }
          else // if u is not a local maximum
          {
            if (dmin1 <= dmin2)
            {
              P.push_back(u);
              idtree[u->id] = k;
            }
            else
            {
              inN[i] = true;
              N.push_back(u);
            }
          }
        }
      }
    }

    // Keep the point in N and redo the loop with remining points
    std::vector<PointXYZ*> temp;
    temp.reserve(N.size()-1);

    for (unsigned int i = 0 ; i < n ; i++)
    {
      if(inN[i])
        temp.push_back(U[i]);
      else
        delete U[i];
    }

    U.swap(temp);
    n = U.size();
    k++;                        // Increase current tree id
  }

  delete dummy;

  return idtree;
}