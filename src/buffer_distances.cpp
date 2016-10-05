#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;


//' Calculating distances between points subject to a buffer-zone edge correction.
//'
//' Calculates pairwise distances between points subject to a buffer-zone edge correction.
//'
//' @inheritParams fit.ns
//'
//' @export
// [[Rcpp::export]]
NumericVector buffer_distances(const NumericMatrix& points, const NumericMatrix& lims, const double& R){
  int n_dims = points.ncol();
  int n_points = points.nrow();
  int i, j, k, m;
  // Working out internal points.
  LogicalVector is_internal(n_points);
  int n_internal = 0;
  for (i = 0; i < n_points; i++){
    is_internal(i) = true;
    for (j = 0; j < n_dims; j++){
      if (points(i, j) < (lims(j, 0) + R) || points(i, j) > (lims(j, 1) - R)){
	is_internal(i) = false;
      }
    }
    if (is_internal(i)){
      n_internal++;
    }
  }
  // Calculating distances from each internal point to all others.
  NumericVector out(n_points*n_internal - n_internal);
  double dist_sq, diff;
  m = 0;
  for (i = 0; i < n_points; i++){
    // Only bother continuing if ith point is internal.
    if (is_internal(i)){
      for (j = 0; j < n_points; j++){
	// Skip if i == j.
	if (i != j){
	  // Initialising variable for sum of squared differences.
	  dist_sq = 0;
	  for (k = 0; k < n_dims; k++){
	    // Calculating squared difference on the kth dimension.
	    dist_sq += pow(points(i, k) - points(j, k), 2);
	  }
	  // Saving distance.
	  out(m) = sqrt(dist_sq);
	  m++;
	}
      }
    }
  }
  return out;
}
