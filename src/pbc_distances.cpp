#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector pbc_distances(const NumericMatrix& points, const NumericMatrix& lims){
  int n_dims = points.ncol();
  int n_points = points.nrow();
  NumericVector out((pow(n_points, 2) - n_points)/2);
  double diff, dist_sq;
  int i, j, k, m;
  // Difference between limits for each dimension.
  NumericVector lim_diffs(n_dims);
  for (i = 0; i < n_dims; i++){
    lim_diffs(i) = lims(i, 1) - lims(i, 0);
  }
  m = 0;
  for (i = 0; i < (n_points - 1); i++){
    for (j = i + 1; j < n_points; j++){
      // Initialising variable for sum of squared differences.
      dist_sq = 0;
      for (k = 0; k < n_dims; k++){
	// Calculating difference on the kth dimension.
	diff = abs(points(i, k) - points(j, k));
	// Adjustment for periodic boundary constraints.
	if (diff > 0.5*lim_diffs(k)){
	  diff = lim_diffs(k) - diff;
	}
	// Squared difference on the kth dimension.
	dist_sq += pow(diff, 2);
      }
      // Saving distances (duplication for consistency with method of
      // Tanaka et al).
      out(m) = sqrt(dist_sq);
      m++;
    }
  }
  return out;
}
