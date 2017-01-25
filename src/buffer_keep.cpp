#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
LogicalMatrix buffer_keep(const NumericMatrix& points, const NumericMatrix& lims, const double& R){
  int n_dims = points.ncol();
  int n_points = points.nrow();
  int i, j;
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
  LogicalMatrix out(n_points, n_points);
  for (i = 0; i < n_points; i++){
    for (j = 0; j < n_points; j++){
      if (i == j || !is_internal(i)){
	out(i, j) = false;
      } else {
	out(i, j) = true;
      }
    }
  }
  return out;
}
