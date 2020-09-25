#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix euc_distances(const NumericVector& x1, const NumericVector& y1, const NumericVector& x2, const NumericVector& y2){
  int n1 = x1.size();
  int n2 = x2.size();
  NumericMatrix out (n1, n2);
  for (int j = 0; j < n2; j++){
    for (int i = 0; i < n1; i++){
      out(i, j) = pow(pow(x1(i) - x2(j), 2) + pow(y1(i) - y2(j), 2), 0.5);
    }
  }
  return out;
}
