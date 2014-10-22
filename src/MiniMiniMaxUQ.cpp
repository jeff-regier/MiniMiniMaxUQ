#include <Rcpp.h>
#include <iostream>
#include <bitset>
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
double lp_distance(NumericVector a, NumericVector b, double p) {
  double ret = 0.;
	int m = a.size();
  if (p < 0) {  // negative values indicate the infinity norm
    for(int i = 0; i < m; ++i) {
  		double diff = abs(a[i] - b[i]);
      if (diff > ret)
          ret = diff;
  	}
  }
  else {
    for(int i = 0; i < m; ++i)
  		ret += pow(abs(a[i] - b[i]), p);
    ret = pow(ret, 1. / p);
  }
	return ret;
}


// [[Rcpp::export]]
double find_K_hat(NumericMatrix xs, NumericVector ys, double p_norm=-1., double min_dist=1e-4) {
  double K_hat = 0.;
  int n = ys.size();
  for (int i = 0; i < (n - 1); ++i) {
    for (int j = i + 1; j < n; ++j) {
      double dist = lp_distance(xs(i, _), xs(j, _), p_norm);
      if (dist < min_dist)
        continue;
      double abs_slope = abs(ys[i] - ys[j]) / dist;
      if (abs_slope > K_hat)
        K_hat = abs_slope;
    }
  }
  return K_hat;
}

// [[Rcpp::export]]
NumericVector f_deviations(NumericMatrix xs, double K_hat, NumericVector x) {
  int n = xs.nrow();
  NumericVector ret(n);
  for (int i = 0; i < n; ++i) {
    ret[i] = K_hat * lp_distance(x, xs(i, _), -42.);
  }
  return ret;
}

// [[Rcpp::export]]
double do_pointwise_uncertainty(NumericVector ys, NumericVector deviations) {
  double envelope_ub = min(ys + deviations);
  double envelope_lb = max(ys - deviations);
  return (envelope_ub - envelope_lb) / 2.;
}


// [[Rcpp::export]]
double pointwise_uncertainty(NumericMatrix xs, NumericVector ys, double K_hat, NumericVector x) {
  NumericVector deviations = f_deviations(xs, K_hat, x);
  return do_pointwise_uncertainty(ys, deviations);
}
