#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector rwC (double sigma, double x0, int N) {
    NumericVector u = runif(N);
    NumericVector x(N);
    x[0] = x0;
    for (int i=1; i <= N-1; i++) {
        NumericVector y = rnorm(1,x[i-1],sigma);
        if (u[i] <= exp(abs(x[i-1])-abs(y[0]))) {
            x[i] = y[0];
        }
        if (u[i] > exp(abs(x[i-1])-abs(y[0]))) {
            x[i] = x[i-1];
        }
    }
    return(x);
}