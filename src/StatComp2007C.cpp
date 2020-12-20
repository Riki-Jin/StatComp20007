#include <Rcpp.h>
using namespace Rcpp;

//' @title A random walk Metropolis sampler using Rcpp
//' @description A random walk Metropolis sampler using Rcpp
//' @param N the number of samples
//' @param sigma the standard deviation used for the proposal distribution
//' @param x0 the initial value
//' @return a random sample of size N
//' @examples
//' \dontrun{
//' sigma <- 2
//' x0 <- 5
//' N <- 2000
//' rnC <- rwC(sigma, x0, N)
//' }
//' @export
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