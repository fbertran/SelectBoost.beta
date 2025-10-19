#include <RcppArmadillo.h>
#include <Rmath.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List beta_work_cpp(const NumericVector& y, const NumericVector& eta, const double phi) {
  int n = y.size();
  NumericVector mu(n), ystar(n), W(n), z(n);

  #ifdef _OPENMP
  #pragma omp parallel for if(n > 5000)
  #endif
  for (int i = 0; i < n; ++i) {
    double et = eta[i];
    double yi = y[i];
    double mui = 1.0 / (1.0 + std::exp(-et));
    if (mui < 1e-6) mui = 1e-6;
    if (mui > 1.0 - 1e-6) mui = 1.0 - 1e-6;
    mu[i] = mui;
    ystar[i] = std::log(yi) - std::log1p(-yi);
  }

  #ifdef _OPENMP
  #pragma omp parallel for if(n > 2000)
  #endif
  for (int i = 0; i < n; ++i) {
    double mui = mu[i];
    double ai = mui * phi;
    double bi = (1.0 - mui) * phi;
    double trig_a = R::trigamma(ai);
    double trig_b = R::trigamma(bi);
    double mustari = R::digamma(ai) - R::digamma(bi);
    double w = std::pow(mui * (1.0 - mui), 2.0) * (phi * phi) * (trig_a + trig_b);
    if (w < 1e-10) w = 1e-10;
    double u = (mui * (1.0 - mui)) * phi * (ystar[i] - mustari);
    W[i] = w;
    z[i] = eta[i] + u / w;
  }

  return List::create(_["mu"] = mu, _["W"] = W, _["z"] = z);
}

// [[Rcpp::export]]
double beta_phi_update_cpp(double phi, const NumericVector& mu, const NumericVector& y, int maxit) {
  double ph = phi;
  int n = y.size();
  for (int it = 0; it < maxit; ++it) {
    double lp = 0.0, lpp = 0.0;
    for (int i = 0; i < n; ++i) {
      double mui = mu[i];
      double ai = mui * ph;
      double bi = (1.0 - mui) * ph;
      lp  += R::digamma(ph) - mui * R::digamma(ai) - (1.0 - mui) * R::digamma(bi)
             + mui * std::log(y[i]) + (1.0 - mui) * std::log1p(-y[i]);
      lpp += R::trigamma(ph) - (mui*mui) * R::trigamma(ai) - ((1.0 - mui)*(1.0 - mui)) * R::trigamma(bi);
    }
    double denom = (std::abs(lpp) < 1e-12) ? (lpp < 0 ? -1e-12 : 1e-12) : lpp;
    double step = lp / denom;
    ph = ph - step;
    if (!R_finite(ph) || ph <= 1e-6) { ph = phi; break; }
    if (std::abs(step) < 1e-8) break;
  }
  if (!R_finite(ph) || ph <= 0) ph = phi;
  return ph;
}
