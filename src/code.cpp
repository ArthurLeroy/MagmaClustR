#include <Rcpp.h>
using namespace Rcpp;

// Function used to speed-up computations of the stationary kernel
// [[Rcpp::export]]
NumericMatrix cpp_dist(NumericMatrix m1, NumericMatrix m2) {
  int nrow1 = m1.nrow();
  int nrow2 = m2.nrow();
  int ncol = m1.ncol();

  if (ncol != m2.ncol()) {
    throw std::runtime_error("Incompatible number of dimensions");
  }

NumericMatrix out(nrow1, nrow2);

  for (int r1 = 0; r1 < nrow1; r1++) {
    for (int r2 = 0; r2 < nrow2; r2++) {
      double total = 0;
      for (int c12 = 0; c12 < ncol; c12++) {
        total += pow(m1(r1, c12) - m2(r2, c12), 2);
      }
      out(r1, r2) = total;
    }
  }
  return out;
}

// Function used to speed-up computations of the periodic kernel
// [[Rcpp::export]]
NumericMatrix cpp_perio(
      NumericMatrix m1,
      NumericMatrix m2,
      double period) {
  int nrow1 = m1.nrow();
  int nrow2 = m2.nrow();
  int ncol = m1.ncol();

  NumericMatrix out(nrow1, nrow2);

  for (int r1 = 0; r1 < nrow1; r1++) {
    for (int r2 = 0; r2 < nrow2; r2++) {
      double total = 0;
      for (int c12 = 0; c12 < ncol; c12++) {
        total += pow(sin(M_PI/exp(period)*abs(m1(r1, c12)-m2(r2, c12))),2);
      }
      out(r1, r2) = total;
    }
  }
  return out;
}

// Function used to speed-up computations of derivatives of the periodic kernel
// [[Rcpp::export]]
NumericMatrix cpp_perio_deriv(
        NumericMatrix m1,
        NumericMatrix m2,
        double period) {
    int nrow1 = m1.nrow();
    int nrow2 = m2.nrow();
    int ncol = m1.ncol();

    NumericMatrix out(nrow1, nrow2);

    for (int r1 = 0; r1 < nrow1; r1++) {
      for (int r2 = 0; r2 < nrow2; r2++) {
        double total = 0;
        for (int c12 = 0; c12 < ncol; c12++) {
          double angle = M_PI/exp(period) * abs(m1(r1, c12) - m2(r2, c12));
          total += 2 * sin(angle) * cos(angle) * angle;
        }
        out(r1, r2) = total;
      }
    }
    return out;
}

// Function used to speed-up computations of the stationary kernel
// [[Rcpp::export]]
NumericMatrix cpp_prod(NumericMatrix m1, NumericMatrix m2) {
        int nrow1 = m1.nrow();
        int nrow2 = m2.nrow();
        int ncol = m1.ncol();

        if (ncol != m2.ncol()) {
          throw std::runtime_error("Incompatible number of dimensions");
        }

        NumericMatrix out(nrow1, nrow2);

        for (int r1 = 0; r1 < nrow1; r1++) {
          for (int r2 = 0; r2 < nrow2; r2++) {
            double total = 0;
            for (int c12 = 0; c12 < ncol; c12++) {
              total += m1(r1, c12) * m2(r2, c12);
            }
            out(r1, r2) = total;
          }
        }
        return out;
}

// Function used to speed-up computations of derivatives of the periodic kernel
// [[Rcpp::export]]
NumericMatrix cpp_noise(
            NumericMatrix m1,
            NumericMatrix m2,
            double noise) {
          int nrow1 = m1.nrow();
          int nrow2 = m2.nrow();
          int ncol = m1.ncol();

          NumericMatrix out(nrow1, nrow2);

          for (int r1 = 0; r1 < nrow1; r1++) {
            for (int r2 = 0; r2 < nrow2; r2++) {
              double diff = 0;
              double value = 0;

              for (int c12 = 0; c12 < ncol; c12++) {
                diff += abs(m1(r1, c12) - m2(r2, c12));
              }

              if(diff < 0.000001){
                value = exp(noise);
              }
              out(r1, r2) = value;
            }
          }
          return out;
}
