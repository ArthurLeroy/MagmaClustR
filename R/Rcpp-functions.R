## Function used to speed-up computations of the stationary kernel
Rcpp::cppFunction('NumericMatrix cpp_dist(NumericMatrix m1, NumericMatrix m2) {
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
}')

## Function used to speed-up computations of the periodic kernel
Rcpp::cppFunction('Rcpp::NumericMatrix cpp_perio(
Rcpp::NumericMatrix m1,
Rcpp::NumericMatrix m2,
double period) {
  int nrow1 = m1.nrow();
  int nrow2 = m2.nrow();
  int ncol = m1.ncol();

  NumericMatrix out(nrow1, nrow2);

  for (int r1 = 0; r1 < nrow1; r1++) {
    for (int r2 = 0; r2 < nrow2; r2++) {
      double total = 0;
      for (int c12 = 0; c12 < ncol; c12++) {
        total += pow(sin(M_PI/std::exp(period)*abs(m1(r1, c12)-m2(r2, c12))),2);
      }
      out(r1, r2) = total;
    }
  }
  return out;
}')

## Function used to speed-up computations of derivatives of the periodic kernel
Rcpp::cppFunction('Rcpp::NumericMatrix cpp_perio_deriv(
Rcpp::NumericMatrix m1,
Rcpp::NumericMatrix m2,
double period) {
  int nrow1 = m1.nrow();
  int nrow2 = m2.nrow();
  int ncol = m1.ncol();

  NumericMatrix out(nrow1, nrow2);

  for (int r1 = 0; r1 < nrow1; r1++) {
    for (int r2 = 0; r2 < nrow2; r2++) {
      double total = 0;
      for (int c12 = 0; c12 < ncol; c12++) {
        double angle = M_PI/std::exp(period) * abs(m1(r1, c12) - m2(r2, c12));
        total += 2 * sin(angle) * cos(angle) * angle;
      }
      out(r1, r2) = total;
    }
  }
  return out;
}')

## Function used to speed-up computations of the stationary kernel
Rcpp::cppFunction('NumericMatrix cpp_prod(
Rcpp::NumericMatrix m1,
Rcpp::NumericMatrix m2) {
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
        total += m1(r1, c12)  * m2(r2, c12);
      }
      out(r1, r2) = total;
    }
  }
  return out;
}')

