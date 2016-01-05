#include <Rcpp.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix Schulze_C(IntegerMatrix Pairs) {
  int nrow = Pairs.nrow();
  IntegerMatrix Schulze(nrow, nrow);
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < nrow; j++) {
      if (i != j) {
        if (Pairs(i, j) > Pairs(j, i)) {
          Schulze(i, j) = Pairs(i, j);
        } else {
          Schulze(i, j) = 0;
        }
      }
    }
  }
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < nrow; j++) {
      if (i != j) {
        for (int k = 0; k < nrow; k++) {
          if ((i != k) && (j != k)) {
            Schulze(j, k) = (std::max)(Schulze(j, k), (std::min)(Schulze(j, i), Schulze(i, k)));
          }
        }
      } else {
        if ((i = j)) {
          Schulze(i, j) = 0;
        }
      }
    }
  }
  return(Schulze);
}