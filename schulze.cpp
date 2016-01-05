#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>
using namespace Rcpp;

// [[Rcpp::export]]
List Schulze_M(IntegerMatrix Ranks ) {
  int n = 1000000000;
  int nrow = Ranks.nrow();
  int ncol = Ranks.ncol();
  IntegerMatrix vec(1,ncol);
  IntegerMatrix mtx(ncol, ncol);
  IntegerMatrix result(ncol, ncol);
  std::fill(vec.begin(), vec.end(), 1);
  Timer timer;
  
    if (nrow == 1) {
    vec = Ranks;
    }  else { 
      for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol-1; j++) { 
          for (int k = j+1; k < ncol; k++)  {
            if(Ranks(i,j)<Ranks(i,k)) mtx(j,k)++;
            if(Ranks(i,k)<Ranks(i,j)) mtx(k,j)++;
          } 
        } 
      }
    }
    for (int i = 0; i < ncol; i++) {
      for (int j = 0; j < ncol; j++) {
         if(i!=j) {
           if (mtx(i,j) > mtx(j,i)) {
             result(i,j) = mtx(i,j);
           } else {
             result(i,j) = 0;
           }
         }}}
    for (int i = 0; i < ncol; i++) {
      for (int j = 0; j < ncol; j++) {
        if (i != j) {
          for (int k = 0; k < ncol; k++) {
            if ((i != k) && (j != k)) {
              result(j, k) = (std::max)(result(j, k), (std::min)(result(j, i), result(i, k)));
            }
          }
        } else {
          if ((i = j)) {
            result(i, j) = 0;
          }
        }
      }
    }
  
  for (int i = 0; i < ncol-1; i++) {
    for (int j = i+1; j < ncol; j++) {            
      if(result(i,j)>result(j,i)) {
        vec(0,j)++;
      } else {
        if(result(i,j)<result(j,i)) vec(0,i)++;
      }
    }}
  Rcpp::colnames(vec) = Rcpp::colnames(Ranks);
  timer.step("elapsed");
  NumericVector res(timer);
  for (int i=0; i<res.size(); i++) {
    res[i] = res[i] / n;
  }

 return Rcpp::List::create(Rcpp::Named("Consensus") = vec,
                           Rcpp::Named("Schulze") = result,
                           Rcpp::Named("Eltime")=res);   
}

