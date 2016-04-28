#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector inRegion(DataFrame df, NumericVector x, NumericVector y) {
  NumericVector x1 = df["x1"];
  NumericVector x2 = df["x2"];
  NumericVector y1 = df["y1"];
  NumericVector y2 = df["y2"];
  IntegerVector Loc2 = df["Loc2"];
  IntegerVector LocRes(x.size(), NA_INTEGER);

  for (int i = 0; i < x.size(); i++) {
    for (int j = 0; j < x1.size(); j++) {
      if ((x[i] >= x1[j]) && (x[i] <= x2[j]) && (y[i] >= y1[j]) && (y[i] <= y2[j])) {
        LocRes[i] = Loc2[j];
        break;
      } else {}
    }
  }

  return LocRes;
}
