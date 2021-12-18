//
//Programmer: Jiayuan Xiao, 12/3/21
//Purpose: 
    

//#include "Sling.h"
#include <Rcpp.h>
#include <numeric> //accumulate
#include <cmath> //sqrt
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double vecSum(NumericVector &vec){  //point iterator
  double sum = 0;
  for (NumericVector::iterator i = vec.begin(); i != vec.end(); i++){
    sum = sum + *i;
  }
  return sum;
}

// [[Rcpp::export]]
double Sum(NumericVector &vec){  //Use R func
  Function f("sum");
  return sum(vec);
}

// [[Rcpp::export]]
double Sum2(NumericVector &vec){  //Rcpp sugar
  return sum(vec);
}

// [[Rcpp::export]]
double Sum3(NumericVector &vec){  //Non point iterator, pass by reference
  int n = vec.size();
  double total = 0;
  for (int i = 0; i < n; i++){
    total += vec[i];
  }
  return total;
}

// [[Rcpp::export]]
double Sum4(NumericVector vec){  //Non point iterator, pass by object
  int n = vec.size();
  double total = 0;
  for (int i = 0; i < n; i++){
    total += vec[i];
  }
  return total;
}

// [[Rcpp::export]]
double Sum5(NumericVector &vec){  //Accumulate
  return accumulate(vec.begin(), vec.end(), 0);
}

// [[Rcpp::export]]
double Inner_Product1(NumericVector &vec1, NumericVector &vec2){  //Faster
  return sum(vec1 * vec2);
}

// [[Rcpp::export]]
double Inner_Product2(NumericVector &vec1, NumericVector &vec2){  //Slower
  return inner_product(vec1.begin(), vec1.end(), vec2.begin(), 0);
}

// [[Rcpp::export]]
NumericMatrix mmult(NumericMatrix &m, NumericMatrix &v)
{
  Environment base("package:base");
  Function mat_Mult = base["%*%"];
  return(mat_Mult(m, v));
}

// [[Rcpp::export]]
NumericMatrix mvmult(NumericMatrix &m, NumericVector &v)
{
  Environment base("package:base");
  Function mv_Mult = base["%*%"];
  return(mv_Mult(m, v));
}

unordered_set<int> convertToSet(NumericVector v)
{
    // Declaring the set
    // using range of vector
    unordered_set<int> s(v.begin(), v.end());
  
    // Return the resultant Set
    return s;
}

// Standard Lasso
double a(NumericMatrix &X, int i){
  return mean(pow(X( _ , i), 2));
}

double z(NumericMatrix &X, NumericVector &y, NumericVector &w,
          int i, int n, double w0){
  double tot = 0;
  for (int j = 0; j < n; j++){
    double yhj = w0 + sum(w * X(j, _ )) - w[i]*X(j, i);
    tot += (y[j] - yhj) * X(j, i);
  }
  return tot/n;
}

bool KKT(NumericMatrix &X, NumericVector &y, NumericVector &w,
         int i, int n, double w0, double l){
  bool res;
  double e = 0.3;
  double g = sum(X( _ , 1) * (y-(mvmult(X, w)+w0))) / (l*n);
  if (w[i] > 0){
    res = ((g>1-e) && (g<1+e));
  }
  else if (w[i] < 0){
    res = ((g>-1-e) && (g<-1+e));
  }
  else{
    res = ((g>-1-e) && (g<1+e));
  }
  return res;
}

// [[Rcpp::export]]
List standardLasso(NumericMatrix &X, NumericVector &y, double l, int max_iter = 100){
  int n = X.nrow();
  int p = X.ncol();
  //double l = 100/(2*n);
  double w0 = mean(y);
  NumericVector w(p);
  w = w+1;
  //DataFrame ws;  //part of result
  NumericVector trainmse;
  NumericVector testmse;
  NumericVector wr(p);
  NumericMatrix v(p);
  
  for (int i = 0; i < p; i++){
    for (int j = 0; j < p; j++){
      v(i, j) = sum(X( _ , i) * X( _ , j));
    }
  }
  int nkkt = p;
  int itr = 0;
  while (nkkt > p/2){
    itr += 1;
    for (int i = 0; i < p; i++){
      NumericVector wr = clone(w);
      double zi = z(X, y, w, i, n, w0);
      double ai = a(X, i);
      if (zi > l){
        w[i] = (zi-l)/ai;
      }
      else if (zi < -l){
        w[i] = (zi+l)/ai;
      }
      else{
        w[i] = 0;
      }
      //NumericVector a = clone(w);
      //ws.push_back(a);    //Super slow
    }
    nkkt = 0;
    for (int i = 0; i < p; i++){
      bool isKKT = KKT(X, y, w, i, n, w0, l);
      if (!isKKT){
        nkkt += 1;
      }
    }
    if (itr >= max_iter){
      
      Rcout << "Break" << endl;  //PRINT
      
      break;
    }
  }
  return List::create(Named("w", w), /*Named("ws", ws),*/ Named("itr", itr));
}


// Sling python ver.
double z_up(int i, int n, NumericVector &w, NumericVector &wr,
            NumericMatrix &v, NumericVector &zr){
  return w[i] - wr[i] + sum(pow(v(i, _ ), 2)) * sum(pow(w-wr, 2)) / n + zr[i];
}

double z_low(int i, int n, NumericVector &w, NumericVector &wr,
             NumericMatrix &v, NumericVector &zr){
  return w[i] - wr[i] - sum(pow(v(i, _ ), 2)) * sum(pow(w-wr, 2)) / n + zr[i];
}

void update(int i, int n, double l, NumericMatrix &X, NumericVector &y,
              NumericVector &w, double w0){
  double zi = z(X, y, w, i, n, w0);
  double ai = a(X, i);
  if (zi > l){
    w[i] = (zi - l) / ai;
  }
  else if (zi < -l){
    w[i] = (zi + l) / ai;
  }
  else{
    w[i] = 0;
  }
}

// [[Rcpp::export]]
List sling(NumericMatrix &X, NumericVector &y, double l, int max_iter = 100){
  int n = X.nrow();
  int p = X.ncol();
  //double l = 100/(2*n);
  double w0 = mean(y);
  NumericVector w(p);
  //DataFrame ws;  //part of result
  NumericVector trainmse;
  NumericVector testmse;
  NumericVector wr(p);
  NumericMatrix v(p);
  for (int i = 0; i < p; i++){
    for (int j = 0; j < p; j++){
      v(i, j) = sum(X( _ , i) * X( _ , j));
    }
  }
  
  NumericVector U;
  bool loop = true;
  int itr = 0;
  int nkkt = p;
  
  while (nkkt > p/2){
    itr += 1;
    loop = false;
    wr = w;
    NumericVector zr;
    NumericMatrix X_T = transpose(X);
    zr = wr + mvmult(X_T, y)/n;
    double err = 1e-4;
    NumericVector prev(p);
    prev = prev + R_PosInf;
    
    while (is_true(all(abs(w-prev)>err))){
      prev = w;  //Need to clone?
      for (int j = 0; j < U.length(); j++){
        int i = U[j];
        double z_low_i = z_low(i, n, w, wr, v, zr);
        double z_up_i = z_up(i, n, w, wr, v, zr);
        
        if ((z_low_i > l) || (z_up_i < -l)){
          update(i, n, l, X, y, w, w0);
          //NumericVector a = clone(w);
          //ws.push_back(a);    //Super slow
        }
      }
    }
    
    wr = w;
    zr = wr + mvmult(X_T, y)/n;
    NumericVector prev2(p);
    prev2 = prev2 + R_PosInf;
    
    while (is_true(all(abs(w-prev2)>err))){
      prev2 = w;  //Need to clone?
      for (int j = 0; j < U.length(); j++){
        int i = U[j];
        double z_low_i = z_low(i, n, w, wr, v, zr);
        double z_up_i = z_up(i, n, w, wr, v, zr);
        
        if ((z_up_i > l) || (z_low_i < -l)){
          update(i, n, l, X, y, w, w0);
          //NumericVector a = clone(w);
          //ws.push_back(a);    //Super slow
        }
        else{
          w[i] = 0;
          //NumericVector a = clone(w);
          //ws.push_back(a);    //Super slow
        }
      }
    }
    
    nkkt = 0;
    for (int i = 0; i < p; i++){
      bool isKKT = KKT(X, y, w, i, n, w0, l);
      if (!isKKT){
        loop = true;
        nkkt += 1;
        if (is_true(all(U != i))){
          U.push_back(i);
        }
      }
    }
    
    if (itr >= max_iter){
      
      Rcout << "Break" << endl;  //PRINT
      
      break;
    }
  }
  return List::create(Named("w", w), /*Named("ws", ws),*/ Named("itr", itr));
}

