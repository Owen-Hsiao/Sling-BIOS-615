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

NumericVector seqStrongRule(NumericMatrix &X, NumericVector &y,
                            NumericVector &w_tilda_lambda,
                            NumericVector &lambdaSeq, int k){
  NumericVector res;
  int n = X.nrow();
  int p = X.ncol();
  if (k > 0){  //k>1 (When k=1, discard all predictors. Tibshirani 2012)
    for (int i = 0; i < p; i++){
      double lhs = abs(sum(X( _ , i) * (y - mvmult(X, w_tilda_lambda)))) / n;
      double rhs = 2 * lambdaSeq[k] - lambdaSeq[k-1];
      if (lhs >= rhs){
        res.push_back(i);
      }
    }
  }
  return res;  //Return the index of predictor survive the Seq Strong check
}

bool isViolateKKT(NumericMatrix &X, NumericVector &y,
                  NumericVector &w_tilda, double lambda, int i){
  bool res;
  int n = X.nrow();
  NumericVector w_tilda_tmp = w_tilda;
  w_tilda_tmp[i] = 0;
  double lhs = abs(sum(X( _ , i) * (y -  mvmult(X, w_tilda)))) / n;
  double rhs = lambda;
  if (lhs <= rhs){
    res = false;  //Did not violate KKT
  }
  else{
    res = true;  //Did violate KKT
  }
  return res;
}

// [[Rcpp::export]]
List sling(NumericMatrix &X, NumericVector &y, NumericVector &lambdas,
           double tol = 0.00001, int max_iter = 10){
  int numLambdas = lambdas.length();
  int n = X.nrow();
  int p = X.ncol();
  NumericMatrix res(p, numLambdas);  //Store results of w_tilde
  Rcout << res.nrow() << " by " << res.ncol() << endl;
  bool isSkViolateKKT = true;
  bool isPViolateKKT = true;
  bool isConverge = false;
  NumericVector P_k;  //P_k-1 is the index of a set of predictors that have nonzero weights                                       //as the solution for tuning parameter λk-1
  for (int k = 0; k < numLambdas; k++){
    Rcout << "loop 1" << endl;
    double lambda = lambdas[k];
    unordered_set<int> U;  //U is the index of a SET of predictors updated in the iterations
    
    //NumericVector P;  //P is the index of a set of all p predictors
    
    if (k != 0){  //k != 1 in Sling line 5
      U = convertToSet(P_k);  //From Numeric Vector to Set
      P_k = {};  //Reset P_k after assign it to P_k-1
    }
    
    NumericVector w(p);
    NumericVector w_tilde(p);
    NumericVector w_r(p);
    NumericVector z(p);
    NumericVector z_r(p);
    Rcout << lambda << endl;
    
    //Initialize w_tilde, if k=1, w_tilde is 0 vector.
    if (k == 1){  //k=2 in equation 11
      w_tilde = res( _ , k - 1);
    }
    else if (k >= 2){  //k=3 in equation 11
      w_tilde = 2 * res( _ , k - 1) - res( _ , k - 2);
    }
    Rcout << k << " " << w_tilde << endl;
    
    int loopCount1 = 0;
    while (isPViolateKKT){
      Rcout << "UU: " << U.size() << endl;
      loopCount1 += 1;
      if (loopCount1 == max_iter){
        Rcout << "Break at loop 1" << endl;
        break;
      }
      Rcout << "loop2" << endl;
      int loopCount2 = 0;
      while (isSkViolateKKT){
        if (loopCount2 == max_iter){
          Rcout << "Break at loop 2" << endl;
          break;
        }
        Rcout << "loop3" << endl;
        w_r = clone(w_tilde);
        //Calculate z_r
        for (unordered_set<int>::iterator iter = U.begin(); iter != U.end(); iter++){//For pi in U
          double tmp = 0;
          for (int j = 0; j < p; j++){
            if (w_r[j] != 0){
              tmp += sum(X( _ , *iter) * X( _ , j)) * w_r[j];
            }
          }
          z_r[*iter] = w_r[*iter] + (sum(X( _ , *iter) * y) - tmp) / n;
        }
        
        double w_norm2 = 0;  //Because w_r = w_tilde
        while (!isConverge){
          NumericVector w_tmp = clone(w_tilde);
          for (unordered_set<int>::iterator iter = U.begin(); iter != U.end(); iter++){//For pi in U
            double v_norm2 = 0;
            NumericVector upperBound(p);
            NumericVector lowerBound(p);
            //Calculate ||v_i||_2
            for (int j = 0; j < p; j++){
              v_norm2 += pow(sum(X( _ , *iter) * X( _ , j)),2);
            }
            v_norm2 = sqrt(v_norm2);
            //Calculate upper and lower bound
            upperBound[*iter] = w_tilde[*iter] - w_r[*iter] + z_r[*iter] + v_norm2 * w_norm2/n;
            lowerBound[*iter] = w_tilde[*iter] - w_r[*iter] + z_r[*iter] - v_norm2 * w_norm2/n;
            if (lowerBound[*iter] > lambda || upperBound[*iter] < -lambda){
              double tmp2 = 0;
              double w_tilde_tmp = w_tilde[*iter];  //w_tilde before update. For w_norm2 update
              //Update w_tilde[i] by equation 5
              for (int j = 0; j < p; j++){
                if (w_tilde[j] != 0){
                  tmp2 += sum(X( _ , *iter) * X( _ , j)) * w_tilde[j];
                }
              }
              z[*iter] = w_tilde[*iter] + (sum(X( _ , *iter) * y) - tmp2) / n;
              if (z[*iter] > 0){
                w_tilde[*iter] = z[*iter] - lambda;
              }
              else{
                w_tilde[*iter] = z[*iter] + lambda;
              }
              //Update w_norm2 by equation 10
              w_norm2 = sqrt(pow(w_norm2, 2) - pow(w_tilde_tmp, 2) + pow(w_tilde[*iter], 2));
            }
          }
          Rcout << "Lhs " << sum((w_tmp - w_tilde) * (w_tmp - w_tilde)) << endl;
          if (sum((w_tmp - w_tilde) * (w_tmp - w_tilde)) < tol){  //Need to adjusted!!!
            isConverge = true;
          }
        }
        
        isConverge = false;  //Reset the flag for next while loop
        w_r = clone(w_tilde);
        //Calculate z_r
        for (unordered_set<int>::iterator iter = U.begin(); iter != U.end(); iter++){//For pi in U
          double tmp = 0;
          for (int j = 0; j < p; j++){
            if (w_r[j] != 0){
              tmp += sum(X( _ , *iter) * X( _ , j)) * w_r[j];
            }
          }
          z_r[*iter] = w_r[*iter] + (sum(X( _ , *iter) * y) - tmp) / n;
        }
        
        w_norm2 = 0;  //Reset w_norm2 since w_r = w_tilde again
        while (!isConverge){
          NumericVector w_tmp = w_tilde;
          for (unordered_set<int>::iterator iter = U.begin(); iter != U.end(); iter++){//For pi in U
            double v_norm2 = 0;
            NumericVector upperBound(p);
            NumericVector lowerBound(p);
            //Calculate ||v_i||_2
            for (int j = 0; j < p; j++){
              v_norm2 += pow(sum(X( _ , *iter) * X( _ , j)),2);
            }
            v_norm2 = sqrt(v_norm2);
            //Calculate upper and lower bound
            upperBound[*iter] = w_tilde[*iter] - w_r[*iter] + z_r[*iter] + v_norm2 * w_norm2/n;
            lowerBound[*iter] = w_tilde[*iter] - w_r[*iter] + z_r[*iter] - v_norm2 * w_norm2/n;
  
            double w_tilde_tmp = w_tilde[*iter];
            if (upperBound[*iter] > lambda || lowerBound[*iter] < - lambda){
              double tmp2 = 0;
              //Update w_tilde[i] by equation 5
              for (int j = 0; j < p; j++){
                if (w_tilde[j] != 0){
                  tmp2 += sum(X( _ , *iter) * X( _ , j)) * w_tilde[j];
                }
              }
              z[*iter] = w_tilde[*iter] + (sum(X( _ , *iter) * y) - tmp2) / n;
              if (abs(z[*iter]) > lambda){
                if (z[*iter] > 0){
                  w_tilde[*iter] = z[*iter] - lambda;
                }
                else{
                  w_tilde[*iter] = z[*iter] + lambda;
                }
              }
              else{
                w_tilde[*iter] = 0;
              }
            }
            else{
              w_tilde[*iter] = 0;
            }
            //Update w_norm2 by equation 10
            w_norm2 = sqrt(pow(w_norm2, 2) - pow(w_tilde_tmp, 2) + pow(w_tilde[*iter], 2));
          }
          if (sum((w_tmp - w_tilde) * (w_tmp - w_tilde)) < tol){  //Need to adjusted!!!
            isConverge = true;
          }
        }
        
        isConverge = false;  //Reset the flag for next while loop
        //Calculate Sk
        NumericVector Sk;
        if (k != 0){
          NumericVector w_tilde_previous = res( _ , k-1);
          Sk = seqStrongRule(X, y, w_tilde_previous, lambdas, k);
        }
        
        isSkViolateKKT = false;  //Reset the end loop flag
        for (NumericVector::iterator iter = Sk.begin(); iter != Sk.end(); iter++){//For pi in S
          if (isViolateKKT(X, y, w_tilde, lambda, *iter)){
            U.insert(*iter);
            isSkViolateKKT = true;  //Cannot end the loop
          }
        }
      }
      
      isPViolateKKT = false;  //Reset the end loop flag
      for (int l = 0; l < p; l++){//For pi in P
        if (isViolateKKT(X, y, w_tilde, lambda, l)){
          U.insert(l);
          isPViolateKKT = true;
        }
      }
    }
    
    res( _ , k) = clone(w_tilde);
    for (int m = 0; m < p; m++){
      if (w_tilde[m] != 0){
        P_k.push_back(m);
      }
    }
  }
  
  NumericVector b = {1,2,3};
  return List::create(Named("beta") = b, Named("res", res));
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
List sling2(NumericMatrix &X, NumericVector &y, double l){
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
    prev = prev*0 + R_PosInf;
    while (is_true(all(abs(w-prev)>err))){
      prev = w;  //Need to clone?
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
  }
  return List::create(Named("w", w), /*Named("ws", ws),*/ Named("itr", itr));
}

