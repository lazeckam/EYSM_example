#include <Rcpp.h> 
#include <string.h>
#include <math.h> 
using namespace Rcpp;


// [[Rcpp::export]]
double entropy(Rcpp::NumericVector p) {
  
  int n = p.length();
  double result = 0;
  
  for(int i = 0; i < n; i++){
    if(p[i] > 0){
      result -= p[i]*std::log(p[i]);
    }
  }
  
  return result;
}

// [[Rcpp::export]]
double conditional_mutual_information(Rcpp::NumericVector p_xyz) {
  
  int n = p_xyz.length();
  Rcpp::NumericVector p_xz(n/2), p_yz(n/2), p_z(n/4);
  double result = 0;
  
  result -= entropy(p_xyz);
  
  for(int i = 0; i < n/2; i++){
    p_xz[i] = 0;
    p_yz[i] = 0;
  }
  for(int i = 0; i < n/4; i++){
    p_z[i] = 0;
  }
  
  for(int i = 0; i < n; i++){
    p_xz[i/2] += p_xyz[i];
    p_yz[i % 2 + 2*(i / 4)] += p_xyz[i];
  }
  
  result += entropy(p_xz);
  result += entropy(p_yz);
  
  for(int i = 0; i < n/2; i++){
    p_z[i/2] += p_xz[i];
  }
  
  result -= entropy(p_z);
  
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector projection_CI (Rcpp::NumericVector p_xyz) {
  
  int n = p_xyz.length();
  Rcpp::NumericVector p_xz(n/2), p_yz(n/2), p_z(n/4);
  Rcpp::NumericVector p_ci(n);
  
  for(int i = 0; i < n/2; i++){
    p_xz[i] = 0;
    p_yz[i] = 0;
  }
  for(int i = 0; i < n/4; i++){
    p_z[i] = 0;
  }
  
  for(int i = 0; i < n; i++){
    p_xz[i/2] += p_xyz[i];
    p_yz[i % 2 + 2*(i / 4)] += p_xyz[i];
  }
  
  for(int i = 0; i < n/2; i++){
    p_z[i/2] += p_xz[i];
  }
  
  for(int i = 0; i < n; i++){
    if(p_z[i / 4]!= 0){
      p_ci[i] = p_xz[i/2]*p_yz[i % 2 + i / 4]/p_z[i / 4];
    }
  }
  
  return p_ci;
}


// [[Rcpp::export]]
Rcpp::NumericVector compute_p_yz(Rcpp::NumericVector p_xyz) {
  
  int n = p_xyz.length();
  Rcpp::NumericVector p_xz(n/2);
  
  for(int i = 0; i < n/2; i++){
    p_xz[i] = 0;
  }
  
  for(int i = 0; i < n; i++){
    p_xz[i/2] += p_xyz[i];
  }
  
  return p_xz;
}

// [[Rcpp::export]]
Rcpp::NumericVector compute_p_xz(Rcpp::NumericVector p_xyz) {
  
  int n = p_xyz.length();
  Rcpp::NumericVector p_yz(n/2);
  
  for(int i = 0; i < n/2; i++){
    p_yz[i] = 0;
  }
  
  for(int i = 0; i < n; i++){
    p_yz[i % 2 + 2*(i / 4)] += p_xyz[i];
  }
  
  return p_yz;
}

// [[Rcpp::export]]
Rcpp::NumericVector compute_p_z(Rcpp::NumericVector p_xz) {
  
  int n = 2*p_xz.length();
  Rcpp::NumericVector p_z(n/4);
  
  for(int i = 0; i < n/4; i++){
    p_z[i] = 0;
  }
  
  for(int i = 0; i < n/2; i++){
    p_z[i/2] += p_xz[i];
  }
  
  return p_z;
}

// [[Rcpp::export]]
Rcpp::NumericVector compute_p_x_z(Rcpp::NumericVector p_xyz) {
  
  int n = p_xyz.length();
  Rcpp::NumericVector p_xz(n/2), p_z(n/4), p_x_z(n/2);
  
  for(int i = 0; i < n/2; i++){
    p_xz[i] = 0;
  }
  for(int i = 0; i < n/4; i++){
    p_z[i] = 0;
  }
  
  for(int i = 0; i < n; i++){
    p_xz[i % 2 + 2*(i / 4)] += p_xyz[i];
  }
  
  for(int i = 0; i < n/2; i++){
    p_z[i/2] += p_xz[i];
  }
  
  for(int i = 0; i < n/2; i++){
    if(p_z[i / 2]!= 0){
      p_x_z[i] = p_xz[i]/p_z[i/2];
    }
  }
  
  return p_x_z;
}
