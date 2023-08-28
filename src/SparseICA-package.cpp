#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double soft_thresh(double x, double alpha){
  if (x > alpha){
    return(x-alpha);
  }else if (x < -alpha){
    return(x + alpha);
  }else{
    return(0);
  }
}


// [[Rcpp::export]]
arma::mat procrustes(arma::mat& X, arma::mat& V){
  arma::mat U;

  arma::vec s;
  arma::mat Q, R;

  arma::svd_econ(Q, s, R, X * V);
  U = Q * R.t();
  return(U);
}


// [[Rcpp::export]]
Rcpp::List relax_laplace(const arma::mat& xData,const arma::mat& newV, double nu, double lambda, double maxit, double eps = 0.000001){
  // initialize
  int p = newV.n_rows;
  int r = newV.n_cols;
  int t = xData.n_cols;
  arma::mat newU(t,r);
  arma::vec converge(maxit);
  arma::mat XTV(t,r);

  // Intermediate V and XU
  arma::mat lagU = arma::eye(t,t);
  arma::mat newV2 = newV;
  arma::mat XU(p,r);
  arma::mat IdenMat=arma::eye(r,r);

  for (int i = 0; i < maxit; i++)
  {
    // update U, use irina's function
    XTV = xData.t()*newV2;
    newU = procrustes(XTV,IdenMat);

    // update V
    XU = xData*newU;
    for (int j = 0; j < p; j++){
      for (int k = 0; k < r; k++){
        newV2(j, k) = soft_thresh(XU(j, k), nu/lambda);
      }
    }

    // check convergence
    arma::mat temp = newU*lagU.t();
    converge[i] = max(abs(abs(temp.diag()) - 1));
    //converge[i] = arma::norm(lagU-newU,"fro");
    //converge[i] = accu(square(lagV-newV2))/(p*r*r);
    if (converge[i] < eps){
      break;
    }
    lagU = newU;
  }
  return Rcpp::List::create(Rcpp::Named("newU") = newU, Rcpp::Named("newV") = newV2, Rcpp::Named("converge") = converge);
}


// [[Rcpp::export]]
Rcpp::List runs_relax_laplace(const arma::mat& xData, Rcpp::List& W,double runs, double r,double nu, double lambda, double maxit, double eps = 0.000001){
  //arma::mat myW = W[0];
  //arma::mat newV = xData*myW;

  int p = xData.n_rows;
  double my_norm;
  arma::mat newV(p,r);
  arma::vec loglik(runs);
  double max_location;

  for (int m = 0; m < runs; m++){
    arma::mat myW = W[m];
    newV = xData*myW;
    Rcpp::List Rcpp_res=relax_laplace(xData,newV,nu,lambda,maxit,eps);

    arma::mat Rcpp_newV = Rcpp_res["newV"];
    arma::mat Rcpp_newU = Rcpp_res["newU"];
    my_norm = arma::norm(Rcpp_newV-xData*Rcpp_newU,"fro");
    loglik[m]=accu(-abs(Rcpp_newV)/lambda-log(2*lambda))-1/(2*nu)*my_norm*my_norm;
  }

  max_location = std::max_element(loglik.begin(),loglik.end())-loglik.begin();
  arma::mat myW = W[max_location];
  newV = xData*myW;
  Rcpp::List Rcpp_res=relax_laplace(xData,newV,nu,lambda,maxit,eps);

  return Rcpp::List::create(Rcpp::Named("newU") = Rcpp_res["newU"], Rcpp::Named("newV") = Rcpp_res["newV"], Rcpp::Named("converge") = Rcpp_res["converge"],Rcpp::Named("loglik") =loglik.max());
}

