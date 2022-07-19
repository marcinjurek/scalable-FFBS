#define ARMA_64BIT_WORD 1
#define ARMA_DONT_PRINT_ERRORS

#include <iostream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
double multiply(arma::uvec p, arma::uvec inds, arma::vec x, uword i, uword j){
  
  double prod = 0;
  
  arma::uvec iinds = inds.subvec( p(i), p(i + 1) - 1 );
  arma::uvec jinds = inds.subvec( p(j), p(j + 1) - 1 );
  
  uword k = 0;
  uword l = 0;
  while( k < iinds.size() & l < jinds.size() ){
    if( iinds(k) < jinds(l) ){
      k++;
    } else if( iinds(k) > jinds(l) ) {
      l++;
    } else {
      prod += x(p(i)+k) * x(p(j)+l);
      k++;
    }
  }
  return prod;
}



// [[Rcpp::export]]
arma::mat getMatCovFromFactorCpp(S4 F, arma::Mat<arma::uword> revNNarray){
  
  
  
  // obtain dim, i, p. x from S4 object
  IntegerVector dims = F.slot("Dim");
  arma::uvec inds = Rcpp::as<arma::urowvec>(F.slot("i")).t();
  arma::uvec p = Rcpp::as<arma::urowvec>(F.slot("p")).t();
  arma::vec x     = Rcpp::as<arma::vec>(F.slot("x"));
  
  arma::mat sigSel = arma::zeros<arma::mat>(revNNarray.n_rows, revNNarray.n_cols);
  
  
  
  
  
  for( uword j = 0; j < revNNarray.n_rows; j++ ){
    for( uword i = 0; i < revNNarray.n_cols; i++ ){
      if( revNNarray(i, j) > 0 ){
        sigSel(i, j) = multiply(p, inds, x, i, revNNarray(i,j) - 1);
      }
      
    }
  }
  
  return sigSel;
}


// [[Rcpp::export]]
arma::mat getMatCovFromFactorCppOld(const arma::sp_mat F, const arma::umat revNNarray){
  
  arma::mat sigSel = arma::zeros<arma::mat>(revNNarray.n_rows, revNNarray.n_cols);
  
  for(uword i=0; i < revNNarray.n_rows; i++) {
    
    arma::uvec r = revNNarray.row( i ).t();
    arma::uvec inds = find( r );
    arma::uvec cols = r.elem( inds ) - 1;
    arma::sp_mat thisCol = F.col( i ).t();
    
    for(uword colnum=0; colnum<cols.n_rows; colnum++){
      
      arma::sp_mat cl = F.col( cols(colnum) ); 
      arma::sp_mat val = thisCol * cl;
      sigSel( i, inds(colnum) ) = val( 0, 0 );
    }
    
    
    
  }
  return sigSel;
}
