#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

const double SINGULAR_EPS = 10e-10; // criteria for matrix singularity

//' @description Updates parameter estimates of a finite mixture of
//' Gaussian regressions by the EM algorithm.
//' @export
//' @title cppRegmixMLE_homo
//' @name cppRegmixMLE_homo
//' @param bs (m-1) + (q+1)m + 1 + p by ninits matrix of initial values of (alpha,mu,beta,sigma,gamma).
//' @param ys n by 1 vector of data for y.
//' @param xs x n by q matrix of data for x.
//' @param zs n by p matrix of regressor associated with gamma.
//' @param m number of components in the mixture.
//' @param p dimension of the regressor associated with gamma.
//' @param maxit maximum number of iterations.
//' @param ninits number of initial values.
//' @param tol Convergence is declared when the penalized log-likelihood increases by less than \code{tol}.
//' @return  A list with items:
//' \item{loglikset}{vector of the maximized value of the log-likelihood.}
//' \item{notcg}{vector that records whether EM steps converged or not for each initial value.}
//' \item{post}{n*m by ininits matrix of posterior probabilities for observations.}
//'
// [[Rcpp::export]]
List cppRegmixMLE_homo(NumericMatrix bs,
                      NumericVector ys,
                      NumericMatrix xs,
                      NumericMatrix zs,
                      int m,
                      int p,
                      int maxit = 2000,
                      int ninits = 10,
                      double tol = 1e-8) {
  int n = ys.size();
  int q = xs.ncol();
  int q1 = q + 1;
  arma::mat b(bs.begin(), bs.nrow(), bs.ncol(), false);
  arma::vec y(ys.begin(), ys.size(), false);
  arma::mat x(xs.begin(), xs.nrow(), xs.ncol(), false);
  arma::mat z(zs.begin(), zs.nrow(), zs.ncol(), false);
  arma::vec b_jn(bs.nrow());
  arma::vec alpha(m), ssr(m), alp_sig(m);
  arma::mat mubeta(q1,m);
  arma::vec r(m), l_j(m);
  arma::mat w(m,n);
  arma::mat post(m*n,ninits);
  arma::vec notcg(ninits), loglikset(ninits);
  arma::vec gamma(p);
  arma::vec ytilde(n);
  arma::vec wtilde(n);
  arma::mat xtilde(n,q1);
  arma::mat ztilde(n,p);
  arma::mat zz(p,p);
  arma::mat ze(p,1);
  int emit, sing;
  double oldloglik, diff, minr, w_j, sum_l_j, sigma;
  arma::mat x1(n,q1);
  double loglik = 0; // force initialization
  double ll = 0; // force initialization
  notcg.zeros();  // initialization

  x1.col(0) = arma::ones(n);
  for (int i=1; i < q1; ++i) {
    x1.col(i) = x.col(i-1);
  }

  /* iteration over ninits initial values of b */
  for (int jn=0; jn<ninits; jn++) {

    /* initialize EM iteration */
    b_jn = b.col(jn);
    for (int j=0; j < m; ++j){
      alpha(j) = b_jn(j);
      for (int i=0; i<q1; ++i){
        mubeta(i,j) = b_jn(m+q1*j+i);
      }
    }
    sigma = b_jn((q1+1)*m);
    if (p>0) {
      for (int j=0; j < p; j++){
        gamma(j) = b_jn((q1+1)*m+1+j);
      }
    }
    oldloglik = R_NegInf;
    emit = 0;
    diff = 1.0;
    sing = 0;

    /* EM loop begins */
    for (int iter = 0; iter < maxit; iter++) {
      ll = - (double)n * M_LN_SQRT_2PI; /* n/2 times log(2pi) */
      alp_sig = alpha/sigma;

      if (p==0) {
        ytilde = y;
      } else {
        ytilde = y - z*gamma;
      }

      for (int i = 0; i < n; i++) {
        /* standardized squared residual */
        r = (1.0/sigma) * ( ytilde(i) - trans(x1.row(i)*mubeta) );
        r = 0.5 * (r % r); /* This is faster than r = pow( r, 2.0 ) */
        minr = min(r);
        /* posterior for i */
        /* normalizing with minr avoids the problem of dividing by zero */
        l_j =  alp_sig % exp( minr-r );
        sum_l_j = sum( l_j );
        w.col(i) = l_j/sum_l_j; /* w(j,i) = alp_j*l_j / sum_j (alp_j*l_j) */
        /* loglikelihood*/
        ll +=  log(sum_l_j) - minr; /* subtract back minr */
      } /* end for (i=0; i<n; i++) loop */

      loglik = ll;
      diff = loglik - oldloglik;
      oldloglik = loglik;

      /* Normal exit */
      if (diff < tol || emit>=maxit){
        break;
      }
       emit++;

      /* update alpha, mu, and sigma */
      for (int j = 0; j < m; j++) {
        w_j = sum( w.row(j) ); /* w_j(j) = sum_i w(i,j) */
        alpha(j) = w_j / n;
        wtilde = trans(w.row(j));
        for (int ii = 0; ii < q1; ii++) {
          xtilde.col(ii) = wtilde % x1.col(ii);
        }

        arma::mat design_matrix = trans(xtilde) * x1;
        if (rcond(design_matrix) < SINGULAR_EPS)
        {
          sing = 1;
          break;
        }

        mubeta.col(j) = solve( design_matrix , trans(xtilde) * ytilde );
        //mubeta.col(j) = solve( trans(xtilde) * x1 , trans(xtilde) * ytilde );
        ssr(j) = sum( trans(w.row(j)) % pow(  ytilde - x1*mubeta.col(j) , 2 ) );
        alpha(j) = fmin( fmax(alpha(j),0.01), 0.99);
      }
      sigma = sqrt( sum(ssr) / n );

      if (p>0) { /* update gamma */
        zz.zeros();
        ze.zeros();
        for (int j = 0; j < m; j++) {
          wtilde = trans(w.row(j)) ;
          for (int ii = 0; ii < p; ii++) {
            ztilde.col(ii) = wtilde % z.col(ii);
          }
          zz = zz + trans(ztilde) * z;
          ze = ze + trans(ztilde) * (y - x1*mubeta.col(j));
        }

        if (rcond(zz) < SINGULAR_EPS)
        {
          sing = 1;
          break;
        }

        gamma = solve(zz,ze);
      }

      /* Check singularity */
      for (int j=0; j<m; j++) {
        if (alpha(j) < 1e-8 || std::isnan(alpha(j)) || sigma < 1e-8){
          sing = 1;
        }
      }

      /* Exit from the loop if singular */
      if (sing) {
        notcg(jn) = 1;
        break;
      }

    }/* EM loop ends */

    loglikset(jn) = ll;
    for (int j=0; j < m; j++){
      b_jn(j) = alpha(j);
      for (int i=0; i<q1; ++i){
        b_jn(m+q1*j+i) = mubeta(i,j);
      }
    }
    b_jn((q1+1)*m) = sigma;
    if (p>0) {
      for (int j=0; j < p; j++){
        b_jn((q1+1)*m+1+j) = gamma(j);
      }
    }
    b.col(jn) = b_jn; /* b is updated */
    post.col(jn) = vectorise(trans(w));

  } /* end for (jn=0; jn<ninits; jn++) loop */

  return Rcpp::List::create(Named("loglikset") = wrap(loglikset),
                            Named("notcg") = wrap(notcg),
                              Named("post") = wrap(post)
  );
}

