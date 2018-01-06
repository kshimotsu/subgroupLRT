#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

const double SINGULAR_EPS = 10e-10; // criteria for matrix singularity

//' @description Updates parameter estimates of a finite mixture of
//' Gaussian regressions by the EM algorithm.
//' @export
//' @title cppSubgroupMLE_homo
//' @name cppSubgroupMLE_homo
//' @param bs q2 + (q+1)m + 1 + p by ninits matrix of initial values of (tau,mu,beta,sigma,gamma).
//' @param ys n by 1 vector of data for y.
//' @param xs x n by q1-1 matrix of data for x.
//' @param vs v n by q2-1 matrix of data for v
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
List cppSubgroupMLE_homo(NumericMatrix bs,
                      NumericVector ys,
                      NumericMatrix xs,
                      NumericMatrix vs,
                      NumericMatrix zs,
                      int m,
                      int p,
                      int maxit = 2000,
                      int ninits = 10,
                      double tol = 1e-8) {
  int n = ys.size();
  int q = xs.ncol();
  int qv = vs.ncol();
  int q1 = q + 1;
  int q2 = qv + 1;
  arma::mat b(bs.begin(), bs.nrow(), bs.ncol(), false);
  arma::vec y(ys.begin(), ys.size(), false);
  arma::mat x(xs.begin(), xs.nrow(), xs.ncol(), false);
  arma::mat v(vs.begin(), vs.nrow(), vs.ncol(), false);
  arma::mat z(zs.begin(), zs.nrow(), zs.ncol(), false);
  arma::vec b_jn(bs.nrow());
  arma::vec alpha(m), ssr(m), alp_sig(m);
  arma::mat mubeta(q1,m);
  arma::mat tau(q2,1);
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
  int emit, sing, tauit;
  double oldloglik, oldlogliktau, diff, difftau, minr, w_j, sum_l_j, sigma;
  arma::mat x1(n,q1);
  arma::mat v1(n,q2);
  double vtau =0;
  double asum;
  double loglik = 0; // force initialization
  double ll = 0; // force initialization
  notcg.zeros();  // initialization
  arma::vec pi(n);
  double logliktau = 0;
  arma::mat grad(q2,1);
  arma::mat mhess(q2,q2);
  arma::mat Imat(q2,q2);
  arma::vec lambda(n);
  arma::vec restau(n);

  x1.col(0) = arma::ones(n);
  for (int i=1; i < q1; ++i) {
    x1.col(i) = x.col(i-1);
  }

  v1.col(0) = arma::ones(n);
  for (int i=1; i < q2; ++i) {
    v1.col(i) = v.col(i-1);
  }


  /* iteration over ninits initial values of b */
  for (int jn=0; jn<ninits; jn++) {

    /* initialize EM iteration */
    b_jn = b.col(jn);
    for (int j=0; j < q2; ++j){
      tau(j) = b_jn(j);
    }

    for (int j=0; j < m; ++j){
      for (int i=0; i<q1; ++i){
        mubeta(i,j) = b_jn(q2+q1*j+i);
      }
    }
  sigma = b_jn(q2+q1*m);
  if (p>0) {
    for (int j=0; j < p; j++){
      gamma(j) = b_jn(q2+q1*m+1+j);
    }
  }
  oldloglik = R_NegInf;
  emit = 0;
  diff = 1.0;
  sing = 0;

    /* EM loop begins */
    for (int iter = 0; iter < maxit; iter++) {
      ll = - (double)n * M_LN_SQRT_2PI; /* n/2 times log(2pi) */
  //     alp_sig = alpha/sigma;

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
        /* alpha */
        vtau = as_scalar(v1.row(i)*tau);
        // alpha(0) = exp( vtau ) / ( 1 + exp ( vtau ));
        asum = 0;
        for (int j = 0; j < m-1; j++) {
          alpha(j) = exp( vtau ) / ( 1 + exp( vtau ) );
          asum += alpha(j);
        }
        alpha(m-1) = 1 - asum;
        alp_sig = alpha/sigma;
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
      //Rprintf("diff = %f \n", diff);

      /* update pi vector */
      pi = w.row(0).t();

      /* Normal exit */
      if (diff < tol || emit>=maxit){
        break;
      }
      emit++;

      /* update tau, mu, and sigma */
      for (int j = 0; j < m; j++) {
        w_j = sum( w.row(j) ); /* w_j(j) = sum_i w(i,j) */
  //       alpha(j) = w_j / n;
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
  //
        mubeta.col(j) = solve( design_matrix , trans(xtilde) * ytilde );
        //mubeta.col(j) = solve( trans(xtilde) * x1 , trans(xtilde) * ytilde );
        ssr(j) = sum( trans(w.row(j)) % pow(  ytilde - x1*mubeta.col(j) , 2 ) );
  //       alpha(j) = fmin( fmax(alpha(j),0.01), 0.99);
      }
      sigma = sqrt( sum(ssr) / n );

      /* Update tau */
      tauit = 0;
      difftau = 0;
      oldlogliktau = R_NegInf;
      grad = arma::zeros(q2,1);
      mhess = arma::zeros(q2,q2);
      for (int itertau = 0; itertau < maxit; itertau++) {
        for (int i = 0; i < n; i++) {
          vtau = as_scalar(v1.row(i)*tau);
          lambda(i) = exp(vtau) / ( 1 + exp(vtau));
          restau(i) = pi(i) - lambda(i);
          grad += trans( restau(i) * v1.row(i) );
          mhess -= lambda(i)*(1-lambda(i))*trans( v1.row(i) ) * v1.row(i); /* H */
          // mhess += lambda(i)*restau(i)*trans( v1.row(i) ) * v1.row(i); /* -H */
        }
        if (rcond(mhess) > SINGULAR_EPS) {
          tau = tau - inv_sympd(mhess) * grad; /* update tau */
        } else {
          tau = tau - 1.0*Imat.eye() * grad;
          /* This one is slow but aviods errors.
           * In the future, we should use BFGS. */
        }

        logliktau = 0;
        for (int i = 0; i < n; i++) {
          vtau = as_scalar(v1.row(i)*tau);
          lambda(i) = exp(vtau) / ( 1 + exp(vtau));
          logliktau += pi(i) * log( lambda(i) ) + ( 1 - pi(i) ) * log( 1 - lambda(i) ); /* compute the value of log likelihood */
        }
        difftau = logliktau - oldlogliktau;
        oldlogliktau = logliktau;

        /* Normal exit */
        if (difftau < tol || tauit>=maxit){
          break;
        }

        /* Check singularity. Exit from the loop if singular. */
          // if (alpha(j) < 1e-8 || std::isnan(alpha(j))|| std::isnan(logliktau)){
        if (std::isnan(logliktau)){
          sing = 1;
          break;
        }

        tauit++;
      } /* end of itertau loop */


  //
  //     /* test */
  //     if (p>0) { /* update gamma */
  //       zz.zeros();
  //       ze.zeros();
  //       for (int j = 0; j < m; j++) {
  //         wtilde = trans(w.row(j)) ;
  //         for (int ii = 0; ii < p; ii++) {
  //           ztilde.col(ii) = wtilde % z.col(ii);
  //         }
  //         zz = zz + trans(ztilde) * z;
  //         ze = ze + trans(ztilde) * (y - x1*mubeta.col(j));
  //       }
  //
  //       if (rcond(zz) < SINGULAR_EPS)
  //       {
  //         sing = 1;
  //         break;
  //       }
  //
  //       gamma = solve(zz,ze);
  //     }
  //
  //     /* Check singularity */
  //     for (int j=0; j<m; j++) {
  //       if (alpha(j) < 1e-8 || std::isnan(alpha(j)) || sigma < 1e-8){
  //         sing = 1;
  //       }
  //     }
  //
      /* Exit from the loop if singular */
      if (sing) {
        notcg(jn) = 1;
        break;
      }

     }/* EM loop ends */

     loglikset(jn) = ll;
     for (int j = 0; j<q2; j++) {
       b_jn(j) = tau(j);
     }
     for (int j=0; j < m; j++){
      for (int i=0; i<q1; ++i){
        b_jn(q2+q1*j+i) = mubeta(i,j);
      }
     }
     b_jn(q2+q1*m) = sigma;
     if (p>0) {
       for (int j=0; j < p; j++){
         b_jn(q2 + q1*m+1+j) = gamma(j);
       }
     }
     b.col(jn) = b_jn; /* b is updated */
     post.col(jn) = vectorise(trans(w));

  } /* end for (jn=0; jn<ninits; jn++) loop */
  //
  // return Rcpp::List::create(Named("loglikset") = wrap(loglikset),
  //                           Named("notcg") = wrap(notcg),
  //                             Named("post") = wrap(post)
  // );
  return Rcpp::List::create(Named("loglikset") = wrap(loglikset),
                            Named("post") = wrap(post),
                            Named("notcg") = wrap(notcg),
                            Named("logliktau") = wrap(logliktau)
                            );
}

