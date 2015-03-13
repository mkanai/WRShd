#include <math.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Remove missing numeric values in objects
// [[Rcpp::export(".na.omit")]]
NumericVector na_omit(NumericVector x, int cores = 1) {
    #ifdef _OPENMP
    if (cores > 0) {
        omp_set_num_threads(cores);
    }
    #endif
    int n = x.size();
    std::vector<double> ret(n);

    int k = 0;
    #pragma omp parallel for reduction(+:k)    
    for (int i = 0; i < n; i++) {
        if (x[i] != NA_REAL) {ret[k++] = x[i];}
    }
    #pragma omp barrier
    
    ret.resize(k);
    return Rcpp::wrap(ret);
}

NumericVector sample(NumericVector x, int size, int replace = FALSE, NumericVector prob = NumericVector::create()) {
    return Rcpp::RcppArmadillo::sample(x, size, replace, prob);
}

//' Compute the Harrell-Davis estimate of the qth quantile
// [[Rcpp::export(".hd")]]
double hd(NumericVector x, double q = 0.5, int na_rm = TRUE, int cores = 1) {
    #ifdef _OPENMP
    if (cores > 0) {
        omp_set_num_threads(cores);
    }
    #endif

    NumericVector xx = na_rm ? na_omit(x, cores) : x;
    int n = xx.size();
    double m1 = (n + 1.0) * q, m2 = (n + 1.0) * (1.0 - q);
    double output = 0.0;
    arma::vec xsort = arma::sort(as<arma::vec>(xx));
    
    #pragma omp parallel for reduction(+:output)
    for(int i = 1; i <= n; i++){
        output += (R::pbeta(i * 1.0 / n, m1, m2, 1, 0) - 
                   R::pbeta((i - 1.0) / n, m1, m2, 1, 0)) * 
                   xsort(i - 1);
    }

    return(output);
}

double hd(arma::vec x, double q, int cores = 1) {
    #ifdef _OPENMP
    if (cores > 0) {
        omp_set_num_threads(cores);
    }
    #endif

    int n = x.n_elem;
    double m1 = (n + 1.0) * q, m2 = (n + 1.0) * (1.0 - q);
    double output = 0.0;
    arma::vec xsort = arma::sort(x);

    #pragma omp parallel for reduction(+:output)
    for(int i = 1; i <= n; i++) {
        output += (R::pbeta(i * 1.0 / n, m1, m2, 1, 0) - 
                   R::pbeta((i - 1.0) / n, m1, m2, 1, 0)) * 
                   xsort(i - 1);
    }

    return output;
}

//' Compute a bootstrap standard error of the Harrell-Davis estimate of the qth quantile
// [[Rcpp::export(".hdseb")]]
double hdseb(NumericVector x, double q = 0.5, int nboot = 100, int cores = 1) {
    #ifdef _OPENMP
    if (cores > 0) {
        omp_set_num_threads(cores);
    }
    #endif
    
    int n = x.size();
    NumericVector bsample = sample(x, n * nboot, TRUE);
    arma::mat data = arma::mat(bsample.begin(), nboot, n, FALSE);
    
    arma::vec bvec(nboot);    
    #pragma omp parallel
    for (int i = 0; i < nboot; i++) {
        bvec(i) = hd(data.row(i).t(), q);
    }
    #pragma omp barrier
    
    return sqrt(var(bvec));
}

//' Compute a bootstrap confidence interval for the Harrell-Davis estimate of the qth quantile
// [[Rcpp::export(".hdpb")]]
List hdpb(NumericVector x, double q = 0.5, double alpha = 0.05, int nboot = 2000, double nv = 0, int cores = 1) {
    #ifdef _OPENMP
    if (cores > 0) {
        omp_set_num_threads(cores);
    }
    #endif

    NumericVector xx = na_omit(x, cores);
    int n = xx.size();
    NumericVector bsample = sample(xx, n * nboot, TRUE);
    arma::mat data = arma::mat(bsample.begin(), nboot, n, FALSE);

    int nv_ex = 0, nv_eq = 0;
    arma::vec bvec(nboot);
    
    #pragma omp parallel for reduction(+:nv_ex, nv_eq)
    for (int i = 0; i < nboot; i++) {
        bvec(i) = hd(data.row(i).t(), q);
        if (bvec(i) > nv) {nv_ex++;}
        else if (bvec(i) == nv) {nv_eq++;}
    }
    #pragma omp barrier
    
    arma::vec bvecsort = arma::sort(bvec);
    int low = round((alpha / 2) * nboot);
    int up = nboot - low - 1;
    double pv = 1.0 * nv_ex / nboot + 0.5 * nv_eq / nboot;
    pv = 2 * std::min(pv, 1 - pv);
    double estimate = hd(xx, q);
    return List::create(Named("ci") = NumericVector::create(bvecsort[low], bvecsort[up]),
                        Named("n") = n,
                        Named("estimate") = estimate,
                        Named("p.value") = pv
                       );
}


// [[Rcpp::export(".qcom.sub")]]
List qcom_sub(NumericVector x, NumericVector y, double q, double alpha = 0.05, int nboot = 2000, int cores = 1) {
    #ifdef _OPENMP
    if (cores > 0) {
        omp_set_num_threads(cores);
    }
    #endif
    
    NumericVector xx = na_omit(x, cores);
    NumericVector yy = na_omit(y, cores);
    int nx = xx.size();
    int ny = yy.size();
    NumericVector bsamplex = sample(xx, nx * nboot, TRUE);
    NumericVector bsampley = sample(yy, ny * nboot, TRUE);
    arma::mat datax = arma::mat(bsamplex.begin(), nboot, nx, FALSE);
    arma::mat datay = arma::mat(bsampley.begin(), nboot, ny, FALSE);

    int z_ex = 0, z_eq = 0;    
    arma::vec bvec(nboot);

    #pragma omp parallel for reduction(+:z_ex, z_eq) 
    for (int i = 0; i < nboot; i++) {
        double bvecx = hd(datax.row(i).t(), q);
        double bvecy = hd(datay.row(i).t(), q);
        bvec(i) = bvecx - bvecy;
        if (bvec(i) > 0) {z_ex++;}
        else if (bvec(i) == 0) {z_eq++;}
    }
    #pragma omp barrier

    arma::vec bvecsort = arma::sort(bvec);
    int low = round((alpha / 2) * nboot);
    int up = nboot - low - 1;
    double pv = 1.0 * z_ex / nboot + 0.5 * z_eq / nboot;
    pv = 2 * std::min(pv, 1 - pv);
    double se = var(bvec);
    return List::create(Named("est.1") = hd(xx, q),
                        Named("est.2") = hd(yy, q),
                        Named("ci") = NumericVector::create(bvec[low], bvec[up]),
                        Named("p.value") = pv,
                        Named("sq.se") = se,
                        Named("n1") = nx,
                        Named("n2") = ny
                       );
}
