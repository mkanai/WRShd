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
NumericVector na_omit(NumericVector x) {
    int n = x.size();
    std::vector<double> ret(n);

    int k = 0;
    for (int i = 0; i < n; i++) {
        if (x[i] != NA_REAL) {ret[k++] = x[i];}
    }

    ret.resize(k);
    return Rcpp::wrap(ret);
}

NumericVector sample(NumericVector x, int size, bool replace = true, NumericVector prob = NumericVector::create()) {
    return Rcpp::RcppArmadillo::sample(x, size, replace, prob);
}

//' Compute the Harrell-Davis estimate of the qth quantile
// [[Rcpp::export(".hd")]]
double hd(NumericVector x, double q = 0.5, bool na_rm = true, int cores = 1) {
    #ifdef _OPENMP
    if (cores > 0) {
        omp_set_num_threads(cores);
    }
    #endif

    NumericVector xx = na_rm ? na_omit(x) : x;
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

double hd(arma::vec x, double q, bool na_rm = true, int cores = 1) {
    int n = x.n_elem;
    double m1 = (n + 1.0) * q, m2 = (n + 1.0) * (1.0 - q);
    double output = 0.0;
    arma::vec xsort = arma::sort(x);

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
    NumericVector bsample = sample(x, n * nboot, true);
    arma::mat data = arma::mat(bsample.begin(), nboot, n, false);
    
    arma::vec bvec(nboot);    
    #pragma omp parallel for
    for (int i = 0; i < nboot; i++) {
        bvec(i) = hd(data.row(i).t(), q);
    }
    #pragma omp barrier
    
    return sqrt(var(bvec));
}


//' Compute a 1-alpha confidence for the Harrell-Davis estimate of the qth quantile
// [[Rcpp::export(".hdci")]]
List hdci(NumericVector x, double q = 0.5, int nboot = 100, bool pr = true, int cores = 1) {
    NumericVector xx = na_omit(x);
    int n = xx.size();
    if (pr && sum(duplicated(xx)) > 0) {Rprintf("Duplicate values detected; use hdpb\n");}
    double se = hdseb(x, q, nboot, cores);
    
    double crit = 0.5064 / pow(n, 0.25) + 1.96;;
    if ((q <= 0.2 || q >= 0.8) && n <= 20) {
        crit = (-6.23) / n + 5.01;
    }
    if ((q <= 0.1 || q >= 0.9) && n <= 40) {
        crit = 36.2 / n + 1.31;
    }
    if (n <= 10){
        Rprintf("The number of observations is less than 11.");
        Rprintf("Accurate critical values have not been determined for this case.");
    }
    
    double low = hd(xx, q, cores) - crit*se;
    double hi = hd(xx, q, cores) + crit*se;
    return List::create(Named("ci") = NumericVector::create(low, hi),
                        Named("crit") = crit,
                        Named("se") = se
                       );
}

//' Compute a bootstrap 1-alpha confidence for the Harrell-Davis estimate of the qth quantile
// [[Rcpp::export(".hdpb")]]
List hdpb(NumericVector x, double q = 0.5, double alpha = 0.05, int nboot = 2000, double nv = 0, int cores = 1) {
    #ifdef _OPENMP
    if (cores > 0) {
        omp_set_num_threads(cores);
    }
    #endif

    NumericVector xx = na_omit(x);
    int n = xx.size();
    NumericVector bsample = sample(xx, n * nboot, true);
    arma::mat data = arma::mat(bsample.begin(), nboot, n, false);

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
    double estimate = hd(xx, q, cores);
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
    
    NumericVector xx = na_omit(x);
    NumericVector yy = na_omit(y);
    int nx = xx.size();
    int ny = yy.size();
    NumericVector bsamplex = sample(xx, nx * nboot, true);
    NumericVector bsampley = sample(yy, ny * nboot, true);
    arma::mat datax = arma::mat(bsamplex.begin(), nboot, nx, false);
    arma::mat datay = arma::mat(bsampley.begin(), nboot, ny, false);

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
    int up = nboot - low - 2;
    double pv = 1.0 * z_ex / nboot + 0.5 * z_eq / nboot;
    pv = 2 * std::min(pv, 1 - pv);
    double se = var(bvec);
    return List::create(Named("est.1") = hd(xx, q, cores),
                        Named("est.2") = hd(yy, q, cores),
                        Named("ci") = NumericVector::create(bvecsort[low], bvecsort[up]),
                        Named("p.value") = pv,
                        Named("sq.se") = se,
                        Named("n1") = nx,
                        Named("n2") = ny
                       );
}
