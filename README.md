# WRShd [![Build Status](https://travis-ci.org/mkanai/WRShd.svg?branch=master)](https://travis-ci.org/mkanai/WRShd)
An Experimental Refactor of the Harrell-Davis Estimate Functions from R. R. Wilcox' Robust Statistics Package `WRS` https://github.com/nicebread/WRS

## Motivation
* R. R. Wilcox' robust statistics functions and its R package `WRS` are definitely great work, yet they are a bit too large to maintain, depending on lots of packages.
* Bootstrapping requires substantial iterations, which can be speeded up by C++ sub-routines.
* Current way of loading `WRScpp` package is truly platform-specific; The repository `WRScpp` only provides a Mac binary, and other platform users should install the appropriate binaries provided in other repositories (e.g. `WRScppLin64` by [Joe Johnston](https://github.com/JoeJohnston/WRScppLin64)).


## Installation
```{r}
library(devtools)
install_github("mkanai/WRShd")
```

## Examples
```{r}
> library(WRShd)
> set.seed(2)
> x <- rnorm(100)
> y <- rnorm(100)
> hd(x, q = .75)
[1] 0.8353369
> hdci(x, q = .75)
$ci
[1] 0.4135639 1.2571100

$crit
[1] 2.120138

$se
[1] 0.1989366

> hdpb(x, q = .75)
$ci
[1] 0.4556359 1.2127681

$n
[1] 100

$estimate
[1] 0.8353369

$p.value
[1] 0

> hdseb(x, q = .75)
[1] 0.1989366
> qcomhd(x, y, cores = 4)
     q  n1  n2      est.1      est.2 est.1_minus_est.2      ci.low     ci.up     p_crit p.value signif
1 0.10 100 100 -1.5281745 -1.2540107       -0.27416376 -0.71975393 0.1630491 0.01250000   0.215     NO
2 0.25 100 100 -0.8650431 -0.7889496       -0.07609355 -0.41089472 0.2110202 0.02500000   0.593     NO
3 0.50 100 100 -0.1392823  0.1286732       -0.26795546 -0.67572195 0.2126300 0.01666667   0.280     NO
4 0.75 100 100  0.8353369  0.7689238        0.06641318 -0.39215194 0.4737672 0.05000000   0.818     NO
5 0.90 100 100  1.7388228  1.2116368        0.52718604 -0.04679927 0.8050822 0.01000000   0.059     NO

```

## Benchmark

### setup
```{r}
> library(rbenchmark)
> library(WRShd)
> set.seed(2)
> options(mc.cores = 4)
> x <- rnorm(100)
> y <- rnorm(100)
```

### `hd`
```{r}
> benchmark(WRS::hd(x), WRShd::hd(x), WRShd::hd(x, cores = 4), order = "relative")
                     test replications elapsed relative user.self sys.self user.child sys.child
3 WRShd::hd(x, cores = 4)          100   0.009    1.000     0.033        0          0         0
2            WRShd::hd(x)          100   0.017    1.889     0.017        0          0         0
1              WRS::hd(x)          100   0.052    5.778     0.052        0          0         0
```

### `hdci`
```{r}
> benchmark(WRS::hdci(x), WRShd::hdci(x), WRShd::hdci(x, cores = 4), order = "relative")
                       test replications elapsed relative user.self sys.self user.child sys.child
3 WRShd::hdci(x, cores = 4)          100   0.353    1.000     1.406    0.000          0         0
2            WRShd::hdci(x)          100   1.247    3.533     1.241    0.005          0         0
1              WRS::hdci(x)          100   5.333   15.108     5.323    0.002          0         0
```

### `hdpb`
```{r}
> benchmark(WRS::hdpb(x), WRShd::hdpb(x), WRShd::hdpb(x, cores = 4), order = "relative")
                       test replications elapsed relative user.self sys.self user.child sys.child
3 WRShd::hdpb(x, cores = 4)          100   6.434    1.000    25.615    0.012          0         0
2            WRShd::hdpb(x)          100  24.248    3.769    24.193    0.017          0         0
1              WRS::hdpb(x)          100 109.325   16.992   109.078    0.082          0         0
```

### `hdseb`
```{r}
> benchmark(WRS::hdseb(x), WRShd::hdseb(x), WRShd::hdseb(x, cores = 4), order = "relative")
                        test replications elapsed relative user.self sys.self user.child sys.child
3 WRShd::hdseb(x, cores = 4)          100   0.327    1.000     1.297    0.003          0         0
2            WRShd::hdseb(x)          100   1.220    3.731     1.215    0.003          0         0
1              WRS::hdseb(x)          100   5.214   15.945     5.205    0.001          0         0
```

### `qcomhd` and `qcomhdMC`
```{r}
> benchmark(WRS::qcomhd(x, y), WRS::qcomhdMC(x, y), WRShd::qcomhd(x, y), WRShd::qcomhd(x, y, cores = 4), order = "relative")
                            test replications  elapsed relative user.self sys.self user.child sys.child
4 WRShd::qcomhd(x, y, cores = 4)          100   90.841    1.000   359.966    0.260      0.007     0.021
3            WRShd::qcomhd(x, y)          100  331.000    3.644   328.542    0.907      0.014     0.042
2            WRS::qcomhdMC(x, y)          100  442.429    4.870   116.586   22.796   1107.782   144.736
1              WRS::qcomhd(x, y)          100 1182.963   13.022  1178.741    0.709      0.001     0.005
```


## Credits
* The original R script containing over 1,100 robust statistics functions are written by Rand R. Wilcox. http://dornsife.usc.edu/labs/rwilcox/software/
* The R package `WRS` is maintained by Felix Schönbrodt. https://github.com/nicebread/WRS
* The R package `WRScpp` which provides several C++ sub-routines for `WRS` is written by Xiao He. https://github.com/mrxiaohe/WRScpp

