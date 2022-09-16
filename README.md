Introduction to CoxTEL
================
Chih-Yuan Hsu
Sept 16, 2022

## How to use CoxTEL

![Figure 1: (Demo, simulated data) Reported KM curves and HR with 95% CI](https://https://github.com/cyhsuTN/CoxTEL/blob/master/Demo_Fig.png)

### First step: manually extract probabilities from KM curves

``` r
extract_prob # extracted probabilities
```

    ##                6mo 12mo 18mo 24mo 30mo 36mo 42mo 48mo cured prop.
    ## Arm1 (blue, %)  34   23   19   14   13   12   11   11          11
    ## Arm0 (red, %)   41   22   12    8    6    4    3    3           3

### Second step: input the extracted probabilities and the reported Cox HR as well as CI

``` r
library(CoxTEL)
s1mix.chosen <- c(34, 23, 19, 14, 13, 12, 11, 11)/100 # extracted prob. in treatment arm
s0mix.chosen <- c(41, 22, 12,  8,  6,  4,  3,  3)/100 # extracted prob. in control arm
pi1.est <- 1 - 0.11 # uncured proportion in treatment arm
pi0.est <- 1 - 0.03 # uncured proportion in control arm
HR_cox <- 0.92 # reported Cox HR
HR_cox_CI <- c(0.77, 1.10) # reported Cox HR CI

adjustment(HR_cox, HR_cox_CI, s1mix.chosen, s0mix.chosen, pi1.est, pi0.est)
```

    ## $Adj.before.after
    ##        Cox_HR    Cox_HR_CIL    Cox_HR_CIU     CoxTEL_HR CoxTEL_HR_CIL 
    ##         0.920         0.770         1.100         1.185         0.992 
    ## CoxTEL_HR_CIU     CoxTEL_DP CoxTEL_DP_CIL CoxTEL_DP_CIU 
    ##         1.417         0.080         0.038         0.123 
    ## 
    ## $Proportion.LTS
    ##     Arm0 Arm0_CIL Arm0_CIU     Arm1 Arm1_CIL Arm1_CIU 
    ##    0.030    0.010    0.042    0.110    0.080    0.133 
    ## 
    ## $s1mix.chosen
    ## [1] 0.34 0.23 0.19 0.14 0.13 0.12 0.11 0.11
    ## 
    ## $s0mix.chosen
    ## [1] 0.41 0.22 0.12 0.08 0.06 0.04 0.03 0.03

CoxTEL HR: 1.185 (0.992-1.417); CoxTEL DP: 0.080 (0.038-0.123)


