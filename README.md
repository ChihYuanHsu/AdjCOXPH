CoxTEL
================
Chih-Yuan Hsu

Aug 15, 2023

## How to use CoxTEL
### Installation
Download CoxTEL_1.2.0.tar.gz and locally install it, or execute the following code:
``` r
library(devtools)
install_github("cyhsuTN/CoxTEL")
```

![Figure 1: (Demo, simulated data) Reported KM curves and HR with 95% CI](https://github.com/cyhsuTN/CoxTEL/blob/master/Demo_Fig.png)

Figure 1: (Demo, simulated data) Reported KM curves and HR with 95% CI

### First step: manually extract probabilities from KM curves

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

CoxTEL HR: 1.185 (0.992-1.417) for short-term survivors; CoxTEL DP: 0.080 (0.038-0.123) for long-term survivors.

Proportion.LTS: the estimates for the proportions of long-term survivors in both arms.

### References
C.Y. Hsu, E.P.Y. Lin, Y. Shyr. (2021). Development and Evaluation of a Method to Correct Misinterpretation of Clinical Trial Results With Long-term Survival. JAMA Oncol. doi:10.1001/jamaoncol.2021.0289.

E.P.Y. Lin, C.Y. Hsu, J.F. Chiou et al. (2022). Cox Proportional Hazard Ratios Overestimate Survival Benefit of Immune Checkpoint Inhibitors (ICI): Cox-TEL Adjustment and Meta-analyses of PD-L1 Expression and ICI Survival Benefit. J Thorac Oncol. 17, 1365-1374. doi:10.1016/j.jtho.2022.08.010.

E.P.Y. Lin, C.Y. Hsu, L. Berry, P. Bunn, Y. Shyr (2022). Analysis of Cancer Survival associated with Immune Checkpoint Inhibitors after Statistical Adjustment: A Systematic Review and Meta-analyses. JAMA Network Open. 5(8):e2227211. doi:10.1001/jamanetworkopen.2022.27211.

