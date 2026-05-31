# LongSubTest
Subgroup testing in change-plane models for longitudinal data.

Change-plane testing methods including the supremum and the weighted average of score test are proposed for longitudinal data. 

# Installation

    #install Rtools 3.5 (http://cran.r-project.org/bin/windows/Rtools)
    #install.packages("devtools")
    #install.packages("Rcpp")
    library(devtools)
    install_github("PanpanRen/LongSubTest")

# Usage

- [x] [LongSubTest-manual](https://github.com/PanpanRen/LongSubTest/tree/main/inst/LongSubTest-manual.pdf) ------------ Details of the usage of the package.

# Example

    library(Matrix)
    library(mnormt)
    library(LongSubTest)
    
    n = 200
    m = 4
    p1 = 2
    p2 = 1
    p3 = 3
    alpha = rep(1, p1)
    beta  = rep(1, p2)/2
    gamma = c(1, seq(-1,1,length.out = p3-1)) 
    rho = 0
    set.seed(100)
    data = generate_data(n, m, alpha, beta, gamma, rho)
    fit <- LongSubTestST(data)
    fit$pval
    fit <- longSubTestWAST(data)
    fit$pval


# References

Subgroup testing for longitudinal data. Manuscript.

# Development
The R-package is developed by Panpan Ren (panpanren@stu.sufe.edu.cn).



