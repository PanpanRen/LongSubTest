generate_data <- function(n, m, alpha, beta, gamma, rho){
    p = length(alpha)-1
    q = length(beta)
    r = length(gamma)-1

    Y = matrix(NA,n,m)
    U = matrix(NA,n,r)

    rv = matrix(rnorm(n*m*p),ncol=p)
    X = matrix(rv[,1:p],n*m,p)
    Z = matrix(rv[,1:q],n*m,q)
    U = matrix(rnorm(n * r), ncol = r)
    etaU = U%*%gamma[-1]
    threshold = quantile(etaU, probs = 0.6)
    etaU = (etaU>threshold)

    varcov = diag(m)
    eps = rmnorm(n,mean = rep(0,m),varcov = outer(1:m,1:m,FUN = function(x,y) rho^(abs(x - y))))  # n*m
    # print(solve(cov(eps)))
    # kkk

    Y = alpha[1] + t(matrix(X%*%alpha[-1],m,n)) + t(matrix(Z%*%beta,m,n)%*%diag(as.vector(etaU))) + eps

    # for(i in 1:n){
    #     Y[i,] = alpha[1] + X[(i*m-m+1):(i*m),]%*%alpha[-1] + Z[(i*m-m+1):(i*m),]*beta*etaU[i]  + eps[i,]
    # }

    return(list(Y=Y,X=X,Z=Z,U=U))
}