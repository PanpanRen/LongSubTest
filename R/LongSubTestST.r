LongSubTestST <- function(data, Gamma, B = 1000, K = 1000, qlb = 0.1, seed = 1) {
    tic = proc.time()

    Y = data$Y  # n*m
    X = data$X
    Z = data$Z
    U = data$U

    n = nrow(Y)
    m = ncol(Y)
    p = ncol(X)+1
    q = ncol(Z)
    r = ncol(U)+1

    if (missing(Gamma)) {
        set.seed(seed)
        cols = apply(U, 2, var) != 0
        Gamma = gam.init(K, U[,cols], lb.quantile=qlb, ub.quantile=1-qlb, ss=1)
        rm(cols)
    }
    ind = cbind(1,U)%*%t(Gamma)>=0
    set.seed(seed+1)
    xi = matrix(rnorm(n*B),n,B)
    
    fit = .Call(
        "_ST1",
        as.integer(n),
        as.integer(m),
        as.integer(p),
        as.integer(q),
        as.numeric(cbind(1,X)),
        as.numeric(t(Y)),
        as.numeric(t(Z)),
        as.numeric(ind),
        as.integer(B),
        as.integer(K),
        as.numeric(xi))
    
    toc = proc.time()
    fit$time = toc[3] - tic[3]
    
    return(fit)
}