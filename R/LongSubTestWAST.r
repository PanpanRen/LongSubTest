LongSubTestWAST <- function(data, family = 1, B = 1000, seed = 1) {
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

    set.seed(seed+1)
    V = matrix(rnorm(n*B), n, B)
    
    fit = .Call(
        "_WAST1",
        as.integer(n),
        as.integer(m),
        as.integer(p),
        as.integer(q),
        as.integer(r),
        as.numeric(cbind(1,X)),
        as.numeric(t(Y)),
        as.numeric(t(Z)),
        as.numeric(cbind(1,U)),
        as.integer(family),
        as.integer(B),
        as.numeric(V))
    
    toc = proc.time()
    fit$time = toc[3] - tic[3]
    
    return(fit)
}