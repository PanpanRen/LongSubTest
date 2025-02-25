gam.init.sub = function(q,n.initials){
    out = matrix(rnorm(n.initials*q), n.initials, q)
    dd 	= apply(out^2,1,sum)
    out	= out/sqrt(dd)
    return(out)
}

gam.init = function(n.initials, Z, lb.quantile, ub.quantile, ss=1){
    q = ifelse(is.matrix(Z), ncol(Z), 1)
    if(q==1){
        gamma.initials = matrix(1,n.initials,q+1)
        gamma.initials[,1] = -quantile(Z,seq(lb.quantile,ub.quantile,length=n.initials))
    }else{
        gamma.initials = gam.init.sub(q, n.initials/ss)
        Z.gamma.initials = Z %*% t(gamma.initials)
        ll=round(n.initials/ss)
        qtile = sample(seq(lb.quantile,ub.quantile,length=n.initials),n.initials)
        gamma.initials.1 = sapply(1:n.initials,function(x)return(
                            -quantile(Z.gamma.initials[,x-floor((x-0.1)/ll)*ll],qtile[x])
                        ))

        gamma.initials.1=ifelse(gamma.initials.1==(-1)*apply(Z.gamma.initials,2,min),gamma.initials.1-0.001,gamma.initials.1)
        gamma.initials.1=ifelse(gamma.initials.1==(-1)*apply(Z.gamma.initials,2,max),gamma.initials.1+0.001,gamma.initials.1)
        gamma.initials.aug=do.call("rbind", rep(list(gamma.initials), ss))
        gamma.initials = cbind(gamma.initials.1, gamma.initials.aug)
    }
    return(gamma.initials)
}