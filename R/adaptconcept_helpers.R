is.in.lims <- function(xx,lims) {
  all(xx >= lims[,1], xx <= lims[,2])
}
msfunc <- function(func1,lims,pow=1L,batch=F, n=1e3) {#browser()
  # Find mean square of function over limits using grid sample
  #X1 <- simple.grid(10,nrow(lims),scaledto=lims) # Too big in high dimensions, switching to just random points
  d <- nrow(lims)
  X1 <- simple.random(n=n, d=d, scaledto=lims)
  if(batch) {return(mean(func1(X1)^pow))}
  mean(apply(X1,1,func1)^pow)
}
maxgridfunc <- function(func1,lims,batch=F) {
  # Find max of function over limits using grid sample
  X1 <- simple.grid(10,nrow(lims),scaledto=lims)
  if(batch) {return(max(func1(X1)))}
  max(apply(X1,1,func1))
}

#' @importFrom stats runif
msecalc <- function(truefunc, guessfunc,lims, n=500) {
  #X1 <- simple.grid(20,nrow(lims),scaledto=lims)
  #X1 <- lhs::maximinLHS(n, nrow(lims))
  d <- nrow(lims)
  X1 <- matrix(runif(n*d), n, d)
  mean((apply(X1,1,function(xx){truefunc(xx) - guessfunc(xx)}))^2)
}
outer.inttoind <- function(i,a) {
  i <- i-1
  1+sapply(1:length(a),function(j){ifelse(j==1,i%%a[j],i%/%prod(a[1:(j-1)])%%a[j])})
}
outer.d1n <- function(...,func) {
  a <- c(...)
  b <- array(1:prod(a),dim=a)
  apply(b,1:length(a),function(xx){func(outer.inttoind(xx,a))})
}
maxN <- function(x, N=2,all.indices=F){
  # Find second max
  # all.indices will give order of N top indices
  len <- length(x)
  if(N>len){
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  if(all.indices) {return(order(x)[len:(len-N+1)])}
  sort(x,partial=len-N+1)[len-N+1]
}


#' @importFrom stats runif
simple.random <- function(n,d,scaledto=NULL) {
  m <- matrix(runif(n*d), ncol=d, nrow=n)
  if(!is.null(scaledto)) {#browser()
    m <- m * matrix(scaledto[,2]-scaledto[,1],nrow=nrow(m),ncol=ncol(m),byrow=T) + matrix(scaledto[,1],nrow=nrow(m),ncol=ncol(m),byrow=T)
  }
  m
}


#' @importFrom stats runif
simple.LHS <- function(n,d,scaled=TRUE,centered=FALSE) {
  m <- matrix(rep(1:n,d),n,d)
  m <- apply(m,2,function(xx){sample(xx)})
  if(scaled) m <- (m - runif(n*d) ) / n
  if(centered) m <- m - ifelse(scaled,.5,n/2+.5)
  m
}
is.LHS <- function(m,scaled=TRUE) {
  if(is.matrix(m)) {
    return(apply(m,2,function(mm){is.LHS(m=mm,scaled=scaled)}))
  } else {
    if(scaled) m <- ceiling(m*length(m))
    return(all(sort(m) == 1:length(m)))
  }
  stop('Error 0239357')
}
is.OA <- function(m,strength=2) {

}

#' @importFrom stats setNames
simple.grid <- function(n,d,scaled=TRUE,random=TRUE,centered=FALSE,scaledto=NULL) {
  m <- matrix(1:n,ncol=1)
  if (d>1) {
    for(di in 2:d) {
      m1 <- cbind(m,1)
      for(j in 2:n) {
        m1 <- rbind(m1,cbind(m,setNames(j,NULL)))
      }
      m <- m1
    }
  }
  if(random) m <- m - runif(n^d*d)
  if(scaled) m <- (m - ifelse(random,0,.5)) / n
  if(!is.null(scaledto)) {#browser()
    m <- m * matrix(scaledto[,2]-scaledto[,1],nrow=nrow(m),ncol=ncol(m),byrow=T) + matrix(scaledto[,1],nrow=nrow(m),ncol=ncol(m),byrow=T)
  }
  if(centered) m <- m - ifelse(scaled,.5,n/2+.5)
  m
}
