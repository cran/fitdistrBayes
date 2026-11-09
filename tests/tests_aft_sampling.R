library(fitdistrBayes)
E <- asNamespace("fitdistrBayes")
assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)

# Validate the centering against independent gamma quantiles and tail limits.
for (prob in c(.05, .1, .3, .5)) {
  for (q in c(-5, -1, -.01, .01, 1, 5)) {
    a <- q^-2
    ref <- (log(qgamma(if (q > 0) prob else 1-prob, a))-log(a))/(q*sqrt(1+q*q))
    assert(abs(E$.aft_quantile_offset(q,prob)-ref)<1e-9,"Gamma quantile centering")
  }
  for(q in c(-1e-5,0,1e-5))
    assert(abs(E$.aft_quantile_offset(q,prob)-qnorm(prob))<2e-5,"Lognormal centering limit")
  assert(abs(E$.aft_quantile_offset(1e100,prob)-log(prob))<1e-12,"Positive shape centering tail")
  assert(abs(E$.aft_quantile_offset(-1e100,prob)+log1p(-prob))<1e-12,"Negative shape centering tail")
}

# Intercept translation is bijective with unit determinant. The remaining
# log-scale/shape Jacobian is sigma*cosh(zeta), as in the original target.
v <- c(.4,-.7,.2,log(1.3),asinh(-2))
prob <- .1; p <- 3L; d <- p+2L
map <- function(v) E$.aft_physical(E$.aft_uncenter(v,prob),p)
J <- vapply(seq_len(d),function(j){
  h <- rep(0,d); h[j] <- 1e-5
  (map(v+h)-map(v-h))/(2e-5)
},numeric(d))
assert(abs(det(J)/(map(v)[p+1L]*cosh(v[d]))-1)<1e-7,"Centering Jacobian")

# Evaluate the transformed posterior against a direct physical-scale formula.
y<-log(c(.4,1,2,3,4,5));status<-c(1L,1L,0L,1L,0L,1L)
X<-cbind(1,c(-1,0,1,-.5,.5,2),c(0,1,0,1,0,1))
th<-map(v)
for(pr in c("reference","simple")){
  prior<-E$.aft_prior(pr)
  lp<-if(pr=="reference")E$.aft_logomega(th[d]) else -log1p(th[d]^2)
  expected<-sum(E$.aft_loglik(th,y,status,X))+lp-log(th[p+1L])+log(abs(det(J)))
  actual<-E$.aft_logtarget(E$.aft_uncenter(v,prob),y,status,X,prior)
  assert(abs(expected-actual)<1e-7,"Unchanged posterior including all Jacobians")
}

# Independent analytic target for the full mixed-basis transition.
mu<-c(.2,-.3,.4); V<-matrix(c(.5,.25,.1,.25,.8,-.15,.1,-.15,.6),3,3)
H<-solve(V);target<-function(v){x<-v-mu;-.5*sum(x*as.vector(H%*%x))}
mode<-list(center=mu,V=V)
set.seed(824)
a<-E$.aft_sample(target,mode,1L,1200L,200L,1L,4L,E$.aft_control(list()),
                c("beta:(Intercept)","sigma","Q"),.3)
z<-do.call(rbind,a$chains)
z[,3]<-asinh(z[,3]);z[,2]<-log(z[,2])+log(cosh(z[,3]))
z[,1]<-z[,1]+exp(z[,2])*vapply(z[,3],function(x)E$.aft_quantile_offset(sinh(x),.3),0.)
assert(max(abs(colMeans(z)-mu))<.08,"Slice means versus Gaussian target")
assert(max(abs(cov(z)-V))<.1,"Slice covariance versus Gaussian target")
assert(all(a$movement>0) && all(a$evaluations>1200*3),"Sampling diagnostics record work")
assert(length(a$starts)==4L,"Independent chain starts are retained")
cat("AFT centering, Jacobians and mixed-basis Gaussian-target checks passed.\n")

# Independent Gaussian limit for the quantile MCSE approximation.
set.seed(436)
normal_chains<-lapply(1:4,function(k)matrix(rnorm(3000),ncol=1,dimnames=list(NULL,"x")))
prec<-E$.aft_quantile_precision(normal_chains)
analytic<-sqrt(c(.025,.5,.975)*c(.975,.5,.025)/12000)/dnorm(qnorm(c(.025,.5,.975)))
estimated<-unlist(prec[,c("mcse_q025","mcse_median","mcse_q975")])
assert(all(estimated/analytic>.6&estimated/analytic<1.5),"Gaussian quantile MCSE versus analytic limit")
stuck<-normal_chains;stuck[[1]][,1]<-0
assert(is.na(E$.aft_quantile_precision(stuck)$mcse_q025),"Stuck chain does not certify quantile precision")
assert(all(is.finite(E$.aft_quantile_precision(a$chains)$mcse_median)),"Sampler quantile precision")
