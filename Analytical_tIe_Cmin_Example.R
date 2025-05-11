
library(mvtnorm)

## parameter settings
alpha <- 0.025
beta <- 0.1
tau <- 0.5
rho_xy <- 0.7
rho_xz <- 0.5
t_S <- 0.5
t_F <- 0.5
nmaxovern_F <- 2
nmaxovern_S <- 2

W_S <- qnorm(1-alpha)*sqrt(t_S)+qnorm(1-beta)*sqrt(t_S*(1-t_S))
W_F <- qnorm(1-alpha)*sqrt(t_F)+qnorm(1-beta)*sqrt(t_F*(1-t_F))

## a function to calculate analytical type I error given c and Delta
## the intermediate steps have been labelled correspondingly based on the equation numbers in the manuscrip
calc_tIe <- function(c, Delta){
  d <- c+Delta
  y <- c(c,d)
  
  # \psi function
  phi_bvn <- function(s,f){
    dmvnorm(c(s,f), mean=c(0,0), sigma=matrix(c(1,sqrt(tau),sqrt(tau),1), nrow=2, byrow=T))
  }
  phi_bvn_vec <- Vectorize(phi_bvn, vectorize.args = c("s","f"))
  # Equation (11)
  f1_F <- function(z){
    pnorm((z*sqrt(t_F)-qnorm(1-alpha))/sqrt(1-t_F))
  }
  f1_S <- function(z){
    pnorm((z*sqrt(t_S)-qnorm(1-alpha))/sqrt(1-t_S))
  }
  # Equation (12)
  f2_F <- function(z){
    min1 <- ifelse(nmaxovern_F < t_F+t_F/z^2*(W_F-z*t_F)^2/(t_F*(1-t_F)),
                   nmaxovern_F,
                   t_F+t_F/z^2*(W_F-z*t_F)^2/(t_F*(1-t_F)))
    min2 <- ifelse(sqrt(nmaxovern_F-t_F) < abs(sqrt(t_F)*(W_F-z*t_F)/(z*sqrt(t_F*(1-t_F)))),
                   sqrt(nmaxovern_F-t_F),
                   abs(sqrt(t_F)*(W_F-z*t_F)/(z*sqrt(t_F*(1-t_F)))))
    pnorm( (z*sqrt(t_F)-qnorm(1-alpha)*sqrt(min1))/(min2) )
  }
  f2_S <- function(z){
    min1 <- ifelse(nmaxovern_S < t_S+t_S/z^2*(W_S-z*t_S)^2/(t_S*(1-t_S)),
                   nmaxovern_S,
                   t_S+t_S/z^2*(W_S-z*t_S)^2/(t_S*(1-t_S)))
    min2 <- ifelse(sqrt(nmaxovern_S-t_S) < abs(sqrt(t_S)*(W_S-z*t_S)/(z*sqrt(t_S*(1-t_S)))),
                   sqrt(nmaxovern_S-t_S),
                   abs(sqrt(t_S)*(W_S-z*t_S)/(z*sqrt(t_S*(1-t_S)))))
    pnorm( (z*sqrt(t_S)-qnorm(1-alpha)*sqrt(min1))/(min2) )
  }
  # Equation (13)
  g_F <- function(z, y){
    as.numeric(pmvnorm(lower=c(-Inf,y[1]), upper=c(y[2],Inf),
                       mean=c(sqrt(tau)*rho_xz*z, rho_xz*z),
                       sigma=matrix(c(1-tau*rho_xz^2, sqrt(tau)*(1-rho_xz^2), sqrt(tau)*(1-rho_xz^2), 1-rho_xz^2), nrow=2, byrow=T)))
  }
  g_F_vec <- Vectorize(g_F, vectorize.args = "z")
  # Equation (14)
  g_S <- function(z, y){
    as.numeric(pmvnorm(lower=c(y[2],y[1]), upper=c(Inf,Inf),
                       mean=c(rho_xz*z, sqrt(tau)*rho_xz*z),
                       sigma=matrix(c(1-rho_xz^2, sqrt(tau)*(1-rho_xz^2), sqrt(tau)*(1-rho_xz^2), 1-tau*rho_xz^2), nrow=2, byrow=T)))
  }
  g_S_vec <- Vectorize(g_S, vectorize.args = "z")
  # Equation (5)
  part_Y_S <- adaptIntegrate(function(x) {
    pnorm((rho_xy*(d+x[1]/(1-x[1]))-qnorm(1-alpha))/sqrt(1-rho_xy^2))*phi_bvn(d+x[1]/(1-x[1]), c-(1-x[2])/x[2])*1/(1-x[1])^2/x[2]^2
  }, c(0,0), c(1,1))$integral
  # Equation (6)
  part_Y_F <- adaptIntegrate(function(x) {
    pnorm((rho_xy*(c-(1-x[2])/x[2])-qnorm(1-alpha))/sqrt(1-rho_xy^2))*phi_bvn(d-(1-x[1])/x[1], c-(1-x[2])/x[2])*1/x[1]^2/x[2]^2
  }, c(0,0), c(1,1))$integral
  # Equation (9)
  part_Z2_F <- integrate(function(p){
    f1_F(W_F+p/(1-p))*g_F_vec(z=W_F+p/(1-p), y=y)*dnorm(W_F+p/(1-p))/(1-p)^2
  },0,1)$value + integrate(function(p){                           
    f2_F(W_F-(1-p)/p)*g_F_vec(z=W_F-(1-p)/p, y=y)*dnorm(W_F-(1-p)/p)/p^2
  }, 0, 1)$value
  # Equation (10)
  part_Z2_S <- integrate(function(p){
    f1_S(W_S+p/(1-p))*g_S_vec(z=W_S+p/(1-p), y=y)*dnorm(W_S+p/(1-p))/(1-p)^2
  }, 0, 1)$value + integrate(function(p){
    f2_S(W_S-(1-p)/p)*g_S_vec(z=W_S-(1-p)/p, y=y)*dnorm(W_S-(1-p)/p)/p^2
  }, 0, 1)$value
  # Equation (15)
  tIe <- part_Y_S + part_Y_F + part_Z2_S + part_Z2_F
  tIe
}

## Example
#> calc_tIe(c=1.2, Delta=0.3)
#[1] 0.0204397
  
## Given Delta and solve Cmin
Delta <- 0
Cmin <- uniroot(function(x) {calc_tIe(c=x, Delta=Delta)-alpha},interval=c(-3,4),extendInt="yes")$root

#> Cmin
#[1] -0.04606817

