
library(mvtnorm)
library(MASS)
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


nsims <- 10^6
set.seed(2023)
temp <- mvrnorm(n=nsims, mu=c(0,0,0,0), Sigma=matrix(c(1, sqrt(tau), rho_xz, sqrt(tau)*rho_xz,
                                                       sqrt(tau), 1, sqrt(tau)*rho_xz, rho_xz,
                                                       rho_xz, sqrt(tau)*rho_xz, 1, sqrt(tau),
                                                       sqrt(tau)*rho_xz, rho_xz, sqrt(tau), 1),nrow=4, byrow=T))

temp <- as.data.frame(temp)
names(temp) <- c("X_S","X_F","Z1_S","Z1_F")

mean_Y_F <- rho_xy*temp$X_F
mean_Y_S <- rho_xy*temp$X_S

generate_Y <- function(mean_F, mean_S){
  rmvnorm(1, mean=c(mean_F, mean_S), sigma=matrix(c(1,sqrt(tau),sqrt(tau), 1), nrow=2, byrow=T)*(1-rho_xy^2))
}
generate_Y_vec <- Vectorize(generate_Y, vectorize.args = c("mean_F", "mean_S"))

Y <- generate_Y_vec(mean_Y_F, mean_Y_S)
temp$Y_F <- Y[1,]
temp$Y_S <- Y[2,]

temp$Z2tilde_S <- rnorm(nsims, 0, 1)
temp$Z2tilde_F <- rnorm(nsims, 0, 1)

W_S <- qnorm(1-alpha)*sqrt(t_S)+qnorm(1-beta)*sqrt(t_S*(1-t_S))
W_F <- qnorm(1-alpha)*sqrt(t_F)+qnorm(1-beta)*sqrt(t_F*(1-t_F))

tau_S <- sapply(temp$Z1_S, function(x) x^(-2)*((qnorm(1-alpha)-x*sqrt(t_S))/sqrt(1-t_S)+qnorm(1-beta))^2) # tau->n2*/n1
temp$tau_S <- ifelse(temp$Z1_S<=W_S, tau_S, (1-t_S)/t_S) #n2*/n1
tau_F <- sapply(temp$Z1_F, function(x) x^(-2)*((qnorm(1-alpha)-x*sqrt(t_F))/sqrt(1-t_F)+qnorm(1-beta))^2) # tau->n2*/n1
temp$tau_F <- ifelse(temp$Z1_F<=W_F, tau_F, (1-t_F)/t_F)

if (!is.infinite(nmaxovern_S)){
  taumax_S <- nmaxovern_S/t_S-1
  temp$tau_S <- ifelse(temp$tau_S<taumax_S, temp$tau_S, taumax_S)
}

if (!is.infinite(nmaxovern_F)){
  taumax_F <- nmaxovern_F/t_F-1
  temp$tau_F <- ifelse(temp$tau_F<taumax_F, temp$tau_F, taumax_F)
}

temp$Z2_S <- temp$Z1_S*sqrt(1/(1+temp$tau_S))+temp$Z2tilde_S*sqrt(temp$tau_S/(1+temp$tau_S))
temp$Z2_F <- temp$Z1_F*sqrt(1/(1+temp$tau_F))+temp$Z2tilde_F*sqrt(temp$tau_F/(1+temp$tau_F))


empirical_tIe <- function(c, Delta, temp){
  d <- c+Delta
  (sum(temp$X_S>d & temp$X_F<c & temp$Y_S>qnorm(1-alpha)) + 
     sum(temp$X_S<d & temp$X_F<c & temp$Y_F>qnorm(1-alpha)) + 
     sum(temp$X_S>d & temp$X_F>c & temp$Z2_S>qnorm(1-alpha)) + 
     sum(temp$X_S<d & temp$X_F>c & temp$Z2_F>qnorm(1-alpha)))/
    dim(temp)[1]
}

## Example
#> empirical_tIe(c=1.2, Delta=0.3, temp=temp)
#[1] 0.020426

## Given Delta and solve Cmin
Delta <- 0
Cmin <- uniroot(function(x) {empirical_tIe(c=x, Delta=Delta, temp=temp)-alpha},interval=c(-3,2),extendInt="yes")$root
#> Cmin
#[1] -0.07405802

