## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(lattice)
n <- seq(5, 45, 5)
x <- rnorm(sum(n))
y <- factor(rep(n, n), labels=paste("n =", n))
densityplot(~ x | y,
panel = function(x, ...) {
panel.densityplot(x, col="DarkOliveGreen", ...)
panel.mathdensity(dmath=dnorm,
args=list(mean=mean(x), sd=sd(x)),
col="darkblue")
})


## -----------------------------------------------------------------------------
stdv = function(x) {  if (NROW(x)>1)   sigma <- round(sd(x),2) 
                                   else   sigma <- 0
                                   sigma }
smry <- function(X){     mu <- apply(X,2,mean)  
sigma <- apply(X,2,stdv)
c(mu=mu,sigma=sigma)}

patients <- c("0071021198307012008001400400150",
                       "0071201198307213009002000500200",
                      "0090903198306611007013700300000",
                      "0050705198307414008201300900000",
                     "0050115198208018009601402001500",
                     "0050618198207017008401400800400",
                     "0050703198306414008401400800200")
id <- substr(patients,1,3)    # First 3 digits signify the patient ID
date <- as.Date(substr(patients,4,11), format = "%m%d%Y")    # 4-11 for treatment date
hr <- as.numeric(substr(patients,12,14))    # Heart rate
sbp <- as.numeric(substr(patients,15,17))    # Systolic blood pressure
dbp <- as.numeric(substr(patients,18,20))    # Diastolic blood pressure
dx <- substr(patients,21,23)
docfee <- as.numeric(substr(patients,24,27))  
labfee <- as.numeric(substr(patients,28,31))

tapply(hr, id, mean)
tapply(hr, id, stdv)

# Show these results in a more compact way
PATIENTS <- data.frame(id, hr, sbp, dbp, docfee, labfee)
str(PATIENTS)
smry(PATIENTS[id=='005',2:6])
smry(PATIENTS[id=='007',2:6])
smry(PATIENTS[id=='009',2:6])
by(PATIENTS[2:6], id, smry) 
by(PATIENTS[2:6], id, summary)

# Calculate the difference between the first and the last observations of HR, SBP and DBP
HrSbpDbp  <- data.frame(id, date, hr, sbp, dbp)
# Sort by ID first, then sort by treatment date
HrSbpDbpSorted <- HrSbpDbp[order(HrSbpDbp$id, HrSbpDbp$date), ] 
HrSbpDbpSorted

## -----------------------------------------------------------------------------

n <- 10000
a <- 2
b <- 2
unifrv <- runif(n)
paretorv <- b/(1-unifrv)^(1/a)    # inverse transformation
hist(paretorv[paretorv>0 & paretorv<20], freq = FALSE, breaks = seq(0,20,0.5), main = "Histogram of the Pareto sample with the true density curve",xlab = "Sample value")    # graph the density histogram
# I found that F(20)=0.99, which is very close to 1. For the sake of the tidiness and better description for the feature of the histogram, I made a truncation on x=20, and only consider the variables between 0 and 20.
f <- function(x) {a*b^a/x^(a+1)}    # true pdf
curve(f, 2, 20, col = 2, lwd = 3, add = TRUE)    # add the true density curve
legend(12,0.6,"true density", col = 2, lwd = 3)    # add a legend


## -----------------------------------------------------------------------------
u1 <- runif(10000, min = -1, max = 1)    # consider a sample number of 10000
u2 <- runif(10000, min = -1, max = 1)
u3 <- runif(10000, min = -1, max = 1)
u <- ifelse((abs(u3)>abs(u2) & abs(u3)>abs(u1)), u2, u3)
hist(u, freq = FALSE, breaks = seq(-1,1,0.02), main = "Histogram with the true density curve", xlab = "Sample value")
f <- function(x) {3/4*(1-x^2)}
curve(f, -1, 1, col = 2, lwd = 3, add = TRUE)    # add the true density curve
legend(0.5,0.85,"true density", col = 2, lwd = 3, cex=0.6)    # add a legend


## -----------------------------------------------------------------------------
set.seed(3333)
monte_carlo_rep <- 100000
uniform_x <- runif(monte_carlo_rep,0,pi/3)
single_estimate <- sin(uniform_x)*pi/3
final_estimate <- mean(single_estimate)
final_estimate
#=0.5006204, very close to the exact value 0.5.

## -----------------------------------------------------------------------------
set.seed(2222)
monte_carlo_r <- 100
var_reduc_store <- rep(0,monte_carlo_r)
for (i in 1:monte_carlo_r){
  u <- runif(10000)
  simple_mc <- exp(1)^u
  antithetic <- (exp(1)^u+exp(1)^(1-u))/2
  var_reduc_store[i] <- 1- var(antithetic) / var(simple_mc)
}
mean(var_reduc_store)
# empirical estimate of the percent reduction is 98.38%, very close to the theoretical value 98.39%.

## -----------------------------------------------------------------------------
set.seed(3333)
g <- function(x) x^2*exp(-x^2/2)/sqrt(2*pi)
f1 <- function(x) sqrt(exp(1))/2*exp(-x/2)
f2 <- function(x) exp(-(x-1))
plot(g,1,5,ylim=c(0,1),ylab='y')
par(new=TRUE)
plot(f1,1,5,xlab='',ylab='',main='',col='blue',ylim=c(0,1))
par(new=TRUE)
plot(f2,1,5,xlab='',ylab='',main='',col='red',ylim=c(0,1))
legend(4,0.9,c("g","f1","f2"), col = c('black','blue','red'),lwd = c(1,1,1))

## -----------------------------------------------------------------------------
set.seed(3333)
m <- 100000
se <- rep(0,2)
theta.hat <- rep(0,2)
g <- function(x) x^2*exp(-x^2/2)/sqrt(2*pi)
f1 <- function(x) sqrt(exp(1))/2*exp(-x/2)
f2 <- function(x) exp(-(x-1))

u1 <- runif(m,0,0.5)    # to confirm that all the x1>1
x1 <- -2*log(u1 * 2 / sqrt(exp(1)))    # f1, using inverse transform method
f1g <- g(x1)/f1(x1)
theta.hat[1] <- mean(f1g)
se[1] <- sd(f1g)

x2 <- rexp(m)+1
f2g <- g(x2)/f2(x2)
theta.hat[2] <- mean(f2g)
se[2] <- sd(f2g)

theta.hat
# 0.4014080 0.4002661, both are close to the real value.
se
# 0.302011 0.157805, f2 has smaller standard error, as I expected.

## -----------------------------------------------------------------------------
set.seed(3333)
m <- 10000
g <- function(x) exp(-x - log(1+x^2))
f1 <- function(x) exp(-x)/(1-exp(-1/5))
f2 <- function(x) exp(-x)/(exp(-1/5)-exp(-2/5))
f3 <- function(x) exp(-x)/(exp(-2/5)-exp(-3/5))
f4 <- function(x) exp(-x)/(exp(-3/5)-exp(-4/5))
f5 <- function(x) exp(-x)/(exp(-4/5)-exp(-1))

u <- runif(m)
x1 <- -log(1-u*((1-exp(-1/5))))
x2 <- -log(exp(-1/5)-u*((exp(-1/5)-exp(-2/5))))
x3 <- -log(exp(-2/5)-u*((exp(-2/5)-exp(-3/5))))
x4 <- -log(exp(-3/5)-u*((exp(-3/5)-exp(-4/5))))
x5 <- -log(exp(-4/5)-u*((exp(-4/5)-exp(-1))))

fg1 <- g(x1)/f1(x1)
fg2 <- g(x2)/f2(x2)
fg3 <- g(x3)/f3(x3)
fg4 <- g(x4)/f4(x4)
fg5 <- g(x5)/f5(x5)

fg <- fg1+fg2+fg3+fg4+fg5

mean(fg)
# 0.5248477, very close to the result in Example 5.10

sd(fg)
# 0.01694816, much smaller than the result in Example 5.10


## -----------------------------------------------------------------------------
library(stats)
set.seed(3333)
n <- 100
m <- 100
l_bound <- rep(0,m)
u_bound <- rep(0,m)
for (i in 1:100) {
  x <- rlnorm(n,5,1)
  mu_hat <- sum(log(x))/n
  q <- qt(1-0.05/2,n-1)
  sigmasqu_hat <- sum((log(x)-mu_hat)^2)/n
  se_hat <- sigmasqu_hat/sqrt(n)
  l_bound[i] <- mu_hat - q*se_hat
  u_bound[i] <- mu_hat + q*se_hat
}
lb <- mean(l_bound)
ub <- mean(u_bound)
c(lb,ub)
# an empirical estimate of the confidence level is [4.808141, 5.199664], while the true value is 5.

## -----------------------------------------------------------------------------
set.seed(3333)
n <- 20
m <- 10000
lb_non <- rep(0,m)
ub_non <- rep(0,m)
lb_nor <- rep(0,m)
ub_nor <- rep(0,m)
for (i in 1:m) {
  non_normal_x <- rchisq(n,2)
  lb_non[i] <- mean(non_normal_x) - sd(non_normal_x)*qt(1-0.05/2,n-1)/sqrt(n)
  ub_non[i] <- mean(non_normal_x) + sd(non_normal_x)*qt(1-0.05/2,n-1)/sqrt(n)
  normal_x <- rnorm(n,mean=0,sd=2)
  lb_nor[i] <- mean(normal_x) - sd(normal_x)*qt(1-0.05/2,n-1)/sqrt(n)
  ub_nor[i] <- mean(normal_x) + sd(normal_x)*qt(1-0.05/2,n-1)/sqrt(n)
}
cov_prob_non <- mean(lb_non<2 & ub_non>2)
cov_prob_nor <- mean(lb_nor<0 & ub_nor>0)
c(cov_prob_non,cov_prob_nor)
# The normal samples have a coverage probability of 94.9%, which is effective. However, the non-normal samples only have a coverage probability of 91.73%, which suggests that the probability that the confidence interval covers the mean is not necessarily equal to 0.95 for non-normal samples.

