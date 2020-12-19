## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.width=7, fig.height=7) 

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

## -----------------------------------------------------------------------------
set.seed(3333)
sk <- function(x) {
  # skewness computation
  m2 <- mean((x-mean(x))^2)
  m3 <- mean((x-mean(x))^3)
  m3/(m2*sqrt(m2))
}

n <- 40
m <- 2000
alpha_test <- 0.1    # confidence level alpha
alpha <- c(seq(0,1,0.1),seq(1,10,0.5))    # parameter alpha for symmetric Beta distribution
N <- length(alpha)
powerrate <- rep(0,N)
critical_value <- qnorm(1-alpha_test/2, 0, sqrt((n-2)*6 / ((n+3)*(n+1))))

for (i in 1:N) {
  test_rej <- rep(0,m)
  for (j in 1:m) {
    x <- rbeta(n,alpha[i],alpha[i])
    test_rej[j] <- (critical_value <= abs(sk(x)))
  }
  powerrate[i] <- mean(test_rej)
}

plot(alpha,powerrate,xlab=bquote(alpha),type='b',ylim=c(0,1))
abline(h = alpha_test, lty = 3)
se <- sqrt(powerrate*(1-powerrate)/m)    # add standard errors
lines(alpha, powerrate+se, lty = 3)
lines(alpha, powerrate-se, lty = 3)

nu <- c(seq(0.1,1,0.1),seq(2,20,1))    # paramer nu for t distribution
N_t <- length(nu)
powerrate_t <- rep(0,N_t)

for (i in 1:N_t) {
  test_rej <- rep(0,m)
  for (j in 1:m) {
    x <- rt(n,nu[i])
    test_rej[j] <- (critical_value <= abs(sk(x)))
  }
  powerrate_t[i] <- mean(test_rej)
}

plot(nu,powerrate_t,xlab=bquote(nu),type='b',ylim=c(0,1))
abline(h = alpha_test, lty = 3)
se_t <- sqrt(powerrate_t*(1-powerrate_t)/m)    # add standard errors
lines(nu, powerrate_t+se_t, lty = 3)
lines(nu, powerrate_t-se_t, lty = 3)

# The result of heavy-tailed symmetric distribution (such as t-distribution) and light-tailed symmetric distribution (such as beta-distribution with alpha>1) are significantly different! Heavy-tailed symmetric distribution has higher power, while light-tailed symmetric distribution has very low power.

## -----------------------------------------------------------------------------
set.seed(3333)

count5test <- function(x, y) {
  x_central <- x - mean(x)
  y_central <- y - mean(y)
  x_outlier <- sum(x_central > max(y_central)) + sum(x_central < min(y_central))
  y_outlier <- sum(y_central > max(x_central)) + sum(y_central < min(x_central))
  # return 1 (reject) or 0 (do not reject H0)
  as.integer(max(c(x_outlier,y_outlier)) > 5)
}

n_small <- 20
n_medium <- 80
n_large <- 200
# small, medium and large sample sizes

m <- 2000
sigma1 <- 1
sigma2 <- 1.5
# generate samples under H1 to estimate power

power_s_c5 <- mean(replicate(m, expr={
x <- rnorm(n_small, 0, sigma1)
y <- rnorm(n_small, 0, sigma2)
count5test(x, y)
}))

power_m_c5 <- mean(replicate(m, expr={
x <- rnorm(n_medium, 0, sigma1)
y <- rnorm(n_medium, 0, sigma2)
count5test(x, y)
}))

power_l_c5 <- mean(replicate(m, expr={
x <- rnorm(n_large, 0, sigma1)
y <- rnorm(n_large, 0, sigma2)
count5test(x, y)
}))

power_s_f <- mean(replicate(m, expr={
x <- rnorm(n_small, 0, sigma1)
y <- rnorm(n_small, 0, sigma2)
as.integer(var.test(x,y)$p.value < 0.055)
}))

power_m_f <- mean(replicate(m, expr={
x <- rnorm(n_medium, 0, sigma1)
y <- rnorm(n_medium, 0, sigma2)
as.integer(var.test(x,y)$p.value < 0.055)
}))

power_l_f <- mean(replicate(m, expr={
x <- rnorm(n_large, 0, sigma1)
y <- rnorm(n_large, 0, sigma2)
as.integer(var.test(x,y)$p.value < 0.055)
}))

c(power_s_c5, power_s_f)
# For small sample sizes, Count Five Test power=0.3115, F-Test power=0.4220
c(power_m_c5, power_m_f)
# For medium sample sizes, Count Five Test power=0.8060, F-Test power=0.9500
c(power_l_c5, power_l_f)
# For large sample sizes, Count Five Test power=0.9550, F-Test power=1

## -----------------------------------------------------------------------------
set.seed(3333)
library(MASS)
mul_sk <- function(x) {
  # computes the Mardia multivariate sample skewness coeff.
  n <- nrow(x)
  sigma <- var(x)
  sigmainv <- solve(sigma)
  mu <- apply(x,2,mean)
  x_cen <- x-mu
  sk <- 0
  for (i in 1:n)
    for (j in 1:n)
      sk <- sk + (t(as.matrix(x_cen[i,])) %*% sigmainv %*% as.matrix((x_cen[j,])))^3
  sk/n^2
}

d <- 2    # We suppose the dimension is 2
n <- c(5, 10, 20, 30)    # sample sizes
cv <- qchisq(0.975,d*(d+1)*(d+2)/6)

# n is a vector of sample sizes
# we are doing length(n) different simulations
p.reject <- numeric(length(n))    # to store sim. results
m <- 100    # num. repl. each sim.  It should be 1000, but to shorten the vignettes creation time, I manually switch it to 100.
for (i in 1:length(n)) {
  sktests <- numeric(m)    # test decisions
  for (j in 1:m) {
    x <- mvrnorm(n[i],c(0,0),matrix(c(1,0,0,1),2,2))    # test decision is 1 (reject) or 0
    sktests[j] <- as.integer(n[i]*abs(mul_sk(x))/6 >= cv)
  }
  p.reject[i] <- mean(sktests)    # proportion rejected
}
p.reject
# when n=5,10,20,30, p.reject=0.222,0.146,0.114,0.107
# This is a repetition of Example 6.8 for multivariate case.

alpha <- .1
n <- 30
m <- 30    # It should be 100, but to shorten the vignettes creation time, I manually switch it to 30.
epsilon <- c(seq(0, .15, .05), seq(.2, 1, .2))
N <- length(epsilon)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qchisq(0.975,d*(d+1)*(d+2)/6)
for (j in 1:N) {    # for each epsilon
  e <- epsilon[j]
  sktests <- numeric(m)
  for (i in 1:m) {    # for each replicate
    x <- (1-e)*mvrnorm(n, c(0,0), matrix(c(1,0,0,1),2,2)) + e*mvrnorm(n, c(0,0), matrix(c(100,0,0,100),2,2))
    sktests[i] <- as.integer(n*abs(mul_sk(x))/6 >= cv)
  }
  pwr[j] <- mean(sktests)
}
# plot power vs epsilon
plot(epsilon, pwr, type = "b", xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m)    # add standard errors
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

# This is a repetition of Example 6.10 for multivariate case.

## -----------------------------------------------------------------------------
set.seed(3333)
library(bootstrap)
law <- as.matrix(law)
corr <- cor(law[,1],law[,2])
n <- nrow(law)
jackknife_store <- rep(0,n)
for (i in 1:n) {
  jacksample <- law[-i,]
  jackknife_store[i] <- cor(jacksample[,1],jacksample[,2])
}
jackknife_bias <- (n-1)*(mean(jackknife_store)-corr)
jackknife_se <- (n-1)*sd(jackknife_store)/sqrt(n)
c(jackknife_bias, jackknife_se)
# Bias=-0.006473623, SE=0.1425186

## -----------------------------------------------------------------------------
set.seed(3333)
library(boot)
aircondit <- as.matrix(aircondit)
boot.mean <- function(x,i) mean(x[i])
boot.obj <- boot(aircondit, statistic=boot.mean, R=2000)
print(boot.ci(boot.obj, type = c("norm","basic","perc","bca")))

## -----------------------------------------------------------------------------
set.seed(3333)
library(bootstrap)
scor <- as.matrix(scor)
n <- nrow(scor)
lambda <- eigen(cov(scor))$values
theta <- lambda[1] / (lambda[1]+lambda[2]+lambda[3]+lambda[4]+lambda[5])
jackknife_theta <- rep(0,n)
for (i in 1:n) {
  jacksample <- scor[-i,]
  lambda_j <- eigen(cov(jacksample))$values
  jackknife_theta[i] <- lambda_j[1] / (lambda_j[1]+lambda_j[2]+lambda_j[3]+lambda_j[4]+lambda_j[5])
}
jackknife_bias <- (n-1)*(mean(jackknife_theta)-theta)
jackknife_se <- (n-1)*sd(jackknife_theta)/sqrt(n)
c(jackknife_bias,jackknife_se)
# Bias=0.001069139, SE=0.049552307

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
n <- length(magnetic)
e1 <- numeric(n*(n-1)/2)
e2 <- numeric(n*(n-1)/2)
e3 <- numeric(n*(n-1)/2)
e4 <- numeric(n*(n-1)/2)
count <- 0
for (i in 1:(n-1))
  for (j in (i+1):n) {
    count <- count+1
    y <- magnetic[-c(i,j)]
    x <- chemical[-c(i,j)]
    
    P1 <- lm(y~x)
    y1_1 <- chemical[i]*P1$coef[2] + P1$coef[1]
    y1_2 <- chemical[j]*P1$coef[2] + P1$coef[1]
    e1[count] <- (magnetic[i]-y1_1)^2+(magnetic[j]-y1_2)^2
    
    P2 <- lm(y~x+I(x^2))
    y2_1 <- P2$coef[1] + P2$coef[2] * chemical[i] + P2$coef[3] * chemical[i]^2
    y2_2 <- P2$coef[1] + P2$coef[2] * chemical[j] + P2$coef[3] * chemical[j]^2
    e2[count] <- (magnetic[i]-y2_1)^2+(magnetic[j]-y2_2)^2
    
    P3 <- lm(log(y)~x)
    y3_1 <- exp(P3$coef[1] + P3$coef[2] * chemical[i])
    y3_2 <- exp(P3$coef[1] + P3$coef[2] * chemical[j])
    e3[count] <- (magnetic[i]-y3_1)^2+(magnetic[j]-y3_2)^2
    
    P4 <- lm(log(y)~log(x))
    y4_1 <- exp(P4$coef[1] + P4$coef[2] * log(chemical[i]))
    y4_2 <- exp(P4$coef[1] + P4$coef[2] * log(chemical[j]))
    e4[count] <- (magnetic[i]-y4_1)^2+(magnetic[j]-y4_2)^2
  }

c(mean(e1)/2,mean(e2)/2,mean(e3)/2,mean(e4)/2)
# 19.57227 17.87018 18.45491 20.46718
# According to the prediction error criterion, Model 2, the quadratic model, would again be the best fit for the data.

## -----------------------------------------------------------------------------
# Case without Permutation:

set.seed(123)

count5test <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
return(as.integer(max(c(outx, outy)) > 5))
}

n1 <- 20
n2 <- 30
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
m <- 1000
alphahat <- mean(replicate(m, expr={
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
x <- x - mean(x) #centered by sample mean
y <- y - mean(y)
count5test(x, y)
}))
print(alphahat)

# 0.104, suggests that the “Count Five” criterion does not necessarily control Type I error at alpha <= 0.0625 when the sample sizes are unequal.

# Now we use Permutation Method to deal with the same data:

per_countfivetest <- function(z,B) {
  n <- length(z)
  rst <- rep(0,B)
  for (i in 1:B) {
    z <- sample(z)
    rst[i] <- count5test(z[1:(n/2)],z[-(1:(n/2))])
  }
  mean(rst)
}

alphahat2 <- mean(replicate(m, expr={
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
x <- x - mean(x) #centered by sample mean
y <- y - mean(y)
per_countfivetest(c(x,y),50)
}))
print(alphahat2)

# 0.06154, this is a good Type-I error, which is much closer to 0.0625, and suggests that this is an available permutation test for equal variance that applies when sample sizes are not necessarily equal.

## -----------------------------------------------------------------------------
set.seed(333)
library(Ball)
library(energy)
library(boot)
library(RANN)

m <- 100
k <- 3
p <- 2
mu <- 0.2
R <- 99
n <- n1+n2
N <- c(n1,n2)
alpha <- 0.1
sample_size <- c(10,20,30,50,80,100)

Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) # what's the first column?
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n)
}

eqdist.nn <- function(z,sizes,k){
boot.obj <- boot(data=z,statistic=Tn,R=R,
sim = "permutation", sizes = sizes,k=k)
ts <- c(boot.obj$t0,boot.obj$t)
p.value <- mean(ts>=ts[1])
list(statistic=ts[1],p.value=p.value)
}

compare_1 <- function(sample_size) {
  n1 <- n2 <- sample_size
  n <- n1+n2
  N <- c(n1,n2)
  p.values <- matrix(NA,m,3)
  for(i in 1:m){
  x <- matrix(rnorm(n1*p,0,1.7),ncol=p);
  y <- cbind(rnorm(n2,0,1.7),rnorm(n2));
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=99,seed=i)$p.value
  }
  colMeans(p.values<alpha)
}

case1 <- matrix(0,nrow=6,ncol=3)
for (i in 1:6)
  case1[i,] <- compare_1(sample_size[i])

plot(x=sample_size, y=case1[,1], type='b', xlab='Sample size', ylab='Power', main='Unequal variances and equal expectations',col=1,lwd=2,xlim=c(5,105),ylim=c(0,1))
lines(x=sample_size, y=case1[,2], type="b", lwd=2, col=2)
lines(x=sample_size, y=case1[,3], type="b", lwd=2, col=3)
legend(x=5,y=1,legend=c('NN','Energy','Ball'),lwd=2,col=c(1,2,3))

# Unequal variances and equal expectations

## ----eval=FALSE---------------------------------------------------------------
#  
#  # To save the time of create vignette, we don't run the remaining three cases.
#  compare_2 <- function(sample_size) {
#    n1 <- n2 <- sample_size
#    n <- n1+n2
#    N <- c(n1,n2)
#    p.values <- matrix(NA,m,3)
#    for(i in 1:m){
#    x <- matrix(rnorm(n1*p,0,1.7),ncol=p);
#    y <- cbind(rnorm(n2,0,1.7),rnorm(n2,mean=mu));
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x=x,y=y,num.permutations=99,seed=i)$p.value
#    }
#    colMeans(p.values<alpha)
#  }
#  
#  case2 <- matrix(0,nrow=6,ncol=3)
#  for (i in 1:6)
#    case2[i,] <- compare_2(sample_size[i])
#  
#  plot(x=sample_size, y=case2[,1], type='b', xlab='Sample size', ylab='Power', main='Unequal variances and unequal expectations',col=1,lwd=2,xlim=c(5,105),ylim=c(0,1))
#  lines(x=sample_size, y=case2[,2], type="b", lwd=2, col=2)
#  lines(x=sample_size, y=case2[,3], type="b", lwd=2, col=3)
#  legend(x=70,y=0.4,legend=c('NN','Energy','Ball'),lwd=2,col=c(1,2,3))
#  
#  
#  # Unequal variances and unequal expectations.
#  
#  compare_3 <- function(sample_size) {
#    n1 <- n2 <- sample_size
#    n <- n1+n2
#    N <- c(n1,n2)
#    p.values <- matrix(NA,m,3)
#    for(i in 1:m){
#    x <- matrix(rt(n1*p,1),ncol=p);
#    y <- cbind(rt(n2,1),rt(n2,100));
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x=x,y=y,num.permutations=99,seed=i)$p.value
#    }
#    colMeans(p.values<alpha)
#  }
#  
#  case3 <- matrix(0,nrow=6,ncol=3)
#  for (i in 1:6)
#    case3[i,] <- compare_3(sample_size[i])
#  
#  plot(x=sample_size, y=case3[,1], type='b', xlab='Sample size', ylab='Power', main=' t-distribution',col=1,lwd=2,xlim=c(5,105),ylim=c(0,1))
#  lines(x=sample_size, y=case3[,2], type="b", lwd=2, col=2)
#  lines(x=sample_size, y=case3[,3], type="b", lwd=2, col=3)
#  legend(x=5,y=1,legend=c('NN','Energy','Ball'),lwd=2,col=c(1,2,3))
#  
#  # Non-normal distributions: t distribution with 1 df.
#  
#  compare_4 <- function(sample_size) {
#    n1 <- n2 <- sample_size
#    n <- n1+n2
#    N <- c(n1,n2)
#    p.values <- matrix(NA,m,3)
#    for(i in 1:m){
#    x <- matrix(rnorm(n1*p,(sample(10)<3)*1,1),ncol=p);
#    y <- cbind(rnorm(n2,(sample(10)<3)*1,1),rnorm(n2,(sample(10)>=3)*1,1));
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x=x,y=y,num.permutations=99,seed=i)$p.value
#    }
#    colMeans(p.values<alpha)
#  }
#  
#  case4 <- matrix(0,nrow=6,ncol=3)
#  for (i in 1:6)
#    case4[i,] <- compare_4(sample_size[i])
#  
#  plot(x=sample_size, y=case4[,1], type='b', xlab='Sample size', ylab='Power', main='mixture normal distribution',col=1,lwd=2,xlim=c(5,105),ylim=c(0,1))
#  lines(x=sample_size, y=case4[,2], type="b", lwd=2, col=2)
#  lines(x=sample_size, y=case4[,3], type="b", lwd=2, col=3)
#  legend(x=5,y=1,legend=c('NN','Energy','Ball'),lwd=2,col=c(1,2,3))
#  
#  # Non-normal distributions: mixture of two normal distributions.
#  
#  compare_5 <- function(sample_size) {
#    n1 <- n2 <- sample_size
#    n <- n1+n2
#    N <- c(n1,n2)
#    p.values <- matrix(NA,m,3)
#    for(i in 1:m){
#    x <- matrix(rnorm(n1*(p+3),0,1.7),ncol=p+3);
#    y <- cbind(rnorm(n2,0,1.7),rnorm(n2),rnorm(n2),rnorm(n2),rnorm(n2));
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x=x,y=y,num.permutations=99,seed=i)$p.value
#    }
#    colMeans(p.values<alpha)
#  }
#  
#  case5 <- matrix(0,nrow=6,ncol=3)
#  for (i in 1:6)
#    case5[i,] <- compare_5(sample_size[i])
#  
#  plot(x=sample_size, y=case5[,1], type='b', xlab='Sample size', ylab='Power', main='Unbalanced samples with the Case 1',col=1,lwd=2,xlim=c(5,105),ylim=c(0,1))
#  lines(x=sample_size, y=case5[,2], type="b", lwd=2, col=2)
#  lines(x=sample_size, y=case5[,3], type="b", lwd=2, col=3)
#  legend(x=70,y=0.4,legend=c('NN','Energy','Ball'),lwd=2,col=c(1,2,3),cex=0.8)
#  
#  # Unbalanced samples with the Case 1. (1:4)
#  
#  # From the comparison, we can find that Ball test is usually the most powerful method of these three. When dealing with mixture normal distribution, Energy test behaves better.

