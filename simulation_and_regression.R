## ------------------------------------------------------------------------
#Variables which govern the size of the simulation
# And the size of our causal effects
nclass <- 5
nstudent <- 25
Eff <- 5
EffSD <- 3
# Simulate data
set.seed(1977)
Yr1ClassType <- rep(c(1,0),nclass*nstudent)
Yr2ClassType <- sample(Yr1ClassType,replace=FALSE)
Yr1Score <- rnorm(2*nclass*nstudent,76+Yr1ClassType*5,9)
# Fixed margins randomization
Trt  <- sample(Yr1ClassType,replace=FALSE)

## ------------------------------------------------------------------------
#There is an independent effect of class type each year
# Variance is different across class types in year 2
CtlOutcome <- rnorm(2*nclass*nstudent,Yr1Score+
                      Yr2ClassType*3,9-Yr2ClassType*4)
# Treatment effect is random, but with expectation Eff
Yr2Obs <- CtlOutcome + 
  Trt * rnorm(2*nclass*nstudent,Eff,EffSD)

summary(lm(Yr2Obs~Trt))$coefficients[2,]
summary(lm(Yr2Obs~Trt+Yr1Score))$coefficients[2,]

## ------------------------------------------------------------------------
pdf("scatter.pdf", 3.5,3.5)
plot(jitter(Trt),Yr2Obs,axes=F,xlab="Treatment",ylab="Test Result (Yr 2)",col="grey")
axis(2)
axis(1,at=c(0,1))
# Calculate quantities for plotting CIs
mns <- tapply(Yr2Obs,Trt,mean)
# SEs could also be pulled from the linear models we fit above with:
ses <- tapply(Yr2Obs,Trt,function(x) sd(x)/sqrt(length(x)))
points(c(0,1),mns,col="red",pch=19)
# Note the loop so that I only write this code once
for(tr in unique(Trt)) {
  for(q in c(.25,.025)) {
    upr<-mns[as.character(tr)]+qnorm(1-q)*ses[as.character(tr)]
    lwr <- mns[as.character(tr)]-qnorm(1-q)*ses[as.character(tr)]
    segments(tr,upr,tr,lwr,lwd=(-4/log(q)))
  }
}
dev.off()

## ----2-resid-it, fig.cap='',fig.width=12,fig.height=3.75-----------------
OutcomeRes <- residuals(lm(Yr2Obs~Yr1Score+0))
TrtRes <- residuals(lm(Trt~Yr1Score+0))
# Diagnostics
par(mfrow=c(1,4))
plot(Yr1Score,Yr2Obs)
plot(Yr1Score,jitter(Trt))
hist(OutcomeRes)
hist(TrtRes)

## ------------------------------------------------------------------------
pdf("residualized.pdf", 3.5,3.5)
par(mfrow=c(1,1))
plot(jitter(TrtRes),OutcomeRes,axes=F,xlab="Treatment (residuals)",ylab="Test Result (residuals)",col="grey")
axis(2)
axis(1)
# Pull information from the new bivariate model
mns<-coef(lm(OutcomeRes~TrtRes))
ses<-summary(lm(OutcomeRes~TrtRes))$coefficients[,2]
TrtResMns<-tapply(TrtRes,Trt,mean)
names(ses)<-names(mns)<-names(TrtResMns)
points(TrtResMns,mns,col="red",pch=19)
for(tr in names(TrtResMns)) {
  for(q in c(.25,.025)) {
    upr<-mns[tr]+qnorm(1-q)*ses[tr]
    lwr <- mns[tr]-qnorm(1-q)*ses[tr]
    segments(TrtResMns[tr],upr,TrtResMns[tr],lwr,lwd=(-4/log(q)))
  }
}
dev.off()

## ----2-coef-plot-code----------------------------------------------------
ests <- c(coef(lm(Yr2Obs~Trt))[2],coef(lm(Yr2Obs~Trt+Yr1Score))[2])
ses <- c(summary(lm(Yr2Obs~Trt))$coefficients[2,2],summary(lm(Yr2Obs~Trt+Yr1Score))$coefficients[2,2])
var.names <- c("Unadjusted","Adjusted")

## ----2-coef-plots,eval=FALSE---------------------------------------------
## par(
##   family = "serif",
##   oma = c(0,0,0,0),
##   mar = c(5,10,4,2)
## )
## 
## plot(NULL,
##   xlim = c(-0.2, 8),
##   ylim = c(.7, length(ests) + .3),
##   axes = F, xlab = NA, ylab = NA)
## 
## for (i in 1:length(ests)) {
##   points(ests[i], i, pch = 19, cex = .5)
##   lines(c(ests[i] + 1.64*ses[i], ests[i] - 1.64*ses[i]), c(i, i))
##   lines(c(ests[i] + .67*ses[i], ests[i] - .67*ses[i]), c(i, i), lwd = 3)
##   text(-1.1, i, var.names[i], xpd = T, cex = .8,pos=2)
## }
## 
## axis(side = 1)
## abline(v = 0, lty = 3, col = "black")
## mtext(side = 1, "Estimated Effect", line = 3)
## mtext(side = 3, "Adjusted vs Unadjusted Regression", line = 1)
## box()

## ----2-coef-plots-show,echo=FALSE,fig.cap='',fig.width=5,fig.height=4----
par(
  family = "serif",
  oma = c(0,0,0,0),
  mar = c(5,10,4,2)
)
  
plot(NULL,
  xlim = c(-0.2, 8),
  ylim = c(.7, length(ests) + .3),
  axes = F, xlab = NA, ylab = NA)
  
for (i in 1:length(ests)) {
  points(ests[i], i, pch = 19, cex = .5)
  lines(c(ests[i] + 1.64*ses[i], ests[i] - 1.64*ses[i]), c(i, i))
  lines(c(ests[i] + .67*ses[i], ests[i] - .67*ses[i]), c(i, i), lwd = 3)
  text(-1.1, i, var.names[i], xpd = T, cex = .8,pos=2)
}

axis(side = 1)
abline(v = 0, lty = 3, col = "black")
abline(v=Eff, lty = 5, col = "blue")
mtext(side = 1, "Estimated Effect", line = 3)
mtext(side = 3, "Adjusted vs Unadjusted Regression", line = 1)
box()    

## ------------------------------------------------------------------------
rm(list=ls())

set.seed(20140714)
N <- 2000
N.treated <- 1000
Replications <- 10000

true.treatment.effect <- 1

# Create pre-treatment covariates
owns.id.card <- rbinom(n = N, size = 1, prob = .18)
has.formal.schooling <- rbinom(n = N, size = 1, prob = .6)
age <- round(rnorm(n = N, mean = 37, sd = 16))
age[age<18] <- 18
age[age>65] <- 65
TV.access <- rbinom(n = N, size = 1, prob = .7)
epsilon <- rnorm(n = N, mean = 0, sd = 2)

## ------------------------------------------------------------------------
# Create potential outcomes correlated 
#with pre-treatment covariates
Y0 <- round(owns.id.card + 2*has.formal.schooling + 
             3*TV.access + log(age) + epsilon)
Y1 <- Y0 + true.treatment.effect

# Assign treatment repeatedly
Z.mat <- replicate(Replications, 
                  ifelse(1:N %in% 
                           sample(1:N, N.treated), 1, 0))

# %in% is a logical vector which indicates whether a 
#match was located for vector1 in vector2





## ------------------------------------------------------------------------
# Generate observed outcomes
Y.mat <- Y1 * Z.mat + Y0 * (1 - Z.mat)

diff.in.means <- function(Y, Z) {
  coef(lm(Y ~ Z))[2]
}
ols.adjust <- function(Y, Z) {
  coef(lm(Y ~ Z + owns.id.card + 
            has.formal.schooling + age + TV.access))[2]
}
unadjusted.estimates <- rep(NA, Replications)
adjusted.estimates   <- rep(NA, Replications)

## ------------------------------------------------------------------------
for (i in 1:Replications) {
  unadjusted.estimates[i]  =  
    diff.in.means(Y.mat[,i], Z.mat[,i])
  adjusted.estimates[i]    =  
    ols.adjust(Y.mat[,i], Z.mat[,i])
}
# Estimated variability 
sd.of.unadj <- sd(unadjusted.estimates)
sd.of.unadj
sd.of.adj   <- sd(adjusted.estimates)
sd.of.adj

## ------------------------------------------------------------------------

# Estimated bias of each estimator
mean(unadjusted.estimates) - true.treatment.effect
mean(adjusted.estimates) - true.treatment.effect

# Confidence interval for the bias
1.96 * sd.of.unadj / sqrt(Replications)
1.96 * sd.of.adj   / sqrt(Replications)

