##################################################################################################
## This is an example of running the functions in ILCM to Simulated Data 


## Set working Directory and run the R code containing the function, along with dependent packages
setwd("C:/Renogram Project/Project Long2")
source("ILCM.R")

if ( !require(plotrix) )
{
  install.packages("plotrix")
  library(plotrix)
}

if ( !require(fda) )
{
  install.packages("fda")
  library(fda)
}





####### Load Simulated Data ###########
load("simulated_data.Rdata")

## Training Data
y1_train <- simulated_data$y1_train
y2_train <- simulated_data$y2_train
x_train <- simulated_data$x_train
x1_train <- simulated_data$x1_train
x2_train <- simulated_data$x2_train
xs_train <- simulated_data$xs_train
w_train <- simulated_data$w_train
C_train <- simulated_data$C_train
O_train <- simulated_data$O_train
sbj_train <- simulated_data$sbj_train
t1 <- simulated_data$t1
t2 <- simulated_data$t2

## testing data
y1_test <- simulated_data$y1_test
y2_test <- simulated_data$y2_test
x_test <- simulated_data$x_test
x1_test <- simulated_data$x1_test
x2_test <- simulated_data$x2_test
xs_test <- simulated_data$xs_test
w_test <- simulated_data$w_test
C_test <- simulated_data$C_test
O_test <- simulated_data$O_test
sbj_test <- simulated_data$sbj_test




######  Plot of Renogram Curves - Figure 1 ######
par(mfrow=c(1,2))
plot(t1,y1_train[,1],ylim=c(0,25),type='l',ylab="MAG3/1000",xlab="Time Since Injection (secs)",
     main="Baseline Renogram")
for (i in 2:180) lines(t1,y1_train[,i],col=i)
plot(t2,y2_train[,1],ylim=c(0,40),type='l',ylab="MAG3/1000",xlab="Time Since Injection (secs)",
     main="Post-furosemide Renogram")
for (i in 2:180) lines(t2,y2_train[,i],col=i)






############################# TRAIN the Proposed Algorithm ################################
##Set up hyperparameters
nsample <- 30000
p <- 15
nu <- 2
psi <- 1/p
q <- c(6,6)
ups <- 5
a <- c(5,5)
a_sigma <- 2
b_sigma <- 2
a_zeta <- 2
b_zeta <- 2
a_iota <- 2
b_iota <- 2
a_xi <- 5 #same as tau to control variability in expert ratings
b_xi <- 1
tau <- 5
kappa <- 5


# Setting theta based on expert ratings
n_train <- dim(w_train)[2]
sum <- apply(w_train,2,sum)
mj <- rep(0, n_train)
for (i in 1:n_train) {
  if (sum[i] < 4) {
    mj[i] = 0
  }
  else {
    mj[i] = 1
  }
}
ep <- sum(mj)/n_train
theta <- log(-(ep-1)/ep)



## Run the Code!
set.seed(227920)
result_train <- Nuclear(nsample,sbj=sbj_train,y1=y1_train,y2=y2_train,X1=x1_train,X2=x2_train,X3=xs_train,W=w_train,
                        t1,t2,p,q,nu,psi,ups,a,a_sigma,b_sigma,a_zeta,b_zeta,a_iota,b_iota,a_xi,b_xi,tau,kappa,theta)










############################# PREDICT THE DISEASE STATUS ################################
n_test <- dim(w_test)[2] # number of kidneys in testing data
n_sbj_test <- length(O_test) # number of subjects in testing data
burn <- 1:20000  # burn-in iterations
marginal_prob <- matrix(NA,nrow=n_test,ncol=nsample-length(burn)) # marginal probability MCMC iterations
for (i in 1:n_sbj_test) {
  prob <- Nuclear_predict(fit=result_train,bi=length(burn),y1=y1_test[,c(2*i-1,2*i)],y2=y2_test[,c(2*i-1,2*i)],
                          x1=x1_test[,i],x2=x2_test[,c(2*i-1,2*i)],x3=xs_test[,c(2*i-1,2*i)])
  marginal_prob[2*i-1,] <- prob[,2,1]+prob[,2,2] #left kidney Posterior predictive prob. of obstruction
  marginal_prob[2*i,] <- prob[,1,2]+prob[,2,2]   #right kidney Posterior predictive prob. of obstruction
  print(i)
}
marginal_prob_mean <- apply(marginal_prob,1,mean) # Predictive prob of obstruction of each kidney (marginal prob)
pred <- ifelse(marginal_prob_mean > 0.5, 1, 0) # Predicted obstruction status based on the 0.5 cutoff point

## Evaluating the predictive performance
pred_table <- table(C_test,pred) # Confusion Matrix
CR <- (pred_table[1,1]+pred_table[2,2])/sum(pred_table) # Concordance Rate
sen <- (pred_table[2,2])/(pred_table[2,1]+pred_table[2,2]) # Sensitivity
spe <- (pred_table[1,1])/(pred_table[1,1]+pred_table[1,2]) # Specificity
ppv <- (pred_table[2,2])/(pred_table[1,2]+pred_table[2,2]) # Posterior Predictive Value
npv <- (pred_table[1,1])/(pred_table[1,1]+pred_table[2,1]) # Negative Predictive Value
brier <- mean((marginal_prob_mean-C_test)^2) # Brier Score




##### Probability Plot - Figure 2 ######
par(mfrow=c(1,1))
marginal_prob_up <- apply(marginal_prob,1,quantile, 0.975) # upper limit of the interval
marginal_prob_low <- apply(marginal_prob,1,quantile, 0.025) # lower limit of the interval
ti <- 1:n_test

ppdata <- cbind(marginal_prob_mean,marginal_prob_up,marginal_prob_low)
plotCI(ti[which(C_test==0)],ppdata[which(C_test==0),1],ui=ppdata[which(C_test==0),2],li=ppdata[which(C_test==0),3],
       sfrac=0.003,cex=0.7, ylab="Predictive Probability of Obstruction",xlab="Kidney Index",xaxt = "n",xaxt = "n",ylim=c(0,1),xlim=c(1,n_test),xaxt = "n") # CI plot
plotCI(ti[which(C_test==1)],ppdata[which(C_test==1),1],ui=ppdata[which(C_test==1),2],li=ppdata[which(C_test==1),3],col='green',
       sfrac=0.003,cex=0.7, ylab="Predictive Probability of Obstruction",xlab="Kidney Index",xaxt = "n",add=TRUE) # CI plot
axis(1,at=1:n_test,cex.axis = 1)
abline(h=0.5,lty=2,lwd=1,col="red")
text(2,0.55,"Cutoff = 0.5",col="red")

