## Right censored data with a cured fraction

# Write a function to generate survival data with a cured fraction
generate_cured_data <- function(n, b_true, g_true,
                                censoringtimemax=3){
  
  #create X and Z
  X <- as.matrix(rbinom(n, 1, 0.5))
  Z <- as.matrix(runif(n, 0, 1))
  Z.mod <- cbind(1, Z)
  
  #generate results from logistic regression
  pi_i <- exp(Z.mod %*% b_true)/(1 + exp(Z.mod %*% b_true))
  
  #generate cure indicator
  U_C <- runif(n)
  cure_i <- as.numeric(U_C < pi_i) #1 is not cured
  
  #generate observation times
  TL_i <- TR_i <- Y_i <- censor <- rep(0,n)
  
  for(i in 1:n){
    #generate right censoring/"dropout" time
    max.followup <- runif(1, 0.5, censoringtimemax)
    if(cure_i[i] == 0){ #if they are cured
      Y_i[i] <- Inf #they have no event time
      #they are right censored their dropout time
      TL_i[i] <- max.followup
      TR_i[i] <- Inf
    }else{ #if they are not cured
      #generate an event time
      U_Y <- runif(1)
      Y_i[i] <- (-log(U_Y)/exp(X[i,]%*%g_true))^(1/3)
      if(Y_i[i] < max.followup){
        #if their event time is before their 
        #dropout/censoring time
        TL_i[i] <- Y_i[i]
        TR_i[i] <- Y_i[i]
        #they have an event time
      }else{
        #they are right censored
        TL_i[i] <- max.followup
        TR_i[i] <- Inf
      }
    }
    censor[i] <- as.numeric(TR_i[i]!=Inf)
  }
  data <- data.frame(Y_i,TL_i,TR_i,censor,cure_i,Z,X)
  return(data)
}

dat <- generate_cured_data(n = 500, b_true = c(1, -0.5), 
                           g_true = 1, censoringtimemax=1.5)

# Check simulated data
sum(dat$cure_i)/500 #should be ~70%
sum(dat$censor[dat$cure_i == 1])/sum(dat$cure_i == 1) 
#should be ~80%

# Fit model using smcure package
library(smcure)
smc.fit <- smcure(Surv(TL_i, censor) ~ X, cureform = ~ Z, 
                  data = dat, model = "ph", link = "logit")

# Compute pi values for each i
pi_i <- exp(as.matrix(dat[,6:7]) %*% as.matrix(smc.fit$b))/
  (1 + exp(as.matrix(dat[,6:7]) %*% as.matrix(smc.fit$b)))

hist(pi_i)

# Compute and plot the baseline survival function
pred.smc <- predictsmcure(smc.fit, newX = c(0,0), 
                          newZ = c(0), model = "ph")
baseSurv.smc <- (pred.smc$prediction[,1])[order(pred.smc$prediction[,3])]
times.smc <- (pred.smc$prediction[,3])[order(pred.smc$prediction[,3])]
plot(baseSurv.smc ~ times.smc, type = "s", ylim = c(0,1), 
     main = "Predicted baseline survival",
     xlab = "Time", ylab = "Survival probability")
