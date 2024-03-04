## Melanoma example with a cured fraction

# Write a function to generate data similar to the melanoma study
generate_melanoma_example <- function(n, b_true, g_true,event=FALSE,event_adj=0,censoringtimemax=3){
  
  # create X and Z
  BreslowThickness <- rbinom(n,1,0.5)
  SexMale <- rbinom(n,1,0.5)
  Mitosis <- rbinom(n,1,0.5)
  Z <- (cbind(BreslowThickness, SexMale, Mitosis))
  
  Ulceration <- rbinom(n,1,0.5)
  X <- (cbind(Ulceration, SexMale))
  
  # generate results from logistic regression
  pi_i <- rep(0,n)
  for(i in 1:n){
    pi_i[i] <- exp(b_true[1]+sum(b_true[-c(1)]*Z[i,]))/(1+exp(b_true[1]+sum(b_true[-c(1)]*Z[i,])))
  }
  
  # generate cure indicator
  U_C <- runif(n)
  cure_i <- as.numeric(U_C < pi_i) #1 is not cured
  uncure.n <- length(cure_i[cure_i]) 
  #number of uncured in sample
  cure.n <- length(cure_i[!cure_i]) 
  #number of cured in sample
  
  # generate observation times
  TL_i <- TR_i <- rep(0,n)
  Y_i <- rep(0,n)
  
  for(i in 1:n){ 
    #generate censoring time 
    max.followup <- as.numeric(cure_i[i] == 1) * runif(1,min=1,max=censoringtimemax) +
      as.numeric(cure_i[i] == 0) * runif(1,min=1,max=2*censoringtimemax) 
    
    # simulate event monitoring process
    exam.num <- as.numeric(cure_i[i] == 1) * rpois(1,8) + 1 +
      as.numeric(cure_i[i] == 0) * rpois(1,16) + 1
    exam.int <- runif(exam.num, min=0.3, max=0.7)
    examtimes <- cumsum(exam.int)
    
    if(cure_i[i]==0){ #if they are cured
      Y_i[i] <- Inf #they have no event time
      max.ind <- max(which(examtimes<max.followup)) 
      #they are right censored at the last 
      #exam time before they exit the study
      TL_i[i] <- examtimes[max.ind]
      TR_i[i] <- Inf
      
    }else{ #if they are not cured
      # generate an event time
      U_Y <- runif(1)
      Y_i[i]=(-log(U_Y)/exp(X[i,]%*%g_true))^(1/3) 
      #F^-1 method to sample from weibull(1, 3)
      
      if(max.followup < min(examtimes)){ #if they exit the study before the first exam
        TL_i[i] <- max.followup 
        #record their censoring time, they are right censored
        TR_i[i] <- Inf
        
      }else if(Y_i[i] < max.followup & max.followup > min(examtimes)){ #if they have at least one exam AND have their event before exiting the study
        last.exam.ind <- max(which(examtimes < max.followup))
        #they might be left censored
        if(Y_i[i] < min(examtimes)){
          TL_i[i] <- -Inf
          TR_i[i] <- min(examtimes)
          
        }else if(examtimes[last.exam.ind] < Y_i[i]){
          #they might be right censored
          TL_i[i] <- examtimes[last.exam.ind]
          TR_i[i] <- Inf
          
        }else{
          #they might be interval censored
          left <- max(which(examtimes<Y_i[i]))
          right <- min(which(examtimes>Y_i[i]))
          TL_i[i] <- examtimes[left]
          TR_i[i] <- examtimes[right]
          
          #in some cases
          #they might have an event time
          if(event==TRUE & (TR_i[i]-TL_i[i]) < event_adj){
            TL_i[i]=Y_i[i]
            TR_i[i]=Y_i[i]
          }
        }
      }else if(Y_i[i] > max.followup & max.followup > min(examtimes)){ 
        #if they have at least one exam but do not have their event before exiting
        last.exam.ind <- max(which(examtimes<max.followup))
        TL_i[i] <- examtimes[last.exam.ind]
        TR_i[i] <- Inf
      }
    }
  }
  
  data  <- data.frame(Y_i,TL_i,TR_i,cure_i,Z,X)
  return(data)
  
}

# Generate data set

# use below commented out code to generate a new dataset
#mel.dat <- generate_melanoma_example(n=2000, b_true = c(-1.5, 0.5, 2, 1), g_true = c(0.5, -2), censoringtimemax = 1.5)

# Read in data set used in the book
mel.dat <- read.csv("mel_dat_save.csv")

mel.surv <- Surv(time=mel.dat$TL, time2=mel.dat$TR, type = "interval2")

# Fit mixture cure model using MPL method
source("phmc_mpl.R")
source("phmc_mpl_summary.R")
source("penalty_phmc.R")
source("knots_phmc.R")
source("basis_phmc.R")

m1 <- phmc_mpl(mel.surv ~ mel.dat$BreslowThickness + mel.dat$SexMale + mel.dat$Mitosis + mel.dat$Ulceration,
               pi.formula = ~mel.dat$BreslowThickness 
               + mel.dat$SexMale + mel.dat$Mitosis 
               + mel.dat$Ulceration,
               data = mel.dat,
               control = phmc_mpl.control(smooth=NULL,new_constraint = 1,new_criteria = 1,new_knots = 1,
                                          maxIter = c(10,4000,10000),n.knots = 3,conv_limit=1e-5))


# Compare with non mixture cure model
control.mel <- coxph_mpl.control(n.obs=NA, basis="mspline", smooth=NULL, n.knots=c(3,0), max.iter = c(100,2000,4000))
m2 <- coxph_mpl(mel.surv ~ mel.dat$BreslowThickness + mel.dat$SexMale + mel.dat$Mitosis + mel.dat$Ulceration, data = mel.dat, control = coxph_mpl.control(n.obs=NA, basis="mspline", smooth=NULL, n.knots=c(5,0), max.iter = c(100,2000,4000)))

