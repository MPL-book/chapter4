
# Discrimination of the logistic model for incidence

# Calculate individual pi and 1 - pi from the fitted model
prob_uncured <- m1$pi_estimate
prob_cured <- 1 - m1$pi_estimate

# Calculate the survival function value at t 
times <- mel.dat$TL_i
X <- m1$data$X
Z <- m1$data$Z

Psi <- mSpline(times, knots = c(m1$knots$Alpha[2:4]), 
               Boundary.knots = c(m1$knots$Alpha[1], m1$knots$Alpha[5]), 
               integral = TRUE)
H0_t <- Psi %*% m1$theta
H_t <- H0_t * exp(X %*% m1$gamma)
S_t <- exp(-H_t)
S_t[which(is.nan(S_t))] <- 0

# ASANO ET AL. METHOD
#clinical cure tau = 1.5
# whether cure status is known exactly
v_i = as.numeric(m1$data$censoring > 0) 
v_i[which(m1$data$censoring == 0 & 
            mel.dat$TL_i > 1.5)] = 1

# cure status value
u_i = as.numeric(m1$data$censoring == 0) 
u_i[which(m1$data$censoring == 0 & mel.dat$TL_i < 1.5)] = 
  prob_cured[which(m1$data$censoring == 0 & 
                     mel.dat$TL_i < 1.5)]

#sensitivity
sum(as.numeric(prob_uncured < 0.25) * 
      (v_i * u_i + (1-v_i) * (1-prob_cured)))/
  sum((v_i * u_i + (1-v_i) * (1-prob_cured)))

#specificity
1 - sum(as.numeric(prob_uncured < 0.25) * 
          (v_i * (1 - u_i) + (1-v_i)*prob_cured))/
  sum(v_i * (1 - u_i) + (1-v_i)*prob_cured)

#clinical cure tau = 2
## whether cure status is known exactly
v_i = as.numeric(m1$data$censoring > 0) 
v_i[which(m1$data$censoring == 0 & mel.dat$TL_i > 2)] = 1

## cure status value
u_i = as.numeric(m1$data$censoring == 0) 
u_i[which(m1$data$censoring == 0 & 
            mel.dat$TL_i < 2)] = 
  prob_cured[which(m1$data$censoring == 0 
                   & mel.dat$TL_i < 2)]

#sensitivity
sum(as.numeric(prob_uncured < 0.25) * 
      (v_i * u_i + (1-v_i) * (1-prob_cured)))/
  sum((v_i * u_i + (1-v_i) * (1-prob_cured)))

#specificity
1 - sum(as.numeric(prob_cured < 0.75) * 
          (v_i * (1 - u_i) + (1-v_i)*prob_cured))/
  sum(v_i * (1 - u_i) + (1-v_i)*prob_cured)

## AMICO ET AL. METHOD
# clinical cure tau = 1.5
w_i0 = as.numeric(mel.dat$TL_i > 1.5) + 
  (m1$data$censoring == 0) * 
  (as.numeric(mel.dat$TL_i < 1.5) * 
     prob_cured/(prob_cured + prob_uncured*S_t))
w_i1 = 1 - w_i0

#sensitivity
1/((2000-sum(w_i0))) * sum(w_i1 * 
                             as.numeric(prob_uncured < 0.25))

#specificity
1 - (1/sum(w_i0)) * sum(w_i0 * 
                          as.numeric(prob_uncured < 0.25))

# clinical cure tau = 2
w_i0 = as.numeric(mel.dat$TL_i > 2) + 
  (m1$data$censoring == 0) * 
  (as.numeric(mel.dat$TL_i < 2) * 
     prob_cured/(prob_cured + prob_uncured*S_t))
w_i1 = 1 - w_i0

#sensitivity
1 - (1/sum(w_i0)) * 
  sum(w_i0 * as.numeric(prob_cured < 0.75))

#specificity
1/((2000-sum(w_i0))) * 
  sum(w_i1 * as.numeric(prob_cured < 0.75))

## ZHANG METHOD
# no clinical cure
w_i = as.numeric(m1$data$censoring != 0) + 
  (m1$data$censoring == 0) * 
  (prob_uncured*S_t)/(prob_cured + prob_uncured*S_t)

#sensitivity
sum(w_i * as.numeric(prob_uncured < 0.25))/sum(w_i)

#specificity
sum((1-w_i) * as.numeric(prob_uncured > 0.25))/sum(1-w_i)

