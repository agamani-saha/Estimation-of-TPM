###################################### GET DATASET ######################################################################
library(GSODR)
library(bit64)
library(dplyr)
library(lubridate)
library(MCMCprecision)
library(markovchain)

heathrow_rain <- read.csv("HeathrowRainMonthly1995To31Aug2025.csv")

monthly_rain <- heathrow_rain

heathrow_rain$states <- round(monthly_rain$PRCP/10)
plot(monthly_rain$states)
###################################### SPLIT DATA ##############################################
full_data <- monthly_rain$states
train_ratio <- 0.5
train_size <- floor(train_ratio * nrow(monthly_rain))
train_data <- monthly_rain[1:train_size, ]

cat("Full data:", nrow(full_data), "rows\n")
cat("Training data:", nrow(train_data), "rows\n")

###################################### FUNCTION TO COMPUTE TPM ##############################################
compute_tpm <- function(states, prob) {
  uniq_states <- 0:max(monthly_rain$states)
  #uniq_states <- 0:max(r)
  d <- length(uniq_states)
  
  # Initialize transition matrix
  tpm <- matrix(0, nrow = d, ncol = d, dimnames = list(uniq_states, uniq_states))
  
  # Count transitions
  for (i in 1:(length(states) - 1)) {
    tpm[as.character(states[i]), as.character(states[i + 1])] <-
      tpm[as.character(states[i]), as.character(states[i + 1])] + 1
  }
  
  # Convert counts to probabilities
  if(prob == TRUE){
    tpm <- tpm / rowSums(tpm)
  } 
  tpm[is.na(tpm)] <- 0
  return(tpm)
}

###################################### COMPUTE TPM ##############################################
TPM_full <- compute_tpm(monthly_rain$states, TRUE)
write.csv(TPM_full, "TPM.csv")
TPM_train <- compute_tpm(train_data$states, TRUE)
TPM_count_train <- compute_tpm(train_data$states, FALSE)
###################################### OUTPUT TPM ##############################################
cat("\n--- Full Dataset TPM ---\n")
print(TPM_full)

cat("\n--- Training Dataset TPM ---\n")
print(TPM_train)

###################################### PLOT RAIN GRAPH ##############################################
monthly_rain$Date <- as.Date(paste0(monthly_rain$YearMonth, "-01"))
train_data$Date   <- as.Date(paste0(train_data$YearMonth, "-01"))
library(ggplot2)
ggplot() +
  geom_line(data = monthly_rain, aes(x = Date, y = states), color = "brown") +
  geom_line(data = train_data, aes(x = Date, y = states), color = "blue") +
  labs(title = "",
       x = "Year", y = "Monthly Rainfall (mm)") +
  theme_minimal()
###################################### SPECIFY INPUTS ##############################################################
library(ggplot2)
set.seed(1)
alpha = 4 
beta = 3  
b0 = 10
N = 2
Burn.in = 1000
M = 1000
d = nrow(TPM_full)

############################# Sampling from Cover Function ###############################################################

fkTilde <- function(x, J, A, B) x^(A - 1)*exp(-B*x)/(gamma(x))^J

SamplingFromCoverFunction <- function(N, J, A, B){
  ###initialize (2N+2) length vector of knots
  m <- rep(0, 2*N + 2)
  
  ### initialize (2N+2) length vectors of intercepts and slopes of the tangent lines
  a <- rep(0, 2*N + 2)
  lambda <- rep(0, 2*N + 2)
  
  ### initialize (2N+2) length vector of points of intersection between the tangent lines
  ### concatenate a 0 in front later on to make q[0] = 0
  q <- rep(0, 2*N + 2)
  
  ### log of the k-th density function
  hkTilde <- function(x) {
    retVal <- -J*log(gamma(x)) + (A - 1)*log(x) - B*x
    return(retVal)
  }
  
  ### and its derivative
  hkPrimeTilde <- function(x) {
    retVal <- -J*psigamma(x, deriv = 0) + (A - 1)/x - B
    return(retVal)
  }
  
  ### central knot
  M <- 1.5*(B > 0) + exp(1 - B/J)*(B < 0)
  m[N + 1] <- optimise(hkTilde, interval = c(0, M), maximum = TRUE, 
                       tol = 0.000001)$maximum
  
  ### last knot 
  m[2*N + 2] <- (m[N + 1] + 1.5)*(B > 0) + exp(1 - B/J)*(B < 0)
  
  ### first knot
  m[1] <- 0.5*m[N + 1]
  
  ### last-but-one knot
  m[2*N + 1] <- mean(c(m[N + 1], m[2*N + 2]))
  
  ### rest of the notes filled up in an equi-spaced manner
  if (N > 1){
    m[2:N] <- seq(m[1], m[N + 1], length.out = N + 1)[-c(1, N + 1)]
    m[(N + 2):(2*N)] <- seq(m[N + 1], m[2*N + 1], length.out = N + 1)[-c(1, N + 1)]
  }
  
  for (i in 1:(2*N + 2)) {
    a[i] <- hkTilde(m[i])
    lambda[i] <- hkPrimeTilde(m[i])
  }
  
  for (i in 1:(2*N + 1)) q[i] <- (a[i + 1] - a[i] + m[i]*lambda[i] - m[i + 1]*
                                    lambda[i + 1])/(lambda[i]- lambda[i + 1])
  q[2*N + 2] <- Inf
  q <- c(0, q)
  
  ### integrating constants
  Cgki <- rep(0, 2*N + 2)
  Cgki[N + 1] <- exp(a[N+1])*(q[N + 2] - q[N + 1])
  
  for (i in setdiff(1:(2*N + 2), N + 1)) Cgki[i] <- exp(a[i] - m[i]*lambda[i])*
    (exp(q[i+1]*lambda[i]) - exp(q[i]*lambda[i]))/lambda[i]
  Cgk <- sum(Cgki)
  
  
  ii <- sample(1:(2*N+2), size = 1, prob = Cgki/Cgk)
  ### now sample from gki
  if (ii == N + 1){
    sim <- runif(1, q[ii], q[ii + 1])
  } else {
    u <- runif(1)
    sim <- log(u*exp(lambda[ii]*q[ii + 1]) + (1 - u)*exp(lambda[ii]*q[ii]))/lambda[ii]
  }
  
  for (i in 1:(2*N + 2)) if (q[i] <= sim && sim < q[i + 1]) gkTildeSim <- 
    exp(a[i] + (sim - m[i])*lambda[i])
  
  return(list(sim = sim, gkTildeSim = gkTildeSim))
}

##### samples from fkTilde using the cover function and rejection sampling method
RejectionSampling <- function(N, J, A, B){
  retVal <- NA
  while (is.na(retVal)){
    ss <- SamplingFromCoverFunction(N, J, A, B)
    gkTildeSim <- ss$gkTildeSim 
    u <- runif(1)
    if (u <= fkTilde(ss$sim, J, A, B)/gkTildeSim) retVal <- ss$sim
  }
  return(retVal)
}

#################### BLOCKED GIBBS SAMPLER #####################################
# Function to draw samples from a dirichlet distribution.
# alpha : concentration parameters

update_Pi = function(n, t){
  
  Pi =  matrix(NA, nrow = d, ncol = d)
  
  for(i in 1:d){
    
    # n_j = (n_j1, ..., n_jL), n_jk = \sum_i I(Z_ji = k)
    n_i = n[i, ]
    
    # Pi |... follows Dirichlet (n_j + t)
    Pi.draw = rdirichlet(1, n_i + t)
    
    # setting an lower bound of 10^(-10) for \pi_jk's to avoid numerical issues
    Pi.ind = which(Pi.draw < 1e-10)
    
    if(length(Pi.ind) > 0){
      
      Pi.draw[Pi.ind] = 1e-10
      
      excess.P = sum(Pi.draw) - 1
      
      ind.max = which.max(Pi.draw)
      Pi.draw[ind.max] = Pi.draw[ind.max] - excess.P
    }
    
    Pi[i, ] = Pi.draw
  }
  return(Pi)
}

update_u = function(d, sum_t){
  u = rgamma(d, shape = sum_t, scale = 1)
  return(u)
}

update_w = function(d, alpha, t){
  w = rep(NA, d)
  for(j in 1:d){
    sum_t_j = sum(t[j:(length(t))])
    w[j] = rgamma(1, shape = alpha, scale = sum_t_j)
  }
  return(w)
}

update_Beta = function(n, N, Pi, alpha, beta, B, u, w){
  
  d <- ncol(n)
  delta <- c(rep(alpha, d - 1), beta)
  
  log_u = sum(log(u))
  log.p = colSums(log(Pi))
  
  # log.p tends to -\infty for unoccupied clusters
  # setting a threshold of -10^10 for log.p to avoid numerical issues
  log.p = sapply(1:d, function(x) ifelse(var(Pi[,x]) == 0, -1e+10, log.p))
  
  # draw t_j from f_j(.), j = 1, 2,..., d.
  t = rep(NA, d)
  for(i in 1:d){
    A = delta[i]
    sum_w = sum(w[1:i])
    B = b0 + sum_w - log.p[i] - log_u
    t[i] = RejectionSampling(N = N, J = d, A = A, B = B)
  }
  sum_t = sum(t)
  Beta = t/sum_t
  
  return(list("Beta" = Beta, "sum_t" = sum_t, "t" = t))
}

blocked_gibbs = function(N, n, Burn.in, M, d, alpha, beta, b0){
  
  # set initial values for running the Gibbs sampler
  Pi = matrix(1/d, nrow = d, ncol = d)
  
  t = rep(1/d, d)
  u = rgamma(n = d, shape = 1, rate = 1)
  w = rgamma(n = d, shape = 1, rate = 1)
  
  # list to store the posterior samples
  Iterates = vector(mode = "list", length = M)
  
  for(m in 1:(M + Burn.in)){
    
    # time at the beginning
    T1 = Sys.time()
    
    # update Pi
    Pi = update_Pi(n = n, t = t)
    
    # update t and Beta 
    res = update_Beta(n = n, N = N, Pi = Pi, alpha = alpha, beta = beta, 
                      B = B, u = u, w = w)
    Beta = res$Beta
    sum_t = res$sum_t
    t = res$t
    
    # update u
    u = update_u(sum_t = sum_t, d = d)
    
    # update w
    w = update_w(d = d, alpha = alpha, t = t)
    
    # time at the end of all updates
    T2 = Sys.time()
    Tdiff =  difftime(T2, T1, units = "secs")
    
    
    # print every 200th iteration
    if(m %% 200 == 0){
      print(paste("iteration :", m))
    }
    
    # store samples after Burn in
    if(m > Burn.in){
      Iterates[[m-Burn.in]] = list("Pi" = Pi, "Beta" = Beta, 
                                   "t" = t, "w" = w,
                                   "u" = u, "time" = Tdiff)
    }
  }
  return(Iterates)
  
}
########################## Calculate MLE and GSBP #################################
  n <- TPM_count_train
  
  bg <- blocked_gibbs(N = N, n = n, Burn.in = Burn.in, M = M, d = d, alpha = alpha, 
                      beta = beta, b0 = b0)
  
  ################################# calculate posterior mean ###############################
  Pi_sum =  matrix(rep(0, d*d), nrow = d)
  Pi_error = rep(0, d)
  for(i in 1:M){
    Pi_sum = Pi_sum + bg[[i]]$Pi
    Pi_error[i] = sqrt(mean((bg[[i]]$Pi - TPM_full)^2))
  }
  Pi_mean = Pi_sum/M
  ################################################################################
  TPM_count_train_Extended = n/rowSums(n)
  TPM_count_train_Extended[is.nan(TPM_count_train_Extended)] <- 0
  TrueVsMLE = TPM_full - TPM_count_train_Extended
  SampVsEst = TPM_count_train_Extended - Pi_mean
  TrueVsEst = TPM_full - Pi_mean
  
  RMSE_MLE<- sqrt(mean(TrueVsMLE^2)) * 100
  RMSE_GSBP <- sqrt(mean(TrueVsEst^2)) * 100
  RMSE_MLE
  RMSE_GSBP
################################## ##############################################
 ##################################### GRAPH of TEST #####################################
  states <- as.character(0:(nrow(Pi_mean)-1))
  rownames(TPM_count_train_Extended) <- states
  colnames(TPM_count_train_Extended) <- states
  rownames(Pi_mean) <- states
  colnames(Pi_mean) <- states
  ###################################### SPLIT TEST DATA ##############################################
  test_data <- monthly_rain[(train_size + 1):nrow(monthly_rain),]
  
  cat("Test data:", nrow(test_data), "rows\n")
  
  ###################################### SIMULATE TEST PERIOD BASED ON TPM ##############################################
  set.seed(12345)
  simulate_from_tpm <- function(initial_state, tpm, steps) {
    states <- numeric(steps)
    states[1] <- initial_state
    uniq_states <- as.numeric(rownames(tpm))
    
    for (i in 2:steps) {
      current_state <- as.character(states[i - 1])
      probs <- tpm[current_state, ]
      if (sum(probs) == 0) {
        # If no outgoing transition, stay in same state
        states[i] <- states[i - 1]
      } else {
        states[i] <- sample(uniq_states, 1, prob = probs)
      }
    }
    return(states)
  }
  
  # Starting state = last state of training data
  initial_state <- tail(train_data$states, 1)
  steps <- nrow(test_data)
  
  # Simulate using MLE TPM
  pred_states_MLE <- simulate_from_tpm(initial_state, TPM_count_train_Extended, steps)
  
  # Simulate using GHSBP TPM
  pred_states_GHSBP <- simulate_from_tpm(initial_state, Pi_mean, steps)
  
  # Convert states to approximate volumes (states represent floored volume)
  pred_vol_MLE <- pred_states_MLE
  pred_vol_GHSBP <- pred_states_GHSBP
  
  ###################################### CREATE COMBINED DATA FOR PLOT ##############################################
  plot_data <- data.frame(
    Date = c(train_data$Date, test_data$Date, test_data$Date, test_data$Date, test_data$Date),
    Rain = c(train_data$states,
               test_data$states,
               pred_vol_MLE,
               pred_vol_HSBP,
               pred_vol_GHSBP),
    Type = c(rep("Training", nrow(train_data)),
             rep("Actual", nrow(test_data)),
             rep("MLE", nrow(test_data)),
             rep("HSBP", nrow(test_data)),
             rep("GHSBP", nrow(test_data)))
  )
###################################### FILTER TEST PERIOD DATA ##############################################
  test_plot_data <- subset(plot_data, Type != "Training")  # Remove training
  test_states <- test_data$states 
###################################### PLOT TEST PERIOD ONLY ##############################################
  library(ggplot2)
  
  ggplot(test_plot_data, aes(x = Date, y = Rain, color = Type)) +
    geom_line(linewidth = 0.5) +
    facet_wrap(~Type, ncol = 1, scales = "fixed") +
    labs(title = "",
         x = "Date", y = "states") +
    theme_minimal() +
    theme(legend.position = "none")
  
  

