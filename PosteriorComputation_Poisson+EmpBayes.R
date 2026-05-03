###########################################################################################
library(MCMCprecision)
library(markovchain)
set.seed(1)
###########################################################################################
  
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
    m[N + 1] <- optimise(hkTilde, interval = c(0, M), maximum = TRUE, tol = 0.000001)$maximum
    
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
    
    for (i in 1:(2*N + 1)) q[i] <- (a[i + 1] - a[i] + m[i]*lambda[i] - m[i + 1]*lambda[i + 1])/(lambda[i]- lambda[i + 1])
    q[2*N + 2] <- Inf
    q <- c(0, q)
    
    ### integrating constants
    Cgki <- rep(0, 2*N + 2)
    Cgki[N + 1] <- exp(a[N+1])*(q[N + 2] - q[N + 1])
    
    for (i in setdiff(1:(2*N + 2), N + 1)) Cgki[i] <- exp(a[i] - m[i]*lambda[i])*(exp(q[i+1]*lambda[i]) - exp(q[i]*lambda[i]))/lambda[i]
    Cgk <- sum(Cgki)
    
    
    ii <- sample(1:(2*N+2), size = 1, prob = Cgki/Cgk)
    ### now sample from gki
    if (ii == N + 1){
      sim <- runif(1, q[ii], q[ii + 1])
    } else {
      u <- runif(1)
      sim <- log(u*exp(lambda[ii]*q[ii + 1]) + (1 - u)*exp(lambda[ii]*q[ii]))/lambda[ii]
    }
    
    for (i in 1:(2*N + 2)) if (q[i] <= sim && sim < q[i + 1]) gkTildeSim <- exp(a[i] + (sim - m[i])*lambda[i])
    
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
  ######## Validation of the rejection sampling ###################################
  ### inputs
  J <- 5
  A <- 0.4
  B <- 1
  N <- 2
  
  val <- rep(0, 10000)
  for (i in 1:10000) val[i] <- SamplingFromCoverFunction(N, J, A, B)$sim
  
  par(mfrow = c(1, 2))
  xSeq <- seq(0.1, 3, length.out = 1000)
  plotVal1 <- sapply(1:1000, FUN = function(i) fkTilde(xSeq[i], J, A, B))
  plot(xSeq, plotVal1, type = 'l', main = "true distribution")
  hist(val, 50, main = "empirical distribution")
  
  ### function to generate n samples from Dirichlet Distribution with parameter alpha
  rdirichlet <- function(n, alpha){
    k <- length(alpha)
    g <- matrix(0, n, k)
    for (i in 1:k){
      g[, i] <- rgamma(n, shape = alpha[i])
    }
    g <- g/matrix(rowSums(g), nrow = n, ncol = k, byrow = FALSE)
    return(g)
  }
  
  
  #### returns a list with d components
  #### j-th component is pi[[j]]
  samplingPI <- function(gammaTilde, nTilde, alpha0){
    pi <- vector(mode = "list", length = d)
    pi[[j]] <- rdirichlet(nTilde[[j]] + alpha0*gammaTilde)
    return(pi)
  }
  
  
  #### returns a vector of length d
  #### j-th component follows gamma(t_j) independently
  samplingUTilde <- function(tTilde){
    u <- rep(0, d)
    for (j in 1:d) u[j] <- rgamma(1, shape = tTilde[j], scale = 1)
    return(u)
  }
  
  #### returns a vector of length d
  #### j-th component follows independently gamma(alpha, sum_{k=j}^d t[j])
  samplingWTilde <- function(alpha, tTilde){
    w <- rep(0, d)
    for (j in 1:d){
      w[j] <- rgamma(n = 1, shape = alpha, rate = sum(tTilde[j:d]))
    }
    return(w)
  }
  
  #### beeta has been used to avoid conflict with the R base function beta
  samplingTTilde <- function(PI, UTilde, WTilde){
    deltaj <- alpha
    for (j in 1:d){
      if (j == d) deltaj <- beeta
      Bj <- b0 + sum(WTilde[1:j]) - sum(PI[, j]) - sum(log(UTilde))
      tTilde[j] <- rTiltedGamma(1, d, deltaj, Bj, alpha0)
    }
  }
  
#################################################################################
  sample.size = 10^6
  random.start = 0
  sim.data = rep(NA, sample.size)
  sim.data[1] = random.start
  for(i in 2:sample.size){
    sim.data[i] = rpois(n = 1, lambda = log(sim.data[i-1] + 10))
  }
  trans.matrix <- function(X, prob = F)
  {
    tt <- table(c(X[-length(X)]), c(X[-1]))
    if(prob == T) tt <- tt / rowSums(tt)
    tt
  }
  trans.matrix.sample = trans.matrix(sim.data)
  fullmat_dim <- 1 + max(as.numeric(c(rownames(trans.matrix.sample), 
                                           colnames(trans.matrix.sample))))
  fullnames <- as.character(0:(fullmat_dim - 1))
  fullmat <- matrix(0, nrow = fullmat_dim, ncol = fullmat_dim, 
                    dimnames = list(fullnames, fullnames))
  fullmat[rownames(trans.matrix.sample), colnames(trans.matrix.sample)] <- 
    trans.matrix.sample
  trans.matrix.sample <- fullmat
  sample.prob.matrix <- trans.matrix.sample/rowSums(trans.matrix.sample)
  sample.prob.matrix[is.nan(sample.prob.matrix)] <- 0
  
########################## Put the values of the hyperparameters #########################
d = nrow(trans.matrix.sample)
n <- trans.matrix.sample
alpha = 1
beta = 1
b0 = 10
N = 2
Burn.in = 1000
M = 2000
bg <- blocked_gibbs(N = N, n = n, Burn.in = Burn.in, M = M, d = d, alpha = alpha, 
                    beta = beta, b0 = b0)
  
################################# calculate sample mean ########################
Pi_sum =  matrix(rep(0, d*d), nrow = d)
for(i in 1:M){
  Pi_sum = Pi_sum + bg[[i]]$Pi
}
Pi_mean = Pi_sum/M
Pi_mean
########################### True TPM ###########################################
trans.matrix.true = matrix(rep(NA, d*d), nrow = d)
for(i in 0:(d-1)){
  for(j in 0:(d-1)){
    lambda = (log(i+10))
    trans.matrix.true[i+1,j+1] = exp(-lambda)*(lambda^j)/factorial(j)
  }
}
trans.matrix.true
############################## Empirical Bayes #######################################################
empirical_bayes_markov <- function(Y, tol = 1e-6, max_iter = 100) {
  # Y: k x k matrix of transition proportions (y_ij = n_ij / n_i.)
  
  k <- ncol(Y)
  
  # ---- Step 1: Compute geometric means G_j ----
  G <- apply(Y, 2, function(col) {
    exp(mean(log(col)))
  })
  
  # ---- Step 2: Initial alpha estimates ----
  sumG <- sum(G)
  
  alpha_dot <- ((k - 1) / 2) / (1 - sumG) + 1/2
  
  alpha_vec <- 1/2 + (alpha_dot - 1/2) * G
  
  # ---- Step 3: Newton-Raphson Iteration ----
  
  for (iter in 1:max_iter) {
    
    alpha_old <- alpha_vec
    alpha_dot <- sum(alpha_vec)
    
    # ---- Score vector ----
    S <- numeric(k)
    for (j in 1:k) {
      S[j] <- k * digamma(alpha_dot) - k * digamma(alpha_vec[j]) + sum(log(Y[, j]))
    }
    
    # ---- Information matrix ----
    trig_alpha <- trigamma(alpha_vec)
    trig_sum <- trigamma(alpha_dot)
    
    D <- diag(trig_alpha)
    I <- k * (D - trig_sum * matrix(1, k, k))
    
    # ---- Update step ----
    delta <- delta <- MASS::ginv(I) %*% S
    alpha_vec <- alpha_vec + delta
    
    # ---- Convergence check ----
    if (max(abs(alpha_vec - alpha_old)) < tol) {
      cat("Converged in", iter, "iterations\n")
      break
    }
  }
  
  return(alpha_vec)
}
sample.prob.matrix[sample.prob.matrix == 0] <- 1e-10
alpha_est <- empirical_bayes_markov(sample.prob.matrix)

print(alpha_est)
##
estimate_TPM_EB <- function(n, alpha_est) {
  
  d <- nrow(n)
  alpha_sum <- sum(alpha_est)
  
  P_EB <- matrix(0, d, d)
  
  for (i in 1:d) {
    row_sum <- sum(n[i, ])
    
    if (row_sum == 0) {
      # fallback: use prior
      P_EB[i, ] <- alpha_est / alpha_sum
    } else {
      P_EB[i, ] <- (n[i, ] + alpha_est) / (row_sum + alpha_sum)
    }
  }
  
  return(P_EB)
}
P_EB <- estimate_TPM_EB(n, alpha_est)

round(P_EB, 3)
############################## Empirical Bayes Hyperparameters different #######################################################
empirical_bayes_markov_unequal <- function(n, tol = 1e-6, max_iter = 100) {
  # n: k x k transition count matrix
  
  k <- nrow(n)
  
  alpha_vec2 <- matrix(1, k, k)  # initial values
  
  for (i in 1:k) {
    
    ni <- n[i, ]
    ni_dot <- sum(ni)
    
    alpha_i <- rep(1, k)
    
    for (iter in 1:max_iter) {
      
      alpha_old <- alpha_i
      alpha_dot <- sum(alpha_i)
      
      # ---- Score vector ----
      S <- numeric(k)
      for (j in 1:k) {
        S[j] <- digamma(alpha_dot) - digamma(alpha_i[j]) +
          digamma(ni[j] + alpha_i[j]) - digamma(ni_dot + alpha_dot)
      }
      
      # ---- Information matrix ----
      trig_alpha <- trigamma(alpha_i)
      trig_nalpha <- trigamma(ni + alpha_i)
      
      trig_sum <- trigamma(alpha_dot)
      trig_nsum <- trigamma(ni_dot + alpha_dot)
      
      d_vec <- trig_alpha - trig_nalpha
      c_val <- trig_sum - trig_nsum
      
      D <- diag(d_vec)
      I <- D - c_val * matrix(1, k, k)
      
      # ---- Update step ----
      delta <- MASS::ginv(I) %*% S
      alpha_i <- alpha_i + as.vector(delta)
      
      # ---- Stability ----
      alpha_i[alpha_i <= 1e-10] <- 1e-10
      
      # ---- Convergence ----
      if (max(abs(alpha_i - alpha_old)) < tol) {
        cat("Row", i, "converged in", iter, "iterations\n")
        break
      }
    }
    
    alpha_vec2[i, ] <- alpha_i
  }
  
  return(alpha_vec2)
}
alpha_vec2 <- empirical_bayes_markov_unequal(n)

print(alpha_vec2)
estimate_TPM_EB_unequal <- function(n, alpha_vec2) {
  
  d <- nrow(n)
  P_EB2 <- matrix(0, d, d)
  
  for (i in 1:d) {
    
    row_sum <- sum(n[i, ])
    alpha_sum <- sum(alpha_vec2[i, ])
    
    if (row_sum == 0) {
      P_EB2[i, ] <- alpha_vec2[i, ] / alpha_sum
    } else {
      P_EB2[i, ] <- (n[i, ] + alpha_vec2[i, ]) / (row_sum + alpha_sum)
    }
  }
  
  return(P_EB2)
}
P_EB2 <- estimate_TPM_EB_unequal(n, alpha_vec2)

round(P_EB2, 3)
################################################################################
round(trans.matrix.true[1:d,1:d], 3)
TrueVsMLE = trans.matrix.true[1:d,1:d] - sample.prob.matrix
SampVsEst = sample.prob.matrix - Pi_mean
TrueVsEst = trans.matrix.true[1:d,1:d] - Pi_mean
TrueVsPB = trans.matrix.true[1:d,1:d] - P_EB 
TrueVsPB2 = trans.matrix.true[1:d,1:d] - P_EB2 


round(trans.matrix.true,3)
round(Pi_mean,3)
round(mean(abs(TrueVsMLE))*100, 3)
round(mean(abs(TrueVsEst))*100, 3)
round(mean(abs(TrueVsPB))*100, 3)
round(mean(abs(TrueVsPB2))*100, 3)
RMSE_MLE <- sqrt(mean(TrueVsMLE^2)) * 100
RMSE_GSBP <- sqrt(mean(TrueVsEst^2)) * 100
RMSE_PB <- sqrt(mean(TrueVsPB^2)) * 100
RMSE_PB2 <- sqrt(mean(TrueVsPB2^2)) * 100
round(RMSE_MLE, 3)
round(RMSE_GSBP, 3)
round(RMSE_PB, 3)
round(RMSE_PB2, 3)

################################################################################
