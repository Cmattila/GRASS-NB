#' @title Function for the NB variable selection (HS or SS) with group updates (with or without space)
#'
#' @param K Design matrix including an intercept column (first column all 1s)#'
#' @param y is the vector of the count outcome
#' @param group_ind is the vector of feature grouping indicators
#' @param NeighborhoodList is an object of class \code{nb}, typically created with \code{\link[spdep]{poly2nb}}. Each element of the list corresponds to a spatial unit and contains the indices of its neighboring units, as defined by the chosen contiguity rile (e.g queen or rook). Used to encode the spatial adjacency structure for analysis.
#' @param which.prior denotes the model to be used, "HS" or "SS", horseshoe and spike-and-slab respectively
#' @param niter is the number of MCMC iterations
#' @param verbose is TRUE or FALSE, to show or suppress progression updates
#' @param pop_col is a vector of the population of each spatial unit for the model offset, if modeling with space
#' @param scale is an integer of the desired scale for outcome rates (e.g 100000 for per 100,000 rates)

#' @return It returns the estimated beta's, r, shrinkage parameters, phi, Delta and pDelta (if applicable )
#' @export

VS_Group_offset <- function(K, y, group_ind, NeighborhoodList, which.prior, niter, verbose, pop_col, scale){

  # currently all algorithms update r with Chinese restaurant table (CRT) and have a burn-in of half

  # variable definitions
  N <- nrow(K)
  G <- length(unique(group_ind)) # number of groups
  Mg_tab <- table(group_ind)         # Mg is a vector with the number of variables in each group
  Mg <- as.integer(Mg_tab)
  cumulative_M <- cumsum(Mg)
  scale <- as.numeric(scale)
  which.prior <- which.prior

  p = ncol(K)               #includes intercept
  if (sum(Mg) != (p - 1)) {
    stop("Group indicator vector does not match the number of covariates!")
  }

  # MCMC features
  niter<-niter	              # Number of MCMC Iterations
  thin<-1				            # Thinning interval
  burn<-niter/2		          # Burnin
  lastit<-(niter-burn)/thin	# Last stored value

  # Spatial variables
  if(!is.null(NeighborhoodList)){
    neighbor <- NeighborhoodList
    nspace <- length(neighbor) # Number of spatial units
    A <- spdep::nb2mat(neighbor, style = "B")
    neighbor_count <- rowSums(A)  # M (diagonal matrix) neighbor count for each spatial unit (notation from Mutiso, Neelon 2022)
    D <- diag(neighbor_count)  # Create diagonal matrix... Degree Matrix
    Q <- D - A     # CAR precision matrix

    #spatial parameter inits
    phi_init <- c(spam::rmvnorm(1, sigma=diag(.01, nspace)))	          # Random effects
    phi <- phi_init-mean(phi_init) # centered phi
    s2phi <- var(phi)              # var phi
    sphi <- sqrt(s2phi)            # std phi
    tauphi <- 1/s2phi              # precision of phi
  }

  # Storage objects
  Beta <- matrix(0,lastit,p)              # store betas
  Delta <- pDelta<-matrix(0,lastit,p)     # store SS inclusion indicator and probability of inclusion
  Lambda <- matrix(0,lastit,p-1)            # store local-local shrinkage parameters
  Tau2 <- matrix(0,lastit,G)              # store global-local shrinkage parameter
  Zeta2<-matrix(0,lastit,1)               # store global-global shrinkage parameters
  R <- matrix(0,lastit,1)                 # store over-dispersion parameter
  Sphi.sq <- matrix(0,lastit,1)           # store variance of spatial random effect, phi
  phi.mat <- matrix(0,lastit,length(y))     # store spatial random effect

  # Initialize
  beta <- rep(0, p)
  tau <- rep(0.1, p)
  r <- 1
  adjLambda2 <- Lambda2 <- rep(0.1, p-1)
  a_eta<-eta2<-rep(0.1, p-1)
  gamma <- rep(0.1, p-1)
  tau2 <- rep(0.1, G)
  omega <- rep(0, G)
  lambda2 <- 5
  del1<-1
  del2<-1                  # lambda hyper parameters
  nugget<-0.05
  T0<-0.1                 # prior on features
  s<-0.01
  epsilon<-1
  c2<-1
  c2v<-4
  c2s<-2
  wa0<-1
  wb0<-1
  delta <- pdelta <- c(1, rep(0, p-1))
  lambda0 <- rep(1, p)
  a.p0 <- 1
  b.p0 <- 1
  tau0 <- 1
  a0prior <- rep(0, p)
  iA0 <- diag(T0, p, p)
  psi <- rep(1, p - 1)
  psi.nu <- 5               # hyper-parameter of slab variance
  V <- 5                    # slab variance
  psi.Q <- V*(psi.nu - 1)

  epsilon <- rep(1, G)
  zeta2 <- 1
  xi <- 1


  l <- rep(0, N)      # latent vector for CRT-based update of r
  a <- b <- 0.01      # Gamma hyperparms for r for CRT
  Acc <- 0            # counter for MH based update for r, unused

  # Starting VS algorithms
  if(verbose == "TRUE" & which.prior == "HS" & !is.null(NeighborhoodList)){print(paste0("MUSC: NB model fitting has started with group HS prior for spatial count data with offset. "))}
  if(verbose == "TRUE" & which.prior == "SS" & !is.null(NeighborhoodList)){print(paste0("MUSC: NB model fitting has started with group SS prior for spatial count data with offset."))}
  if(verbose == "TRUE" & which.prior == "HS" & is.null(NeighborhoodList)){print(paste0("MUSC: NB model fitting has started with group HS prior for aspatial count data with offset."))}
  if(verbose == "TRUE" & which.prior == "SS" & is.null(NeighborhoodList)){print(paste0("MUSC: NB model fitting has started with group SS prior for aspatial count data with offset."))}

  if(which.prior == "HS" & !is.null(NeighborhoodList)){
    #algorithm for group HS VS for spatial count data (with spatial random effect)----

    for (i in 1:niter){
      #============================================================================#
      ## Algorithm 1st step: update the predictor function and probability vector q ----
      #============================================================================#
      eta <- K%*%beta + rep(phi,1) + log(pop_col) - log(scale)  #modeling rates with offset for population and scale

      q<-1/(1+exp(eta)) # dnbinom fn uses q=1-psi
      q<-ifelse(q<=0, 10^-5,q)


      #============================================================================#
      ## Algorithm 2nd step: Chinese restaurant table-based r update ----
      #============================================================================#
      for(j in 1:N) l[j]<-sum(rbinom(y[j],1,round(r/(r+1:y[j]-1),6)))
      r<-rgamma(1,a+sum(l),b-sum(log(q)))


      #============================================================================#
      ## Algorithm 3rd step: Polya-gamma weights and latent response update ----
      #============================================================================#
      w<-BayesLogit::rpg(N,y+r,eta)                               # Polya weights
      z <- ((y-r)/(2*w)) - log(pop_col) +log(scale)     # latent response


      #============================================================================#
      ## Algorithm 4th step: updating posterior covariance of beta and sample betas ----
      #============================================================================#
      # based on priors, first coefficient i.e., intercept is ignored from penalization
      K_w <- sqrt(w)*K
      Sigma_inv <- diag(c(T0, 1/Lambda2/rep(tau2, times = Mg)/zeta2)) + crossprod(K_w) + diag(nugget, p)
      M <- t(K_w)%*%(sqrt(w)*z)
      beta <- c(spam::rmvnorm.canonical(1, M, as.matrix(Sigma_inv)))
      beta_cov <- beta[-1]


      #============================================================================#
      ## Algorithm 5th step: update shrinkage parameters ----
      #============================================================================#
      summand_vec <- NULL
      for(g in 1:G) {
        for(j in 1:Mg[g]) {
          if(g > 1){                                                    # g = 1 corresponds to the first group, cumulative_M captures M1, M1+ M2, M1+M2+M3,...
            Lambda2[j + cumulative_M[g - 1]] <- MCMCpack::rinvgamma(1, 1, scale =
                                                                      1/gamma[j + cumulative_M[g - 1]] + beta_cov[j + cumulative_M[g - 1]]^2/2/tau2[g]/zeta2)
            gamma[j + cumulative_M[g - 1]] <- MCMCpack::rinvgamma(1, 1,
                                                                  scale = 1 + 1/Lambda2[j + cumulative_M[g - 1]])
          }else{
            Lambda2[j] <- MCMCpack::rinvgamma(1, 1, scale =               # the if else statement could potentially be merged if you define cumulative_M to store 0 as well, i.e., 0, M1, M1 + M2,...
                                                1/ gamma[j + Mg[1]] + beta_cov[j]^2/tau2[1]/zeta2/2)
            gamma[j] <- MCMCpack::rinvgamma(1, 1, scale =  1 + 1/Lambda2[j])
          }
        }
        if(g > 1) {
          summand <- sum(beta_cov[(cumulative_M[g - 1] + 1):cumulative_M[g]]^2/      # carefully note which sum_j beta_j^2's are being used in each case
                           Lambda2[(cumulative_M[g - 1] + 1):cumulative_M[g]])
        }else{
          summand <- sum(beta_cov[1:Mg[1]]^2/Lambda2[1:Mg[1]])
        }
        summand_vec <- c(summand_vec, summand)
        tau2[g] <- MCMCpack::rinvgamma(1, (Mg[g] + 1)/2, scale = 1/epsilon + summand/zeta2/2)
        tau2[g] <- ifelse(tau2[g] > exp(5), exp(5),
                          ifelse(tau2[g] < exp(-5), exp(-5), tau2[g]))    #thresholding
        epsilon[g] <- MCMCpack::rinvgamma(1, 1, scale = 1 + 1/tau2[g])
      }
      zeta2 <- MCMCpack::rinvgamma(1, p/2,                           # all sum_j beta_j^2 is being used for Zeta update
                                   scale = 1/xi + sum(summand_vec/tau2)/2)
      zeta2 <- ifelse(zeta2 > exp(5), exp(5),
                      ifelse(zeta2 < exp(-5), exp(-5), zeta2))    #thresholding
      xi <-  MCMCpack::rinvgamma(1, 1, scale = 1 + 1/zeta2)


      #============================================================================#
      #============================================================================#
      priorprec<-1/(sphi^2)*Q                    # Prior precision of phi1
      priormean<- 0
      prec<-priorprec+spam::as.spam(diag(w, nspace, nspace))
      M<-c(w*(z-K%*%beta))                           # note priormean is 0 and only data contributes to M
      phi<-spam::rmvnorm.canonical(1, M, prec)

      # center phi and update tauphi
      phi <- phi-mean(phi)

      tauphi<-rgamma(1, .1+(nspace-1)/2, .1+(phi)%*%Q%*%t(phi)/2)      # n-1 since rank of Q is n-1
      sphi <- sqrt(1/tauphi)


      #============================================================================#
      # Store assigning ----
      #============================================================================#
      if (i> burn & i%%thin==0) {
        j<-(i-burn)/thin
        Beta[j,]<-beta
        Delta[j, ]<-delta
        pDelta[j, ]<-pdelta
        Lambda[j, ]<-Lambda2
        Zeta2[j] <- zeta2
        Tau2[j,] <-tau2
        R[j]<-r
        phi.mat[j,]<-phi
        Sphi.sq[j] <- sphi^2                       # spatial random effect variance
      }

      if(verbose == "TRUE"){
        svMisc::progress(i, niter, progress.bar = FALSE)
      }else{if(i%%500 == 0){print(paste0(i, " / ", niter))}
      }
    }


  }else if(which.prior == "SS" & !is.null(NeighborhoodList)){
    #algorithm for group SS VS for spatial count data (with spatial random effect)----
    for (i in 1:niter){

      #============================================================================#
      ## Algorithm 1st step: update the predictor function and probability vector q ----
      #============================================================================#
      eta <- K%*%beta + rep(phi,1) + log(pop_col) - log(scale)  #modeling rates with offset for population and scale

      q<-1/(1+exp(eta)) # dnbinom fn uses q=1-psi
      q<-ifelse(q<=0, 10^-5,q)  #prevents overly small q values for stability


      #============================================================================#
      ## Algorithm 2nd step: Chinese restaurant table-based r update ----
      #============================================================================#
      for(j in 1:N) l[j]<-sum(rbinom(y[j],1,round(r/(r+1:y[j]-1),6))) # Could try to avoid loop
      r<-rgamma(1,a+sum(l),b-sum(log(q)))


      #============================================================================#
      ## Algorithm 3rd step: Polya-gamma weights and latent response update ----
      #============================================================================#
      w<-BayesLogit::rpg(N,y+r,eta)                               # Polya weight
      z <- ((y-r)/(2*w)) - log(pop_col) +log(scale)


      #============================================================================#
      ## Algorithm 4th step: updating posterior covariance of beta and sample betas ----
      #============================================================================#
      # based on priors, first coefficient i.e., intercept is ignored from penalization

      # this is a fusion between SS and HS, group-level HS and inside group simply the earlier SS
      # we are looping over groups to generate group-level omega for SS
      incfix <- NULL
      for(g in 1:G){
        if(g == 1){
          delta_sub <- delta[2:(Mg[1] + 1)]
        }else{
          delta_sub <- delta[(cumulative_M[g - 1] + 2):(cumulative_M[g] + 1)]
        }
        incfix <- c(incfix, sum(delta_sub == 1))                         # how many delta_j's are nonzero at this iteration, first element is always 1 as it corresponds to the intercept
        omega[g]  <- rbeta(1, wa0 + incfix[g], wb0 + Mg[g] - incfix[g])               # mixture weight from Eq. 9, page 56 from Dvorzak
      }

      invA0 <- diag(c(T0, 1/psi), nrow = p)            # diagonal precision matrix with first element being fixed as it corresponds the intercept

      del_up <- update_delta_g(delta, omega, invA0, z, w, K, p, G, Mg, cumulative_M, a0prior) # delta vector update now loops over each group to update each delta_gi
      delta <- del_up[[1]]                                             # first element of the list correspond to delta_j's
      pdelta <- del_up[[2]]                                            # second element of the list correspond to p(delta_j = 1)'s
      Mstar <- del_up[[3]]                                             # important for nu update

      index <- which(delta == 1)                         # indices of the non-zero delta_j's, first element always selected as it corresponds to the intercept
      invA0 <- invA0[index, index, drop = FALSE]         # we only simulate the non-zero betas, thus only those rows/columns are selected
      Xsel  <- K[, index, drop = FALSE]*sqrt(w) 		     # select the columns of covariate matrix K for which beta_j's are nonzero, adjustment by Polya weights w
      yc <- sqrt(w)*z                                    # adjusting z by Polya weights w, recall that we had the term: (z - K beta)^T diag(w) (z - K beta) inside the exponent,
      # these weight adjustments simplifies it as  (yc - Xsel beta)^T  (yc - Xsel beta), getting rid of the diag(w), note that Xsel is
      # just a transformed version of K, to be specific, a few columns of K for which beta_j's are non-zero

      k_sel <- ncol(Xsel)
      Sigma_inv <- invA0 + t(Xsel)%*%Xsel                            # posterior precision matrix
      sim_beta <- spam::rmvnorm.canonical(1, (t(Xsel)%*%(yc)),       # simulated beta_j's (for only non-zero delta_j's)
                                          as.matrix(Sigma_inv))
      beta[index] <- sim_beta                                        # store the non-zero betas at appropriate indices
      beta[-index] <- 0                                              # rest are simply 0
      beta_cov <- beta[-1]


      #============================================================================#
      ## Algorithm 5th step: update shrinkage variance parameter ----
      #============================================================================#
      # this is a fusion between SS and HS, group-level HS and inside group simply the earlier SS
      summand_vec <- NULL
      for(g in 1:G) {
        if(g > 1) {
          summand <- sum(beta_cov[(cumulative_M[g - 1] + 1):cumulative_M[g]]^2)
        }else{
          summand <- sum(beta_cov[1:Mg[1]]^2)
        }
        summand_vec <- c(summand_vec, summand)
        tau2[g] <- MCMCpack::rinvgamma(1, (Mstar[g] + 1)/2, scale = 1/epsilon + summand/zeta2/2)
        tau2[g] <- ifelse(tau2[g] > exp(5), exp(5),
                          ifelse(tau2[g] < exp(-5), exp(-5), tau2[g]))    #thresholding
        epsilon[g] <- MCMCpack::rinvgamma(1, 1, scale = 1 + 1/tau2[g])
      }
      zeta2 <- MCMCpack::rinvgamma(1, sum(Mstar)/2,
                                   scale = 1/xi + sum(summand_vec/tau2)/2)
      zeta2 <- ifelse(zeta2 > exp(5), exp(5),
                      ifelse(zeta2 < exp(-5), exp(-5), zeta2))    #thresholding
      xi <-  MCMCpack::rinvgamma(1, 1, scale = 1 + 1/zeta2)
      psi <- t(rep(tau2, times = Mg)*zeta2)


      #============================================================================#
      ## Algorithm 7th step: update phi (icar prior) ----
      #============================================================================#
      priorprec<-1/(sphi^2)*Q                    # Prior precision of phi1
      priormean<- 0
      prec<-priorprec+spam::as.spam(diag(w, nspace, nspace))
      M<-c(w*(z-K%*%beta))                           # note priormean is 0 and only data contributes to M
      phi<-spam::rmvnorm.canonical(1, M, prec)

      # center phi and update tauphi
      phi <- phi-mean(phi)

      tauphi<-rgamma(1, .1+(nspace-1)/2, .1+(phi)%*%Q%*%t(phi)/2)      # n-1 since rank of Q is n-1
      sphi <- sqrt(1/tauphi) #save this and see if it is behaving properly should be 1/nu


      #============================================================================#
      # Store assigning ----
      #============================================================================#
      if (i> burn & i%%thin==0) {
        j<-(i-burn)/thin
        Beta[j,]<-beta
        Delta[j, ]<-delta
        pDelta[j, ]<-pdelta
        Zeta2[j] <- zeta2
        Tau2[j,] <-tau2
        R[j]<-r
        phi.mat[j,]<-phi
        Sphi.sq[j] <- sphi^2                       # spatial random effect variance
      }

      if(verbose == "TRUE"){
        svMisc::progress(i, niter, progress.bar = FALSE)
      }else{if(i%%500 == 0){print(paste0(i, " / ", niter))}
      }
    }

  }else if(which.prior == "HS" & is.null(NeighborhoodList)){
    #algorithm for group HS VS for aspatial count data (without spatial random effect)----

    for (i in 1:niter){
      #============================================================================#
      ## Algorithm 1st step: update the predictor function and probability vector q ----
      #============================================================================#
      eta<-K%*%beta + log(pop_col) - log(scale) # modeling rates

      q<-1/(1+exp(eta)) # dnbinom fn uses q=1-psi
      q<-ifelse(q<=0, 10^-5,q)


      #============================================================================#
      ## Algorithm 2nd step: Chinese restaurant table-based r update ----
      #============================================================================#
      for(j in 1:N) l[j]<-sum(rbinom(y[j],1,round(r/(r+1:y[j]-1),6)))
      r<-rgamma(1,a+sum(l),b-sum(log(q)))


      #============================================================================#
      ## Algorithm 3rd step: Polya-gamma weights and latent response update ----
      #============================================================================#
      w<-BayesLogit::rpg(N,y+r,eta)                               # Polya weights
        z<-((y-r)/(2*w)) - log(pop_col) +log(scale)


      #============================================================================#
      ## Algorithm 4th step: updating posterior covariance of beta and sample betas ----
      #============================================================================#
      # based on priors, first coefficient i.e., intercept is ignored from penalization
      K_w <- sqrt(w)*K
      Sigma_inv <- diag(c(T0, 1/Lambda2/rep(tau2, times = Mg)/zeta2)) + crossprod(K_w) + diag(nugget, p)   # notice the rep command
      M <- t(K_w)%*%(sqrt(w)*z)
      beta <- c(spam::rmvnorm.canonical(1, M, as.matrix(Sigma_inv)))
      beta_cov <- beta[-1]


      #============================================================================#
      ## Algorithm 5th step: update shrinkage parameters ----
      #============================================================================#
      summand_vec <- NULL
      for(g in 1:G) {
        for(j in 1:Mg[g]) {
          if(g > 1){                                                    # g = 1 corresponds to the first group, cumulative_M captures M1, M1+ M2, M1+M2+M3,...
            Lambda2[j + cumulative_M[g - 1]] <- MCMCpack::rinvgamma(1, 1, scale =
                                                                      1/gamma[j + cumulative_M[g - 1]] + beta_cov[j + cumulative_M[g - 1]]^2/2/tau2[g]/zeta2)
            gamma[j + cumulative_M[g - 1]] <- MCMCpack::rinvgamma(1, 1,
                                                                  scale = 1 + 1/Lambda2[j + cumulative_M[g - 1]])
          }else{
            Lambda2[j] <- MCMCpack::rinvgamma(1, 1, scale =               # the if else statement could potentially be merged if you define cumulative_M to store 0 as well, i.e., 0, M1, M1 + M2,...
                                                1/ gamma[j + Mg[1]] + beta_cov[j]^2/tau2[1]/zeta2/2)
            gamma[j] <- MCMCpack::rinvgamma(1, 1, scale =  1 + 1/Lambda2[j])
          }
        }
        if(g > 1) {
          summand <- sum(beta_cov[(cumulative_M[g - 1] + 1):cumulative_M[g]]^2/      # carefully note which sum_j beta_j^2's are being used in each case
                           Lambda2[(cumulative_M[g - 1] + 1):cumulative_M[g]])
        }else{
          summand <- sum(beta_cov[1:Mg[1]]^2/Lambda2[1:Mg[1]])
        }
        summand_vec <- c(summand_vec, summand)
        tau2[g] <- MCMCpack::rinvgamma(1, (Mg[g] + 1)/2, scale = 1/epsilon + summand/zeta2/2)
        tau2[g] <- ifelse(tau2[g] > exp(5), exp(5),
                          ifelse(tau2[g] < exp(-5), exp(-5), tau2[g]))    #thresholding
        epsilon[g] <- MCMCpack::rinvgamma(1, 1, scale = 1 + 1/tau2[g])
      }
      zeta2 <- MCMCpack::rinvgamma(1, p/2,                           # all sum_j beta_j^2 is being used for Zeta update
                                   scale = 1/xi + sum(summand_vec/tau2)/2)
      zeta2 <- ifelse(zeta2 > exp(5), exp(5),
                      ifelse(zeta2 < exp(-5), exp(-5), zeta2))    #thresholding
      xi <-  MCMCpack::rinvgamma(1, 1, scale = 1 + 1/zeta2)


      #============================================================================#
      # Store assigning ----
      #============================================================================#
      if (i> burn & i%%thin==0) {
        j<-(i-burn)/thin
        Beta[j,]<-beta
        Delta[j, ]<-delta
        pDelta[j, ]<-pdelta
        Lambda[j, ]<-Lambda2
        Zeta2[j] <- zeta2
        Tau2[j,] <-tau2
        R[j]<-r
      }

      if(verbose == "TRUE"){
        svMisc::progress(i, niter, progress.bar = FALSE)
      }else{if(i%%500 == 0){print(paste0(i, " / ", niter))}
      }
    }


  }else if(which.prior == "SS" & is.null(NeighborhoodList)){
    #algorithm for group SS VS for aspatial count data (without spatial random effect)----
    for (i in 1:niter){

      #============================================================================#
      ## Algorithm 1st step: update the predictor function and probability vector q ----
      #============================================================================#
      eta<-K%*%beta + log(pop_col) - log(scale) # modeling rates

      q<-1/(1+exp(eta)) # dnbinom fn uses q=1-psi
      q<-ifelse(q<=0, 10^-5,q)  #prevents overly small q values for stability


      #============================================================================#
      ## Algorithm 2nd step: Chinese restaurant table-based r update ----
      #============================================================================#
      for(j in 1:N) l[j]<-sum(rbinom(y[j],1,round(r/(r+1:y[j]-1),6))) # Could try to avoid loop
      r<-rgamma(1,a+sum(l),b-sum(log(q)))


      #============================================================================#
      ## Algorithm 3rd step: Polya-gamma weights and latent response update ----
      #============================================================================#
      w<-BayesLogit::rpg(N,y+r,eta)                               # Polya weights
      z <- ((y-r)/(2*w)) - log(pop_col) +log(scale)


      #============================================================================#
      ## Algorithm 4th step: updating posterior covariance of beta and sample betas ----
      #============================================================================#
      # based on priors, first coefficient i.e., intercept is ignored from penalization

      # this is a fusion between SS and HS, group-level HS and inside group simply the earlier SS
      # we are looping over groups to generate group-level omega for SS
      incfix <- NULL
      for(g in 1:G){
        if(g == 1){
          delta_sub <- delta[2:(Mg[1] + 1)]
        }else{
          delta_sub <- delta[(cumulative_M[g - 1] + 2):(cumulative_M[g] + 1)]
        }
        incfix <- c(incfix, sum(delta_sub == 1))                         # how many delta_j's are nonzero at this iteration, first element is always 1 as it corresponds to the intercept
        omega[g]  <- rbeta(1, wa0 + incfix[g], wb0 + Mg[g] - incfix[g])               # mixture weight from Eq. 9, page 56 from Dvorzak
      }

      invA0 <- diag(c(T0, 1/psi), nrow = p)            # diagonal precision matrix with first element being fixed as it corresponds the intercept

      del_up <- update_delta_g(delta, omega, invA0, z, w, K, p, G, Mg, cumulative_M, a0prior) # delta vector update now loops over each group to update each delta_gi
      delta <- del_up[[1]]                                             # first element of the list correspond to delta_j's
      pdelta <- del_up[[2]]                                            # second element of the list correspond to p(delta_j = 1)'s
      Mstar <- del_up[[3]]                                             # important for nu update

      index <- which(delta == 1)                         # indices of the non-zero delta_j's, first element always selected as it corresponds to the intercept
      invA0 <- invA0[index, index, drop = FALSE]         # we only simulate the non-zero betas, thus only those rows/columns are selected
      Xsel  <- K[, index, drop = FALSE]*sqrt(w) 		     # select the columns of covariate matrix K for which beta_j's are nonzero, adjustment by Polya weights w
      yc <- sqrt(w)*z                                    # adjusting z by Polya weights w, recall that we had the term: (z - K beta)^T diag(w) (z - K beta) inside the exponent,
      # these weight adjustments simplifies it as  (yc - Xsel beta)^T  (yc - Xsel beta), getting rid of the diag(w), note that Xsel is
      # just a transformed version of K, to be specific, a few columns of K for which beta_j's are non-zero

      k_sel <- ncol(Xsel)
      Sigma_inv <- invA0 + t(Xsel)%*%Xsel + diag(nugget, k_sel, k_sel)                            # posterior precision matrix
      sim_beta <- spam::rmvnorm.canonical(1, (t(Xsel)%*%(yc)),       # simulated beta_j's (for only non-zero delta_j's)
                                          as.matrix(Sigma_inv))
      beta[index] <- sim_beta                                        # store the non-zero betas at appropriate indices
      beta[-index] <- 0                                              # rest are simply 0
      beta_cov <- beta[-1]


      #============================================================================#
      ## Algorithm 5th step: update shrinkage variance parameter ----
      #============================================================================#
      # this is a fusion between SS and HS, group-level HS and inside group simply the earlier SS
      summand_vec <- NULL
      for(g in 1:G) {
        if(g > 1) {
          summand <- sum(beta_cov[(cumulative_M[g - 1] + 1):cumulative_M[g]]^2)
        }else{
          summand <- sum(beta_cov[1:Mg[1]]^2)
        }
        summand_vec <- c(summand_vec, summand)
        tau2[g] <- MCMCpack::rinvgamma(1, (Mstar[g] + 1)/2, scale = 1/epsilon + summand/zeta2/2)
        tau2[g] <- ifelse(tau2[g] > exp(5), exp(5),
                          ifelse(tau2[g] < exp(-5), exp(-5), tau2[g]))    #thresholding
        epsilon[g] <- MCMCpack::rinvgamma(1, 1, scale = 1 + 1/tau2[g])
      }
      zeta2 <- MCMCpack::rinvgamma(1, sum(Mstar)/2,
                                   scale = 1/xi + sum(summand_vec/tau2)/2)
      zeta2 <- ifelse(zeta2 > exp(5), exp(5),
                      ifelse(zeta2 < exp(-5), exp(-5), zeta2))    #thresholding
      xi <-  MCMCpack::rinvgamma(1, 1, scale = 1 + 1/zeta2)
      psi <- t(rep(tau2, times = Mg)*zeta2)


      #============================================================================#
      # Store assigning ----
      #============================================================================#
      if (i> burn & i%%thin==0) {
        j<-(i-burn)/thin
        Beta[j,]<-beta
        Delta[j, ]<-delta
        pDelta[j, ]<-pdelta
        Zeta2[j] <- zeta2
        Tau2[j,] <-tau2
        R[j]<-r
      }

      if(verbose == "TRUE"){
        svMisc::progress(i, niter, progress.bar = FALSE)
      }else{if(i%%500 == 0){print(paste0(i, " / ", niter))}
      }
    }

  }

  # Ending VS algorithms
  if(verbose == "TRUE" & which.prior == "HS" & !is.null(NeighborhoodList)){print(paste0("MUSC: NB model fitting with group HS prior for spatial count data with offset completed!"))}
  if(verbose == "TRUE" & which.prior == "SS" & !is.null(NeighborhoodList)){print(paste0("MUSC: NB model fitting with group SS prior for spatial count data with offset completed!"))}
  if(verbose == "TRUE" & which.prior == "HS" & is.null(NeighborhoodList)){print(paste0("MUSC: NB model fitting with group HS prior for aspatial count data with offset completed!"))}
  if(verbose == "TRUE" & which.prior == "SS" & is.null(NeighborhoodList)){print(paste0("MUSC: NB model fitting with group SS prior for aspatial count data with offset completed!"))}

  return(NB = list(Beta = Beta,  r = R, Delta = Delta, pDelta = pDelta,phi = phi.mat, sphiSq=Sphi.sq, Zeta2 =Zeta2, Tau2 =Tau2, Lambda=Lambda))

}



#============================================================================#
# Additional Spike and Slab functions ----
#============================================================================#

#============================================================================#
## Update delta and pdelta ----
#============================================================================#
# delta   		current delta
# omega   		current omega
# invA0			  inverse prior variance (or precision) of beta_j's
# z   		pseudo response vector
# w   		polya weights
# K      	design matrix
# p       number of covariates (including intercept)
# G       number of groups
# Mg      vector of M1, M2, M3, ...
# cumulative_M   vector of M1, M1 + M2, M1 + M2 + M3, ...
# a0prior	prior mean of beta_j's, kept to be a vector of 0's

update_delta_g <- function(delta, omega, invA0, z, w, K, p, G, Mg, cumulative_M, a0prior){
  pdelta <- matrix(NA, 1, p)                                   # initiating a matrix to store p(delta_j = 1)'s
  pdelta[1] <- 1
  Mstar <- NULL
  # first element is 1, as it corresponds to the intercept
  # number of covariates except intercept
  for(g in 1:G) {
    if(g == 1){
      iDel <- 2:(1 + Mg[1])
    }else{iDel <- (cumulative_M[g - 1] + 2):(cumulative_M[g] + 1)}                                                  # first element of delta is always 1, so not perturbed

    ranOrdDelta <- sample(Mg[g])                                 # randomized update order of delta_j's, each delta_j is updated individually
    which_delta_1 <- 0
    for(i in 1:Mg[g]) {
      j         <- iDel[ranOrdDelta[i]]                          # which delta_j is getting updated
      delta.new <- delta                                         # new delta vector is set to be the old delta vector, new delta vector will be perturbed element-wise and
      # MH will be used to accept/reject the new value for every element

      lp <- matrix(0, 2, 1)                                      # store the likelihood value for both cases: 1) new value of the selected delta_j = 0, and
      # 2) new value of the selected delta_j = 1
      for (ii in 0:1){
        delta.new[j] <- ii
        lprior       <- ii*log(omega[g]) + (1 - ii)*log(1 - omega[g]) # prior bernoulli likelihood
        llik <- lmarglik(z, w, K, p, delta.new, a0prior, invA0)       # marginal likelihood
        lp[ii+1] <- llik + lprior
      }
      # maxL  <- max(lp)
      # expl  <- exp(lp - maxL)                                   # maxL adjustment stabilizes the exponent calculation, no other impact on MH
      # lprob <- expl/sum(expl)

      maxL <- max(lp)
      expl <- exp(lp - maxL)
      denom <- sum(expl)
      if (!is.finite(denom) || denom <= 0) {
        lprob <- c(0.5, 0.5)  # safe fallback
      } else {
        lprob <- expl / denom
      }


      deltaj <- runif(1) > lprob[1]                             # first element of lprob is basically p(delta_j = 0)
      if (deltaj != delta[j]) delta[j] <- deltaj                # delta_j update is accepted
      pdelta[j] <- lprob[2]                                     # the second element of lprob is p(delta_j = 1) and is always updated
      which_delta_1 <- which_delta_1 + delta[j]
    }
    Mstar <- c(Mstar, which_delta_1)
  }
  return(list(deltanew = delta, pdeltanew = pdelta, Mstar = Mstar))
}


#============================================================================#
## Marginal likelihood update from Walli and Wagner (also used in pogit package) ----
#============================================================================#

# z   		pseudo response vector
# w   		polya weights
# K      	design matrix
# delta.new is the proposed delta vector
# a0prior	prior mean of beta_j's, kept to be a vector of 0's
# invA0			inverse prior variance (or precision) of beta_j's


lmarglik <- function(z, w, K, p, delta.new, a0prior, invA0){

  index <- which(delta.new == 1)                           # select the rows/columns corresponding to non-zero delta_j's only
  invA0 <- invA0[index, index, drop = FALSE]
  a0    <- a0prior[index]
  Xsel  <- sqrt(w)*K[, index, drop = FALSE]
  yc <- sqrt(w)*z

  cholApost <- chol(invA0 + t(Xsel)%*%Xsel)                # cholesky of posterior precision matrix
  Apost <- Matrix::chol2inv(cholApost)                     # posterior variance matrix
  apost <- (t(Xsel)%*%(yc))                                # posterior mean without the pre-multiplication of Apost on the right,
  # i.e., Apost %*% apost is actually the posterior mean

  # conditional (log) marginal likelihood
  h <- - 2*sum(log(diag(cholApost))) + sum(log(diag(invA0))) # - (-log(det(invA0)))  # difference between posterior and prior log-determinant of variance matrices
  # using the property: log(det(chol(A))) = sum(log(diag(chol(A))))
  Q <- t(yc)%*%yc - t(apost)%*%Apost%*%apost               # difference between log of exponents, notice how not pre-multiplying apost by Apost simplifies the calculation
  lml <- 0.5*(h - Q)                                       # these updates follow Eq. 4 from W and W

  return(lml)
}






