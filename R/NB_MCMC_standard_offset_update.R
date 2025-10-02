#' @title Function for the NB variable selection (HS or SS) with standard updates (with or without space)
#'
#' @param K Design matrix including an intercept column (first column all 1s)
#' @param y is the vector of the count outcome
#' @param NeighborhoodList is an object of class \code{nb}, typically created with \code{\link[spdep]{poly2nb}}. Each element of the list corresponds to a spatial unit and contains the indices of its neighboring units, as defined by the chosen contiguity rile (e.g queen or rook). Used to encode the spatial adjacency structure for analysis.
#' @param which.prior denotes the model to be used, "HS" or "SS", horseshoe and spike-and-slab respectively
#' @param niter is the number of MCMC iterations
#' @param verbose is TRUE or FALSE, to show or suppress progression updates
#' @param pop_col is a vector of the population of each spatial unit for the model offset, if modeling with space
#' @param scale is an integer of the desired scale for outcome rates (e.g 100000 for per 100,000 rates)

#' @return It returns the estimated beta's, r, shrinkage parameters, phi, Delta and pDelta (if applicable )
#' @export


VS_Standard_offset <- function(K, y, NeighborhoodList, which.prior,niter, verbose, pop_col, scale){
  # currently all algorithms update r with Chinese restaurant table (CRT) and have a burn-in of half

  # variable definitions
  N = nrow(K)                                        # number of observations
  p = ncol(K)                                       # number of features (includes intercept)

  # Spatial variables
  if(!is.null(NeighborhoodList)){
    neighbor <- NeighborhoodList
    nspace <- length(neighbor) # Number of spatial units
    A <- spdep::nb2mat(neighbor, style = "B")
    neighbor_count <- rowSums(A)  # M (diagonal matrix) neighbor count for each spatial unit (notation from Mutiso, Neelon 2022)
    D <- diag(neighbor_count)  # Create diagonal matrix... Degree Matrix
    Q <- D - A     # CAR precision matrix

    #spatial inits
    #spatial parameter inits
    phi_init <- c(spam::rmvnorm(1, sigma=diag(.01, nspace)))	          # Random effects
    phi <- phi_init-mean(phi_init) # centered phi
    s2phi <- var(phi)              # var phi
    sphi <- sqrt(s2phi)            # std phi
    tauphi <- 1/s2phi              # precision of phi
  }

  # MCMC features
  thin <- 1				            # Thinning interval
  burn <- niter/2		          # Burnin
  lastit <- (niter-burn)/thin	# Last stored value

  # Storage objects
  Beta <- matrix(0,lastit,p)              # store betas
  Delta <- pDelta<-matrix(0,lastit,p)     # store SS inclusion indicator and probability of inclusion
  Lambda <- matrix(0,lastit,p-1)            # store local shrinkage parameters
  Tau2 <- matrix(0,lastit,1)                # store global shrinkage parameter
  R <- matrix(0,lastit,1)                 # store over-dispersion parameter
  Sphi.sq <- matrix(0,lastit,1)           # store variance of spatial random effect, phi
  phi.mat <- matrix(0,lastit,length(y))     # store spatial random effect

  # Initialize
  beta <- rep(0, p)
  r <- 1
  Lambda2 <- rep(0.1, p-1)
  gamma <- rep(0.1, p-1)
  tau2 = 0.1
  del1 <- 1
  del2 <- 1                  # lambda hyper parameters
  nugget <- 0.05              # small diagonal element to be added to the covariance
  T0 <- 0.1                # prior on features
  s <- 0.01
  epsilon<-1
  wa0<-1
  wb0<-1
  delta <- c(1, rep(0, p-1))
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


  l <- rep(0, N)      # latent vector for CRT-based update of r
  a <- b <- 0.01      # Gamma hyperparms for r for CRT
  Acc <- 0            # counter for MH based update for r, unused

  # Starting VS algorithms
  if(verbose == "TRUE" & which.prior == "HS" & !is.null(NeighborhoodList)){print(paste0("MUSC: NB model fitting has started with standard HS prior for spatial count data with offset."))}
  if(verbose == "TRUE" & which.prior == "SS" & !is.null(NeighborhoodList)){print(paste0("MUSC: NB model fitting has started with standard SS prior for spatial count data with offset."))}
  if(verbose == "TRUE" & which.prior == "HS" & is.null(NeighborhoodList)){print(paste0("MUSC: NB model fitting has started with standard HS prior for aspatial count data with offset."))}
  if(verbose == "TRUE" & which.prior == "SS" & is.null(NeighborhoodList)){print(paste0("MUSC: NB model fitting has started with standard SS prior for aspatial count data with offset."))}

  if(which.prior == "HS" & !is.null(NeighborhoodList)){
    # Algorithm for standard HS VS for spatial count data (with spatial random effect)----

    for (i in 1:niter){
      #============================================================================#
      ## Algorithm 1st step: update the predictor function and probability vector q ----
      #============================================================================#
      eta <- K%*%beta + rep(phi,1) + log(pop_col) - log(scale)  #modeling rates with offset for population and scale

      q<-1/(1+exp(eta)) # dnbinom fn uses q=1-psi
      q<-ifelse(q<=0, 10^-5,q) # prevents overly small q values for stability


      #============================================================================#
      ## Algorithm 2nd step: Chinese restaurant table-based r update ----
      #============================================================================#
      for(j in 1:N) l[j] <- sum(rbinom(y[j],1,round(r/(r+1:y[j]-1),6)))
      r <- rgamma(1,a+sum(l),b-sum(log(q)))


      #============================================================================#
      ## Algorithm 3rd step: Polya-gamma weights and latent response update ----
      #============================================================================#
      w <- BayesLogit::rpg(N,y+r,eta)                               # Polya weights
      z <- ((y-r)/(2*w)) - log(pop_col) +log(scale)     # latent response


      #============================================================================#
      ## Algorithm 4th step: updating posterior covariance of beta and sample betas ----
      #============================================================================#
      # based on priors, first coefficient i.e., intercept is ignored from penalization
      Sigma = diag(c(T0, 1/Lambda2/tau2)) +
        crossprod(K*sqrt(w)) +
        diag(nugget, p)

      beta = c(spam::rmvnorm.canonical(1, t(sqrt(w)*K)%*%(sqrt(w)*(z-rep(phi, 1))), as.matrix(Sigma)))


      #============================================================================#
      ## Algorithm 5th step: update local shrinkage parameter ----
      #============================================================================#
      for(j in 2:p) {
        Lambda2[j-1]<-MCMCpack::rinvgamma(1, 1, scale =
                                            1/gamma[j-1] + beta[j]^2/2/tau2)
        gamma[j-1]<-MCMCpack::rinvgamma(1, 1, scale =
                                          1 + 1/Lambda2[j-1])
      }

      #============================================================================#
      ## Algorithm 6th step: update global shrinkage parameter ----
      #============================================================================#
      tau2 <- MCMCpack::rinvgamma(1, p/2,
                                  scale = 1/epsilon + sum(beta[-1]^2/Lambda2)/2)
      tau2 <- ifelse(tau2 > exp(5), exp(5),
                        ifelse(tau2 < exp(-5), exp(-5), tau2))    #thresholding
      epsilon <-MCMCpack::rinvgamma(1, 1, scale = 1 + 1/tau2)


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
      sphi <- sqrt(1/tauphi)

      #============================================================================#
      # Store assigning ----
      #============================================================================#
      if (i> burn & i%%thin==0 ) {
        j <- (i-burn)/thin
        Beta[j,] <- beta
        R[j] <- r
        Tau2[j] <-tau2
        Lambda[j,] <- Lambda2
        phi.mat[j,] <- phi
        Sphi.sq[j] <- sphi^2                       # spatial random effect variance
      }

      if(verbose == "TRUE"){
        svMisc::progress(i, niter, progress.bar = FALSE)
      }else{if(i%%500 == 0){print(paste0(i, " / ", niter))}
      }
    }


  }else if(which.prior == "SS" & !is.null(NeighborhoodList)){
    #algorithm for standard SS VS for spatial count data (with spatial random effect) ----

    for (i in 1:niter){
      #============================================================================#
      ## Algorithm 1st step: update the predictor function and probability vector q ----
      #============================================================================#
      eta <- K%*%beta + rep(phi,1) + log(pop_col) - log(scale)  #modeling rates with offset for population and scale

      q<-1/(1+exp(eta)) # dnbinom fn uses q=1-psi
      q<-ifelse(q<=0, 10^-5,q) # prevents overly small q values for stability


      #============================================================================#
      ## Algorithm 2nd step: Chinese restaurant table-based r update ----
      #============================================================================#
      for(j in 1:N) l[j] <- sum(rbinom(y[j],1,round(r/(r+1:y[j]-1),6)))
      r <- rgamma(1,a+sum(l),b-sum(log(q)))


      #============================================================================#
      ## Algorithm 3rd step: Polya-gamma weights and latent response update ----
      #============================================================================#
      w <- BayesLogit::rpg(N,y+r,eta)                               # Polya weights
      z <- ((y-r)/(2*w)) - log(pop_col) +log(scale)     # latent response


      #============================================================================#
      ## Algorithm 4th step: updating posterior covariance of beta and sample betas ----
      #============================================================================#
      # based on priors, first coefficient i.e., intercept is ignored from penalization
      incfix <- sum(delta == 1)                          # how many delta_j's are nonzero at this iteration, first element is always 1 as it corresponds to the intercept
      omega  <- rbeta(1, wa0 + incfix, wb0 + p - incfix) # mixture weight from Eq. 9, page 56 from Dvorzak
      invA0 <- diag(c(1/T0, 1/psi), nrow = p)            # diagonal precision matrix with first element being fixed as it corresponds the intercept

      del_up <- update_delta(delta, omega, invA0, z, w, K, p, a0prior) # delta vector update
      delta <- del_up[[1]]                                             # first element of the list correspond to delta_j's
      pdelta <- del_up[[2]]                                            # second element of the list correspond to p(delta_j = 1)'s

      index <- which(delta == 1)                         # indices of the non-zero delta_j's, first element always selected as it corresponds to the intercept
      invA0 <- invA0[index, index, drop = FALSE]         # we only simulate the non-zero betas, thus only those rows/columns are selected
      Xsel  <- sqrt(w)*K[, index, drop = FALSE] 		     # select the columns of covariate matrix K for which beta_j's are nonzero, adjustment by Polya weights w
      yc <- sqrt(w)*z                                    # adjusting z by Polya weights w, recall that we had the term: (z - K beta)^T diag(w) (z - K beta) inside the exponent,
      # these weight adjustments simplifies it as  (yc - Xsel beta)^T  (yc - Xsel beta), getting rid of the diag(w), note that Xsel is
      # just a transformed version of K, to be specific, a few columns of K for which beta_j's are non-zero
      Sigma_inv <- invA0 + t(Xsel)%*%Xsel                            # posterior precision matrix
      sim_beta <- spam::rmvnorm.canonical(1, (t(Xsel)%*%(yc)),       # simulated beta_j's (for only non-zero delta_j's)
                                          as.matrix(Sigma_inv))
      beta[index] <- sim_beta                                        # store the non-zero betas at appropriate indices
      beta[-index] <- 0                                              # rest are simply 0


      #============================================================================#
      ## Algorithm 5th step: update psi ----
      #============================================================================#
      psiv <- 1/rgamma(p - 1, shape = psi.nu + t(delta[-1])/2,
                       rate = (psi.Q + 0.5*beta[-1]^2))                   # simulate variance terms psi_j's from Eq.9 from Dvorzak, (except index 1 which corresponds to the intercept)
      psi <- t(psiv)


      #============================================================================#
      ## Algorithm 6th step: update phi (icar prior) ----
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
      if (i> burn & i%%thin==0 ) {
        j <- (i-burn)/thin
        Beta[j,] <- beta
        R[j] <- r
        Tau2[j] <-tau2
        Lambda[j,] <- Lambda2
        Delta[j, ]<-delta
        pDelta[j, ]<-pdelta
        phi.mat[j,] <- phi
        Sphi.sq[j] <- sphi^2                       # spatial random effect variance
      }

      if(verbose == "TRUE"){
        svMisc::progress(i, niter, progress.bar = FALSE)
      }else{if(i%%500 == 0){print(paste0(i, " / ", niter))}
      }
    }


  }else if(which.prior == "HS" & is.null(NeighborhoodList)){
    # Algorithm for standard HS VS for aspatial count data (without spatial random effect)----

    for (i in 1:niter){
      #============================================================================#
      ## Algorithm 1st step: update the predictor function and probability vector q ----
      #============================================================================#
      eta <- K%*%beta + log(pop_col) - log(scale)  #modeling rates with offset for population and scale

      q<-1/(1+exp(eta)) # dnbinom fn uses q=1-psi
      q<-ifelse(q<=0, 10^-5,q) # prevents overly small q values for stability


      #============================================================================#
      ## Algorithm 2nd step: Chinese restaurant table-based r update ----
      #============================================================================#
      for(j in 1:N) l[j] <- sum(rbinom(y[j],1,round(r/(r+1:y[j]-1),6)))
      r <- rgamma(1,a+sum(l),b-sum(log(q)))


      #============================================================================#
      ## Algorithm 3rd step: Polya-gamma weights and latent response update ----
      #============================================================================#
      w <- BayesLogit::rpg(N,y+r,eta)                               # Polya weights
      z <- ((y-r)/(2*w)) - log(pop_col) +log(scale)     # latent response


      #============================================================================#
      ## Algorithm 4th step: updating posterior covariance of beta and sample betas ----
      #============================================================================#
      # based on priors, first coefficient i.e., intercept is ignored from penalization
      Sigma = diag(c(T0, 1/Lambda2/tau2)) +
        crossprod(K*sqrt(w)) +
        diag(nugget, p)

      beta = c(spam::rmvnorm.canonical(1, t(sqrt(w)*K)%*%(sqrt(w)*(z-rep(phi, 1))), as.matrix(Sigma)))


      #============================================================================#
      ## Algorithm 5th step: update local shrinkage parameter ----
      #============================================================================#
      for(j in 2:p) {
        Lambda2[j-1]<-MCMCpack::rinvgamma(1, 1, scale =
                                            1/gamma[j-1] + beta[j]^2/2/tau2)
        gamma[j-1]<-MCMCpack::rinvgamma(1, 1, scale =
                                          1 + 1/Lambda2[j-1])
      }

      #============================================================================#
      ## Algorithm 6th step: update global shrinkage parameter ----
      #============================================================================#
      tau2 <- MCMCpack::rinvgamma(1, p/2,
                                  scale = 1/epsilon + sum(beta[-1]^2/Lambda2)/2)
      tau2 <- ifelse(tau2 > exp(5), exp(5),
                     ifelse(tau2 < exp(-5), exp(-5), tau2))    #thresholding
      epsilon <-MCMCpack::rinvgamma(1, 1, scale = 1 + 1/tau2)


      #============================================================================#
      # Store assigning ----
      #============================================================================#
      if (i> burn & i%%thin==0 ) {
        j <- (i-burn)/thin
        Beta[j,] <- beta
        R[j] <- r
        Tau2[j] <-tau2
        Lambda[j,] <- Lambda2
      }

      if(verbose == "TRUE"){
        svMisc::progress(i, niter, progress.bar = FALSE)
      }else{if(i%%500 == 0){print(paste0(i, " / ", niter))}
      }
    }


  }else if(which.prior == "SS" & is.null(NeighborhoodList)){
    #algorithm for standard SS VS for aspatial count data (without spatial random effect) ----

    for (i in 1:niter){
      #============================================================================#
      ## Algorithm 1st step: update the predictor function and probability vector q ----
      #============================================================================#
      eta <- K%*%beta + log(pop_col) - log(scale)  #modeling rates with offset for population and scale

      q<-1/(1+exp(eta)) # dnbinom fn uses q=1-psi
      q<-ifelse(q<=0, 10^-5,q) # prevents overly small q values for stability


      #============================================================================#
      ## Algorithm 2nd step: Chinese restaurant table-based r update ----
      #============================================================================#
      for(j in 1:N) l[j] <- sum(rbinom(y[j],1,round(r/(r+1:y[j]-1),6)))
      r <- rgamma(1,a+sum(l),b-sum(log(q)))


      #============================================================================#
      ## Algorithm 3rd step: Polya-gamma weights and latent response update ----
      #============================================================================#
      w <- BayesLogit::rpg(N,y+r,eta)                               # Polya weights
      z <- ((y-r)/(2*w)) - log(pop_col) +log(scale)     # latent response


      #============================================================================#
      ## Algorithm 4th step: updating posterior covariance of beta and sample betas ----
      #============================================================================#
      # based on priors, first coefficient i.e., intercept is ignored from penalization
      incfix <- sum(delta == 1)                          # how many delta_j's are nonzero at this iteration, first element is always 1 as it corresponds to the intercept
      omega  <- rbeta(1, wa0 + incfix, wb0 + p - incfix) # mixture weight from Eq. 9, page 56 from Dvorzak
      invA0 <- diag(c(1/T0, 1/psi), nrow = p)            # diagonal precision matrix with first element being fixed as it corresponds the intercept

      del_up <- update_delta(delta, omega, invA0, z, w, K, p, a0prior) # delta vector update
      delta <- del_up[[1]]                                             # first element of the list correspond to delta_j's
      pdelta <- del_up[[2]]                                            # second element of the list correspond to p(delta_j = 1)'s

      index <- which(delta == 1)                         # indices of the non-zero delta_j's, first element always selected as it corresponds to the intercept
      invA0 <- invA0[index, index, drop = FALSE]         # we only simulate the non-zero betas, thus only those rows/columns are selected
      Xsel  <- sqrt(w)*K[, index, drop = FALSE] 		     # select the columns of covariate matrix K for which beta_j's are nonzero, adjustment by Polya weights w
      yc <- sqrt(w)*z                                    # adjusting z by Polya weights w, recall that we had the term: (z - K beta)^T diag(w) (z - K beta) inside the exponent,
      # these weight adjustments simplifies it as  (yc - Xsel beta)^T  (yc - Xsel beta), getting rid of the diag(w), note that Xsel is
      # just a transformed version of K, to be specific, a few columns of K for which beta_j's are non-zero

      Sigma_inv <- invA0 + t(Xsel)%*%Xsel                            # posterior precision matrix
      sim_beta <- spam::rmvnorm.canonical(1, (t(Xsel)%*%(yc)),       # simulated beta_j's (for only non-zero delta_j's)
                                          as.matrix(Sigma_inv))
      beta[index] <- sim_beta                                        # store the non-zero betas at appropriate indices
      beta[-index] <- 0                                              # rest are simply 0


      #============================================================================#
      ## Algorithm 5th step: update psi ----
      #============================================================================#
      psiv <- 1/rgamma(p - 1, shape = psi.nu + t(delta[-1])/2,
                       rate = (psi.Q + 0.5*beta[-1]^2))                   # simulate variance terms psi_j's from Eq.9 from Dvorzak, (except index 1 which corresponds to the intercept)
      psi <- t(psiv)


      #============================================================================#
      # Store assigning ----
      #============================================================================#
      if (i> burn & i%%thin==0 ) {
        j <- (i-burn)/thin
        Beta[j,] <- beta
        R[j] <- r
        Tau2[j] <-tau2
        Lambda[j,] <- Lambda2
        Delta[j, ]<-delta
        pDelta[j, ]<-pdelta
      }

      if(verbose == "TRUE"){
        svMisc::progress(i, niter, progress.bar = FALSE)
      }else{if(i%%500 == 0){print(paste0(i, " / ", niter))}
      }
    }
  }

  # Ending VS algorithms
  if(verbose == "TRUE" & which.prior == "HS" & !is.null(NeighborhoodList)){print(paste0("MUSC: NB model fitting with standard HS prior for spatial count data with offset completed!"))}
  if(verbose == "TRUE" & which.prior == "SS" & !is.null(NeighborhoodList)){print(paste0("MUSC: NB model fitting with standard SS prior for spatial count data with offset completed!"))}
  if(verbose == "TRUE" & which.prior == "HS" & is.null(NeighborhoodList)){print(paste0("MUSC: NB model fitting with standard HS prior for aspatial count data with offset completed!"))}
  if(verbose == "TRUE" & which.prior == "SS" & is.null(NeighborhoodList)){print(paste0("MUSC: NB model fitting with standard SS prior for aspatial count data with offset completed!"))}

  return(NB = list(Beta = Beta,  r = R, Delta = Delta, pDelta = pDelta,phi = phi.mat, sphiSq=Sphi.sq, Tau2 =Tau2, Lambda=Lambda))

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
# a0prior	prior mean of beta_j's, kept to be a vector of 0's

update_delta<- function(delta, omega, invA0, z, w, K, p, a0prior){
  nDelta <- p - 1                                                # number of covariates except intercept
  ## update indicators for regression effects
  if (nDelta > 0){                                               # trivially true in our cases
    iDel <- 2:p                                                  # first element of delta is always 1, so not perturbed
    ranOrdDelta <- sample(p - 1)                                 # randomized update order of delta_j's, each delta_j is updated individually
    pdelta <- matrix(NA, 1, p)                                   # initiating a matrix to store p(delta_j = 1)'s
    pdelta[1] <- 1                                               # first element is 1, as it corresponds to the intercept

    for (i in 1:nDelta){
      j         <- iDel[ranOrdDelta[i]]                          # which delta_j is getting updated
      delta.new <- delta                                         # new delta vector is set to be the old delta vector, new delta vector will be perturbed element-wise and
      # MH will be used to accept/reject the new value for every element

      lp <- matrix(0, 2, 1)                                      # store the likelihood value for both cases: 1) new value of the selected delta_j = 0, and
      # 2) new value of the selected delta_j = 1
      for (ii in 0:1){
        delta.new[j] <- ii
        lprior       <- ii*log(omega) + (1 - ii)*log(1 - omega) # prior bernoulli likelihood
        llik <- lmarglik(z, w, K, p, delta.new, a0prior, invA0) # marginal likelihood
        lp[ii+1] <- llik + lprior
      }
      maxL  <- max(lp)
      expl  <- exp(lp - maxL)                                   # maxL adjustment stabilizes the exponent calculation, no other impact on MH
      lprob <- expl/sum(expl)

      deltaj <- runif(1) > lprob[1]                             # first element of lprob is basically p(delta_j = 0)
      if (deltaj != delta[j]) delta[j] <- deltaj                # delta_j update is accepted
      pdelta[j] <- lprob[2]                                     # the second element of lprob is p(delta_j = 1) and is always updated

    }
  }else {
    pdelta <- matrix(1, 1, p)
  }
  return(list(deltanew = delta, pdeltanew = pdelta))
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
  h <- - 2*sum(log(diag(cholApost))) - (-log(det(invA0)))  # difference between posterior and prior log-determinant of variance matrices
  # using the property: log(det(chol(A))) = sum(log(diag(chol(A))))
  Q <- t(yc)%*%yc - t(apost)%*%Apost%*%apost               # difference between log of exponents, notice how not pre-multiplying apost by Apost simplifies the calculation
  lml <- 0.5*(h - Q)                                       # these updates follow Eq. 4 from W and W

  return(lml)
}










