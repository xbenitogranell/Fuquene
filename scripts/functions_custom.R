## Functions from Jovanovska et al. on diatom assembly processes in Lake Ohrid https://github.com/thauffe/OhridDiatomAssembly 
# -----------------------------------------------------------------------------

getSesFric <- function (Comm, TraitsScaled, Runs, Ncores) {
  
  sesFricOneComm <- function(i, Comm, TraitsScaled) {
    S <- propAbundResamp(CommMat = Comm, Index = i - 1)
    TraitsScaledTmp <- TraitsScaled[S, ]
    rownames(TraitsScaledTmp) <- colnames(Comm)
    TraitsDistTmp <- daisy(TraitsScaledTmp, metric = "gower")
    FdTmp <- dbFD(x = TraitsDistTmp,
                  # We need two communities to use this function (i.e. a matrix)
                  # A vector or 1-dimensional matrix for our focal community causes an error
                  # (probably because of some matrix operations)
                  # Hence we create a matrix of our focal community in the first row and a fake community in the second.
                  # We only extract the functional richness of the former and dump the latter.
                  a = rbind(Comm[i, ], # Focal community
                            rep(400, 32)), # Fake community
                  m = 2,
                  w.abun = TRUE, # Weight by species' abundancies
                  ord = "podani",
                  corr = "cailliez", # Make non-euclidean distances euclidean using the cailliez correction
                  calc.FRic = TRUE, # Calculate functional richness
                  # Do not calculate other indices of functional diversity to safe some time
                  stand.FRic = FALSE,
                  calc.CWM = FALSE,
                  calc.FDiv = FALSE,
                  messages = FALSE) # Print no messages or progress on screen
    return(FdTmp$FRic[1])
  }
  # Pairwise distance matrix between species according to their differences in traits.
  # We use the gower distance because many traits are categorical/discrete traits (e.g. samll, medium, large)
  TraitsDist <- daisy(TraitsScaled, metric = "gower")
  # Initialize empty vector to store Functional richness values
  SesFric <- rep(NA_real_, nrow(Comm))
  registerDoParallel(Ncores)
  # For each community starting with the 2nd oldest 
  # (the oldest does not have a community before to do the resampling of species)
  for (i in 2:nrow(Comm)) { 
    # Observed Functional richness of the focal community (i.e. at one moment in time)
    FdObs <- dbFD(x = TraitsDist,
                  # We need two communities to use this function (i.e. a matrix)
                  # A vector or 1-dimensional matrix for our focal community causes an error
                  # (probably because of some matrix operations)
                  # Hence we create a matrix of our focal community in the first row and a fake community in the second.
                  # We only extract the functional richness of the former and dump the latter.
                  a = rbind(Comm[i, ], # Focal community
                            rep(400, 32)), # Fake community
                  m = 2,
                  w.abun = TRUE, # Weight by species' abundancies
                  ord = "podani",
                  corr = "cailliez", # Make non-euclidean distances euclidean using the cailliez correction
                  calc.FRic = TRUE, # Calculate functional richness
                  # Do not calculate other indices of functional diversity to safe some time
                  stand.FRic = FALSE, 
                  calc.CWM = FALSE,
                  calc.FDiv = FALSE,
                  messages = FALSE) # Print no messages or progress on screen
    # Get Runs-times (e.g. 999x) a functional richness for a community assembled according to the Null model
    # We use parallel computation on Ncores to speed this up
    FdNull <- foreach(iter = 1:Runs, 
                      .combine = c, # combine results in vector
                      .packages = c("FD", "cluster"), # Packages needed
                      # There is no specific order of the Null communities
                      .inorder = FALSE) %dopar%  sesFricOneComm(i, Comm, TraitsScaled) # Get functional richness for the i-th sample in our community matrix
    # Calculate standardized effect size of functional richness
    # (Observed FR - mean FR of the Null model) / standard deviation of the Null model
    SesFric[i] <- (FdObs$FRic[1] - mean(FdNull, na.rm = TRUE)) / sd(FdNull, na.rm = TRUE)
  }
  stopImplicitCluster()
  return(SesFric)
}
