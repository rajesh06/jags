  #------------------------------------------------------------------------------

  # Jags-Ymet-Xnom1grp-MlogNormal-Script.R
  # April 19, 2016. John K. Kruschke.
  # Requires auxiliary R scripts that accompany the book,
  # Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition:
  # A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

  # graphics.off()
  # rm(list=ls(all=TRUE))

source("D:/Users/Rajesh/Documents/Bayesian Data Analysis/DBDA2Eprograms/DBDA2E-utilities.R")
fileNameRoot = "Jags-Ymet-Xnom1grp-MlogNormal-Script"
graphFileType = "png"

#------------------------------------------------------------------------------
# THE DATA.

# Generate some random data from known parameter values:

set.seed(47405)
trueLogM = 5.0 # true mean of log(y)
trueLogSD = 0.5 # true sd of log(y)

y = rnorm(n = 125)  # standard normal distribution of log-scale values
hist(y)

LogY = (y - mean(y))/sd(y) * trueLogSD + trueLogM
hist(LogY)

y = exp(LogY) # y is exponentiated values of log(y)
hist(y)

# OR JUST PUT IN YOUR ACTUAL DATA HERE:
# y = ...
# Package the data for shipping to JAGS:

dataList = list(
  y = y,
  N = length(y) ,
  meanOfLogY = mean(log(y)) ,
  sdOfLogY = sd(log(y))
)

#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
model {
  for( i in 1 : N ) {
  y[i] ~ dlnorm( muOfLogY , 1/sigmaOfLogY^2 )
  }

  sigmaOfLogY ~ dunif( 0.001*sdOfLogY , 1000*sdOfLogY )
  muOfLogY ~ dunif( 0.001*meanOfLogY , 1000*meanOfLogY )

  muOfY <- exp(muOfLogY + sigmaOfLogY^2/2)
  modeOfY <- exp(muOfLogY - sigmaOfLogY^2)
  sigmaOfY <- sqrt(exp(2 * muOfLogY + sigmaOfLogY^2) * (exp(sigmaOfLogY^2) - 1))
}
" # close quote for modelstring
writeLines(modelstring, con = "model.txt")

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.
# Let JAGS do it

#------------------------------------------------------------------------------
# RUN THE CHAINS

require(rjags)

parameters = c("muOfLogY" , "sigmaOfLogY" , "muOfY" , "modeOfY" , "sigmaOfY" )
adaptSteps = 1000         # Number of steps to "tune" the samplers.
burnInSteps = 1000        # Number of steps to "burn-in" the samplers.
nChains = 3               # Number of chains to run.
numSavedSteps = 20000     # Total number of steps in chains to save.
thinSteps = 1             # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling((numSavedSteps * thinSteps) / nChains) # Steps per chain.
# Create, initialize, and adapt the model:

jagsModel = jags.model( "model.txt" , data = dataList ,
  n.chains = nChains, n.adapt = adaptSteps )

# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update(jagsModel, n.iter = burnInSteps )

# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
mcmcCoda = coda.samples(model = jagsModel, variable.names = parameters ,
  n.iter = nPerChain, thin = thinSteps )

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names
for (parName in parameterNames) {
  diagMCMC(codaObject = mcmcCoda , parName = parName,
    saveName = fileNameRoot, saveType = graphFileType )
}


# Convert coda-object codaSamples to matrix object for easier handling.
mcmcChain = as.matrix(mcmcCoda)
chainLength = NROW(mcmcChain)

openGraph(width = 10, height = 6)
layout(matrix(1:6, nrow = 2, byrow = TRUE))
# posterior predictive
hist(dataList$y , xlab = "y" , main = "Data w. Post. Pred.", breaks = 30,
  col = "pink" , border = "white", prob = TRUE , cex.lab = 1.5)

pltIdx = floor(seq(1, chainLength, length = 20))

xComb = seq(min(dataList$y), max(dataList$y) , length=501 )

for (chnIdx in pltIdx) {
  lines(xComb ,
    dlnorm( xComb, mcmcChain[chnIdx,"muOfLogY"], mcmcChain[chnIdx,"sigmaOfLogY"] ),
    col = "skyblue" )
}
# param's of log(y)
postInfo = plotPost(mcmcChain[,"muOfLogY"] , xlab = "mu of log(y)" )
postInfo = plotPost(mcmcChain[,"sigmaOfLogY"] , xlab = "sigma of log(y)" )
# param's of y
postInfo = plotPost(mcmcChain[,"modeOfY"] , xlab = "mode of y" )
postInfo = plotPost(mcmcChain[,"muOfY"] , xlab = "mu of y" )
postInfo = plotPost( mcmcChain[,"sigmaOfY"] , xlab="sigma of y" , cenTend="mode")

saveGraph(file = fileNameRoot, type=graphFileType)

#-------------------------------------------------------------------------------
