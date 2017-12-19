# functions to assist in the psuedo-crossvalidation exercise

#Libs
library(fda)

#GlobalEnv Context
data("growth")

### functions

# function that takes 1 input number, removes that persons data from the 
# DF, and splices the first half of it onto the end of the DF
# doing this with the data from the males 

dataRearranger <- function(indiv.index, obs.to.drop) {
  indiv.data <- growth$hgtm[ ,indiv.index]
  temp.mat.1 <- growth$hgtm[, -indiv.index]
  temp.df <- data.frame(x = growth$age, y = c(as.numeric(temp.mat.1), indiv.data))
  temp.df <- temp.df[1:(length(temp.df$x) - obs.to.drop),]
  return(temp.df)
}

# generates a person vector / position vector for the obsertvations 
posVecGenerator <- function(obs.to.drop) {
  n.indiv <- dim(growth$hgtm)[2]
  n.obs <- dim(growth$hgtm)[1]
  pos.vec.male <- c(1, (1:n.indiv * n.obs) + 1)
  pos.vec.male[length(pos.vec.male)] <- pos.vec.male[length(pos.vec.male)] - obs.to.drop
  return(pos.vec.male)
}

# XNew sequence generator
newXvalGenerator <- function(n.x.new) {
  x.new <- seq(from = min(growth$age), to = max(growth$age), length.out = n.x.new)
  return(x.new)
}

# fit the stan model and return the stan object
fitREFmodel <- function(indiv.index, obs.to.drop, q, mcmc.iter, rand.param.count, n.x.new,
                        prefit = ref.mod.prefit, chains) {
  if (q %% 2 == 0) {
    k.1 <- (q / 2) - 1
    k.2 <- (q / 2) - 1
  } else {
    k.1 <- (q - 1) / 2
    k.2 <- k.1 - 1
  }
  
  # get data
  # this data needs to be "in the loop" while this pos.vec doesnt  
    data.df <- dataRearranger(indiv.index, obs.to.drop)
    pos.vec <- posVecGenerator(obs.to.drop)
    x.new   <- newXvalGenerator(n.x.new)
    
  # fit model
    data.in <- list(N_x = length(data.df$x), N_i = length(pos.vec) - 1,
                    x_vector = data.df$x, y_vector = data.df$y,
                    position_vector = pos.vec, q = q, k_1 = k.1, k_2 = k.2,
                    alpha = 1, rand_params = rand.param.count,
                    N_x_new =  n.x.new, x_new_vector = x.new)
    
    model.fit <- stan(fit = prefit, data = data.in, chains = chains, cores = chains,
                      iter = mcmc.iter, refresh = 250, warmup = 0.975*mcmc.iter, thin = 5, save_warmup = FALSE)
    
    return(model.fit) 
  # what to return? The stanfit object?
  # yea might as well 
  
}

## all of these functions are written with the knowledge that the last person is sans 
## obs.to.drop number of observations

# this function is one of the worst bits of code i've written,
# as it is just hacked from the full case
# but it still works so????
indivPlotAndDataProducer <- function(fit.object, indiv.index, obs.to.drop, num.indiv, q, x.new.vec) { # unsure if I will need the pos- vec
  # this function should print the combo plot to pdf 
  # and write out the csf files of the fits.
  person.vec <- c()
  for (i in 1:num.indiv) {
    person.vec <- c(person.vec, rep(i, 31))
  }
  full.plot.df <- data.frame(x = dataRearranger(indiv.index, 0)$x, 
                             y = dataRearranger(indiv.index, 0)$y,
                             person.vec = as.factor(person.vec))
  
  
  colour.vec <- c(rep("blue", length(dataRearranger(indiv.index, obs.to.drop)$x)), 
                  rep("red", length(full.plot.df$x) - length(dataRearranger(indiv.index, obs.to.drop)$x)))
  
  full.plot.df$colour <- colour.vec
  
  # make the xmatrix because it's easier than using evalpoly
  x.mat <- cbind(1, growth$age)
  
  for (ii in 2:q) {
    x.mat <- cbind(x.mat, growth$age ^ ii)
  }
  
  # extract relevant beta information 
  beta.final.mat.samples <- extract(fit.object, "beta_final")
  beta.final.mat <- apply(beta.final.mat.samples[[1]], c(2:3), mean)
  mu.pred.samples <- extract(fit.object, "mu_new")
  samples.of.interest <- mu.pred.samples[[1]][,39,]
  
  # get the lower and upper credible interval vectors 
  
  fitted.values.temp.1 <- beta.final.mat.samples[[1]][,39,] %*% t(x.mat)
  full.plot.df$lower.interval.temp <- c(rep(0, 38*length(growth$age)) ,apply(fitted.values.temp.1, 2, function(X){quantile(X, 0.025)}))  
  full.plot.df$upper.interval.temp <- c(rep(0, 38*length(growth$age)) ,apply(fitted.values.temp.1, 2, function(X){quantile(X, 0.975)}))
  
  # get the fitted values 
  fitted.vals <- c()
  for (i in 1:num.indiv) {
    fitted.vals <- rbind(fitted.vals, x.mat %*% beta.final.mat[i, ])
  }
  full.plot.df$fitted.vals <- fitted.vals
  
  ## legend key
  key=list(space="top", columns = 2,
           lines=list(col=c("red","black", "blue", "springgreen4"), lty=c(1,1, 2, 2), lwd=3),
           text=list(c("Mean (Beta)","Cred/Pred int (Beta)", "Mean (MuPred)", "Cred/Pred (MuPred)"))
  )
  
  # make the individual plot
  temp.plot.df <- subset(full.plot.df, person.vec == 39)
  temp.plot.df$person.vec <- as.factor(indiv.index) # change the factor for plotting purposes.
  indiv.plot.1 <- xyplot(y ~ x|person.vec, data = temp.plot.df, type = "p", 
                         groups = colour.vec[(length(colour.vec) - length(growth$age) + 1):length(colour.vec)],  
                         xlab = "Age (Years)", ylab = "Height (cm)", strip = T, 
                         pch =  21, cex = 1.2, key = key )
  indiv.plot.2 <- xyplot(fitted.vals ~ x, data = temp.plot.df, type = "l", col = "Red")
  indiv.plot.3 <- xyplot(lower.interval.temp  ~ x, data = temp.plot.df, type = "l", col = "black")
  indiv.plot.4 <- xyplot( upper.interval.temp ~ x, data = temp.plot.df, type = "l", col = "black")
  mu.pred.mean <- apply(samples.of.interest,2,mean)
  mu.pred.lower.bound <- apply(samples.of.interest, 2, function(X){quantile(X, 0.025)})
  mu.pred.upper.bound <- apply(samples.of.interest, 2, function(X){quantile(X, 0.975)})
  
  temp.df.2 <- data.frame(x = x.new.vec, mean = mu.pred.mean, upper = mu.pred.upper.bound, lower = mu.pred.lower.bound)
  temp.plot.obj.1 <- xyplot(mean ~ x, data = temp.df.2, type = "l", col = "blue", lty = 2)
  temp.plot.obj.2 <- xyplot(lower ~ x, data = temp.df.2, type = "l", col = "springgreen4", lty = 2)
  temp.plot.obj.3 <- xyplot(upper ~ x, data = temp.df.2, type = "l", col = "springgreen4", lty = 2)
  
  whole.plot.object <- indiv.plot.1 + indiv.plot.2 + indiv.plot.3 + indiv.plot.4 + temp.plot.obj.1 + temp.plot.obj.2 + temp.plot.obj.3
  
  # pdf some things
  pdf(file = paste("Person", as.character(indiv.index),"obsDropped", as.character(obs.to.drop), ".pdf", sep = ""), width = 10, height = 10)
  print(whole.plot.object)
  dev.off()
  # 
  ## write the data frames to csv
  write.csv(x = temp.plot.df, file = paste("betaFit", as.character(indiv.index), "obsDropped", as.character(obs.to.drop), ".csv", sep = ""))
  write.csv(x = temp.df.2, file = paste("muPredFit", as.character(indiv.index), "obsDropped", as.character(obs.to.drop),".csv", sep = ""))
  
  coverage.vector <- with(temp.plot.df[(length(temp.plot.df$x) - obs.to.drop + 1):(length(temp.plot.df$x)), ], (lower.interval.temp < y) & (y < upper.interval.temp))
  write.csv(x = coverage.vector, file = paste("coverage", as.character(indiv.index),"obsDropped", as.character(obs.to.drop), ".csv", sep = ""))
  
  
  return(list(whole.plot.object, coverage.vector))
  
}

# coverageChecker <- function(fit.object) {
#   
# }
