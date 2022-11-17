# code to simulate and analyze hierarchical occupancy models
# eDNA case study

Sys.setenv(JAGS_HOME='C:/Program Files')
library(jagsUI)


# specify study dimensions
n.spp <- 4 # number of species
n.lakes <- 12 # number of sampled lakes (occupancy fixed at this level)
n.sites <- 10 # unique sample sites per lake
n.times <- 5 # seasonal sampling periods
n.reps <- 3 # number of subsamples from each water sample

################################################
# specify initial occupancy status
 # imagining 1 spp is ubiquitous (e.g. common carp?)
 # other spp variably present
 # for purposes of this simulation, I'm treating as data
 # but could (should??) also be estimated by model
################################################

occupancy <- matrix(c(0, 0, 0, 1,
                0, 0, 0, 1,
                0, 0, 1, 1,
                0, 0, 1, 1,
                0, 1, 1, 1,
                0, 1, 1, 1,
                1, 1, 1, 1,
                1, 1, 1, 1,
                1, 1, 1, 1,
                1, 1, 0, 1,
                1, 1, 0, 1,
                1, 0, 0, 1), nrow = n.lakes, ncol = n.spp, byrow = TRUE)

####################################################
# specify parameters affecting detection probability
# models will use logit-link
####################################################

# postulate species-specific intercepts
 # e.g. some species release more detectable DNA
set.seed(42)
spp.mu <- 0 # on logit scale, so mean 0.5
spp.sd <- 1
spp.effect <- rnorm(n.spp, spp.mu, spp.sd)
print(round(plogis(spp.effect), 2)) # print spp-specific means, real scale

# model lake specific variation
 # e.g. water chemistry, depth profile that favors/disfavors detection
set.seed(8675309)
lake.mu <- 0
lake.sd <- 0.25 # less effect than spp or spp_lake?
lake.effect <- rnorm(n.lakes, lake.mu, lake.sd)

# model lake-by-species interaction
 # e.g. variation in abundance by lake that affects detection
spp_lake.effect <- matrix(NA, nrow = n.lakes, ncol = n.spp)
spp_lake.mu <- 0
spp_lake.sd <- 0.5
set.seed(123)
for (i in 1:n.spp){
  for (j in 1:n.lakes){
    spp_lake.effect[j,i] <- spp.effect[i] + 
      lake.effect[j] + 
      rnorm(1, spp_lake.mu, spp_lake.sd)
  }
}

plot(1:n.lakes, plogis(spp_lake.effect[,1]), pch = 19, type = 'b', 
     ylim = c(0, 1), ylab = "Detection probability", xlab = "Lake number")
for (i in 2:n.spp){
  points(1:n.lakes, plogis(spp_lake.effect[,i]), pch = 19, type = 'b', col = i)
}

# simulate covariates for sampling location effect
 # imagine 1 strong covariate (depth) 
 # 1 weak covariate (dist_shore)
 # plus random effect for unmeasured covariates
 # model covariate effects as if standardized
depth <- dist_shore <- matrix(NA, nrow = n.lakes, ncol = n.sites)
b_depth <- 0.5 # hypothesized strong coefficient
b_dist <- 0.25 # hypothesized weak coefficient
set.seed(456)
for (i in 1:n.lakes){
  for (j in 1:n.sites){
    depth[i,j] <- rnorm(1, 0, 1)
    dist_shore[i,j] <- rnorm(1, 0, 1)
    # could have lake-specific means and SD's but ignoring for now
  }
}

# simulate spp-specific patterns for season effects
 # good ability to add covariates here, but need thoughtful cov models
 # for now I'm just going to imagine spp have different seasonal patterns
 # including no pattern for spp 1
season.mu <- matrix(c(0, -2,  -1,  1,
                      0, -2,   0,  0.5,
                      0,  0,  1.5,  0,
                      0, 1.5, 0.5, -0.5,
                      0, 2.5, -1, -1), 
                    nrow = n.times, ncol = n.spp, byrow = TRUE)

# patterns (these effects will be added to other covariates)
plot(1:n.times, plogis(season.mu[,1]), pch = 19, type = 'b', 
     ylim = c(0,1), ylab = "Detection probability", xlab = "Season")
for (i in 2:n.spp){
  points(1:n.times, plogis(season.mu[,i]), pch = 19, type = 'b', col = i)
}

# simulate long-form data
set.seed(42)
longform.data <- matrix(NA, nrow = n.spp * n.lakes * n.times * n.sites,
                        ncol = 9)
for (i in 1:n.spp){
  for (j in 1:n.lakes){
    for (t in 1:n.times){
      # these terms will overwrite each iteration, create 3-dimensional arrays to save as data
      time.sd <- dnorm(1, 0, 0.3) # large spp-lake-time effect
      for (s in 1:n.sites){
        site.sd <- dnorm(1, 0, 0.1) # small spp-lake-season-site effect (resid)
        # this overwrites each iteration but saved as data
        p <- plogis(spp_lake.effect[j,i] + # includes spp, lake, and spp-by-lake effects (graph 1)
                    b_depth*depth[j,s] + b_dist*dist_shore[j,s] + # site covariate effects
                    season.mu[t,i] + # these are spp-by-season effects (graph 2)
                    time.sd + site.sd) # true detection probability
        longform.data[(i-1)*n.lakes*n.times*n.sites +
                       (j-1)*n.times*n.sites +
                       (t-1)*n.sites + s, 1] <- i # spp.number in col 1
        longform.data[(i-1)*n.lakes*n.times*n.sites +
                       (j-1)*n.times*n.sites +
                       (t-1)*n.sites + s, 2] <- j # lake.number in col 2
        longform.data[(i-1)*n.lakes*n.sites*n.times +
                       (j-1)*n.times*n.sites +
                       (t-1)*n.sites + s, 3] <- t # sampling period in col 3
        longform.data[(i-1)*n.lakes*n.sites*n.times +
                       (j-1)*n.times*n.sites +
                       (t-1)*n.sites + s, 4] <- s # site number in col 4
        longform.data[(i-1)*n.lakes*n.sites*n.times +
                      (j-1)*n.times*n.sites +
                      (t-1)*n.sites + s, 5] <- occupancy[j,i] # true occupancy status
        longform.data[(i-1)*n.lakes*n.sites*n.times +
                      (j-1)*n.times*n.sites +
                      (t-1)*n.sites + s, 6] <- p # true detection probability
        longform.data[(i-1)*n.lakes*n.sites*n.times +
                      (j-1)*n.times*n.sites +
                      (t-1)*n.sites + s, 7:9] <- rbinom(3,1,occupancy[j,i]*p) # detected or not each of 3 subsamples
      }
    }
  }
}
tail(longform.data)


# Bundle data for JAGS
jags.data <- list(y = as.matrix(longform.data[,7:9]), 
                  n.obs = nrow(longform.data),
                  n.spp = n.spp, n.lakes = n.lakes,
                  n.sites = n.sites, n.times = n.times, 
                  n.reps = n.reps, occupancy = occupancy,
                  spp = longform.data[,1], lake = longform.data[,2],
                  time = longform.data[,3], site = longform.data[,4],
                  depth = depth, dist_shore = dist_shore)

########################################################
# Specify multi-level occupancy model in BUGS language
# want species specific estimates of detection at all combinations
########################################################

sink("eDNA_occ.txt")
cat("
model {
# priors
  # regression coefficients
  b.depth ~ dnorm(0, 0.33) # vague prior for regression slope of depth
  b.dist ~ dnorm(0, 0.33) # regression slope for dist_shore

  # SD's for simple random effects (detection: p)
  p.spp.sd ~ dunif(0, 1) # SD for spp effects
  p.spp.tau <- pow(p.spp.sd, -2) # convert to precision
  p.lake.sd ~ dunif(0, 1) # SD for lake effects
  p.lake.tau <- pow(p.lake.sd, -2) 
  p.time.sd ~ dunif(0, 1) # SD for time effects
  p.time.tau <- pow(p.time.sd, -2) 
  # sites are uniquely nested within lakes, no context as simple RE

  # SDs and precision for complex random effects (all presumed spp-specific)
  p.spp.lake.sd ~ dunif(0, 1) # spp-lake effect
  p.spp.lake.tau <- pow(p.spp.lake.sd, -2) 
  p.spp.time.sd ~ dunif(0, 1) # spp-time effect
  p.spp.time.tau <- pow(p.spp.time.sd, -2) 
  p.spp.lake.site.sd ~ dunif(0, 1) # SD for site effects within lakes within spp
  p.spp.lake.site.tau <- pow(p.spp.lake.site.sd, -2) 
  p.spp.lake.time.sd ~ dunif(0, 1) # SD for time effects within lakes & spp
  p.spp.lake.time.tau <- pow(p.spp.lake.time.sd, -2) 

# species specific priors 
  for (i in 1:n.spp){
    psi.spp[i] ~ dunif(0, 1) # species-specific occupancy
    p.spp[i] ~ dnorm(0, p.spp.tau) # spp mean detection (logit scale)
    # species-lake effect
    for (j in 1:n.lakes){
      p.spp.lake[i,j] ~ dnorm(0, p.spp.lake.tau) # lake-by-spp adjustments
      for (t in 1:n.times){
        p.spp.lake.time[i,j,t] ~ dnorm(0, p.spp.lake.time.tau) # spp-lake-time
      }
      for (s in 1:n.sites){
        p.spp.lake.site[i,j,s] ~ dnorm(0, p.spp.lake.site.tau) # spp-lake-time
      }
    }
    # species-time effect
    for (t in 1:n.times){
      p.spp.time[i,t] ~ dnorm(0, p.spp.time.tau) # time-by-spp adjustment
    }
  }

# Ecological model for true detection
  for (i in 1:n.spp) {
    for (j in 1:n.lakes){
      for (t in 1:n.times){
        for (s in 1:n.sites){
          logit(p[i,j,t,s]) <- p.spp[i] + p.spp.lake[i,j] + p.spp.time[i,t] + 
            p.spp.lake.time[i,j,t] + p.spp.lake.site[i,j,s] +
            b.depth * depth[j,s] + b.dist * dist_shore[j,s]
        }
      } # end time loop (t)
    } # end spp loop (j)
  } # end site loop (i)

# Observation model for replicated detection/nondetection observations
    for (n in 1:n.obs) {
      z[n] <- occupancy[lake[n],spp[n]]
      for (k in 1:n.reps){
        y[n,k] ~ dbern(z[n]*p[spp[n],lake[n],time[n],site[n]]) 
      } # close k
    } # close n
} # end jags model
",fill = TRUE)  # back to R code
sink()

# Parameters monitored
# can monitor p, but there are 2400 of them
params <- c("b.depth", "b.dist", "p.spp", "p.lake.sd", "p.time.sd", 
            "p.spp.lake.sd", "p.spp.time.sd",
            "p.spp.lake.time.sd", "p.spp.lake.site.sd",
            "p.spp.lake", "p.spp.time")

# MCMC settings
na <- 1000; ni <- 5000; nt <- 1; nb <- 1000; nc <- 3 # short for test run
#na <- 1000; ni <- 60000; nt <- 10; nb <- 10000; nc <- 3

# Call JAGS from R (~44 min at 60,000 ni)
eDNA1.out <- jags(jags.data, inits = NULL, params, "eDNA_occ.txt", 
                  n.adapt = na, n.chains = nc, n.thin = nt, 
                  n.iter = ni, n.burnin = nb, parallel = TRUE)
print(eDNA1.out, dig = 2)
max.rhat <- function(x){max(max_rhat <- c(sapply(x$Rhat, max, na.rm=TRUE)))} 
max.rhat(eDNA1.out) # 

par(mfrow = c(1,2))
plot(density(eDNA1.out$sims.list$b.depth), xlim = c(0,1), 
     xlab = "Parameter estimate", main = "Depth effect")
abline(v = b_depth, col = "red")

plot(density(eDNA1.out$sims.list$b.dist), xlim = c(0,1), 
     xlab = "Parameter estimate", main = "Dist_shore effect")
abline(v = b_dist, col = "red")

par(mfrow = c(1,2))
plot(density(eDNA1.out$sims.list$p.spp.lake.sd), xlim = c(0,1), 
     xlab = "Parameter estimate", main = "Species-lake effect")
abline(v = spp_lake.sd, col = "red")

plot(density(eDNA1.out$sims.list$p.lake.sd), xlim = c(0,1), 
     xlab = "Parameter estimate", main = "Dist_shore effect")
abline(v = lake.sd, col = "red")

par(mfrow = c(2,1))
plot(1:n.lakes, plogis(spp_lake.effect[,1]), pch = 19, type = 'b', 
     ylim = c(0, 1), ylab = "Detection probability", xlab = "", main = "Simulated means")
for (i in 2:n.spp){
  points(1:n.lakes, plogis(spp_lake.effect[,i]), pch = 19, type = 'b', col = i)
}
plot(1:12, plogis(eDNA1.out$mean$p.spp[1] + eDNA1.out$mean$p.spp.lake[1,]), 
     pch = 19, type = 'b', ylim = c(0,1), ylab = "Detection probability", 
     xlab = "Lake number", main = "Model estimates")
for (i in 2:n.spp){
  points(1:12, plogis(eDNA1.out$mean$p.spp[i] + eDNA1.out$mean$p.spp.lake[i,]), 
       pch = 19, type = 'b', col = i) 
}

