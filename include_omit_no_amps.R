# How does dropping no-amps from qPCR replicates influence estimated copy number?

require("scales") # needed for visualization

# function for simulation
estimate_qPCR <- function(x){
  copies <- rpois(10, x)                                                        # ten replicates, Poisson distribution
  drop_zeros <- ifelse(copies == 0, NA, copies)                                 # no copies = no amplification
  mean_copies <- mean(copies)                                                   # retain no-amps as zero
  mean_drop_copies <- mean(drop_zeros, na.rm = T)                               # drop no-amps from average
  mean_drop_copies <- ifelse(is.na(mean_drop_copies) == T, 0, mean_drop_copies)
  print(c(mean_copies, mean_drop_copies))                                       # print both estimates
}

# run 1K simulations at each mean copy number level
onecopy <- matrix(, ncol = 10, nrow = 1000)
for(i in 1:1000){
  onecopy[i,] <- estimate_qPCR(1)
}

twocopy <- matrix(, ncol = 10, nrow = 1000)
for(i in 1:1000){
  twocopy[i,] <- estimate_qPCR(2)
}

threecopy <- matrix(, ncol = 10, nrow = 1000)
for(i in 1:1000){
  threecopy[i,] <- estimate_qPCR(3)
}

fourcopy <- matrix(, ncol = 10, nrow = 1000)
for(i in 1:1000){
  fourcopy[i,] <- estimate_qPCR(4)
}

# plot
plot(NULL, NULL,
     xlim = c(0.8,4.2), ylim = c(0,7),
     xlab = "True copies (mean/rxn)", ylab = "Estimated copies",
     axes = F)
axis(1, at = 1:4)
axis(2)

points(onecopy[,1] ~ c(jitter(rep(1, 1000)) - 0.1), pch = 19, col = alpha("black", 0.2))
points(onecopy[,2] ~ c(jitter(rep(1, 1000)) + 0.1), pch = 19, col = alpha("darkblue", 0.2))

points(twocopy[,1] ~ c(jitter(rep(1, 1000)) + 0.9), pch = 19, col = alpha("black", 0.2))
points(twocopy[,2] ~ c(jitter(rep(1, 1000)) + 1.1), pch = 19, col = alpha("darkblue", 0.2))

points(threecopy[,1] ~ c(jitter(rep(1, 1000)) + 1.9), pch = 19, col = alpha("black", 0.2))
points(threecopy[,2] ~ c(jitter(rep(1, 1000)) + 2.1), pch = 19, col = alpha("darkblue", 0.2))

points(fourcopy[,1] ~ c(jitter(rep(1, 1000)) + 2.9), pch = 19, col = alpha("black", 0.2))
points(fourcopy[,2] ~ c(jitter(rep(1, 1000)) + 3.1), pch = 19, col = alpha("darkblue", 0.2))

# lines showing true mean
segments(x0 = c(0.8 + 0:3),
         x1 = c(1.2 + 0:3),
         y0 = c(1:4),
         y1 = c(1:4),
         lwd = 3, col = "red")

legend(x = 1, y = 6, bty = "n", 
       pch = 19, col = c(alpha("black", 0.7),alpha("darkblue", 0.7)),
       legend = c("Include no-amps", "Drop no-amps"))


