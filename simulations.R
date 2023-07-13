# significance level simulations ----

source("simulation_function.R")
source("model.R")

set.seed(123)

s <- 4
frac.seq <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1:10)
n.seq <- round(2^(s+2)*frac.seq)
N <- 500

p <- model_Y_to_XZ(s, 0.5, 1)$probability
p <- projection_CI(p)

results <- NULL

for(n.ind in 1:length(n.seq)){
  
  n <- n.seq[n.ind]
  
  for(i in 1:N){
    
    p.sample <- sample_from_distribution(p, n)
    
    result.tmp <- simulations_for_sample(p.sample, p, resampling.scenario.seq = c("BCI", "CR", "BX", "CP"), B=50, n=n) 
    result.tmp <- as.data.frame.table(result.tmp)
    colnames(result.tmp) <- c("resampling", "test", "pvalue")
    result.tmp <- cbind(result.tmp, n=n, s=s, frac=frac.seq[n.ind])
    
    results <- rbind(result.tmp, results)
  }
}

write.csv(results, "significance.csv")

# power ----

source("simulation_function.R")
source("model.R")

set.seed(123)

frac.seq <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1:5)
s <- 4
n.seq <- round(2^(s+2)*frac.seq)
N <- 500

p <- model_Y_to_XZ(s, 0.5, 1)$probability

results <- NULL

for(n.ind in 1:length(n.seq)){
  
  n <- n.seq[n.ind]

  for(i in 1:N){
    
    p.sample <- sample_from_distribution(p, n)
    
    result.tmp <- simulations_for_sample(p.sample, p, resampling.scenario.seq = c("BCI", "CR", "BX", "CP"), B=50, n=n)
    result.tmp <- as.data.frame.table(result.tmp)
    colnames(result.tmp) <- c("resampling", "test", "pvalue")
    result.tmp <- cbind(result.tmp, n=n, s=s, frac=frac.seq[n.ind])
    
    results <- rbind(result.tmp, results)
  }
}

write.csv(results, "power.csv")