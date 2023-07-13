source("resampling.R")
library("Rcpp")
sourceCpp("cpp_functions.cpp")

simulations_for_sample <- function(p.sample, p, 
                                   resampling.scenario.seq, 
                                   B, n){
  # p.sample and p:
  #    vectors of probabilities p(x,y,z) in the following order: (0,0,z1), (1,0,z1), (0,1,z1), (1,1,z1); (0,0,z2), ... and so on
  #    only for binary X, Y and Z_i for i=1,2,...,s
  #    p.sample = ^p, p=p
  
  # compute proabilities and entropies for sample and joint distribution, that will be used later
  p.sample_xz <- compute_p_xz(p.sample)
  p.sample_yz <- compute_p_yz(p.sample)
  p.sample_z <- compute_p_z(p.sample_xz)
  p.sample.CI <- projection_CI(p.sample)
  
  p_x_z <- compute_p_x_z(p)
  p.sample_x_z <- compute_p_x_z(p.sample)
  
  entropy.sample_z <- entropy(p.sample_z)
  entropy.sample_yz <- entropy(p.sample_yz)
  entropy.sample_xz <- entropy(p.sample_xz)
  entropy.sample_xyz <- entropy(p.sample)
  cmi.sample <- -entropy.sample_xyz + entropy.sample_xz + entropy.sample_yz - entropy.sample_z
  
  s <- round(log(length(p), base=2)) - 2
  
  results <- matrix(NA, length(resampling.scenario.seq), 3)
  colnames(results) <- c("asymptotic", "exact", "df_estimation")
  rownames(results) <- resampling.scenario.seq
  
  for(resampling.scenario.ind in 1:length(resampling.scenario.seq)){
    resampling.scenario <- resampling.scenario.seq[resampling.scenario.ind]
    
    cmi.results <- matrix(NA, B)
    cmi.results <- do.call(paste0("CMI_resampling_", resampling.scenario),
                           list(p.sample=p.sample, p.sample.CI=p.sample.CI, p.sample_x_z=p.sample_x_z, p_x_z=p_x_z,
                                entropy_z=entropy.sample_z, entropy_yz=entropy.sample_yz, entropy_xz=entropy.sample_xz, 
                                B=B, n=n))

    
    test_asymptotic <- 1 - pchisq(2*n*cmi.sample, df=1*1*2^s)
    test_exact <-  (1 + sum(cmi.results >= cmi.sample))/(B + 1)
    test_df_estimation <- ifelse(mean(cmi.results) <= 0, 
                                 ifelse(cmi.sample <= 0, 1, 0),
                                 1 - pchisq(2*n*cmi.sample, df=2*n*mean(cmi.results)))
    
    results[resampling.scenario.ind,] <- c(test_asymptotic, test_exact, test_df_estimation)
  }
  
  return(results)
}
