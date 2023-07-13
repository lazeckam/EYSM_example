
# sampling ----

sample_from_distribution <- function(p, n){
  
  n.sample <- rmultinom(n=1, size=n, prob=p)[,1]
  p.sample <- n.sample/n
  
  return(p.sample)
}

# resampling scenarios ----

resampling_BCI <- function(p.sample.CI, n){
  
  n.resample <- rmultinom(n=1, size=n, prob=p.sample.CI)[,1]
  p.resample <- n.resample/n
  
  return(p.resample)
}

resampling_CR <- function(p.sample, p_x_z, n){
  
  cells <- length(p.sample)
  n.resample <- rep(0, cells)
  n.sample <- round(n*p.sample)
  
  for(i in 1:(cells/2)){
    
    # n(y, z) for y and z fixed
    nyz <- n.sample[2*i - 1] + n.sample[2*i] 
    
    if(nyz != 0){
      # sampling n(x=0, y, z)
      n.resample[2*i-1] <- rbinom(1, size=nyz, prob=c(p_x_z[2*((i - 1) %/% 2) + 1], p_x_z[2*((i - 1) %/% 2) + 2]))
      # computing n(x=1, y, z) = n(y, z) - n(x=0, y, z)
      n.resample[2*i] <- nyz - n.resample[2*i-1]
    }else{
      n.resample[2*i-1] <- 0
      n.resample[2*i] <- 0
    }
    
    
  }
  
  p.resample <- n.resample/n
  
  return(p.resample)
}

resampling_BX <- function(p.sample, p.sample_x_z, n){
  
  p.sample <- resampling_CR(p.sample, p.sample_x_z, n)
  
  return(p.sample)
}

resampling_CP <- function(p.sample, n){
  
  cells <- length(p.sample)
  n.resample <- rep(0, cells)
  n.sample <- n*p.sample
  
  for(i in 1:(cells/4)){
    
    # n(y, z) for y=0,1 and z fixed
    ny0z <- n.sample[4*i - 3] + n.sample[4*i-2] 
    ny1z <- n.sample[4*i - 1] + n.sample[4*i] 
    # n(x, z) for x=0 and z fixed
    nx0z <- n.sample[4*i - 3] + n.sample[4*i-1] 
    
    # sampling n(x=0, y=0, z)
    n.resample[4*i-3] <- rhyper(1, m=ny0z, n=ny1z, k=nx0z)
    # computing n(x=1, y=0, z) = n(y=0, z) - n(x=0, y=0, z)
    n.resample[4*i-2] <- ny0z - n.resample[4*i-3]
    # computing n(x=0, y=1, z) = n(x=0, z) - n(x=0, y=0, z)
    n.resample[4*i-1] <- nx0z - n.resample[4*i-3]
    # computing n(x=1, y=1, z) = n(y=1, z) - n(x=0, y=1, z)
    n.resample[4*i] <- ny1z - n.resample[4*i-1]
    
  }
  
  p.resample <- n.resample/n
  
  return(p.resample)
}

# computing CMI on resampled samples ----

CMI_resampling_BCI <- function(p.sample.CI, B, n,
                               p.sample=NULL, p_x_z=NULL, p.sample_x_z=NULL,
                               entropy_z=NULL, entropy_xz=NULL, entropy_yz=NULL){
  
  cmi.resample <- rep(NA, B)
  
  for(i in 1:B){
    
    p.resample <- resampling_BCI(p.sample.CI, n)
    
    cmi.resample[i] <- conditional_mutual_information(p.resample)
  }
  
  return(cmi.resample)
}

CMI_resampling_CR <- function(p.sample, p_x_z, entropy_z, entropy_yz, B, n,
                              p.sample.CI=NULL, entropy_xz=NULL, p.sample_x_z=NULL){
  
  cmi.resample <- rep(NA, B)
  
  for(i in 1:B){
    
    p.resample <- resampling_CR(p.sample, p_x_z, n)
    
    p_xz <- compute_p_xz(p.resample)
    entropy_xz <- entropy(p_xz)
    entropy_xyz <- entropy(p.resample)
    
    cmi.resample[i] <- -entropy_xyz + entropy_xz + entropy_yz - entropy_z
    
  }
  
  return(cmi.resample)
}

CMI_resampling_BX <- function(p.sample, p.sample_x_z, entropy_z, entropy_yz, B, n,
                              p.sample.CI=NULL, entropy_xz=NULL, p_x_z=NULL){
  
  cmi.resample <- rep(NA, B)
  
  for(i in 1:B){
    
    p.resample <- resampling_BX(p.sample, p.sample_x_z, n)
    
    p_xz <- compute_p_xz(p.resample)
    entropy_xz <- entropy(p_xz)
    entropy_xyz <- entropy(p.resample)
    
    cmi.resample[i] <- -entropy_xyz + entropy_xz + entropy_yz - entropy_z
  }
  
  return(cmi.resample)
  
}

CMI_resampling_CP <- function(p.sample, entropy_z, entropy_yz, entropy_xz, B, n,
                              p.sample.CI=NULL, p.sample_x_z=NULL, p_x_z=NULL){
  
  cmi.resample <- rep(NA, B)
  
  for(i in 1:B){
    
    p.resample <- resampling_CP(p.sample, n)
    
    entropy_xyz <- entropy(p.resample)
    
    cmi.resample[i] <- -entropy_xyz + entropy_xz + entropy_yz - entropy_z
  }
  
  return(cmi.resample)
}
