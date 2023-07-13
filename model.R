model_Y_to_XZ <- function(s, gamma, sigma){
  
  joint.distribution <- expand.grid(split(rep(c(0,1), s+2), rep(1:(s+2), each=2)))
  colnames(joint.distribution) <- c("X", "Y", paste0("Z", 1:s))
  
  probability <- rep(0, nrow(joint.distribution))
  
  Sigma <- diag(rep(sigma, s))
  
  probY <- 0.5
  gamma_vec <- gamma^(1:(s))
  
  
  for(i in 1:nrow(joint.distribution)){
    
    value <- joint.distribution[i,]
    
    
    if(value$Y == 0){
      
      lb <- rep(-Inf, s)
      ub <- rep(Inf, s)
      
      ub[value[,3:(s+2)] == 0] <- (gamma_vec/2)[value[,3:(s+2)] == 0]
      lb[value[,3:(s+2)] == 1] <- (gamma_vec/2)[value[,3:(s+2)] == 1]
      
      probZ <- mvtnorm::pmvnorm(lower=lb, upper=ub, mean=rep(0, s), sigma=Sigma)
      
      
    }else{
      
      lb <- rep(-Inf, s)
      ub <- rep(Inf, s)
      
      ub[value[,3:(s+2)] == 0] <- (gamma_vec/2)[value[,3:(s+2)] == 0]
      lb[value[,3:(s+2)] == 1] <- (gamma_vec/2)[value[,3:(s+2)] == 1]
      
      probZ <- mvtnorm::pmvnorm(lower=lb, upper=ub, mean=gamma_vec, sigma=Sigma)
      
    }
    
    p0.5 <- pnorm(0.5, 0, sigma)
    probX <- ifelse(value[2] == value[1], p0.5, 1 - p0.5)
    
    probability[i] <- probX*probY*probZ
    
  }
  
  joint.distribution <- data.frame(joint.distribution, probability=probability)
  
  return(joint.distribution)
}

