source ('bsplines.R')
library(dplyr)
# beta_old matrix
#y vector
#phi matrix

posterior <- function (y 
                       ,beta_old 
                       ,phimat
                       ,sigma = 1){
    num_group <- ncol(beta_old)
    mean_diff <- rep(0,num_group)
    n <-length(y)
    ratio <- rep(1,num_group)
    
    for (j in 1: num_group){
        mean_diff[j] <- t(y- phimat %*% beta_old[,j]) %*% (y- phimat %*% beta_old[,j])        
    }
    
    baseline <- min(mean_diff)
    
    for (j in 1:num_group ){
      ratio[j] <- exp(-mean_diff[j] + baseline)
    #  cat ('ratio is ', ratio[j], '\n')
    }
    
    for (i in 1:length(ratio)){
      if (ratio[i] < 1*exp(-10)){
          ratio[i] = 0
      }
    }
    
    prob <- ratio/sum(ratio)
  #  cat ('Probability is ', prob,"\n")
    return (prob)
}



update_beta <- function (data
                         ,beta_old
                         ,bb
                         ,GP){
  
  beta_new <- beta_old 
  post <- matrix(0,nrow(data),ncol(beta_old))
  
  for (j in 1:ncol(beta_old)){

    left <- matrix(0,nrow= nrow(beta_old), ncol = nrow(beta_old))
    right <- rep(0,nrow(beta_old))
    
    for (id in 1:n_distinct(data$ID)){
      t <- data %>% filter( ID == id ) %>% arrange(years_seen) %>% select (years_seen)
      t <- t[['years_seen']]
      
#       y <- data %>% filter (ID == id ) %>% arrange(years_seen) %>% select (pfvc_noise)
#       y <- y[['pfvc_noise']]

      y <- data %>% filter (ID == id ) %>% arrange(years_seen) %>% select (pfvc_gp)
      y <- y[['pfvc_gp']]
  
      phimat <- design(t,bb)
      post[id,] <- posterior(y,beta_old,phimat)

      if (!GP){
          left <- left + t(phimat) %*% phimat * post[id,j]
          right <- right + t(phimat) %*% y * post[id,j]
      }else{
        K = FUN(t)+diag(length(t))
        W = solve(K, diag(nrow(K)))
        left <- left + t(phimat) %*% W %*% phimat *post [id, j]
      }
    }
    
    eigen_value <- eigen(left)$values
    
    if (min(eigen_value) < 10^(-5)){
          left <- left +diag(nrow(left)) *10^(-5)      
    }
    
    beta_new[,j] = solve(left,right)
  }

  return (list(beta = beta_new,posterior = post))
  
}





EM <- function(data,beta,bb, ...){ 
  beta_old <-  matrix(0,nrow = nrow(beta),ncol = ncol(beta))
  beta_new <- beta
  count <- 1
  GP <- FALSE
  GP <- list(...)[[1]]
  FUN <- list(...)[[2]]
  
  while (sqrt(sum((beta_new-beta_old)**2)) > 0.5){
    cat ("Iteration number ", count, "\n")
    beta_old <- beta_new
    update <- update_beta(data,beta_old,bb, GP =FALSE)
    beta_new <- update$beta
    posterior <- update$posterior
    loglik = 0
    for (i in 1:n_distinct(data$ID)){
      for (j in 1:ncol(beta)){
     # pf <- data %>% filter( ID == i ) %>% select (pfvc_noise)
     # pfvc <- pf[['pfvc_noise']]
        
      pf <- data %>% filter( ID == i ) %>% select (pfvc_gp)
      pfvc <- pf[['pfvc_gp']]
        
      years <- data %>% filter (ID == i ) %>% select (years_seen)
      years_seen <- years[['years_seen']]
      mu <- design (years_seen, bb) %*% beta_new[,j]
      sigma <- diag(length(years_seen))
      
      if (!GP)  {
          loglik <- -t(pfvc-mu)%*% (pfvc-mu)*posterior[i,j] + loglik             
      }else{
        loglik <- -t(pfvc-mu)%*% (pfvc-mu)*posterior[i,j] - 
                    1/2 *log(det(FUN(years_seen) + 
                                  diag(length(years_seen)))) + loglik  
      }
      
      }
     
    }
    
    cat ("loglikelihood is ", loglik, "\n")
    count <- count +1
    cat ("updated beta is: " ,"\n")
    print (beta_new)
    
    if (count >30)
      return (beta_new)
  }
  
  return (beta_new)

}
