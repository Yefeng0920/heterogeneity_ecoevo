#################################################
# custom functions used in main analysis
#################################################

# function to calculate heterogeneity
h.calc <- function(mod){
  # I2
  # sigma2_v = typical sampling error variance
  sigma2_v <- sum(1 / mod$vi) * (mod$k - 1) /
    (sum(1 / mod$vi)^2 - sum((1 / mod$vi)^2))
  # s^2_t = total variance
  I2_total <- 100 * (sum(mod$sigma2) / (sum(mod$sigma2) + sigma2_v))
  I2_each <- 100 * (mod$sigma2 / (sum(mod$sigma2) + sigma2_v))
  #names(I2_each) <- paste0("I2_", model$s.names)
  #names(I2_total) <- "I2_Total"
  I2s_Shinichi <- c(I2_total, I2_each)
  
# matrix version  
  W <- solve(mod$V)
  X <- model.matrix(mod)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2_total2 <- 100* (sum(mod$sigma2) / (sum(mod$sigma2) + (mod$k - mod$p) / sum(diag(P))))
  I2_each2 <- 100* (mod$sigma2 / (sum(mod$sigma2) + (mod$k - mod$p) / sum(diag(P))))
  #names(I2_each2) <- paste0("I2_", model$s.names)
  #names(I2_total2) <- "I2_Total2"
  I2s_Wolfgang <- c(I2_total2, I2_each2)
  
  
  # CVB
  CV_total <- (sqrt(sum(mod$sigma2)) / abs(mod$beta[1]))
  CV_each <- (sqrt(mod$sigma2) / abs(mod$beta[1]))

  #names(CVB_each) <- paste0("CVB_", mod$s.names)
  #names(CVB_total) <- "CVB_total"
  CVs <- c(CV_total, CV_each)
  
  # M1
  M1_total <- (sum(sqrt(mod$sigma2)) / (sum(sqrt(mod$sigma2)) + abs(mod$beta[1])))
  M1_each <- sqrt(mod$sigma2) / (sum(sqrt(mod$sigma2)) + abs(mod$beta[1]))
  #names(M1_each) <- paste0("CVB_", mod$s.names)
  #names(M1_total) <- "M1_total"
  Ms <- c(M1_total, M1_each)

  hs <- data.frame(I2s_Shinichi,CVs,Ms)
  rownames(hs) <- c("Total", mod$s.names)
  return(hs)

}


# function to estimate typical sampling error variance
sigma2_v <- function(mod){
  sigma2_v <- sum(1 / mod$vi) * (mod$k - 1) /
    (sum(1 / mod$vi)^2 - sum((1 / mod$vi)^2))
  return(sigma2_v)
}
