#devtools::install_github("ZoeGerber/ModeleSIR")
#library(ModeleSIR)



#n0 le nombre de personnes dans la pop a t=0
seir <- function(t, dt,p,beta,gamma, alpha, mu, nu){
  s0 <- p #proportion
  e0 <- (1 -p)/2  #on a supposé que la proportion d'infectés et celle d'exposés étaient égales
  i0 <- (1-p)/2
  r0 <- 0
  N <- s0+e0+i0+r0
  tps  <- floor(t/dt)    #nombre d'itération  : floor(t/dt)
  resS <- rep(s0, tps)
  resE <- rep(e0, tps)
  resI <- rep(i0, tps)
  resR <- rep(r0, tps)
  resN <- rep(N, tps)
  j <- c(1:tps)
  for(i in 2:tps)
  {
    resS[i] <- resS[i-1] + (-beta * resI[i-1] * resS[i-1] - mu * resS[i-1] + nu * resN[i-1]) * dt
    resE[i] <- resE[i-1] + (beta * resI[i-1] * resS[i-1] - alpha * resE[i-1] - mu * resE[i-1] ) * dt
    resI[i] <- resI[i-1] + (alpha * resE[i-1] - gamma * resI[i-1] - mu * resI[i-1]) * dt
    resR[i] <- resR[i-1] + (gamma * resI[i-1] - mu * resR[i-1] ) * dt
    resN[i] <- resS[i-1] + resE[i-1] + resI[i-1] + resR[i-1]
  }
  df <- data.frame(j,resS,resE, resI,resR, resN)
  return(df)
}



tirageBeta <- function(moy,var){
  beta <- sort(abs(rnorm(1,moy,var)))
  return(beta)
}

tirageGamma <- function(m,sd){
  gamma <- sort(abs(rnorm(1,m,sd)))
  return(gamma)
}



picI <- function(beta,gamma, data){
  R0 <- beta/gamma
  pic <- max(data[,3])
  date <- data[which.max(data[,3]),1]
  return(list(PicI = pic, datePicI = date, R0pic = R0))
}
