#devtools::install_github("ZoeGerber/ModeleSIR")
#library(ModeleSIR)

#Initialiser les fractions des sous-populations saines, infectées et guéries
#Au début de l'épidémie on a p personnes saines, et 1-p personnes malades
#Il n'y a aucune personne guéries au début de l'épidémie
#p est donc la proportion de personnes saines
initSir<- function(p){
  s0 <- p           #proportion
  i0 <- 1 - p
  r0 <- 0
  return(c(s0,i0,r0))
}

#Initialiser les taux de transmission et de guérison

tirageBeta <- function(moy,var){
  beta <- abs(rnorm(1,moy,var))
  return(beta)
}

tirageGamma <- function(m,sd){
  gamma <- abs(rnorm(1,m,sd))
  return(gamma)
}


#######################################################


#Au début de l'épidémie on a p personnes saines, et 1-p personnes malades
#Il n'y a aucune personne guéries au début de l'épidémie
#p est donc la proportion de personnes saines
#On doit également définir le temps de la simulation(en jours)


sir <- function(t, dt,p,beta,gamma){
  s0 <- p #proportion
  i0 <- 1 - p
  r0 <- 0
  tps  <- floor(t/dt)    #nombre d'itération  : floor(t/dt)
  resS <- rep(s0, tps)
  resI <- rep(i0, tps)
  resR <- rep(r0, tps)
  j <- c(1:tps)
  for(i in 2:tps)
  {
    resS[i] <- resS[i-1] + (-beta * resI[i-1] * resS[i-1]) * dt
    resI[i] <- resI[i-1] + (beta * resS[i-1] * resI[i-1] - gamma * resI[i-1]) * dt
    resR[i] <- resR[i-1] + (gamma * resI[i-1]) * dt
  }
  df <- data.frame(j,resS,resI,resR)
  return(df)
}

#####Renvoi du pic et de la date du pic

picI <- function(beta,gamma, data){
  R0 <- beta/gamma
  pic <- max(data[,3])
  date <- data[which.max(data[,3]),1]
  return(data.frame(PicI = pic, datePicI = date, R0pic = R0))
}


####################


picI2 <- function(t,dt,p,m1,m2,sd1,sd2){
  beta <- abs(rnorm(1,m1,sd1))
  gamma <- abs(rnorm(1,m2,sd2))
  R0 <- beta/gamma
  s0 <- p #proportion
  i0 <- 1 - p
  r0 <- 0
  tps  <- floor(t/dt)    #nombre d'itération  : floor(t/dt)
  resS <- rep(s0, tps)
  resI <- rep(i0, tps)
  resR <- rep(r0, tps)
  j <- c(1:tps)
  for(i in 2:tps)
  {
    resS[i] <- resS[i-1] + (-beta * resI[i-1] * resS[i-1]) * dt
    resI[i] <- resI[i-1] + (beta * resS[i-1] * resI[i-1] - gamma * resI[i-1]) * dt
    resR[i] <- resR[i-1] + (gamma * resI[i-1]) * dt
  }
  df <- data.frame(j,resS,resI,resR)
  pic <- max(df[,3])
  date <- df[which.max(df[,3]),1]
  return(data.frame(PicI = pic, datePicI = date, R0pic = R0))
}



