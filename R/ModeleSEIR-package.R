#Initialiser les fractions des sous-populations saines, infectées et retirées
#Au début de l'épidémie on a p personnes saines, et 1-p personnes malades
#Il n'y a aucune personne guéries ni mortes au début de l'épidémie
#p est donc la proportion de personnes saines
#on ajoute le compartiment E(t) qui sont les exposés
initSeir<- function(p){
  s0 <- p           #proportion
  e0 <- (1 -p)/2  #on a supposé que la proportion d'infectés et celle d'exposés étaient égales
  i0 <- (1-p)/2
  r0 <- 0
  return(c(s0,e0,i0,r0))
}

####################

#Initialiser les taux de transmission, de guérison et d'incubation
#Tirage d'une valeur aléatoire de beta, gamma et alpha
#par une loi uniforme entre 2 bornes fixées
tirageBeta <- function(min,max){
  beta <- runif(1,min,max)
  return(beta)
}

tirageGamma <- function(max){
  gamma <- runif(1,min,max)
  return(gamma)
}

tirageAlpha <- function(min,max){
  alpha <- runif(1,min,max)
  return(alpha)
}

####################

#N le nombre de personnes dans la pop a t=0
#la pop n'est plus constante
#permet de calculer les valeurs des sous pop
#mise dans un dataframe des valeurs
#permet l'utilisation de ggplot
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


####################

#pic calcul avec le df de la focntion seir()
picISeir <- function(beta,gamma, alpha, data){
  R0 <- beta/gamma
  pic <- max(data[,4])
  date <- data[which.max(data[,4]),1]
  return(list(PicI = pic, datePicI = date, R0pic = R0))
}


####################

#simulation globale
#mu et nu sont fixés
#nombre d'itération  : floor(t/dt)
#p toujours une proportion
picIsimuSeir <- function(t,dt,p,min,max){
  mu <- 0.001
  nu <- 0.009
  beta <- runif(1,min,max)
  gamma <- runif(1,min,max)
  alpha <- runif(1,min,max)
  R0 <- beta/gamma
  s0 <- p
  e0 <- (1 -p)/2
  i0 <- (1 -p)/2
  r0 <- 0
  N <- s0+e0+i0+r0
  tps  <- floor(t/dt)
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
  pic <- max(df[,4])
  date <- df[which.max(df[,4]),1]
  return(data.frame(PicI = pic, datePicI = date, R0pic = R0))
}


