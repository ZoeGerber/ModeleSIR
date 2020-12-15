#Initialiser les fractions des sous-populations saines, infectées et retirées
#Au début de l'épidémie on a p personnes saines, et 1-p personnes malades
#Il n'y a aucune personne guéries ni mortes au début de l'épidémie
#p est donc la proportion de personnes saines
initSir<- function(p){
  s0 <- p           #proportion
  i0 <- 1 - p
  r0 <- 0
  return(c(s0,i0,r0))
}


####################

#Initialiser les taux de transmission et de guérison
#Tirage d'une valeur aléatoire de beta et gamma
#par une loi uniforme entre 2 bornes fixées
tirageBeta <- function(min,max){
  beta <- runif(1,min,max)
  return(beta)
}

tirageGamma <- function(min,max){
  gamma <- runif(1,min,max)
  return(gamma)
}



####################


#Fonction sir() permet de tracer les courbes avec les paramètres fixés
#On connait les proportions initiales des sous populations
##nombre d'itération de la simulation  : floor(t/dt)
#mise des résulats dans un dataframe
sir <- function(t, dt,p,beta,gamma){
  s0 <- p #proportion
  i0 <- 1 - p
  r0 <- 0
  tps  <- floor(t/dt)
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


####################

#Renvoi du pic et de la date du pic
#Eécupère le dataframe de la fonction sir()
#Date du pic, proportion de personne infectées et calcul du R0
#en fonction des beta et gamma de la fonction sir()
picI <- function(beta,gamma, data){
  R0 <- beta/gamma
  pic <- max(data[,3])
  date <- data[which.max(data[,3]),1]
  return(data.frame(PicI = pic, datePicI = date, R0pic = R0))
}


####################

#regroupement des focntion du dessus
#permet de calculer les valeurs des proportions des sous-populations at t+1
#☺on peut tracer les courbes avec ggplot
#beta et gamma aléatoires
#donne la date et le max du pic et le R0
picIsimu <- function(t,dt,p,min,max){
  beta <- runif(1,min,max)
  gamma <- runif(1,min,max)
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



