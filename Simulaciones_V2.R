
#install.packages('vars')
#install.packages("biwavelet")
#install.packages('RTransferEntropy')
#install.packages('mvtnorm')
# tempdir()
# # [1] "C:\Users\XYZ~1\AppData\Local\Temp\Rtmp86bEoJ\Rtxt32dcef24de2"
# dir.create(tempdir())
# Shannon entropy R´enyi
library(biwavelet)
library(RTransferEntropy)
library(vars)
library(mvtnorm)
library(lmtest)
library(dplyr)
library(tidyr)
library(data.table)
############### 1 ###########

Sistema_1 <- function(n, semilla = 12345){
  x1 = c(0)
  x2 = c(0)
  x3 = c(0)
  ind = 0
  media <- c(0,0,0)
  sigma<-diag(3)
  
  
  RuidoBlanco=rmvnorm(1,mean = media,sigma = sigma)
  set.seed(semilla)
  
  while (ind==0) {
    for( i in seq(1,(n-1+500)) ){
      RuidoBlanco=rmvnorm(1,mean = media,sigma = sigma)
      
      x1[i+1] = 3.4*x1[i]*(1-x1[i])^2*exp(-(x1[i]^2))+0.4*RuidoBlanco[1]
      
      x2[i+1] = 3.4*x2[i]*(1-x2[i])^2*exp(-(x2[i]^2))+ 0.5*x1[i]*x2[i] +0.4*RuidoBlanco[2]
      
      x3[i+1] = 3.4*x3[i]*(1-x3[i])^2*exp(-(x3[i]^2))+ 0.3*x2[i] + 0.5*(x1[i]^2) + 0.4*RuidoBlanco[3]
    }
    xM = c(x1,x2,x3)
    xM = matrix(xM,(n+500),3)
    xM = xM[501:(n+500),]
    
    if(is.infinite(sum(sum(xM)))||is.nan(sum(sum(xM)))){
      ind=0
    }else{
      ind=1
    }
  }
  return(xM)
}

grafica_3var <- function(Database, Sistema ,type ="l"){
  plot(Database[,1],ylim=c(min(Database),max(Database)),
       xlab = "Periodo", ylab= "Valor", type = "l",
       main = paste0("Grafica de 3 variables del sistema ",Sistema))
  lines(Database[,2],col="red")
  lines(Database[,3],col="darkblue")
  legend("topright", legend = c("x1","x2","x3"),
         lwd=3, col = c("black","red","darkblue"))
}

grafica_4var <- function(Database, Sistema ,type ="l"){
  plot(Database[,1],ylim=c(min(Database),max(Database)),
       xlab = "Periodo", ylab= "Valor", type = "l",
       main = paste0("Grafica de 4 variables del sistema ",Sistema))
  lines(Database[,2],col="red")
  lines(Database[,3],col="darkblue")
  lines(Database[,3],col="darkgreen")
  legend("topright", legend = c("x1","x2","x3","x4"),
         lwd=3, col = c("black","red","darkblue","darkgreen"))
}


Wavelet_3 <- function(Database,nrands = 200){
  ### Wavelet Coherence
  # Se trata como serie
  A= cbind(rep(1:nrow(Database)),Database[,1])
  B= cbind(rep(1:nrow(Database)),Database[,2])
  C= cbind(rep(1:nrow(Database)),Database[,3])
  
  # nrands Numero de aleatorizaciones montecarlo
  # Se calcula la coherencia de ondeleta para todas las combinaciones posibles
  wtc.AB = wtc(A, B, nrands = nrands)
  wtc.BA = wtc(B, A, nrands = nrands)
  wtc.AC = wtc(A, C, nrands = nrands)
  wtc.CA = wtc(C, A, nrands = nrands)
  wtc.BC = wtc(B, C, nrands = nrands)
  wtc.CB = wtc(C, B, nrands = nrands)
  
  
  # Graficamos
  par(mfrow=c(1,1),oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
  
  plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 1, 
       lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
       plot.cb = TRUE, main = "Coherencia de ondeleta: X1 vs X2")
  n = length(A[, 1])
  abline(v = seq(0, n, n/10), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)
  
  plot(wtc.BA, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
       lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
       plot.cb = TRUE, main = "Coherencia de ondeleta: X2 vs X1")
  abline(v = seq(0, n, n/10), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)
  
  plot(wtc.AC, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
       lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
       plot.cb = TRUE, main = "Coherencia de ondeleta: X1 vs X3")
  abline(v = seq(0, n, n/10), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)
  
  plot(wtc.CA, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
       lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
       plot.cb = TRUE, main = "Coherencia de ondeleta: X3 vs X1")
  abline(v = seq(0, n, n/10), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)
  
  plot(wtc.BC, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
       lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
       plot.cb = TRUE, main = "Coherencia de ondeleta: X2 vs X3")
  abline(v = seq(0, n, n/10), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)
  
  plot(wtc.CB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
       lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
       plot.cb = TRUE, main = "Coherencia de ondeleta: X3 vs X2")
  abline(v = seq(0, n, n/10), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)
}


Wavelet_4 <- function(Database,nrands = 200){
  A= cbind(rep(1:nrow(Database)),Database[,1])
  B= cbind(rep(1:nrow(Database)),Database[,2])
  C= cbind(rep(1:nrow(Database)),Database[,3])
  D= cbind(rep(1:nrow(Database)),Database[,4])
  
  wtc.AB = wtc(A, B, nrands = nrands)
  wtc.BA = wtc(B, A, nrands = nrands)
  wtc.AC = wtc(A, C, nrands = nrands)
  wtc.CA = wtc(C, A, nrands = nrands)
  wtc.AD = wtc(A, D, nrands = nrands)
  wtc.DA = wtc(D, A, nrands = nrands)
  wtc.BC = wtc(B, C, nrands = nrands)
  wtc.CB = wtc(C, B, nrands = nrands)
  wtc.BD = wtc(B, D, nrands = nrands)
  wtc.DB = wtc(D, B, nrands = nrands)
  wtc.CD = wtc(C, D, nrands = nrands)
  wtc.DC = wtc(D, C, nrands = nrands)
  
  
  par(mfrow=c(1,1),oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
  
  plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 1, 
       lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
       plot.cb = TRUE, main = "Wavelet Coherence: X1 vs X2")
  n = length(A[, 1])
  abline(v = seq(0, n, n/10), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)
  
  plot(wtc.BA, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 1, 
       lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
       plot.cb = TRUE, main = "Wavelet Coherence: X2 vs X1")
  abline(v = seq(0, n, n/10), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)
  
  plot(wtc.AC, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 1, 
       lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
       plot.cb = TRUE, main = "Wavelet Coherence: X1 vs X3")
  abline(v = seq(0, n, n/10), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)
  
  plot(wtc.CA, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 1, 
       lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
       plot.cb = TRUE, main = "Wavelet Coherence: X3 vs X1")
  abline(v = seq(0, n, n/10), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)
  
  plot(wtc.AD, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 1, 
       lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
       plot.cb = TRUE, main = "Wavelet Coherence: X1 vs X4")
  abline(v = seq(0, n, n/10), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)
  
  plot(wtc.DA, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 1, 
       lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
       plot.cb = TRUE, main = "Wavelet Coherence: X4 vs X1")
  abline(v = seq(0, n, n/10), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)
  
  plot(wtc.BC, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
       lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
       plot.cb = TRUE, main = "Wavelet Coherence: X2 vs X3")
  abline(v = seq(0, n, n/10), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)
  
  plot(wtc.CB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
       lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
       plot.cb = TRUE, main = "Wavelet Coherence: X3 vs X2")
  abline(v = seq(0, n, n/10), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)
  
  plot(wtc.BD, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
       lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
       plot.cb = TRUE, main = "Wavelet Coherence: X2 vs X4")
  abline(v = seq(0, n, n/10), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)
  
  plot(wtc.DB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
       lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
       plot.cb = TRUE, main = "Wavelet Coherence: X4 vs X2")
  abline(v = seq(0, n, n/10), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)
  
  plot(wtc.CD, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
       lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
       plot.cb = TRUE, main = "Wavelet Coherence: X3 vs X4")
  abline(v = seq(0, n, n/10), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)
  
  plot(wtc.DC, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
       lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
       plot.cb = TRUE, main = "Wavelet Coherence: X4 vs X3")
  abline(v = seq(0, n, n/10), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)
  
}



TE3_boot <- function(Database,entropy = 'Shannon', shuffles = 200, bootstrap = 300, semilla = 12345){
  #### Transferencia de entropia
  x_12 = transfer_entropy(Database[,1], Database[,2],q = 1, entropy = entropy, shuffles = shuffles, type = c('quantiles'),
                          quantiles = c(5, 95), nboot = bootstrap, burn = 60, seed = semilla)
  
  x_13 = transfer_entropy(Database[,1], Database[,3],q = 1, entropy = entropy, shuffles = shuffles, type = c('quantiles'),
                          quantiles = c(5, 95), nboot = bootstrap, burn = 60, seed = semilla)
  
  x_23 = transfer_entropy(Database[,2], Database[,3],q = 1, entropy = entropy, shuffles = shuffles, type = c('quantiles'),
                          quantiles = c(5, 95), nboot = bootstrap, burn = 60, seed = semilla)
  
  TE = matrix(data = c(0, x_12$coef[1,1], x_13$coef[1,1],
                       x_12$coef[2,1], 0 , x_23$coef[1,1],
                       x_13$coef[2,1], x_23$coef[2,1], 0), nrow = 3, ncol = 3, byrow = T)
  pvalue = matrix(data = c(0, x_12$coef[1,4], x_13$coef[1,4],
                           x_12$coef[2,4], 0 , x_23$coef[1,4],
                           x_13$coef[2,4], x_23$coef[2,4], 0), nrow = 3, ncol = 3, byrow = T)
  return(list(TE,pvalue))
}


TE4_boot <- function(Database,entropy = 'Shannon', shuffles = 200, bootstrap = 300, semilla = 12345){
  x_12 = transfer_entropy(Database[,1], Database[,2],q = 1, entropy = entropy, shuffles = shuffles, type = c('quantiles'),
                          quantiles = c(5, 95), nboot = bootstrap, burn = 60, seed = semilla)
  
  x_13 = transfer_entropy(Database[,1], Database[,3],q = 1, entropy = entropy, shuffles = shuffles, type = c('quantiles'),
                          quantiles = c(5, 95), nboot = bootstrap, burn = 60, seed = semilla)
  
  x_14 = transfer_entropy(Database[,1], Database[,4],q = 1, entropy = entropy, shuffles = shuffles, type = c('quantiles'),
                          quantiles = c(5, 95), nboot = bootstrap, burn = 60, seed = semilla)
  
  x_23 = transfer_entropy(Database[,2], Database[,3],q = 1, entropy = entropy, shuffles = shuffles, type = c('quantiles'),
                          quantiles = c(5, 95), nboot = bootstrap, burn = 60, seed = semilla)
  
  x_24 = transfer_entropy(Database[,2], Database[,4],q = 1, entropy = entropy, shuffles = shuffles, type = c('quantiles'),
                          quantiles = c(5, 95), nboot = bootstrap, burn = 60, seed = semilla)
  
  x_34 = transfer_entropy(Database[,3], Database[,4],q = 1, entropy = entropy, shuffles = shuffles, type = c('quantiles'),
                          quantiles = c(5, 95), nboot = bootstrap, burn = 60, seed = semilla)
  
  
  TE = matrix(data = c(0, x_12$coef[1,1], x_13$coef[1,1],x_14$coef[1,1],
                       x_12$coef[2,1], 0 , x_23$coef[1,1], x_24$coef[1,1],
                       x_13$coef[2,1], x_23$coef[2,1], 0, x_34$coef[1,1],
                       x_14$coef[2,1], x_24$coef[2,1] , x_34$coef[2,1], 0), nrow = 4, ncol = 4, byrow = T)
  
  
  pvalue = matrix(data = c(0, x_12$coef[1,4], x_13$coef[1,4],x_14$coef[1,4],
                           x_12$coef[2,4], 0 , x_23$coef[1,4], x_24$coef[1,4],
                           x_13$coef[2,4], x_23$coef[2,4], 0, x_34$coef[1,4],
                           x_14$coef[2,4], x_24$coef[2,4] , x_34$coef[2,4], 0), nrow = 4, ncol = 4, byrow = T)
  
  return(list(TE,pvalue))
  
}


TE3_sbt <- function(Database,entropy = 'Shannon',shuffles = 200, semilla = 12345){
  TE = matrix(data = c(0, calc_te(Database[,1],Database[,2],entropy = entropy,shuffles = shuffles,seed =semilla ), calc_te(Database[,1],Database[,3],entropy = entropy,shuffles = shuffles,seed =semilla),
                       calc_te(Database[,2],Database[,1],entropy = entropy,shuffles = shuffles,seed =semilla), 0 , calc_te(Database[,2],Database[,3],entropy = entropy,shuffles = shuffles,seed =semilla),
                       calc_te(Database[,3],Database[,1],entropy = entropy,shuffles = shuffles,seed =semilla), calc_te(Database[,3],Database[,2],entropy = entropy,shuffles = shuffles,seed =semilla), 0), nrow = 3, ncol = 3, byrow = T)
  
  return(TE)
}


TE4_sbt <- function(Database,entropy = 'Shannon',shuffles = 200, semilla = 12345){
  TE = matrix(data = c(0, calc_te(Database[,1],Database[,2],entropy = entropy,shuffles = shuffles,seed =semilla ), calc_te(Database[,1],Database[,3],entropy = entropy,shuffles = shuffles,seed =semilla),calc_te(Database[,1],Database[,4],entropy = entropy,shuffles = shuffles,seed =semilla),
                       calc_te(Database[,2],Database[,1],entropy = entropy,shuffles = shuffles,seed =semilla), 0 , calc_te(Database[,2],Database[,3],entropy = entropy,shuffles = shuffles,seed =semilla),calc_te(Database[,2],Database[,4],entropy = entropy,shuffles = shuffles,seed =semilla),
                       calc_te(Database[,3],Database[,1],entropy = entropy,shuffles = shuffles,seed =semilla), calc_te(Database[,3],Database[,2],entropy = entropy,shuffles = shuffles,seed =semilla), 0, calc_te(Database[,3],Database[,4],entropy = entropy,shuffles = shuffles,seed =semilla),
                       calc_te(Database[,4],Database[,1],entropy = entropy,shuffles = shuffles,seed =semilla), calc_te(Database[,4],Database[,2],entropy = entropy,shuffles = shuffles,seed =semilla), calc_te(Database[,4],Database[,3],entropy = entropy,shuffles = shuffles,seed =semilla), 0), nrow = 4, ncol = 4, byrow = T)
  
  return(TE)
}



#grangertest(x1,x2)
varfit <- VAR(cbind(x1,x2),p=1,type = 'const')
sum_var =summary(varfit)
sum_var$varresult$x1
VAR(select(as.data.frame(Database),3,1),p=1,type = 'const')

coef.transfer_entropy(A,B)
select(Database,1,2)
Database
# SIMULACION SISTEMA 1
Database <- Sistema_1(2048)
grafica_3var(Database,"1")
Wavelet_3(Database)
TE3_boot(Database)
TE3_sbt(Database)

fwrite(Database,file="F:/Tarea/TESIS/IMAGENES/Sistema_1_2048.csv")

############### 2 ###########

Sistema_2 <- function(n, semilla = 12345){
  set.seed(semilla)
  
  x1 = c(0,runif(1))
  x2 = c(0,runif(1))
  x3 = c(0,runif(1))
  ind = 0
  media <- c(0,0,0)
  sigma<-diag(3)
  
  set.seed(semilla)
  while (ind==0) {
    for( i in seq(2,(n-1+500))  ){
      RuidoBlanco=rmvnorm(1,mean = media,sigma = sigma)
      
      x1[i+1] = 0.7*x1[i] + RuidoBlanco[1]
      x2[i+1] = 0.3*x2[i] + 0.5*x1[i]*x2[i-1] + RuidoBlanco[2]
      x3[i+1] = 0.3*x3[i] + 0.5*x1[i]*(x3[i-1]) + RuidoBlanco[3]
    }
    xM = c(x1,x2,x3)
    xM = matrix(xM,(n+500),3)
    xM = xM[501:(n+500),]
    
    if(is.infinite(sum(sum(xM)))||is.nan(sum(sum(xM)))){
      ind=0
    }else{
      ind=1
    }
  }
  return(xM)
}

#grangertest(x1,x2)
varfit <- VAR(cbind(x1,x2),p=1,type = 'const')
sum_var =summary(varfit)
sum_var$varresult$x1

coef.transfer_entropy(A,B)


# SIMULACION SISTEMA 2
Database <- Sistema_2(2048)
grafica_3var(Database,"2")
Wavelet_3(Database)
TE3_boot(Database)
TE3_sbt(Database)

fwrite(Database,file="F:/Tarea/TESIS/IMAGENES/Sistema_2_2048.csv")

############### 3 ###########

Sistema_3 <- function(n,ac = 0.5 , semilla = 12345){
  set.seed(semilla)
  x1 = c(0,runif(1))
  x2 = c(0,runif(1))
  x3 = c(0,runif(1))
  ind = 0
  
  while (ind==0) {
    for( i in seq(2,(n-1+500)) ){
      
      x1[i+1] = 1.4 - (x1[i]^2) + 0.3*x1[i-1]
      x2[i+1] = 1.4 - ac*x1[i]*x2[i] - (1-ac)*(x2[i]^2) + 0.3*x2[i-1]
      x3[i+1] = 1.4 - ac*x2[i]*x3[i] - (1-ac)*(x3[i]^2) + 0.3*x3[i-1]
    }
    xM = c(x1,x2,x3)
    xM = matrix(xM,(n+500),3)
    xM = xM[501:(n+500),]
    
    if(is.infinite(sum(sum(xM)))||is.nan(sum(sum(xM)))){
      ind=0
    }else{
      ind=1
    }
  }
  return(xM)
}


#grangertest(x1,x2)
varfit <- VAR(cbind(x1,x2),p=1,type = 'const')
sum_var =summary(varfit)
sum_var$varresult$x1

coef.transfer_entropy(A,B)


# SIMULACION SISTEMA 3
Database <- Sistema_3(2048)
grafica_3var(Database,"3")
Wavelet_3(Database)
TE3_boot(Database)
TE3_sbt(Database)

fwrite(Database,file="F:/Tarea/TESIS/IMAGENES/Sistema_3_2048.csv")

############### 4 ###########

Sistema_4 <- function(n,ac = 0.5 , semilla = 12345){
  set.seed(semilla)
  x1 = c(0,runif(1))
  x2 = c(0,runif(1))
  x3 = c(0,runif(1))
  x4 = c(0,runif(1))
  ind = 0
  set.seed(semilla)
  
  while (ind==0) {
    for( i in seq(2,(n-1+500)) ){
      
      x1[i+1] = 1.4 - (x1[i]^2) + 0.3*x1[i-1]
      x4[i+1] = 1.4 - (x4[i]^2) + 0.3*x4[i-1]
      x2[i+1] = 1.4 - (0.5*ac*(x1[i]+x3[i]) + (1-ac)*x2[i])^2 + 0.3*x2[i-1]
      x3[i+1] = 1.4 - (0.5*ac*(x2[i]+x4[i]) + (1-ac)*x3[i])^2 + 0.3*x3[i-1]
    }
    xM = c(x1,x2,x3,x4)
    xM = matrix(xM,(n+500),4)
    xM = xM[501:(n+500),]
    
    if(is.infinite(sum(sum(xM)))||is.nan(sum(sum(xM)))){
      ind=0
    }else{
      ind=1
    }
  }
  return(xM)
}


#grangertest(x1,x2)
varfit <- VAR(cbind(x1,x2),p=1,type = 'const')
sum_var =summary(varfit)
sum_var$varresult$x1

coef.transfer_entropy(A,B)


# SIMULACION SISTEMA 1
Database <- Sistema_4(2048)
grafica_4var(Database,"4")
Wavelet_4(Database)
TE4_boot(Database)
TE4_sbt(Database)


fwrite(Database,file="F:/Tarea/TESIS/IMAGENES/Sistema_4_2048.csv")

############### 5 ########### 

Sistema_5 <- function(n,ac = 0.5 , semilla = 12345){
  ind = 0
  xM = Sistema_3(n, ac = ac, semilla = semilla)
  
  while (ind==0) {
    outlier = sample(n,size = round(n*0.01),replace = F)
    for( i in outlier ){
      xM[i,1] = xM[i,1] + runif(1)
      xM[i,2] = xM[i,2] + runif(1)
      xM[i,3] = xM[i,3] + runif(1)
    }
    if(is.infinite(sum(sum(xM)))||is.nan(sum(sum(xM)))){
      ind=0
    }else{
      ind=1
    }
  }
  return(xM)
}


#grangertest(x1,x2)
varfit <- VAR(cbind(x1,x2),p=1,type = 'const')
sum_var =summary(varfit)
sum_var$varresult$x1

coef.transfer_entropy(A,B)


# SIMULACION SISTEMA 5
Database <- Sistema_5(2048)
grafica_3var(Database,"5")
Wavelet_3(Database)
TE3_boot(Database)
TE3_sbt(Database)

fwrite(Database,file="F:/Tarea/TESIS/IMAGENES/Sistema_5_2048.csv")

############### 6 ###########

Sistema_6 <- function(n,ac = 0.5 , semilla = 12345){
  
  xM = Sistema_3(n,ac,semilla)
  ind = 0 
  n_t1 = c(0)
  n_t2 = c(0)
  n_t3 = c(0)
  while (ind==0) {
    for( i in seq(2,n) ){
      n_t1[i] = (n_t1[i-1] + rnorm(1))
      n_t2[i] = (n_t2[i-1] + rnorm(1))
      n_t3[i] = (n_t3[i-1] + rnorm(1))
      xM[i,1] = xM[i,1] + n_t1[i]
      xM[i,2] = xM[i,2] + n_t2[i]
      xM[i,3] = xM[i,3] + n_t3[i]
    }
    if(is.infinite(sum(sum(xM)))||is.nan(sum(sum(xM)))){
      ind=0
    }else{
      ind=1
    }
  }
  return(xM)
}


#grangertest(x1,x2)
varfit <- VAR(cbind(x1,x2),p=1,type = 'const')
sum_var =summary(varfit)
sum_var$varresult$x1

coef.transfer_entropy(A,B)


# SIMULACION SISTEMA 6
Database <- Sistema_6(2048)
grafica_3var(Database,"6")
Wavelet_3(Database)
TE3_boot(Database)
TE3_sbt(Database)

fwrite(Database,file="F:/Tarea/TESIS/IMAGENES/Sistema_6_2048.csv")

############### 7 ###########

Sistema_7 <- function(n,ac = 0.5 , semilla = 12345){
  xM = Sistema_3(n,ac,semilla)
  ind = 0 
  
  while (ind==0) {
    a1 = rnorm(1,0.01,0.02)
    a2 = rnorm(1,0.01,0.02)
    a3 = rnorm(1,0.01,0.02)
    t = seq(1,n)
    n_t1 = a1*t
    n_t2 = a2*t
    n_t3 = a3*t
    for( i in seq(1,n) ){
      xM[i,1] = xM[i,1] + n_t1[i]
      xM[i,2] = xM[i,2] + n_t2[i]
      xM[i,3] = xM[i,3] + n_t3[i]
    }
    if(is.infinite(sum(sum(xM)))||is.nan(sum(sum(xM)))){
      ind=0
    }else{
      ind=1
    }
  }
  return(xM)
}


#grangertest(x1,x2)
varfit <- VAR(cbind(x1,x2),p=1,type = 'const')
sum_var =summary(varfit)
sum_var$varresult$x1

coef.transfer_entropy(A,B)


# SIMULACION SISTEMA 7
Database <- Sistema_7(2048)
grafica_3var(Database,"7")
Wavelet_3(Database)
TE3_boot(Database)
TE3_sbt(Database)

fwrite(Database,file="F:/Tarea/TESIS/IMAGENES/Sistema_7_2048.csv")

############### 8 ###########

Sistema_8 <- function(n, semilla = 12345,a0 = 0.2,a1 = 0.9,b1 = 0.1,g= 1){
  
  ind = 0
  media <- c(0,0,0)
  sigma<-diag(3)
  
  set.seed(semilla)
  
  # Error = rmvnorm(1,mean = media,sigma = sigma)
  # Error2 = rmvnorm(1,mean = media,sigma = sigma)
  # 
  Error = rnorm(1)
  Error2 = rnorm(1)
  
  error1 = c(Error,Error2)
  error2 = c(Error,Error2)
  error3 = c(Error,Error2)
  
  
  Database = Sistema_2(n,semilla)
  
  sigma2_1 = c(a0,a0+(a1*Error2^2)+b1*a0)
  sigma2_2 = c(a0,a0+(a1*Error2^2)+b1*a0)
  sigma2_3 = c(a0,a0+(a1*Error2^2)+b1*a0)
  
  zt1 = c(sigma2_1[1]*Error,sigma2_1[2]*Error2)
  zt2 = c(sigma2_2[1]*Error,sigma2_2[2]*Error2)
  zt3 = c(sigma2_3[1]*Error,sigma2_3[2]*Error2)
  
  y1 = c(Database[1,1]+g*zt1[1],Database[2,1]+g*zt1[2])
  y2 = c(Database[1,2]+g*zt2[1],Database[2,2]+g*zt2[2])
  y3 = c(Database[1,3]+g*zt3[1],Database[2,3]+g*zt3[2])
  
  
  while (ind==0) {
    for( i in seq(2,(n-1)) ){
      #Error = rmvnorm(1,mean = media,sigma = sigma)
      Error = rnorm(1)
      error1[i+1] = Error
      error2[i+1] = Error
      error3[i+1] = Error
      
      sigma2_1[i+1] = a0 + a1*(error1[i]^2) + b1*(sigma2_1[i]^2)
      sigma2_2[i+1] = a0 + a1*(error2[i]^2) + b1*(sigma2_2[i]^2)
      sigma2_3[i+1] = a0 + a1*(error3[i]^2) + b1*(sigma2_3[i]^2)
      
      zt1[i+1] = sigma2_1[i+1]*error1[i+1]
      zt2[i+1] = sigma2_2[i+1]*error2[i+1]
      zt3[i+1] = sigma2_3[i+1]*error3[i+1]
      
      y1[i+1] = Database[i+1,1]+g*zt1[i+1]
      y2[i+1] = Database[i+1,2]+g*zt2[i+1]
      y3[i+1] = Database[i+1,3]+g*zt3[i+1]
    }
    xM = c(y1,y2,y3)
    xM = matrix(xM,(n),3)
    # xM = xM[501:(n),]
    
    if(is.infinite(sum(sum(xM)))||is.nan(sum(sum(xM)))){
      ind=0
    }else{
      ind=1
    }
  }
  return(xM)
}



#grangertest(x1,x2)
varfit <- VAR(cbind(x1,x2),p=1,type = 'const')
sum_var =summary(varfit)
sum_var$varresult$x1

coef.transfer_entropy(A,B)


# SIMULACION SISTEMA 8
Database <- Sistema_8(2048)
grafica_3var(Database,"8")
Wavelet_3(Database)
TE3_boot(Database)
TE3_sbt(Database)


fwrite(Database,file="F:/Tarea/TESIS/IMAGENES/Sistema_8_2048.csv")


############### 9 ###########

Sistema_9 <- function(n,g=0.2, semilla = 12345,a0 = 0.2,a1 = 0.9,b1 = 0.1){
  ind = 0
  
  media <- c(0,0,0)
  sigma<-diag(3)
  set.seed(semilla)
  x1 = runif(3)
  x2 = runif(3)
  x3 = runif(3)
  
  ERROR =  rnorm(1)
  ERROR2 = rnorm(1)
  ERROR3 = rnorm(1)
  
  sigma2_1 = c(a0,a0+(a1*ERROR^2)+b1*a0, a0+(a1*ERROR2^2)+b1*(a0+(a1*ERROR3^2)+b1*a0))
  sigma2_2 = c(a0,a0+(a1*ERROR^2)+b1*a0, a0+(a1*ERROR2^2)+b1*(a0+(a1*ERROR3^2)+b1*a0))
  sigma2_3 = c(a0,a0+(a1*ERROR^2)+b1*a0, a0+(a1*ERROR2^2)+b1*(a0+(a1*ERROR3^2)+b1*a0))
  
  zt1 = c(sigma2_1[1]*ERROR,sigma2_1[2]*ERROR2,sigma2_1[3]*ERROR3)
  zt2 = c(sigma2_2[1]*ERROR,sigma2_2[2]*ERROR2,sigma2_2[3]*ERROR3)
  zt3 = c(sigma2_3[1]*ERROR,sigma2_3[2]*ERROR2,sigma2_3[3]*ERROR3)
  
  y1 = c(x1[1]+g*zt1[1],x1[2]+g*zt1[2],x1[3]+g*zt1[3])
  y2 = c(x2[1]+g*zt2[1],x2[2]+g*zt2[2],x2[3]+g*zt2[3])
  y3 = c(x3[1]+g*zt3[1],x3[2]+g*zt3[2],x3[3]+g*zt3[3])
  
  error1=c(ERROR,ERROR2,ERROR3)
  error2=c(ERROR,ERROR2,ERROR3)
  error3=c(ERROR,ERROR2,ERROR3)
  
  while (ind==0) {
    for(i in seq(3,(n-1))  ){
      RuidoBlanco=rmvnorm(1,mean = media,sigma = sigma)
      
      x1[i+1] = 0.4*x1[i] + 0.4*x2[i] + 0.5*x3[i] + 0.2*x1[i-1] - 0.2*x2[i-1] - 0.2*x1[i-2] + 0.15*x2[i-2] + 0.1*x3[i-2] + RuidoBlanco[1]
      x2[i+1] = 0.6*x2[i] + 0.2*x2[i-1] + 0.2*x2[i-2] + RuidoBlanco[2]
      x3[i+1] = 0.4*x3[i] + 0.3*x3[i-1] + 0.3*x3[i-2] + RuidoBlanco[3]
      
      Error = rnorm(1)
      error1[i+1] = Error
      error2[i+1] = Error
      error3[i+1] = Error
      
      sigma2_1[i+1] = a0 + a1*(error1[i]^2) + b1*(sigma2_1[i]^2)
      sigma2_2[i+1] = a0 + a1*(error2[i]^2) + b1*(sigma2_2[i]^2)
      sigma2_3[i+1] = a0 + a1*(error3[i]^2) + b1*(sigma2_3[i]^2)
      
      zt1[i+1] = sigma2_1[i+1]*error1[i+1]
      zt2[i+1] = sigma2_2[i+1]*error2[i+1]
      zt3[i+1] = sigma2_3[i+1]*error3[i+1]
      
      y1[i+1] = x1[i+1] +g*zt1[i+1]
      y2[i+1] = x2[i+1] +g*zt2[i+1]
      y3[i+1] = x3[i+1] +g*zt3[i+1]
    }
    xM = c(y1,y2,y3)
    xM = matrix(xM,n,3)
    if(is.infinite(sum(sum(xM)))||is.nan(sum(sum(xM)))){
      ind=0
    }else{
      ind=1
    }
  }
  return(xM)
}


#grangertest(x1,x2)
varfit <- VAR(cbind(x1,x2),p=1,type = 'const')
sum_var =summary(varfit)
sum_var$varresult$x1

coef.transfer_entropy(A,B)
-7.560+.8136

# SIMULACION SISTEMA 9
#Database <- Sistema_9(512)
Database1 <- Sistema_9(1100)
Database2<- Sistema_9(1000)
Database2[,1] = Database2[,1] + tail(Database[,1],1)
Database2[,2] = Database2[,2] + tail(Database[,2],1)
Database2[,3] = Database2[,3] + tail(Database[,3],1)
Database<-rbind(Database,Database2)
Database<-tail(Database,2048)
grafica_3var(Database,"9")
Wavelet_3(Database)
TE3_boot(Database)
TE3_sbt(Database)

fwrite(Database,file="F:/Tarea/TESIS/IMAGENES/Sistema_9_512.csv")






############### Ejemplo 1

Ejemplo_1 <- function(semilla = 12345){
  tiempo = seq(0,50,by=1/12)
  y_t=c()
  set.seed(semilla)
  for (i in 1:length(tiempo)) {
    p1=10
    if (tiempo[i] >= 20 && tiempo[i] <= 30) {
      p2 = 5
      y_t[i]=cos(2*pi*tiempo[i]/p1) + cos(2*pi*tiempo[i]/p2) + runif(1)
    }else{
      p2 = 3
      y_t[i]=cos(2*pi*tiempo[i]/p1) + cos(2*pi*tiempo[i]/p2) + runif(1)
    }
  }
  return(y_t)
}


Database<-Ejemplo_1()


plot(Database,type= "l")
A= cbind(rep(1:length(Database))/12,Database)









############### Ejemplo 2


Ejemplo_2 <- function(semilla = 12345){
  tiempo = seq(0,50,by=1/12)
  
  y_t=c()
  x_t=c()
  set.seed(semilla)
  for (i in 1:length(tiempo)) {
    if(tiempo[i]<=25){
      x_t[i] = 4*sin((2*pi/3) * (tiempo[i]+5/12)) - 3*sin((2*pi/6)*(tiempo[i] -10/12)) + runif(1)
      
      y_t[i] = sin(2*pi*tiempo[i]/3) + 3*sin(2*pi*tiempo[i]/6) + runif(1)
      
    }else{
      x_t[i] = 4*sin((2*pi/3) * (tiempo[i]-5/12)) - 3*sin((2*pi/6)*(tiempo[i] + 10/12)) + runif(1)
      
      y_t[i] = sin(2*pi*tiempo[i]/3) + 3*sin(2*pi*tiempo[i]/6) + runif(1)
      
    }
  }
  return(list(x_t,y_t))
}
Ejemplo_2()


library(tseries)
library(quantmod)
library(fImport)
library(urca)


Database <- read.csv("/data_conjunta_pura.csv")
Database2 <- read.csv("/data_conjunta.csv")

df=data.frame( matrix(ncol = (ncol(Database)-1), nrow = (nrow(Database)-1)))
colnames(df) = c("Date", "Precio", "Volumen", "Vocabulario")
n = nrow(Database)
df[,1] = Database[2:n,2]

#Retorno 
for (i in 3:5) {
  n = length(Database[,i])
  df[,(i-1)] = (Database[2:n,i]-Database[(1:(n-1)),i])/Database[(1:(n-1)),i]
}

plot(ts(df[2:4]), main= "Comportamiento del precio, volumen y vacabulario")

df2=data.frame( matrix(ncol = (ncol(df)), nrow = (nrow(df))))
colnames(df2) = c("Date", "Precio", "Volumen", "Vocabulario")
n = nrow(df)
df2[,1] = df[,1]

#Estandarización
for (i in 2:4) {
  n = length(df[,i])
  media = mean(df[,i])
  desv = sd(df[,i])
  df2[i] = (df[i]-media)/desv
}

# plot(ts(df2[2:4]))
# plot(diff(ts(df2[2:4]), differences = 1, lag=5))


adf.test(ts(Database2[3])) # Prueba del precio en la serie usada
# Augmented Dickey-Fuller Test
# 
# data:  ts(Database2[3])
# Dickey-Fuller = -3.0716, Lag order = 7, p-value = 0.1247
# alternative hypothesis: stationary
# No es estacionaria


adf.test(ts(df2[2])) # Prueba de estacionaridad para el precio
adf.test(ts(df2[3])) # Prueba de estacionaridad para el volumen
adf.test(ts(df2[4])) # Prueba de estacionaridad para el vocabulario

# Augmented Dickey-Fuller Test
# 
# data:  ts(df2[2])
# Dickey-Fuller = -6.7582, Lag order = 7, p-value = 0.01
# alternative hypothesis: stationary
# 
# 
# Augmented Dickey-Fuller Test
# 
# data:  ts(df2[3])
# Dickey-Fuller = -7.7393, Lag order = 7, p-value = 0.01
# alternative hypothesis: stationary
# 
# Augmented Dickey-Fuller Test
# 
# data:  ts(df2[4])
# Dickey-Fuller = -8.0034, Lag order = 7, p-value = 0.01
# alternative hypothesis: stationary
# 
# El valor p sugiere que los datos son muy poco probables dada la hipótesis nula (de integración yt),
# por lo que es más probable la hipótesis alternativa (estacionaria yt)

pp.test(ts(df2[2]), alternative="stationary")
pp.test(ts(df2[3]), alternative="stationary")
pp.test(ts(df2[4]), alternative="stationary")

# Phillips-Perron Unit Root Test
# 
# data:  ts(df2[2])
# Dickey-Fuller Z(alpha) = -503.24, Truncation lag parameter = 5, p-value = 0.01
# alternative hypothesis: stationary
# 
# Phillips-Perron Unit Root Test
# 
# data:  ts(df2[3])
# Dickey-Fuller Z(alpha) = -362.16, Truncation lag parameter = 5, p-value = 0.01
# alternative hypothesis: stationary
# 
# Phillips-Perron Unit Root Test
# 
# data:  ts(df2[4])
# Dickey-Fuller Z(alpha) = -398.7, Truncation lag parameter = 5, p-value = 0.01
# alternative hypothesis: stationary
#
# Sucede lo mismo que el test de Dickey-Fuller, además nos sugiere un lag de 5

df.price <- ur.df(ts(df2[2]), type='none',lags=5, selectlags=c("AIC"))
summary(df.price)

df.vol <- ur.df(ts(df2[3]), type='none',lags=5, selectlags=c("AIC"))
summary(df.vol)

df.voc <- ur.df(ts(df2[4]), type='none',lags=5, selectlags=c("AIC"))
summary(df.voc)

# Residual standard error: 0.9452 on 483 degrees of freedom
# Multiple R-squared:  0.5109,	Adjusted R-squared:  0.5068 
# F-statistic: 126.1 on 4 and 483 DF,  p-value: < 2.2e-16
# 
# 
# Value of test-statistic is: -10.844 
# 
# 
# Residual standard error: 0.9525 on 484 degrees of freedom
# Multiple R-squared:  0.4261,	Adjusted R-squared:  0.4225 
# F-statistic: 119.8 on 3 and 484 DF,  p-value: < 2.2e-16
# 
# 
# Value of test-statistic is: -11.687 
# 
# 
# Residual standard error: 0.9747 on 484 degrees of freedom
# Multiple R-squared:  0.4616,	Adjusted R-squared:  0.4583 
# F-statistic: 138.3 on 3 and 484 DF,  p-value: < 2.2e-16
# 
# 
# Value of test-statistic is: -12.7936 
#
# De igual manera a adf.test() calcula el test Dickey-Fuller para la hipótesis nula que x tiene raíz unitaria
# Se confirma que la serie es estacionaria con las tres pruebas
Database_fin = df2[,2:4]
grafica_3var(Database_fin," Reddit ")
Wavelet_3(Database_fin)
TE3_boot(Database_fin)
TE3_sbt(Database_fin)



