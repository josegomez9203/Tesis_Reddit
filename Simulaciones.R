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

coef.transfer_entropy(A,B)


# SIMULACION SISTEMA 1
Database <- Sistema_1(200)
grafica_3var(Database,"1")
Wavelet_3(Database)
TE3_boot(Database)
TE3_sbt(Database)


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
Database <- Sistema_2(200)
grafica_3var(Database,"2")
Wavelet_3(Database)
TE3_boot(Database)
TE3_sbt(Database)


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
Database <- Sistema_3(200)
grafica_3var(Database,"3")
Wavelet_3(Database)
TE3_boot(Database)
TE3_sbt(Database)



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
Database <- Sistema_4(200)
grafica_4var(Database,"4")
Wavelet_4(Database)
TE4_boot(Database)
TE4_sbt(Database)




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
Database <- Sistema_5(200)
grafica_3var(Database,"5")
Wavelet_3(Database)
TE3_boot(Database)
TE3_sbt(Database)



############### 6 ###########

Sistema_6 <- function(n,ac = 0.5 , semilla = 12345){
  
  xM = Sistema_3(n,ac,semilla)
  ind = 0 
  n_t = c(0)
  while (ind==0) {
    for( i in seq(2,n) ){
      n_t[i] = (n_t[i-1] + rnorm(1))
      xM[i,1] = xM[i,1] + n_t[i]
      xM[i,2] = xM[i,2] + n_t[i]
      xM[i,3] = xM[i,3] + n_t[i]
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
Database <- Sistema_6(200)
grafica_3var(Database,"6")
Wavelet_3(Database)
TE3_boot(Database)
TE3_sbt(Database)


############### 7 ###########

Sistema_7 <- function(n,ac = 0.5 , semilla = 12345){
  xM = Sistema_3(n,ac,semilla)
  ind = 0 
  
  while (ind==0) {
    a1 = rnorm(n,0.01,0.02)
    a2 = rnorm(n,0.01,0.02)
    a3 = rnorm(n,0.01,0.02)
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
Database <- Sistema_7(200)
grafica_3var(Database,"7")
Wavelet_3(Database)
TE3_boot(Database)
TE3_sbt(Database)


############### 8 ###########

Sistema_8 <- function(n, semilla = 12345,a0 = 0.2,a1 = 0.9,b1 = 0.1,g= 1){
  
  ind = 0
  media <- c(0,0,0)
  sigma<-diag(3)
  
  set.seed(semilla)
  
  Error = rmvnorm(1,mean = media,sigma = sigma)
  Error2 = rmvnorm(1,mean = media,sigma = sigma)
  
  error1 = c(Error[1],Error2[1])
  error2 = c(Error[2],Error2[2])
  error3 = c(Error[3],Error2[3])
  
  
  Database = Sistema_2((n),semilla)
  
  sigma2_1 = c(a0,a0+(a1*Error2[1]^2)+b1*a0)
  sigma2_2 = c(a0,a0+(a1*Error2[2]^2)+b1*a0)
  sigma2_3 = c(a0,a0+(a1*Error2[3]^2)+b1*a0)
  
  zt1 = c(sigma2_1[1]*Error[1],sigma2_1[2]*Error2[1])
  zt2 = c(sigma2_2[1]*Error[2],sigma2_2[2]*Error2[2])
  zt3 = c(sigma2_3[1]*Error[3],sigma2_3[2]*Error2[3])
  
  y1 = c(Database[1,1]+g*zt1[1],Database[2,1]+g*zt1[2])
  y2 = c(Database[1,2]+g*zt2[1],Database[2,2]+g*zt2[2])
  y3 = c(Database[1,3]+g*zt3[1],Database[2,3]+g*zt3[2])
  
  
  while (ind==0) {
    for( i in seq(2,(n-1)) ){
      Error = rmvnorm(1,mean = media,sigma = sigma)
      
      error1[i+1] = Error[1]
      error2[i+1] = Error[2]
      error3[i+1] = Error[3]
      
      sigma2_1[i+1] = a0 + a1*(error1[i]^2) + b1*(sigma2_1[i]^2)
      sigma2_2[i+1] = a0 + a1*(error2[i]^2) + b1*(sigma2_2[i]^2)
      sigma2_3[i+1] = a0 + a1*(error3[i]^2) + b1*(sigma2_3[i]^2)
      
      zt1[i+1] = sigma2_1[i+1]*error1[i+1]
      zt2[i+1] = sigma2_2[i+1]*error2[i+1]
      zt3[i+1] = sigma2_3[i+1]*error3[i+1]
      
      y1[i+1] = Database[i,1]+g*zt1[i+1]
      y2[i+1] = Database[i,2]+g*zt2[i+1]
      y3[i+1] = Database[i,3]+g*zt3[i+1]
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
Database <- Sistema_8(500)
grafica_3var(Database,"8")
Wavelet_3(Database)
TE3_boot(Database)
TE3_sbt(Database)








### Wavelet Coherence
# Se trata como serie
A= cbind(rep(1:length(y1)),y1)
B= cbind(rep(1:length(y2)),y2)
C= cbind(rep(1:length(y3)),y3)

nrands = 200 # Numero de aleatorizaciones montecarlo
wtc.AB = wtc(A, B, nrands = nrands)
wtc.AC = wtc(A, C, nrands = nrands)
wtc.BC = wtc(B, C, nrands = nrands)


# Graficamos
par(mfrow=c(1,1),oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)

plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 1, 
     lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
     plot.cb = TRUE, main = "Wavelet Coherence: Y1 vs Y2")
n = length(A[, 1])
abline(v = seq(100, n, 100), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)

plot(wtc.AC, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
     plot.cb = TRUE, main = "Wavelet Coherence: Y1 vs Y3")
abline(v = seq(100, n, 100), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)

plot(wtc.BC, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
     plot.cb = TRUE, main = "Wavelet Coherence: Y2 vs Y3")
abline(v = seq(100, n, 100), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)


#### Transferencia de entropia
transfer_entropy(y1, y2,q = 1, entropy = c('Shannon'), shuffles = 200, type = c('quantiles'),
                 quantiles = c(5, 95), nboot = 300, burn = 60, seed = NULL)

transfer_entropy(y1, y3,q = 1, entropy = c('Shannon'), shuffles = 200, type = c('quantiles'),
                 quantiles = c(5, 95), nboot = 300, burn = 60, seed = NULL)

transfer_entropy(y2, y3,q = 1, entropy = c('Shannon'), shuffles = 200, type = c('quantiles'),
                 quantiles = c(5, 95), nboot = 300, burn = 60, seed = NULL)

#Calculo de transferencia de entropía sin usar bootstrap
#calc_ete(y1,y2) # Transferencia de entropía efectiva
calc_te(y1,y2)
calc_te(y2,y1)

#calc_ete(y1,y3) # Transferencia de entropía efectiva
calc_te(y1,y3)
calc_te(y3,y1)

#calc_ete(y3,y2) # Transferencia de entropía efectiva
calc_te(y3,y2)
calc_te(y2,y3)


############### 10 ###########

Sistema_10 <- function(n,g=0.2, semilla = 12345,a0 = 0.2,a1 = 0.9,b1 = 0.1){
  ind = 0
  
  x1 = runif(3)
  x2 = runif(3)
  x3 = runif(3)
  
  ERROR = rmvnorm(1,mean = media,sigma = sigma)
  ERROR2 = rmvnorm(1,mean = media,sigma = sigma)
  ERROR3 = rmvnorm(1,mean = media,sigma = sigma)
  
  sigma2_1 = c(a0,a0+(a1*ERROR[1]^2)+b1*a0, a0+(a1*ERROR2[1]^2)+b1*(a0+(a1*ERROR[1]^2)+b1*a0))
  sigma2_2 = c(a0,a0+(a1*ERROR[2]^2)+b1*a0, a0+(a1*ERROR2[2]^2)+b1*(a0+(a1*ERROR[2]^2)+b1*a0))
  sigma2_3 = c(a0,a0+(a1*ERROR[3]^2)+b1*a0, a0+(a1*ERROR2[3]^2)+b1*(a0+(a1*ERROR[3]^2)+b1*a0))
  
  zt1 = c(sigma2_1[1]*ERROR[1],sigma2_1[2]*ERROR2[1],sigma2_1[3]*ERROR3[1])
  zt2 = c(sigma2_2[1]*ERROR[2],sigma2_2[2]*ERROR2[2],sigma2_2[3]*ERROR3[2])
  zt3 = c(sigma2_3[1]*ERROR[3],sigma2_3[2]*ERROR2[3],sigma2_3[3]*ERROR3[3])
  
  y1 = c(x1[1]+g*zt1[1],x1[2]+g*zt1[2],x1[3]+g*zt1[3])
  y2 = c(x2[1]+g*zt2[1],x2[2]+g*zt2[2],x2[3]+g*zt2[3])
  y3 = c(x3[1]+g*zt3[1],x3[2]+g*zt3[2],x3[3]+g*zt3[3])
  
  error1=c(ERROR[1],ERROR2[1],ERROR3[1])
  error2=c(ERROR[2],ERROR2[2],ERROR3[2])
  error3=c(ERROR[3],ERROR2[3],ERROR3[3])
  rango = seq(700,from = 3)
  i=3
  
  while (ind==0) {
    for(i in seq(3,(n-1))  ){
      RuidoBlanco=rmvnorm(1,mean = media,sigma = sigma)
      
      x1[i+1] = 0.4*x1[i] + 0.4*x2[i] + 0.5*x3[i] + 0.2*x1[i-1] - 0.2*x2[i-1] - 0.2*x1[i-2] + 0.15*x2[i-2] + 0.1*x3[i-2] + RuidoBlanco[1]
      x2[i+1] = 0.6*x2[i] + 0.2*x2[i-1] + 0.2*x2[i-2] + RuidoBlanco[2]
      x3[i+1] = 0.4*x3[i] + 0.3*x3[i-1] + 0.3*x3[i-2] + RuidoBlanco[3]
      
      Error = rmvnorm(1,mean = media,sigma = sigma)
      error1[i+1] = Error[1]
      error2[i+1] = Error[2]
      error3[i+1] = Error[3]
      
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


# SIMULACION SISTEMA 1
Database <- Sistema_10(500)
grafica_3var(Database,"10")
Wavelet_3(Database)
TE3_boot(Database)
TE3_sbt(Database)

