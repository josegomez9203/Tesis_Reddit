#install.packages('vars')
#install.packages("biwavelet")
#install.packages('RTransferEntropy')
#install.packages('mvtnorm')
# tempdir()
# # [1] "C:\Users\XYZ~1\AppData\Local\Temp\Rtmp86bEoJ\Rtxt32dcef24de2"
# dir.create(tempdir())

library(biwavelet)
library(RTransferEntropy)
library(vars)
library(mvtnorm)
library(lmtest)
############### 1 ###########
x1 = c(0)
x2 = c(0)
x3 = c(0)
media <- c(0,0,0)
sigma<-diag(3)


RuidoBlanco=rmvnorm(1,mean = media,sigma = sigma)

rango = seq(1000)
for( i in rango ){
  RuidoBlanco=rmvnorm(1,mean = media,sigma = sigma)
  x1_1 = 3.4*x1[i]*(1-x1[i])^2*exp(-(x1[i]^2))+0.4*RuidoBlanco[1]
  
  x2_1 = 3.4*x2[i]*(1-x2[i])^2*exp(-(x2[i]^2))+ 0.5*x1[i]*x2[i] +0.4*RuidoBlanco[2]
  
  x3_1 = 3.4*x3[i]*(1-x3[i])^2*exp(-(x3[i]^2))+ 0.3*x2[i] + 0.5*(x1[i]^2) + 0.4*RuidoBlanco[3]
  
  x1 = c(x1,x1_1)
  
  x2 = c(x2,x2_1)
  
  x3 = c(x3,x3_1)
  
}

plot(seq(length(x1)), x1,type="l")
lines(x2,col="red")
lines(x3,col="darkblue")
dev.off()

### Wavelet Coherence
# Se trata como serie
A= cbind(rep(1:length(x1)),x1)
B= cbind(rep(1:length(x2)),x2)
C= cbind(rep(1:length(x3)),x3)

nrands = 200 # Numero de aleatorizaciones montecarlo
wtc.AB = wtc(A, B, nrands = nrands)
wtc.AC = wtc(A, C, nrands = nrands)
wtc.BC = wtc(B, C, nrands = nrands)


# Graficamos
par(mfrow=c(1,1),oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)

plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 1, 
     lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
     plot.cb = TRUE, main = "Wavelet Coherence: X1 vs X2")
n = length(A[, 1])
abline(v = seq(100, n, 100), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)

plot(wtc.AC, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
     plot.cb = TRUE, main = "Wavelet Coherence: X1 vs X3")
abline(v = seq(100, n, 100), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)

plot(wtc.BC, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
     plot.cb = TRUE, main = "Wavelet Coherence: X2 vs X3")
abline(v = seq(100, n, 100), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)



#### Transferencia de entropia
transfer_entropy(x1, x2,q = 1, entropy = c('Shannon'), shuffles = 200, type = c('quantiles'),
                quantiles = c(5, 95), nboot = 300, burn = 60, seed = NULL)

transfer_entropy(x1, x3,q = 1, entropy = c('Shannon'), shuffles = 200, type = c('quantiles'),
                 quantiles = c(5, 95), nboot = 300, burn = 60, seed = NULL)

transfer_entropy(x2, x3,q = 1, entropy = c('Shannon'), shuffles = 200, type = c('quantiles'),
                 quantiles = c(5, 95), nboot = 300, burn = 60, seed = NULL)



#grangertest(x1,x2)
varfit <- VAR(cbind(x3,x1),p=1,type = 'const')
sum_var =summary(varfit)
sum_var$varresult$x3

#Calculo de transferencia de entropía sin usar bootstrap
calc_ete(x3,x1) # Transferencia de entropía efectiva
calc_te(x3,x1)


############### 2 ###########

x1 = c(0,runif(1))
x2 = c(0,runif(1))
x3 = c(0,runif(1))


rango = seq(1000,from = 2)
for( i in rango ){
  RuidoBlanco=rmvnorm(1,mean = media,sigma = sigma)
  x1_1 = 0.7*x1[i] + RuidoBlanco[1]
  
  x2_1 = 0.3*x2[i] + 0.5*x1[i]*x2[i-1] + RuidoBlanco[2]
  
  x3_1 = 0.3*x3[i] + 0.5*x1[i]*(x3[i-1]) + RuidoBlanco[3]
  
  x1 = c(x1,x1_1)
  
  x2 = c(x2,x2_1)
  
  x3 = c(x3,x3_1)
  
}

plot(seq(length(x1)), x1,type="l")
lines(x2,col="red")
lines(x3,col="darkblue")
dev.off()




### Wavelet Coherence
# Se trata como serie
A= cbind(rep(1:length(x1)),x1)
B= cbind(rep(1:length(x2)),x2)
C= cbind(rep(1:length(x3)),x3)

nrands = 200 # Numero de aleatorizaciones montecarlo
wtc.AB = wtc(A, B, nrands = nrands)
wtc.AC = wtc(A, C, nrands = nrands)
wtc.BC = wtc(B, C, nrands = nrands)


# Graficamos
par(mfrow=c(1,1),oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)

plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 1, 
     lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
     plot.cb = TRUE, main = "Wavelet Coherence: X1 vs X2")
n = length(A[, 1])
abline(v = seq(100, n, 100), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)

plot(wtc.AC, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
     plot.cb = TRUE, main = "Wavelet Coherence: X1 vs X3")
abline(v = seq(100, n, 100), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)

plot(wtc.BC, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
     plot.cb = TRUE, main = "Wavelet Coherence: X2 vs X3")
abline(v = seq(100, n, 100), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)



#### Transferencia de entropia
transfer_entropy(x1, x2,q = 1, entropy = c('Shannon'), shuffles = 200, type = c('quantiles'),
                 quantiles = c(5, 95), nboot = 300, burn = 60, seed = NULL)

transfer_entropy(x1, x3,q = 1, entropy = c('Shannon'), shuffles = 200, type = c('quantiles'),
                 quantiles = c(5, 95), nboot = 300, burn = 60, seed = NULL)

transfer_entropy(x2, x3,q = 1, entropy = c('Shannon'), shuffles = 200, type = c('quantiles'),
                 quantiles = c(5, 95), nboot = 300, burn = 60, seed = NULL)



############### 3 ###########


x1 = c(0,runif(1))
x2 = c(0,runif(1))
x3 = c(0,runif(1))
cc = c(0,0.05,0.3,0.5,0.7,0.8,0.9)
ccc = sample(cc,1)
rango = seq(1000,from = 2)

for( i in rango ){
  
  x1_1 = 1.4 - (x1[i]^2) + 0.3*x1[i-1]
  
  x2_1 = 1.4 - ccc*x1[i]*x2[i] - (1-ccc)*(x2[i]^2) + 0.3*x2[i-1]
  
  x3_1 = 1.4 - ccc*x2[i]*x3[i] - (1-ccc)*(x3[i]^2) + 0.3*x3[i-1]
  
  x1 = c(x1,x1_1)
  
  x2 = c(x2,x2_1)
  
  x3 = c(x3,x3_1)
  
}

plot(seq(length(x1)), x1,type="l")
lines(x2,col="red")
lines(x3,col="darkblue")
dev.off()

x1<-x1[-1]
x2<-x2[-1]
x3<-x3[-1]

### Wavelet Coherence
# Se trata como serie
A= cbind(rep(1:length(x1)),x1)
B= cbind(rep(1:length(x2)),x2)
C= cbind(rep(1:length(x3)),x3)

nrands = 200 # Numero de aleatorizaciones montecarlo
wtc.AB = wtc(A, B, nrands = nrands)
wtc.AC = wtc(A, C, nrands = nrands)
wtc.BC = wtc(B, C, nrands = nrands)


# Graficamos
par(mfrow=c(1,1),oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)

plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 1, 
     lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
     plot.cb = TRUE, main = "Wavelet Coherence: X1 vs X2")
n = length(A[, 1])
abline(v = seq(100, n, 100), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)

plot(wtc.AC, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
     plot.cb = TRUE, main = "Wavelet Coherence: X1 vs X3")
abline(v = seq(100, n, 100), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)

plot(wtc.BC, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
     plot.cb = TRUE, main = "Wavelet Coherence: X2 vs X3")
abline(v = seq(100, n, 100), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)



#### Transferencia de entropia
transfer_entropy(x1, x2,q = 1, entropy = c('Shannon'), shuffles = 200, type = c('quantiles'),
                 quantiles = c(5, 95), nboot = 300, burn = 60, seed = NULL)

transfer_entropy(x1, x3,q = 1, entropy = c('Shannon'), shuffles = 200, type = c('quantiles'),
                 quantiles = c(5, 95), nboot = 300, burn = 60, seed = NULL)

transfer_entropy(x2, x3,q = 1, entropy = c('Shannon'), shuffles = 200, type = c('quantiles'),
                 quantiles = c(5, 95), nboot = 300, burn = 60, seed = NULL)



############### 4 ###########


x1 = c(0,runif(1))
x2 = c(0,runif(1))
x3 = c(0,runif(1))
x4 = c(0,runif(1))
cc = c(0,0.2,0.4)
ccc = sample(cc,1)
rango = seq(1000,from = 2)

if(ccc==0){print("Uncoupled case")}else{if(ccc==0.2){print("Weak coupling")}else{print("Strong coupling")}}

for(i in rango ){
  
  x1_1 = 1.4 - (x1[i]^2) + 0.3*x1[i-1]
  
  x4_1 = 1.4 - (x4[i]^2) + 0.3*x4[i-1]
  
  x2_1 = 1.4 - (0.5*ccc*(x1[i]+x3[i]) + (1-ccc)*x2[i])^2 + 0.3*x2[i-1]
  
  x3_1 = 1.4 - (0.5*ccc*(x2[i]+x4[i]) + (1-ccc)*x3[i])^2 + 0.3*x3[i-1]
  
  x1 = c(x1,x1_1)
  
  x2 = c(x2,x2_1)
  
  x3 = c(x3,x3_1)
  
  x4 = c(x4,x4_1)
  
}

plot(seq(length(x1)), x1,type="l")
lines(x2,col="red")
lines(x3,col="darkblue")
lines(x4,col="darkgreen")
dev.off()



x1<-x1[-1]
x2<-x2[-1]
x3<-x3[-1]
x4<-x4[-1]

### Wavelet Coherence
# Se trata como serie
A= cbind(rep(1:length(x1)),x1)
B= cbind(rep(1:length(x2)),x2)
C= cbind(rep(1:length(x3)),x3)
D= cbind(rep(1:length(x4)),x4)

nrands = 200 # Numero de aleatorizaciones montecarlo
wtc.AB = wtc(A, B, nrands = nrands)
wtc.AC = wtc(A, C, nrands = nrands)
wtc.AD = wtc(A, D, nrands = nrands)
wtc.BC = wtc(B, C, nrands = nrands)
wtc.BD = wtc(B, D, nrands = nrands)
wtc.CD = wtc(C, D, nrands = nrands)

# Graficamos
par(mfrow=c(1,1),oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)

plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 1, 
     lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
     plot.cb = TRUE, main = "Wavelet Coherence: X1 vs X2")
n = length(A[, 1])
abline(v = seq(100, n, 100), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)

plot(wtc.AC, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 1, 
     lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
     plot.cb = TRUE, main = "Wavelet Coherence: X1 vs X3")
abline(v = seq(100, n, 100), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)

plot(wtc.AD, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 1, 
     lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
     plot.cb = TRUE, main = "Wavelet Coherence: X1 vs X4")
abline(v = seq(100, n, 100), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)

plot(wtc.BC, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
     plot.cb = TRUE, main = "Wavelet Coherence: X2 vs X3")
abline(v = seq(100, n, 100), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)

plot(wtc.BD, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
     plot.cb = TRUE, main = "Wavelet Coherence: X2 vs X4")
abline(v = seq(100, n, 100), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)

plot(wtc.CD, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
     plot.cb = TRUE, main = "Wavelet Coherence: X3 vs X4")
abline(v = seq(100, n, 100), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)


#### Transferencia de entropia
transfer_entropy(x1, x2,q = 1, entropy = c('Shannon'), shuffles = 200, type = c('quantiles'),
                 quantiles = c(5, 95), nboot = 300, burn = 60, seed = NULL)

transfer_entropy(x1, x3,q = 1, entropy = c('Shannon'), shuffles = 200, type = c('quantiles'),
                 quantiles = c(5, 95), nboot = 300, burn = 60, seed = NULL)

transfer_entropy(x1, x4,q = 1, entropy = c('Shannon'), shuffles = 200, type = c('quantiles'),
                 quantiles = c(5, 95), nboot = 300, burn = 60, seed = NULL)

transfer_entropy(x2, x3,q = 1, entropy = c('Shannon'), shuffles = 200, type = c('quantiles'),
                 quantiles = c(5, 95), nboot = 300, burn = 60, seed = NULL)

transfer_entropy(x2, x4,q = 1, entropy = c('Shannon'), shuffles = 200, type = c('quantiles'),
                 quantiles = c(5, 95), nboot = 300, burn = 60, seed = NULL)

transfer_entropy(x3, x4,q = 1, entropy = c('Shannon'), shuffles = 200, type = c('quantiles'),
                 quantiles = c(5, 95), nboot = 300, burn = 60, seed = NULL)

############### 5 ########### 


x1 = c(0,runif(1))
x2 = c(0,runif(1))
x3 = c(0,runif(1))
cc = c(0,0.05,0.3,0.5,0.7,0.8,0.9)
ccc = sample(cc,1)
ccc = 0.9
n = 1000
rango = seq(n,from = 2)

for( i in rango ){
  
  x1_1 = 1.4 - (x1[i]^2) + 0.3*x1[i-1]
  
  x2_1 = 1.4 - ccc*x1[i]*x2[i] - (1-ccc)*(x2[i]^2) + 0.3*x2[i-1] 
  
  x3_1 = 1.4 - ccc*x2[i]*x3[i] - (1-ccc)*(x3[i]^2) + 0.3*x3[i-1] 
  
  x1 = c(x1,x1_1)
  
  x2 = c(x2,x2_1)
  
  x3 = c(x3,x3_1)
  
}
outlier = sample(1000,size = round(n*0.1),replace = F)
for( i in outlier ){
  
  x1[i] = x1[i] + runif(1)
  
  x2[i] = x2[i] + runif(1)
  
  x3[i] = x3[i] + runif(1)
  
}


plot(seq(length(x1)), x1,type="l")
lines(x2,col="red")
lines(x3,col="darkblue")
dev.off()


x1<-x1[-1]
x2<-x2[-1]
x3<-x3[-1]

### Wavelet Coherence
# Se trata como serie
A= cbind(rep(1:length(x1)),x1)
B= cbind(rep(1:length(x2)),x2)
C= cbind(rep(1:length(x3)),x3)

nrands = 200 # Numero de aleatorizaciones montecarlo
wtc.AB = wtc(A, B, nrands = nrands)
wtc.AC = wtc(A, C, nrands = nrands)
wtc.BC = wtc(B, C, nrands = nrands)


# Graficamos
par(mfrow=c(1,1),oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)

plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 1, 
     lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
     plot.cb = TRUE, main = "Wavelet Coherence: X1 vs X2")
n = length(A[, 1])
abline(v = seq(100, n, 100), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)

plot(wtc.AC, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
     plot.cb = TRUE, main = "Wavelet Coherence: X1 vs X3")
abline(v = seq(100, n, 100), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)

plot(wtc.BC, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.15, ylab = "Escala", xlab = "Periodo", 
     plot.cb = TRUE, main = "Wavelet Coherence: X2 vs X3")
abline(v = seq(100, n, 100), h = 1:16, col = "darkgrey", lty = 1, lwd = 1)


#### Transferencia de entropia
transfer_entropy(x1, x2,q = 1, entropy = c('Shannon'), shuffles = 200, type = c('quantiles'),
                 quantiles = c(5, 95), nboot = 300, burn = 60, seed = NULL)

transfer_entropy(x1, x3,q = 1, entropy = c('Shannon'), shuffles = 200, type = c('quantiles'),
                 quantiles = c(5, 95), nboot = 300, burn = 60, seed = NULL)

transfer_entropy(x2, x3,q = 1, entropy = c('Shannon'), shuffles = 200, type = c('quantiles'),
                 quantiles = c(5, 95), nboot = 300, burn = 60, seed = NULL)



############### 6 ###########


x1 = c(0,runif(1))
x2 = c(0,runif(1))
x3 = c(0,runif(1))
cc = c(0,0.05,0.3,0.5,0.7,0.8,0.9)
ccc = sample(cc,1)
n_t = c(0,rnorm(1))
rango = seq(100,from = 2)
i=17
for( i in rango ){
  n_t_1 = n_t[i] + rnorm(1)
  
  x1_1 = 1.4 - (x1[i]^2) + 0.3*x1[i-1] + n_t_1
  
  x2_1 = 1.4 - ccc*x1[i]*x2[i] - (1-ccc)*(x2[i]^2) + 0.3*x2[i-1] + n_t_1
  
  x3_1 = 1.4 - ccc*x2[i]*x3[i] - (1-ccc)*(x3[i]^2) + 0.3*x3[i-1] + n_t_1
  
  x1 = c(x1,x1_1)
  
  x2 = c(x2,x2_1)
  
  x3 = c(x3,x3_1)
  
  n_t = c(n_t,n_t_1)
  
}


plot(seq(length(x1)), x1,type="l")
lines(x2,col="red")
lines(x3,col="darkblue")
dev.off()


############### 7 ###########


x1 = c(0,runif(1))
x2 = c(0,runif(1))
x3 = c(0,runif(1))
cc = c(0,0.05,0.3,0.5,0.7,0.8,0.9)
ccc = sample(cc,1)
n_t = c(0,0)
rango = seq(100,from = 2)

for(i in rango ){
  a = rnorm(1,0.01,0.02)
  
  n_t_1 = a*(i+1)
  
  x1_1 = 1.4 - (x1[i]^2) + 0.3*x1[i-1] + n_t_1
  
  x2_1 = 1.4 - ccc*x1[i]*x2[i] - (1-ccc)*(x2[i]^2) + 0.3*x2[i-1] + n_t_1
  
  x3_1 = 1.4 - ccc*x2[i]*x3[i] - (1-ccc)*(x3[i]^2) + 0.3*x3[i-1] + n_t_1
  
  x1 = c(x1,x1_1)
  
  x2 = c(x2,x2_1)
  
  x3 = c(x3,x3_1)
  
  n_t = c(n_t,n_t_1)
  
}

plot(seq(length(x1)), x1,type="l")
lines(x2,col="red")
lines(x3,col="darkblue")
dev.off()

############### 8 ###########
Error = rmvnorm(1,mean = media,sigma = sigma)
Error2 = rmvnorm(1,mean = media,sigma = sigma)

x1 = c(0,runif(1))
x2 = c(0,runif(1))
x3 = c(0,runif(1))

a0 = 0.2
a1 = 0.9
b1 = 0.1

sigma2_1 = c(a0,a0+(a1*Error2[1]^2)+b1*a0)
sigma2_2 = c(a0,a0+(a1*Error2[2]^2)+b1*a0)
sigma2_3 = c(a0,a0+(a1*Error2[3]^2)+b1*a0)

zt1 = c(sigma2_1[1]*Error[1],sigma2_1[2]*Error2[1])
zt2 = c(sigma2_2[1]*Error[2],sigma2_2[2]*Error2[2])
zt3 = c(sigma2_3[1]*Error[3],sigma2_3[2]*Error2[3])

g= 0.5
y1 = c(x1[1]+g*zt1[1],x1[2]+g*zt1[2])
y2 = c(x2[1]+g*zt2[1],x2[2]+g*zt2[2])
y3 = c(x3[1]+g*zt3[1],x3[2]+g*zt3[2])

rango = seq(100,from = 2)
for( i in rango ){
  RuidoBlanco=rmvnorm(1,media,sigma)
  x1_1 = 0.7*x1[i] + RuidoBlanco[1]
  
  x2_1 = 0.3*x2[i] + 0.5*x1[i]*x2[i-1] + RuidoBlanco[2]
  
  x3_1 = 0.3*x3[i] + 0.5*x1[i]*(x3[i-1]) + RuidoBlanco[3]
  
  x1 = c(x1,x1_1)
  
  x2 = c(x2,x2_1)
  
  x3 = c(x3,x3_1)
  
  Error = rmvnorm(1,mean = media,sigma = sigma)
  error1 = c(error1,Error[1])
  error2 = c(error2,Error[2])
  error3 = c(error3,Error[3])
  
  sigma21 = a0 + a1*(error1[i]^2) + b1*(sigma2_1[i]^2)
  sigma22 = a0 + a1*(error2[i]^2) + b1*(sigma2_2[i]^2)
  sigma23 = a0 + a1*(error3[i]^2) + b1*(sigma2_3[i]^2)
  
  sigma2_1 = c(sigma2_1,sigma21)
  sigma2_2 = c(sigma2_2,sigma22)
  sigma2_3 = c(sigma2_3,sigma23)
  
  zt1_1 = sigma21*error1[i+1]
  zt2_1 = sigma22*error2[i+1]
  zt3_1 = sigma23*error3[i+1]
  
  zt1 = c(zt1,zt1_1)
  zt2 = c(zt2,zt2_1)
  zt3 = c(zt3,zt3_1)
  
  y1 = c(y1, x1_1+g*zt1_1)
  y2 = c(y2, x2_1+g*zt2_1)
  y3 = c(y3, x3_1+g*zt3_1)
}

y1



plot(seq(length(y1)), y1,type="l")
lines(y2,col="red")
lines(y3,col="darkblue")
dev.off()



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

############### 9 ###########


############### 10 ###########

x1 = runif(3)
x2 = runif(3)
x3 = runif(3)
ERROR = rmvnorm(1,mean = media,sigma = sigma)
ERROR2 = rmvnorm(1,mean = media,sigma = sigma)
ERROR3 = rmvnorm(1,mean = media,sigma = sigma)

a0 = 0.2
a1 = 0.9
b1 = 0.1

sigma2_1 = c(a0,a0+(a1*ERROR[1]^2)+b1*a0, a0+(a1*ERROR2[1]^2)+b1*(a0+(a1*ERROR[1]^2)+b1*a0))
sigma2_2 = c(a0,a0+(a1*ERROR[2]^2)+b1*a0, a0+(a1*ERROR2[2]^2)+b1*(a0+(a1*ERROR[2]^2)+b1*a0))
sigma2_3 = c(a0,a0+(a1*ERROR[3]^2)+b1*a0, a0+(a1*ERROR2[3]^2)+b1*(a0+(a1*ERROR[3]^2)+b1*a0))

zt1 = c(sigma2_1[1]*ERROR[1],sigma2_1[2]*ERROR2[1],sigma2_1[3]*ERROR3[1])
zt2 = c(sigma2_2[1]*ERROR[2],sigma2_2[2]*ERROR2[2],sigma2_2[3]*ERROR3[2])
zt3 = c(sigma2_3[1]*ERROR[3],sigma2_3[2]*ERROR2[3],sigma2_3[3]*ERROR3[3])

g= 0.2
y1 = c(x1[1]+g*zt1[1],x1[2]+g*zt1[2],x1[3]+g*zt1[3])
y2 = c(x2[1]+g*zt2[1],x2[2]+g*zt2[2],x2[3]+g*zt2[3])
y3 = c(x3[1]+g*zt3[1],x3[2]+g*zt3[2],x3[3]+g*zt3[3])


rango = seq(700,from = 3)
i=3
for(i in rango ){
  RuidoBlanco=rmvnorm(1,mean = media,sigma = sigma)
  
  x1_1 = 0.4*x1[i] + 0.4*x2[i] + 0.5*x3[i] + 0.2*x1[i-1] - 0.2*x2[i-1] - 0.2*x1[i-2] + 0.15*x2[i-2] + 0.1*x3[i-2] + RuidoBlanco[1]
  
  x2_1 = 0.6*x2[i] + 0.2*x2[i-1] + 0.2*x2[i-2] + RuidoBlanco[2]
  
  x3_1 = 0.4*x3[i] + 0.3*x3[i-1] + 0.3*x3[i-2] + RuidoBlanco[3]
  
  x1 = c(x1,x1_1)
  
  x2 = c(x2,x2_1)
  
  x3 = c(x3,x3_1)
  
  Error = rmvnorm(1,mean = media,sigma = sigma)
  error1 = c(error1,Error[1])
  error2 = c(error2,Error[2])
  error3 = c(error3,Error[3])
  
  sigma21 = a0 + a1*(error1[i]^2) + b1*(sigma2_1[i]^2)
  sigma22 = a0 + a1*(error2[i]^2) + b1*(sigma2_2[i]^2)
  sigma23 = a0 + a1*(error3[i]^2) + b1*(sigma2_3[i]^2)
  
  sigma2_1 = c(sigma2_1,sigma21)
  sigma2_2 = c(sigma2_2,sigma22)
  sigma2_3 = c(sigma2_3,sigma23)
  
  zt1_1 = sigma21*error1[i+1]
  zt2_1 = sigma22*error2[i+1]
  zt3_1 = sigma23*error3[i+1]
  
  zt1 = c(zt1,zt1_1)
  zt2 = c(zt2,zt2_1)
  zt3 = c(zt3,zt3_1)
  
  y1 = c(y1, x1_1+g*zt1_1)
  y2 = c(y2, x2_1+g*zt2_1)
  y3 = c(y3, x3_1+g*zt3_1)

}

plot(seq(length(y1)), x1,type="l")
lines(y2,col="red")
lines(y3,col="green")
dev.off()


