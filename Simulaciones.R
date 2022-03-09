#install.packages("biwavelet")
library(biwavelet)

A= cbind(rep(1:length(x1)),x1)
B= cbind(rep(1:length(x2)),x2)
nrands = 10
wtc.AB = wtc(A, B, nrands = nrands)

# Plotting a graph
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
plot(wtc.AB, col.coi = "grey", ylab = "Scale", xlab = "Period",  main = "Wavelet Coherence: A vs B")

# Adding grid lines
n = length(t1[, 1])
abline(v = seq(260, n, 260), h = 1:16, col = "brown", lty = 1, lwd = 1)

# Defining x labels
axis(side = 3, at = c(seq(0, n, 260)), labels = c(seq(1999, 2015, 1)))
############### 1 ###########
x1 = c(0)
x2 = c(0)
x3 = c(0)

RuidoBlanco=rnorm(3,0,1)
i=2
rango = seq(10)
for( i in rango ){
  RuidoBlanco=rnorm(3,0,1)
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
############### 2 ###########

x1 = c(0,1)
x2 = c(0,2)
x3 = c(0,3)

rango = seq(10,from = 2)
for( i in rango ){
  RuidoBlanco=rnorm(3,0,1)
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

############### 3 ###########


x1 = c(0,1)
x2 = c(0,2)
x3 = c(0,3)
cc = c(0,0.05,0.3,0.5,0.7,0.8,0.9)
c = sample(cc,1)
rango = seq(10,from = 2)

for( i in rango ){
  
  x1_1 = 1.4 - (x1[i]^2) + 0.3*x1[i-1]
  
  x2_1 = 1.4 - c*x1[i]*x2[i] - (1-c)*(x2[i]^2) + 0.3*x2[i-1]
  
  x3_1 = 1.4 - c*x2[i]*x3[i] - (1-c)*(x3[i]^2) + 0.3*x3[i-1]
  
  x1 = c(x1,x1_1)
  
  x2 = c(x2,x2_1)
  
  x3 = c(x3,x3_1)
  
}

plot(seq(length(x1)), x1,type="l")
lines(x2,col="red")
lines(x3,col="darkblue")
dev.off()



############### 4 ###########


x1 = c(0,1)
x2 = c(0,0.98)
x3 = c(0,1.3)
x4 = c(0,1.5)
cc = c(0,0.2,0.4)
c = sample(cc,1)
rango = seq(10,from = 2)
i=2
if(c==0){print("Uncoupled case")}else{if(c==0.2){print("Weak coupling")}else{print("Strong coupling")}}

for( i in rango ){
  
  x1_1 = 1.4 - (x1[i]^2) + 0.3*x1[i-1]
  
  x4_1 = 1.4 - (x4[i]^2) + 0.3*x4[i-1]
  
  x2_1 = 1.4 - (0.5*c*(x1[i]+x3[i]) + (1-c)*x2[i])^2 + 0.3*x2[i-1]
  
  x3_1 = 1.4 - (0.5*c*(x2[i]+x4[i]) + (1-c)*x3[i])^2 + 0.3*x3[i-1]
  
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



############### 5 ###########


x1 = c(0,1)
x2 = c(0,2)
x3 = c(0,3)
cc = c(0,0.05,0.3,0.5,0.7,0.8,0.9)
c = sample(cc,1)
rango = seq(10,from = 2)

for( i in rango ){
  
  x1_1 = 1.4 - (x1[i]^2) + 0.3*x1[i-1] + sample(c(0,1),size = 1,replace = T,prob = c(0.99,0.01))*runif(1)
  
  x2_1 = 1.4 - c*x1[i]*x2[i] - (1-c)*(x2[i]^2) + 0.3*x2[i-1] + sample(c(0,1),size = 1,replace = T,prob = c(0.99,0.01))*runif(1)
  
  x3_1 = 1.4 - c*x2[i]*x3[i] - (1-c)*(x3[i]^2) + 0.3*x3[i-1] + sample(c(0,1),size = 1,replace = T,prob = c(0.99,0.01))*runif(1)
  
  x1 = c(x1,x1_1)
  
  x2 = c(x2,x2_1)
  
  x3 = c(x3,x3_1)
  
}

plot(seq(length(x1)), x1,type="l")
lines(x2,col="red")
lines(x3,col="darkblue")
dev.off()



############### 6 ###########


x1 = c(0,1)
x2 = c(0,0.2)
x3 = c(0,0.8)
cc = c(0,0.05,0.3,0.5,0.7,0.8,0.9)
c = sample(cc,1)
n_t = c(0,0)
rango = seq(10,from = 2)

for( i in rango ){
  n_t_1 = n_t[i] + rnorm(1)
  
  x1_1 = 1.4 - (x1[i]^2) + 0.3*x1[i-1] + n_t_1
  
  x2_1 = 1.4 - c*x1[i]*x2[i] - (1-c)*(x2[i]^2) + 0.3*x2[i-1] + n_t_1
  
  x3_1 = 1.4 - c*x2[i]*x3[i] - (1-c)*(x3[i]^2) + 0.3*x3[i-1] + n_t_1
  
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


x1 = c(0,1)
x2 = c(0,0.4)
x3 = c(0,0.1)
cc = c(0,0.05,0.3,0.5,0.7,0.8,0.9)
c = sample(cc,1)
n_t = c(0,0)
rango = seq(10,from = 2)

for(i in rango ){
  a = rnorm(1,0.01,0.02)
  
  n_t_1 = a*(i+1)
  
  x1_1 = 1.4 - (x1[i]^2) + 0.3*x1[i-1] + n_t_1
  
  x2_1 = 1.4 - c*x1[i]*x2[i] - (1-c)*(x2[i]^2) + 0.3*x2[i-1] + n_t_1
  
  x3_1 = 1.4 - c*x2[i]*x3[i] - (1-c)*(x3[i]^2) + 0.3*x3[i-1] + n_t_1
  
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

error1 = rnorm(2)
error2 = rnorm(2)
error3 = rnorm(2)

x1 = c(0,1)
x2 = c(0,2)
x3 = c(0,3)

a0 = 0.2
a1 = 0.9
b1 = 0.1

sigma2_1 = c(a0,a0+(a1*error1[1]^2)+b1*a0)
sigma2_2 = c(a0,a0+(a1*error2[1]^2)+b1*a0)
sigma2_3 = c(a0,a0+(a1*error3[1]^2)+b1*a0)

zt1 = c(sigma2[1]*error1[1],sigma2[2]*error1[2])
zt2 = c(sigma2[1]*error2[1],sigma2[2]*error2[2])
zt3 = c(sigma2[1]*error3[1],sigma2[2]*error3[2])

g= 0.5
y1 = c(x1[1]+g*zt1[1],x1[2]+g*zt1[2])
y2 = c(x2[1]+g*zt2[1],x2[2]+g*zt2[2])
y3 = c(x3[1]+g*zt3[1],x3[2]+g*zt3[2])

rango = seq(10,from = 2)
for( i in rango ){
  RuidoBlanco=rnorm(3,0,1)
  x1_1 = 0.7*x1[i] + RuidoBlanco[1]
  
  x2_1 = 0.3*x2[i] + 0.5*x1[i]*x2[i-1] + RuidoBlanco[2]
  
  x3_1 = 0.3*x3[i] + 0.5*x1[i]*(x3[i-1]) + RuidoBlanco[3]
  
  x1 = c(x1,x1_1)
  
  x2 = c(x2,x2_1)
  
  x3 = c(x3,x3_1)
  
  
  error1 = c(error1,rnorm(1))
  error2 = c(error2,rnorm(1))
  error3 = c(error3,rnorm(1))
  
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





plot(seq(length(x1)), x1,type="l")
lines(x2,col="red")
lines(x3,col="darkblue")
dev.off()



############### 9 ###########


############### 10 ###########


x1 = c(0,1)
x2 = c(0,0.2)
x3 = c(0,0.8)
cc = c(0,0.05,0.3,0.5,0.7,0.8,0.9)
c = sample(cc,1)
n_t = c(0,0)
rango = seq(10,from = 2)

for( i in rango ){
  n_t_1 = n_t[i] + rnorm(1)
  
  x1_1 = 1.4 - (x1[i]^2) + 0.3*x1[i-1] + n_t_1
  
  x2_1 = 1.4 - c*x1[i]*x2[i] - (1-c)*(x2[i]^2) + 0.3*x2[i-1] + n_t_1
  
  x3_1 = 1.4 - c*x2[i]*x3[i] - (1-c)*(x3[i]^2) + 0.3*x3[i-1] + n_t_1
  
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


x1 = runif(3)
x2 = runif(3)
x3 = runif(3)

rango = seq(10,from = 3)

for(i in rango ){
  
  x1_1 = 0.4*x1[i] + 0.4*x2[i] + 0.5*x3[i] + 0.2*x1[i-1] - 0.2*x2[i-1] - 0.2*x1[i-2] + 0.15*x2[i-2] + 0.1*x3[i-2] + rnorm(1)
  
  x2_1 = 0.6*x2[i] + 0.2*x2[i-1] + 0.2*x2[i-2] + rnorm(1)
  
  x3_1 = 0.4*x3[i] + 0.3*x3[i-1] + 0.3*x3[i-2] + rnorm(1)
  
  x1 = c(x1,x1_1)
  
  x2 = c(x2,x2_1)
  
  x3 = c(x3,x3_1)

}

plot(seq(length(x1)), x1,type="l")
lines(x2,col="red")
lines(x3,col="green")
dev.off()


