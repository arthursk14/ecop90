### Exercício ###

## Simulação de valores de uma normal trivariada ##

# Usando a decomposição (f(x,y,z)=f(x|y,z)f(y|z)f(z)) #

# Tamanho da amostra
n = 5e3
# Média de x
mx = 10
# Média de y
my = 20
# Média de z
mz = 30
# Desvio-padrão de x
dx = 2
# Desvio-padrão de y
dy = 3
# Desvio-padrão de z
dz = 4
# Covariância de x e y
cxy = 2
# Covariância de x e z
cxz = 3
# Covariância de y e z
cyz = 4

# Coeficientes de correlação 
ro_xy = cxy/(dx*dy)
ro_xz = cxz/(dx*dz)
ro_yz = cyz/(dy*dz)

# Matriz de covariância
sigma = matrix(c(dx^2,cxy,cxz,cxy,dy^2,cyz,cxz,cyz,dz^2),c(3,3))

# Particionando sigma para f(y|Z)
cy_yy = sigma[2,2]
cy_zz = sigma[3,3]
cy_yz = sigma[2,3]
cy_zy = sigma[3,2]

# Particionando sigma para f(x|y,Z)
cx_xx = sigma[1,1]
cx_yz = sigma[2:3,2:3]
cx_xyxz = sigma[1,2:3]
cx_yxyz = sigma[2:3,1]

# Variâncias f(y|z) e f(x|y,z)
dy_z = cy_yy-cy_yz%*%solve(cy_zz)%*%cy_zy
dx_yz = cx_xx-cx_xyxz%*%solve(cx_yz)%*%cx_yxyz

x = matrix(rep(0,n*3),ncol=3)

for (i in 1:n){
  # Distribuição marginal de z
  x[i,3] = rnorm(1,mz,dz)
  # Esperança condicional de y dado z
  my_z = my+cy_yz%*%solve(cy_zz)%*%(x[i,3]-mz)
  # f(y|z)
  x[i,2] = rnorm(1,my_z,dy_z^.5)
  # Esperaça condicional de x dados y e z
  mx_yz = mx+cx_xyxz%*%solve(cx_yz)%*%matrix(c(x[i,2]-my,x[i,3]-mz),c(2,1))
  # f(x|y,z)
  x[i,1] = rnorm(1,mx_yz,dx_yz^.5)
  
}

# Histogramas de x e y
par(mfrow=c(1,2))

for (i in 1:2){
  if (i==1) {
    hist(x[,i],probability=T,main="f(x|y,z)",xlab="x",ylim = c(0,0.2))
    den = density(x[,i])
    lines(den$x,den$y,col="blue")
    lines(den$x,dnorm(den$x,mean(x[,i]),sd(x[,i])),col="red")
    
    }
  else if (i==2){
    hist(x[,i],probability=T,main="f(y|z)",xlab="y",ylim = c(0,0.2))
    den = density(x[,i])
    lines(den$x,den$y,col="blue")
    lines(den$x,dnorm(den$x,mean(x[,i]),sd(x[,i])),col="red")
    
    }
    
}

# Testando a hipótese de normalidade
ks.test(x[,1],"pnorm",mean(x[,1]),sd(x[,1]))
ks.test(x[,2],"pnorm",mean(x[,2]),sd(x[,2]))

par(mfrow=c(2,2))

plot(x[,1],x[,2],main="Dispersão xy",xlab="x",ylab="y")
abline(lm(x[,2]~x[,1]),col="red")
plot(x[,1],x[,3],main="Dispersão xz",xlab="x",ylab="z")
abline(lm(x[,3]~x[,1]),col="red")
plot(x[,2],x[,3],main="Dispersão yz",xlab="y",ylab="z")
abline(lm(x[,3]~x[,2]),col="red")

# Usando Gibbs Sampling #

# Tamanho da amostra
n = 30e3
# Tamanho final
nf = 5e3
# Amostra inicial a ser descartada
bi = 5e3
# Delta
de = 5

# Média de x
mx = 10
# Média de y
my = 20
# Média de z
mz = 30
# Desvio-padrão de x
dx = 2
# Desvio-padrão de y
dy = 3
# Desvio-padrão de z
dz = 4
# Covariância de x e y
cxy = 2
# Covariância de x e z
cxz = 3
# Covariância de y e z
cyz = 4

# Matriz de covariância
sigma = matrix(c(dx^2,cxy,cxz,cxy,dy^2,cyz,cxz,cyz,dz^2),c(3,3))

# Particionando sigma para f(x|y,Z)
cx_xx = sigma[1,1]
cx_yz = sigma[2:3,2:3]
cx_xyxz = sigma[1,2:3]
cx_yxyz = sigma[2:3,1]

# Particionando sigma para f(y|x,Z)
cy_yy = sigma[2,2]
cy_xz = sigma[c(1,3),c(1,3)]
cy_yxyz = sigma[1,c(1,3)]
cy_xyzy = sigma[c(1,3),1]

# Particionando sigma para f(z|x,y)
cz_zz = sigma[3,3]
cz_xy = sigma[1:2,1:2]
cz_zxzy = sigma[1,1:2]
cz_xzyz = sigma[1:2,1]

# Variâncias f(x|y,z), f(y|x,z) e f(z|x,y)
dx_yz = cx_xx-cx_xyxz%*%solve(cx_yz)%*%cx_yxyz
dy_xz = cy_yy-cy_yxyz%*%solve(cy_xz)%*%cy_xyzy
dz_xy = cz_zz-cz_zxzy%*%solve(cz_xy)%*%cz_xzyz

# Criando vetores para as distribuições
x=y=z=rep(0,n)

# Chutando valores iniciais
x[1]=10
y[1]=20
z[1]=30

# Loop para criar a distribuição
for (i in 2:n){
  mx_yz = mx+cx_xyxz%*%solve(cx_yz)%*%matrix(c(y[i-1]-my,z[i-1]-mz),c(2,1))
  # f(x|y,z)
  x[i] = rnorm(1,mx_yz,dx_yz^.5)
  my_xz = my+cy_yxyz%*%solve(cy_xz)%*%matrix(c(x[i]-mx,z[i-1]-mz),c(2,1))
  # f(y|x,z)
  y[i] = rnorm(1,my_xz,dy_xz^.5)
  mz_xy = mz+cz_zxzy%*%solve(cz_xy)%*%matrix(c(x[i]-mx,y[i]-my),c(2,1))
  # f(z|x,y)
  z[i] = rnorm(1,mz_xy,dz_xy^.5)
  
}

# Cria matriz com a distribuição trivariada
xyz = cbind(x,y,z)

# Funções de autocovariância e autocorrelação
par(mfrow=c(2,2))
acf(x)
acf(y)
acf(z)

# Plots das distribuições geradas
par(mfrow=c(2,2))
ts.plot(x)
ts.plot(y)
ts.plot(z)

# Tratamento para excluir as primeiras amostras e pegar valores espaçados,
# de forma a quebrar a possível estrutura de autocorrelação
xc = yc = zc = rep(0,nf)

for (i in 1:nf){
  xc[i] = x[(bi+1)+(i-1)*de]
  yc[i] = y[(bi+1)+(i-1)*de]
  zc[i] = z[(bi+1)+(i-1)*de]
  
}

# Funções de autocovariância e autocorrelação
par(mfrow=c(2,2))
acf(xc)
acf(yc)
acf(zc)

# Plots das distribuições geradas
par(mfrow=c(2,2))
ts.plot(xc)
ts.plot(yc)
ts.plot(zc)

# Histogramas
par(mfrow=c(2,2))

for (i in 1:3){
  if (i==1) {
    hist(xyz[,i],probability=T,main="f(x|y,z)",xlab="xc",ylim = c(0,0.20))
    den = density(xyz[,i])
    lines(den$x,den$y,col="blue")
    lines(den$x,dnorm(den$x,mean(xyz[,i]),sd(xyz[,i])),col="red")
    
  }
  else if (i==2){
    hist(xyz[,i],probability=T,main="f(y|x,z)",xlab="yc",ylim = c(0,0.15))
    den = density(xyz[,i])
    lines(den$x,den$y,col="blue")
    lines(den$x,dnorm(den$x,mean(xyz[,i]),sd(xyz[,i])),col="red")
    
  }
  else if (i==3){
    hist(xyz[,i],probability=T,main="f(z|x,y)",xlab="zc",ylim = c(0,0.12))
    den = density(xyz[,i])
    lines(den$x,den$y,col="blue")
    lines(den$x,dnorm(den$x,mean(xyz[,i]),sd(xyz[,i])),col="red")
    
  }
  
}

# Testando a hipótese de normalidade
ks.test(xyz[,1],"pnorm",mean(xyz[,1]),sd(xyz[,1]))
ks.test(xyz[,2],"pnorm",mean(xyz[,2]),sd(xyz[,2]))
ks.test(xyz[,3],"pnorm",mean(xyz[,3]),sd(xyz[,3]))

par(mfrow=c(2,2))

plot(xyz[,1],xyz[,2],main="Dispersão xy",xlab="x",ylab="y")
abline(lm(xyz[,2]~xyz[,1]),col="red")
plot(xyz[,1],xyz[,3],main="Dispersão xz",xlab="x",ylab="z")
abline(lm(xyz[,3]~xyz[,1]),col="red")
plot(xyz[,2],xyz[,3],main="Dispersão yz",xlab="y",ylab="z")
abline(lm(xyz[,3]~xyz[,2]),col="red")
