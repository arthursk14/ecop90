#### Exercício 1 ###

# Lei dos grandes números #

nmax = 15e3
la = 0.5

xb = rep(0,nmax)
x = rexp(nmax,la)

for (i in 1:nmax){
  xb[i] = mean(x[1:i])
  
}

ts.plot(xb,xlab="Observações", ylab="Média")
linha = rep(1/la,nmax)
lines(linha,col="blue")

# Trocando para distribuição gamma

nmax = 15e3
be = 1/2
al = 2

xb = rep(0,nmax)
x = rgamma(nmax,al,be)

for (i in 1:nmax){
  xb[i] = mean(x[1:i])
  
}

ts.plot(xb,xlab="Observações", ylab="Média")
linha = rep(al/be,nmax)
lines(linha,col="blue")

# Teorema central do limite #

n = c(1,5,20,100,1000)
nn = length(n)
nr =10e3
la = 0.5

xb = matrix(rep(0,nn*nr),ncol=nn)

for (i in 1:nn){
  for (j in 1:nr){
    x = rexp(n[i],la)
    xb[j,i] = mean(x)
    
  }

}

windows()
par(mfrow=c(3,2))

for (i in 1:nn){
  d = density(xb[,i])
  hist(xb[,i],prob=T,main=paste("n =",n[i]),ylim=c(0,max(d$y)))
  lines(d,col="blue")
  
}

ks.test(xb[,1],"pnorm",mean(xb[,1]),sd(xb[,1]))
ks.test(xb[,1],"pexp",1/mean(xb[,1]))

ks.test(xb[,5],"pnorm",mean(xb[,5]),sd(xb[,5]))
ks.test(xb[,5],"pexp",1/mean(xb[,5]))

# Trocando para distribuição gamma

n = c(1,5,20,100,1000)
nn = length(n)
nr =10e3

be = 1/2
al = 2

xb = matrix(rep(0,nn*nr),ncol=nn)

for (i in 1:nn){
  for (j in 1:nr){
    x = rgamma(n[i],al,be)
    xb[j,i] = mean(x)
    
  }
  
}

windows()
par(mfrow=c(3,2))

for (i in 1:nn){
  d = density(xb[,i])
  hist(xb[,i],prob=T,main=paste("n =",n[i]),ylim=c(0,max(d$y)))
  lines(d,col="blue")
  
}

bec = mean(xb[,1])/(sd(xb[,1])**2)
alc = (mean(xb[,1])**2)/(sd(xb[,1])**2)

ks.test(xb[,1],"pnorm",mean(xb[,1]),sd(xb[,1]))
ks.test(xb[,1],"pgamma",alc,bec)

bec = mean(xb[,5])/(sd(xb[,5])**2)
alc = (mean(xb[,5])**2)/(sd(xb[,5])**2)

ks.test(xb[,5],"pnorm",mean(xb[,5]),sd(xb[,5]))
ks.test(xb[,5],"pgamma",alc,bec)

mean(xb[,1])
mean(xb[,5])

sd(xb[,1])
sd(xb[,5])

# Trocando para distribuição uniforme

n = c(1,5,20,100,1000)
nn = length(n)
nr =10e3

min = 10
max = 30

xb = matrix(rep(0,nn*nr),ncol=nn)

for (i in 1:nn){
  for (j in 1:nr){
    x = runif(n[i],min,max)
    xb[j,i] = mean(x)
    
  }
  
}

windows()
par(mfrow=c(3,2))

for (i in 1:nn){
  d = density(xb[,i])
  hist(xb[,i],prob=T,main=paste("n =",n[i]),ylim=c(0,max(d$y)))
  lines(d,col="blue")
  
}

ks.test(xb[,1],"pnorm",mean(xb[,1]),sd(xb[,1]))
ks.test(xb[,5],"pnorm",mean(xb[,5]),sd(xb[,5]))

mean(xb[,1])
mean(xb[,5])

sd(xb[,1])
sd(xb[,5])

### Exercício 2 ###

# Verificar o coeficiente de confiança na construção de Intervalos de Confiança #

nr = 10e3
nv = c(10,100)
nvt = length(nv)
mi = 0
de = 1
ga = 0.90
al2 = (1-ga)/2

ta = rep(0,nvt)

for (n in 1:nvt) {
  for (i in 1:nr){
    x = rnorm(nv[n],mi,de)
    xb = mean(x)
    sn = sd(x)
    
    t = qt(1-al2,nv[n]-1)
    li = xb-t*(sn/(nv[n]^.5))
    ls = xb+t*(sn/(nv[n]^.5))
    
    if (mi>=li & mi<=ls) {
      ta[n]=ta[n]+1
    
    }
    
  }
  
}

print(ta)

# Exercício 3 #

# Teste de hipótese para inferir a média populacional a partir de uma amostra #

# Quantas vezes rejeitamos H0 mesmo ela sendo verdadeira (erro do tipo I)
# Se p-valor < (1-ga), rejeitamos H0
# Ou seja, o resultado esperado é que rejeitemos H0 (1-ga)*nr vezes

nr = 10e3
nv = c(10,20)
nvt = length(nv)
mi = 0
de = 1
ga = 0.90
al2 = (1-ga)/2

ta = rep(0,nvt)

# Hipótese nula -> mi = mic
mic = 0

for (n in 1:nvt) {
  for (i in 1:nr){
    x = rnorm(nv[n],mi,de)
    xb = mean(x)
    sn = sd(x)
    
    # Estatística do teste T
    tobs = (xb-mic)/(sn/(nv[n]^.5))
    t = qt(1-al2,nv[n]-1)

    # Quantas vezes rejeitamos a hipótese nula
    if (tobs>=t | tobs<=(-t)) {
      ta[n]=ta[n]+1
      
    }
    
  }
  
}

print(ta)

# Quantas vezes aceitamos H0 mesmo ela sendo falsa (erro do tipo II)
# Aqui calculamos o poder do teste, que é (1-probabilidade de ocorrência de erro do tipo II), dado mic

nr = 10e3
nv = c(10,20)
nvt = length(nv)
mi = 0
de = 1
ga = 0.90
al2 = (1-ga)/2

ta = rep(0,nvt)

# Hipótese nula -> mi = mic
mic = 1.5

for (n in 1:nvt) {
  for (i in 1:nr){
    x = rnorm(nv[n],mi,de)
    xb = mean(x)
    sn = sd(x)
    
    # Estatística do teste T
    tobs = (xb-mic)/(sn/(nv[n]^.5))
    t = qt(1-al2,nv[n]-1)
    
    # Quantas vezes aceitamos a hipótese nula
    if (tobs<=t & tobs>=(-t)) {
      ta[n]=ta[n]+1
      
    }
    
  }
  
}

print(nr-ta)