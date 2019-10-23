## Exercício opcional ##

# Gibbs Sampling #
# Inferência bayesiana (aplicação popular = DSGE) #


# Seja X~Bernoulli(teta)
# X1, ..., Xn quaisquer (amostra)
# Cada Xi é 0 ou 1

# Item (i): Encontre o estimador de máxima verossimilhança de teta
# Item (ii): Considerando, a priori, que teta ~ Beta(alfa,beta), 
# encontre a distribuição a posteriori de teta

# Fazer simualações deste problema #

# (i) Gerar Bernoulli com vários tetas diferentes
# (ii) Variar alfa e beta nas prioris
# (iii) Fazer gráficos com: priori, verossimilhança e posteriori

## Exercício ##

# Suponha que Xi ~ Poisson(lambda) para i = 1, 2, ..., m-1
# e Xi ~ Poisson(phi) para i = m, m+1, ..., n
# em que Xi, ..., Xm-1 e Xm, ..., Xn representam a.a's.

# O interesse é em estimar os parâmetros desconhecidos (lambda, phi e m)

# Resolução do ponto de vista bayesiano
# Prioris: 
# lambda ~ Gama(alfa,beta)
# phi ~ Gama(gama,ro)
# m ~ Uniforme{1,2,...,n}
# lambda, phi e m são independentes

## Gibbs Sampling para mudança de regime poisson ##

# Gerando a amostra #

# Tamanho da amostra
n = 200

# Ponto de quebra
m = 110

# Parâmetro 1 de Poisson (lambda)
la = 0.1

# Parâmetro 2 de Poisson (phi)
fi = 1.0

# Definido o gerador de números aleatórios
set.seed(988)

x1 = rpois((m-1),la)
x2 = rpois((n-m+1),fi)
x = c(x1,x2)

ts.plot(x)

# Estimando #

# Número de repetições (simulações Gibbs)
nr = 21e3

# Parâmetros da distribuição de lambda ~ Gama(alfa,beta)
al = 0.1
be = 0.1

# Parâmetros da distribuições de phi ~ Gama(gama,ro)
ga = 0.1
ro = 0.05

# Vetores com os parâmetros da poisson a serem estimados
lae = fie = rep(1,nr)

# Vetor com m estimado (assumimos que o salto só pode acontecer a partir de 2, 
# ou seja, não admitimos a hipótese de que não existe troca de regime)
me = rep(2,nr)
pm = rep(0,n)
vm = seq(1,n)

# Somatório dos elementos de x
sox = rep(0,n)
for (i in 1:n) sox[i] = sum(x[1:i])

# Loop
for (i in 2:nr){
  # Gerando lae[i] e fie[i], cuidando para o caso m[i-1]=1
  if(me[i-1]>1){
    p1la = sox[(me[i-1]-1)] + al
    p2la = be + me[i-1] - 1
    
    p1fi = sox[n] - sox[(me[i-1]-1)] + ga
    p2fi = ro + n - me[i-1] + 1
    
  }
  else {
    p1la = al
    p2la = be
    
    p1fi = sox[n] + ga
    p2fi = ro + n
    
  }
  
  lae[i] = rgamma(1,p1la,p2la)
  fie[i] = rgamma(1,p1fi,p2fi)
  
  # Gerando me[i]
  pm[1] = exp(-1*(lae[i]-fie[i]))*fie[i]^sox[n]
  
  for(k in 2:n){
    pm[k] = lae[i]^sox[k-1]*exp(-k*(lae[i]-fie[i]))*fie[i]^(sox[n]-sox[k-1])
    
  }
  
  soma = sum(pm)
  pm = pm/soma
  me[i] = sample(vm,1,prob=pm)
 
}

bi = 1000
de = 5
nf = (nr-bi)/de


laef = fief = mef = rep(0,nf)

for(i in 1:nf){
  laef[i] = lae[(bi+1)+(i-1)*de]
  fief[i] = fie[(bi+1)+(i-1)*de]
  mef[i] = me[(bi+1)+(i-1)*de]
}

ts.plot(laef)
hist(laef,prob=T)
ts.plot(fief)
hist(fief,prob=T)
ts.plot(mef)
hist(mef,prob=T)

acf(me)
acf(mef)
acf(lae)
acf(laef)
acf(fie)
acf(fief)
mean(laef)
median(laef)
mean(fief)
median(fief)
mean(mef)
median(mef)

