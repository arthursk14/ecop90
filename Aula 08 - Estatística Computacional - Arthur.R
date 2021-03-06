### Bootstrap para modelos de regress�o ###

# Set seed
set.seed(999)

# Tamanho da amostra
n = 100
# Beta 0
b0 = 10
# Beta 1
b1 = 1
# Desvio-padr�o do erro
de = 1

# Vetor x
x = seq(0,10,length=n)
# Vetor de erros
e = rnorm(n,0,de)
# Vetor y
y = b0+b1*x+e

# Gr�fico do modelo (x,y)
plot(x,y)

# Cria matriz X
xm = cbind(rep(1,n),x)

# Estimador de MQO da amostra original
beo = solve(t(xm)%*%xm)%*%t(xm)%*%y

## Bootstrap m�todo 1 ##

# Amostrando pares #

# N�mero de repeti��es
nr = 10e3

# Vetor com �ndices i = 1,2,...,n
gri = seq(1,n)

# Vetor de Beta 1 chap�u
b1c = rep(0,nr)

# vetores para y bootstrap e x bootstrap
yb = xb = rep(0,n)

# Loop para o bootstrap
for (i in 1:nr){
  # Loop para a reamostragem (sample para um par de vari�veis)
  k = sample(gri,replace=T) 
  for (j in 1:n){
    xb[j] = x[k[j]]
    yb[j] = y[k[j]]
    
  }
  # Roda o modelo linear com a fun��o do R
  m1 = lm(yb~xb)
  # Pega a estimativa de Beta 1 (Beta 1 chap�u)
  b1c[i] = m1$coef[2]
}

hist(b1c,prob=T)
beo[2]
mean(b1c)

# Fazer IC e Teste de Hip�teses, al�m de comparar com o resultado te�rico 
# supondo normalidade, i.e. b1c ~ N(b1,sigma^2/soma(x_i-xb)^2)

## Bootstrap m�todo 2 ##

# Amostrando res�duos via modelo #

# Regress�o para obter os res�duos (poderia ter calculado por fora, com o vetor beo)
m2 = lm(y~x)
# Res�duos 
r= residuals(m2)
# Res�duos centrados
rp = r-mean(r)

# Beta chap�u original
b0o = m2$coefficients[1]
b1o = m2$coefficients[2]

# Vetor de Beta 1 chap�u
b1c2 = rep(0,nr)

# Loop para o bootstrap
for (i in 1:nr){
  # Reamostragem dos res�duos centrados
  reb = sample(rp,n,replace=T)
  # Encontrar o y bootstrap
  yb = b0o+b1o*x+reb
  # Regress�o de y bootstrap e x original, para estimar beta chap�u bootstrap
  mb = lm(yb~x)
  # Pega a estimativa de Beta 1 (Beta 1 chap�u)
  b1c2[i] = mb$coef[2]
}

hist(b1c2,prob=T)
beo[2]
mean(b1c2)