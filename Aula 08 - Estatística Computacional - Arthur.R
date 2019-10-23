### Bootstrap para modelos de regressão ###

# Set seed
set.seed(999)

# Tamanho da amostra
n = 100
# Beta 0
b0 = 10
# Beta 1
b1 = 1
# Desvio-padrão do erro
de = 1

# Vetor x
x = seq(0,10,length=n)
# Vetor de erros
e = rnorm(n,0,de)
# Vetor y
y = b0+b1*x+e

# Gráfico do modelo (x,y)
plot(x,y)

# Cria matriz X
xm = cbind(rep(1,n),x)

# Estimador de MQO da amostra original
beo = solve(t(xm)%*%xm)%*%t(xm)%*%y

## Bootstrap método 1 ##

# Amostrando pares #

# Número de repetições
nr = 10e3

# Vetor com índices i = 1,2,...,n
gri = seq(1,n)

# Vetor de Beta 1 chapéu
b1c = rep(0,nr)

# vetores para y bootstrap e x bootstrap
yb = xb = rep(0,n)

# Loop para o bootstrap
for (i in 1:nr){
  # Loop para a reamostragem (sample para um par de variáveis)
  k = sample(gri,replace=T) 
  for (j in 1:n){
    xb[j] = x[k[j]]
    yb[j] = y[k[j]]
    
  }
  # Roda o modelo linear com a função do R
  m1 = lm(yb~xb)
  # Pega a estimativa de Beta 1 (Beta 1 chapéu)
  b1c[i] = m1$coef[2]
}

hist(b1c,prob=T)
beo[2]
mean(b1c)

# Fazer IC e Teste de Hipóteses, além de comparar com o resultado teórico 
# supondo normalidade, i.e. b1c ~ N(b1,sigma^2/soma(x_i-xb)^2)

## Bootstrap método 2 ##

# Amostrando resíduos via modelo #

# Regressão para obter os resíduos (poderia ter calculado por fora, com o vetor beo)
m2 = lm(y~x)
# Resíduos 
r= residuals(m2)
# Resíduos centrados
rp = r-mean(r)

# Beta chapéu original
b0o = m2$coefficients[1]
b1o = m2$coefficients[2]

# Vetor de Beta 1 chapéu
b1c2 = rep(0,nr)

# Loop para o bootstrap
for (i in 1:nr){
  # Reamostragem dos resíduos centrados
  reb = sample(rp,n,replace=T)
  # Encontrar o y bootstrap
  yb = b0o+b1o*x+reb
  # Regressão de y bootstrap e x original, para estimar beta chapéu bootstrap
  mb = lm(yb~x)
  # Pega a estimativa de Beta 1 (Beta 1 chapéu)
  b1c2[i] = mb$coef[2]
}

hist(b1c2,prob=T)
beo[2]
mean(b1c2)