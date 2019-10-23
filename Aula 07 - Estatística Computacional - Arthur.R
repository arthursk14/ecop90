## Bootstrap não-paramétrico ##

# Distribuição normal #

# Tamanho da amostra
n = 40
# Média da distribuição
me = 10
# Variância da distribuição
vx = 1
# Desvio padrão da distribuição
de = vx^.5
# Amostra original
x = rnorm(n,me,de) 

# Boostrap do estimador da média amostral #

# Número de repetições 
nr = 10e3
# Cria o vetor para xb
xb = rep(0,nr)
# Média amostral
xbo = mean(x)
# Loop para o bootstrap
for(i in 1:nr) xb[i] = mean(sample(x,n,replace=T))

xbo
mean(xb)

# Variancia de xb simulado (bootstrap)
var(xb)

# Variancia teórica de xb
vxb = (vx/n)
# Desvio padrão teórico de xb
dxb = (vx/n)^.5

# Plot da distribuição teórica de xb (PDF)
gra = seq(9,11,length = 10e3)
den = dnorm(gra,mean=10,sd=dxb)
plot(gra,den,type = "l", ylim=c(0,4))

# Plot da distribuição teórica de xb com a média em xbo (PDF)
gra = seq(9,11,length = 10e3)
dent = dnorm(gra,mean=xbo,sd=dxb)
lines(gra,dent,col = "blue")

# Plot da densidade da distribuição de xb simulada pelo bootstrap
lines(density(xb), col="red")

# Intervalos de confiança #

# Nível de confiança
ga = 0.95
al = 1-ga
al2 = al/2

# Limite inferior (usando a suposição de normalidade)
lit = xbo+qnorm(al2)*de/n^.5
# Limite superior (usando a suposição de normalidade)
lst = xbo-qnorm(al2)*de/n^.5

lit
lst
xbo

# Limites inferior e superior (bootstrap não-paramétrico)
lb = quantile(xb,c(al2,1-al2))

# Teste de hipósteses para bootstrap não paramétrico #
# Devemos impor H0: me = 10 na reamostragem #

# Manipulando a amostra original (população) para impor H0
xh0 = x - xbo + me
mean(xh0)

# Criando a nova distribuição de xb (com média me) via bootstrap
for (i in 1:nr) xb0[i] = mean(sample(xh0,n,replace=T))

# Pontos críticos do teste de hipóteses
quantile(xb0,c(al2,1-al2))
