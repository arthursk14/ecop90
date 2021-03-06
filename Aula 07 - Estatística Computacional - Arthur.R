## Bootstrap n�o-param�trico ##

# Distribui��o normal #

# Tamanho da amostra
n = 40
# M�dia da distribui��o
me = 10
# Vari�ncia da distribui��o
vx = 1
# Desvio padr�o da distribui��o
de = vx^.5
# Amostra original
x = rnorm(n,me,de) 

# Boostrap do estimador da m�dia amostral #

# N�mero de repeti��es 
nr = 10e3
# Cria o vetor para xb
xb = rep(0,nr)
# M�dia amostral
xbo = mean(x)
# Loop para o bootstrap
for(i in 1:nr) xb[i] = mean(sample(x,n,replace=T))

xbo
mean(xb)

# Variancia de xb simulado (bootstrap)
var(xb)

# Variancia te�rica de xb
vxb = (vx/n)
# Desvio padr�o te�rico de xb
dxb = (vx/n)^.5

# Plot da distribui��o te�rica de xb (PDF)
gra = seq(9,11,length = 10e3)
den = dnorm(gra,mean=10,sd=dxb)
plot(gra,den,type = "l", ylim=c(0,4))

# Plot da distribui��o te�rica de xb com a m�dia em xbo (PDF)
gra = seq(9,11,length = 10e3)
dent = dnorm(gra,mean=xbo,sd=dxb)
lines(gra,dent,col = "blue")

# Plot da densidade da distribui��o de xb simulada pelo bootstrap
lines(density(xb), col="red")

# Intervalos de confian�a #

# N�vel de confian�a
ga = 0.95
al = 1-ga
al2 = al/2

# Limite inferior (usando a suposi��o de normalidade)
lit = xbo+qnorm(al2)*de/n^.5
# Limite superior (usando a suposi��o de normalidade)
lst = xbo-qnorm(al2)*de/n^.5

lit
lst
xbo

# Limites inferior e superior (bootstrap n�o-param�trico)
lb = quantile(xb,c(al2,1-al2))

# Teste de hip�steses para bootstrap n�o param�trico #
# Devemos impor H0: me = 10 na reamostragem #

# Manipulando a amostra original (popula��o) para impor H0
xh0 = x - xbo + me
mean(xh0)

# Criando a nova distribui��o de xb (com m�dia me) via bootstrap
for (i in 1:nr) xb0[i] = mean(sample(xh0,n,replace=T))

# Pontos cr�ticos do teste de hip�teses
quantile(xb0,c(al2,1-al2))
