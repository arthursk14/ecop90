## Bootstrap paramétrico ##

# Distribuição normal #

# Tamanho da amostra
n = 20
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
for(i in 1:nr) xb[i] = mean(rnorm(n,xbo,de))

xbo
mean(xb)

# Variancia de xb simulado (bootstrap)
var(xb)

# Variancia teórica de xb
vxb = (vx/n)
# Desvio padrão teórico de xb
dxb = (vx/n)^.5

# Plot da distribuição teórica de xb (PDF)
gra = seq(9,11,length.out = 10e3)
den = dnorm(gra,mean=10,sd=dxb)
plot(gra,den,type = "l")

# Plot da distribuição teórica de xb com a média em xbo (PDF)
gra = seq(9,11,length.out = 10e3)
dent = dnorm(gra,mean=xbo,sd=dxb)
lines(gra,dent,col = "blue")

# Plot da densidade da distribuição de xb simulada pelo bootstrap
lines(density(xb), col="red")
