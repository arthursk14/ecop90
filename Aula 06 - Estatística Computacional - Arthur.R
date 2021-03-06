## Bootstrap param�trico ##

# Distribui��o normal #

# Tamanho da amostra
n = 20
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
for(i in 1:nr) xb[i] = mean(rnorm(n,xbo,de))

xbo
mean(xb)

# Variancia de xb simulado (bootstrap)
var(xb)

# Variancia te�rica de xb
vxb = (vx/n)
# Desvio padr�o te�rico de xb
dxb = (vx/n)^.5

# Plot da distribui��o te�rica de xb (PDF)
gra = seq(9,11,length.out = 10e3)
den = dnorm(gra,mean=10,sd=dxb)
plot(gra,den,type = "l")

# Plot da distribui��o te�rica de xb com a m�dia em xbo (PDF)
gra = seq(9,11,length.out = 10e3)
dent = dnorm(gra,mean=xbo,sd=dxb)
lines(gra,dent,col = "blue")

# Plot da densidade da distribui��o de xb simulada pelo bootstrap
lines(density(xb), col="red")
