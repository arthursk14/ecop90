### Exercício 3 ###
### Yi = Bo + Bi*Xi + E ###

# i)  E ~ N(0,v) (com e sem outliers)
# ii) E ~ t(v) (com e sem outliers)

# nr = 10.000
# n = 20, 100, 500

#Bmq (Mínimos quadrados)
#Bma (Mínimos absolutos)

# Comparar os dois estimadores

# Contaminar os dados com outliers:
# y[50] = 5*y[50]
# y[80] = 10*y[80]

# Comentário: Como o enunciado também pedia para calcular para n = 20, 100 e 500,
# coloquei os outliers da seguinte maneira:
# y[n/2] = 5*y[n/2]
# y[n*3/4] = 10*y[n*3/4]

### Começo do programa ###

# Número de repetições
nr = 10e3

# Tamanho da amostra
n = c(20,100,500)
nn = length(n)

# Beta 0
b0 = 2

# Beta 1
b1 = 1

# Desvio padrão da distribuição normal dos erros
de = 1

# Vetores para os estimadores de b1
b1mq = matrix(rep(0,nn*nr),ncol=nn)
b1ma = matrix(rep(0,nn*nr),ncol=nn)

## Função para calcular os estimadores da regressão de mínimos desvios absolutos 
mabs <- function(x,y){
  
  b = c(1,1)
  
  sda <- function(b,x,y){
    sum(abs(y-b[1]-b[2]*x))
  }
  
  res = optim(b,sda,x=x,y=y)
  res$par
}

# Loop para fazer a estimação
for (i in 1:nn){
  for (j in 1:nr){
    
    # Erros com distribuição normal
    e = rnorm(n[i],0,de)
    
    # Vetor x (variáveis independentes da regressão)
    x = seq(0,10, length=n[i])
    
    # Vetor com a variável dependente
    y = b0 + b1*x + e
    
    # Mínimos quadrados ordinários
    mqo = lm(y~x)
    
    # Mínimos desvios absolutos
    mda = mabs(x,y)
    
    # Coloca os estimadores de b1 no vetor criado anteriormente
    b1mq[j,i] = mqo$coeff[2]
    b1ma[j,i] = mda[2]
  }
}

# Erro da estimação
eb1mq = b1mq - b1
eb1ma = b1ma - b1

# Boxplot do erro
nomes = c("A20 - MQO", "A20 - MDA", "A100 - MQO", "A100 - MDA", "A500 - MQO", "A500 - MDA")
boxplot(eb1mq[,1],eb1ma[,1],eb1mq[,2],eb1ma[,2],eb1mq[,3],eb1ma[,3],outline = F, names = nomes)

## Contaminando os dados com outliers ##

# Vetores para os estimadores de b1 com outliers
ob1mq = matrix(rep(0,nn*nr),ncol=nn)
ob1ma = matrix(rep(0,nn*nr),ncol=nn)

# Loop para fazer a estimação
for (i in 1:nn){
  for (j in 1:nr){
    
    # Erros com distribuição normal
    e = rnorm(n[i],0,de)
    
    # Vetor x (variáveis independentes da regressão)
    x = seq(0,10, length=n[i])
    
    # Vetor com a variável dependente
    y = b0 + b1*x + e
    
    # Colocando os outliers
    y[n[i]/2] = 5*y[n[i]/2]
    y[n[i]*3/4] = 10*y[n[i]*3/4]
    
    # Mínimos quadrados ordinários
    mqo = lm(y~x)
    
    # Mínimos desvios absolutos
    mda = mabs(x,y)
    
    # Coloca os estimadores de b1 no vetor criado anteriormente
    ob1mq[j,i] = mqo$coeff[2]
    ob1ma[j,i] = mda[2]
  }
}

# Erro da estimação com outliers
eob1mq = ob1mq - b1
eob1ma = ob1ma - b1

# Boxplot do erro
nomes = c("A20 - MQO", "A20 - MDA", "A100 - MQO", "A100 - MDA", "A500 - MQO", "A500 - MDA")
boxplot(eob1mq[,1],eob1ma[,1],eob1mq[,2],eob1ma[,2],eob1mq[,3],eob1ma[,3],outline = F, names = nomes)

## Fazendo o erro com distribuição t de student ##

# Graus de liberdade para a distribuição
v = 4

# Vetores para os estimadores de b1
tb1mq = matrix(rep(0,nn*nr),ncol=nn)
tb1ma = matrix(rep(0,nn*nr),ncol=nn)

# Loop para fazer a estimação
for (i in 1:nn){
  for (j in 1:nr){
    
    # Erros com distribuição t de student
    e = rt(n[i],v)
    
    # Vetor x (variáveis independentes da regressão)
    x = seq(0,10, length=n[i])
    
    # Vetor com a variável dependente
    y = b0 + b1*x + e
    
    # Mínimos quadrados ordinários
    mqo = lm(y~x)
    
    # Mínimos desvios absolutos
    mda = mabs(x,y)
    
    # Coloca os estimadores de b1 no vetor criado anteriormente
    tb1mq[j,i] = mqo$coeff[2]
    tb1ma[j,i] = mda[2]
  }
}

# Erro da estimação
etb1mq = tb1mq - b1
etb1ma = tb1ma - b1

# Boxplot do erro
nomes = c("A20 - MQO", "A20 - MDA", "A100 - MQO", "A100 - MDA", "A500 - MQO", "A500 - MDA")
boxplot(etb1mq[,1],etb1ma[,1],etb1mq[,2],etb1ma[,2],etb1mq[,3],etb1ma[,3],outline = F, names = nomes)

## Contaminando os dados com outliers ##

# Vetores para os estimadores de b1 com outliers
otb1mq = matrix(rep(0,nn*nr),ncol=nn)
otb1ma = matrix(rep(0,nn*nr),ncol=nn)

# Loop para fazer a estimação
for (i in 1:nn){
  for (j in 1:nr){
    
    # Erros com distribuição t de student
    e = rt(n[i],v)
    
    # Vetor x (variáveis independentes da regressão)
    x = seq(0,10, length=n[i])
    
    # Vetor com a variável dependente
    y = b0 + b1*x + e
    
    # Colocando os outliers
    y[n[i]/2] = 5*y[n[i]/2]
    y[n[i]*3/4] = 10*y[n[i]*3/4]
    
    # Mínimos quadrados ordinários
    mqo = lm(y~x)
    
    # Mínimos desvios absolutos
    mda = mabs(x,y)
    
    # Coloca os estimadores de b1 no vetor criado anteriormente
    otb1mq[j,i] = mqo$coeff[2]
    otb1ma[j,i] = mda[2]
  }
}

# Erro da estimação com outliers
eotb1mq = otb1mq - b1
eotb1ma = otb1ma - b1

# Boxplot do erro
nomes = c("A20 - MQO", "A20 - MDA", "A100 - MQO", "A100 - MDA", "A500 - MQO", "A500 - MDA")
boxplot(eotb1mq[,1],eotb1ma[,1],eotb1mq[,2],eotb1ma[,2],eotb1mq[,3],eotb1ma[,3],outline = F, names = nomes)

# Enviar por e-mail "trab1comp" "flavioaz@gmail.com"
