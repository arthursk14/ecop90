### Aprendizado de máquinas - Regressões Ridge, Lasso e Elastic Net ###

# Instalação do pacote "glmnet" #
install.packages("glmnet")
require("glmnet")

# Criação dos dados, modelo linear #

# Tamanho da amostra
n = 500
# Número de variáveis
q = 20

# Criando a matriz X, com as variáveis independentes
X = NULL

for (j in 1:q){
  X = cbind(X,rnorm(n,0,1))
}

# Criando o vetor de erros
e = rnorm(n,0,1)

# Coeficientes do modelo linear #

# Beta 0
B0 = 1
# Betas diferentes de 0
B = c(1,2,3,4)
# Vetor Beta para o modelo linear
B = c(B, rep(0,q-length(B)))

# Modelo linear (calcula o y, varíavel independente)
y = B0 + X%*%B + e

# OLS
fit.OLS = lm(formula=y~X)
B.hat.OLS = fit.OLS$coefficients 

# Ridge regression

lambda = 1
fit.RR = glmnet(y=y, x=X, lambda=lambda, alpha=0,)

# Coeficientes do RR
B.hat.RR = fit.RR$beta
B.hat.RR

# Intercepto do RR
B0.hat.RR = fit.RR$a0
B0.hat.RR

#LASSO
fit.LASSO = glmnet(y=y, x=X, lambda=lambda, alpha=1)

#Coeficientes do LASSO
B.hat.LASSO = fit.LASSO$beta
B.hat.LASSO

#Intercepto do LASSO
B0.hat.LASSO = fit.LASSO$a0
B0.hat.LASSO

#Elastic Net
fit.EN = glmnet(y=y, x=X, lambda=lambda, alpha=0.5)

#Coeficientes do Elastic Net
B.hat.EN = fit.EN$beta
B.hat.EN

#Intercepto do Elastic Net
B0.hat.EN = fit.EN$a0
B0.hat.EN

w = (abs(B.hat.LASSO)+(n)^(-1/2))^(-1)

fit.adaLASSO = glmnet(y=y, x=X, lambda=lambda ,alpha=1, penalty.factor=w)

#Coeficientes do adaLASSO
B.hat.adaLASSO = fit.adaLASSO$beta
B.hat.adaLASSO

#Intercepto do adaLASSO
B0.hat.adaLASSO = fit.adaLASSO$a0
B0.hat.adaLASSO

#Matriz com os coeficientes estimados
BETAS = cbind(c(B0,B),
              B.hat.OLS,
              c(B0.hat.RR,as.vector(B.hat.RR)),
              c(B0.hat.LASSO,as.vector(B.hat.LASSO)),
              c(B0.hat.adaLASSO,as.vector(B.hat.adaLASSO)),
              c(B0.hat.EN,as.vector(B.hat.EN))
              )

colnames(BETAS) = c("Betas","OLS","RR","LASSO","adaLASSO","EN")

#Exibe todos os valores de beta (real e estimados)
BETAS

#LASSO com varios valores de lambdas

lambdas = seq(from=0.05, to=0.5, by=0.05)

fit.LASSO2 = glmnet(y=y, x=X, lambda=lambda, alpha=1)

#Coeficientes do LASSO2
B.hat.LASSO2 = fit.LASSO2$beta
B.hat.LASSO2

#Intercepto do adaLASSO
B0.hat.LASSO2 = fit.LASSO2$a0
B0.hat.LASSO2

