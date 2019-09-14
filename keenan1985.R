# Replicando o artigo "A Tukey Nonadditivity-Type Test for Time Series Nonlinearity"
# de "Daniel MacRae Keenan"
# publicado na revista "Biometrika, Vol. 72, No. 1 (Apr., 1985), pp. 39-44" 

# Cria função lag, que será utilizada na função do teste
flag = function(x,k){
  if (k>0){
    return (c(rep(NA,k),x)[1:length(x)]);
  }
  else{
    return (c(x[(-k+1):length(x)],rep(NA,-k)));
  }
}

# Cria a função do teste
keenan_test = function(y,m,...){
  
  if(missing(m)) m=ar(y)$order
  
  # Cria função que faz as regressões necessárias para o teste
  OLS = function(y,x,m){
    X = NULL
    for (i in 1:m) X = cbind(X,flag(x,i))
    X = cbind(y,X)
    X = na.omit(X)
    
    lm(X[,1]~X[,-1])
    
  }
  
  # Parâmetros
  y = as.vector(y)
  n = length(y)
  
  # Passo 1
  reg1 = OLS(y=y,x=y,m=m)
  
  y_hat = reg1$fitted.values
  y_hat_s = y_hat^2
  e_hat = reg1$residuals
  e_hat_s = e_hat^2
  SSR = sum(e_hat_s)
  
  # Passo 2
  reg2 = OLS(y=c(rep(NA,m),y_hat_s),x=y,m=m)
  
  xi_hat = reg2$residuals
  xi_hat_s = xi_hat^2
  
  # Passo 3
  reg3 = lm(e_hat~xi_hat-1)
  
  eta0 = reg3$coefficients
  
  # Encontrando eta
  names(eta0) = NULL
  eta = eta0*(sum(xi_hat_s)^.5)
  
  # Encontrando a estatística de teste F_hat(1,n-2m-2)
  F_hat = (eta^2*(n-2*m-2))/(SSR-eta^2)
  
  # Calculando o p-valor
  p_valor = pf(F_hat, 1, n-2*m-2,lower.tail=FALSE)
  
  # Para replicar a simulação do artigo, só nos interessa a estatística de teste F
  return(F_hat)
}

# Loop para replicar as simulações do autor

# Número de repetições
nr = 350

# Tamanho da amostra (n) e da ordem da aproximação autoregressiva (m)
n = c(70,204)
m = c(4,8)

# Parâmetros para os loops
nn = length(n)
mm = length(m)

# Monta as matrizes que guardam os dados
for (r in 1:nn){
  for (s in 1:mm){
    assign(paste("F_",n[r],"_",m[s],sep = ""), matrix(rep(0,nr*6),ncol = 6))
    
  }
  
}

# Loop
for (l in 1:nr){
  
  for (j in 1:nn){
    
    for (k in 1:mm){
      
    # Ruído branco
    e = rnorm(n[j])
    
    # Criação dos modelos (séries temporais)
    md1 = md2 = md3 = md4 = md5 = md6 = c(0,0)
    
    # Loop para montar os modelos (observações 3:n)
    for(i in 3:n[j]){
      md1[i] = e[i] - 0.4*e[i-1] + 0.3*e[i - 2] 
      md2[i] = e[i] - 0.4*e[i-1] + 0.3*e[i - 2] + 0.5*e[i]*e[i-2]
      md3[i] = e[i] - 0.3*e[i-1] + 0.2*e[i - 2] + 0.4*e[i-1] - 0.25*(e[i-1]^2) 
      md4[i] = e[i] + 0.4*md4[i-1] - 0.3*md4[i - 2] 
      md5[i] = e[i] + 0.4*md5[i-1] - 0.3*md5[i - 2] +0.5*md5[i-1]*e[i-1]
      md6[i] = e[i] + 0.4*md6[i-1] - 0.3*md6[i - 2] +0.5*md6[i-1]*e[i-1] + 0.8*e[i-1]
      
    }
    
    # Declara a matriz utilizada em cada loop
    matriz = paste("F_",n[j],"_",m[k],sep = "")
    
    # Coloca os valores do teste na respectiva matriz
    assign(matriz, `[<-`(get(matriz), cbind(l,1), keenan_test(md1,m[k])))
    assign(matriz, `[<-`(get(matriz), cbind(l,2), keenan_test(md2,m[k])))
    assign(matriz, `[<-`(get(matriz), cbind(l,3), keenan_test(md3,m[k])))
    assign(matriz, `[<-`(get(matriz), cbind(l,4), keenan_test(md4,m[k])))
    assign(matriz, `[<-`(get(matriz), cbind(l,5), keenan_test(md5,m[k])))
    assign(matriz, `[<-`(get(matriz), cbind(l,6), keenan_test(md6,m[k])))
    
    }
    
  }

}

# Cria a matriz final, para comparação com a tabela apresentada no artigo

quantis = c(0.5,0.75,0.9,0.95,0.99)

final = matrix(cbind(
  
              quantile(F_70_4[,3],probs=quantis),
              quantile(F_204_4[,3],probs=quantis),
              quantile(F_70_8[,3],probs=quantis),
              quantile(F_204_8[,3],probs=quantis),
              
              quantile(F_70_4[,2],probs=quantis),
              quantile(F_204_4[,2],probs=quantis),
              quantile(F_70_8[,2],probs=quantis),
              quantile(F_204_8[,2],probs=quantis),
              
              quantile(F_70_4[,1],probs=quantis),
              quantile(F_204_4[,1],probs=quantis),
              quantile(F_70_8[,1],probs=quantis),
              quantile(F_204_8[,1],probs=quantis),
              
              qchisq(quantis,1),
              qchisq(quantis,1),
              
              quantile(F_70_4[,4],probs=quantis),
              quantile(F_204_4[,4],probs=quantis),
              quantile(F_70_8[,4],probs=quantis),
              quantile(F_204_8[,4],probs=quantis),
              
              quantile(F_70_4[,5],probs=quantis),
              quantile(F_204_4[,5],probs=quantis),
              quantile(F_70_8[,5],probs=quantis),
              quantile(F_204_8[,5],probs=quantis),
              
              quantile(F_70_4[,6],probs=quantis),
              quantile(F_204_4[,6],probs=quantis),
              quantile(F_70_8[,6],probs=quantis),
              quantile(F_204_8[,6],probs=quantis)
            
              ), ncol = 10, byrow = T)

# Exporta a tabela final para Latex
require("xtable")
print(xtable(final, type = "latex"), file = "tabela.tex")
