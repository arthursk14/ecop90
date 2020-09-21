
# MCMC - Metropolis Hastings para estimar os parâmetros de uma regressão linear #

  # Seed
    
    r = sample(seq.int(from=1, to=20000, by=1),1)
    set.seed(19152)
  
  # Processo gerador de dados
  
    # Tamanho da amostra
      n = 10e1
    # Beta 0
      b0 = 5
    # Beta 1
      b1 = 2
    # Desvio-padrão do erro
      de = 1
    
    # Vetor x
      x = seq(0,10,length=n)
    # Vetor de erros
      e = rnorm(n,0,de)
    # Vetor y
      y = b0+b1*x+e
    
    # Gráfico dos dados gerados (x,y)
      plot(x,y)
      
    # Cria matriz X
      xm = cbind(rep(1,n),x)
      
  # Estimador de MQO da amostra original
      
    b_hat = solve(t(xm)%*%xm)%*%t(xm)%*%y
    b_hat_r = lm(y~x)
    
  # Plot do modelo com a reta da regressão
 
    plot(y~x)
    abline(b_hat_r, col="blue")
    
  # Função de log-verossimilhança
    
    vero <- function(p){
      b0_hat = p[1]
      b1_hat = p[2]
      de_hat = p[3]
      
      y_hat = b0_hat+b1_hat*x
      
      v = dnorm(y, mean= y_hat, sd= de_hat, log= T)
      soma_v = sum(v)
      
      return(soma_v)
      
    }
  
  # Plot da função de log-verossimilhança (apenas para ilustração)
    
    b1valores <- function(a){
      return(vero(c(b0,a,de)))
      
    }
    
    b1vero = lapply(seq(0,4, by=.05), b1valores)
    
    plot(seq(0,4, by=.05), b1vero, type="l",
         xlab = "valores do parametro b1", ylab = "log verossimilhança")
    
  # Priori
    
    priori <- function(p){
      b0_hat = p[1]
      b1_hat = p[2]
      de_hat = p[3]
      
      b0_prior = dnorm(b0_hat, sd=1, log= T)
      b1_prior = dnorm(b1_hat, sd=1, log= T)   
      de_prior = dnorm(de_hat, sd=1, log= T)
      
      return(b0_prior + b1_prior + de_prior)
      
    }
  
  # Nominador da posteriori
    
    nposteriori <- function(p){
    
      return(vero(p) + priori(p))
      
    }
  
  # Implementação do algoritmo #
  
    # Desvios padrões das distribuições utilizadas para amostrar (padrão do algoritmo)
    MCMC_de = c(0.1,0.1,0.1)
    
    # Função proposta utilizada para gerar as amostras i+1 que serão comparadas com os valores de i
      
      funcao_proposta <- function(p){
        
        aux = rnorm(3,mean= p, sd= MCMC_de)
        
        return(c(aux[1],aux[2],abs(aux[3])))
        
      }
    
    # Loop para implementar o algoritmo
      
    metropolis_MCMC <- function(valor_inicial, it){
      
      # Cria a matriz da cadeia 
      cadeia = array(dim = c(it+1,3))
      
      # Atribui os valores iniciais
      cadeia[1,] = valor_inicial
      
      # Loop para gerar as amostras e comparar
      for (i in 1:it){
        proposta = funcao_proposta(cadeia[i,])
        
        prob = exp(nposteriori(proposta) - nposteriori(cadeia[i,]))
        
        # Aceita a nova amostra se prob for maior que uma probabilidade aleatória
        if (runif(1) < prob){
          
          cadeia[i+1,] = proposta
          
        }else{
          
          cadeia[i+1,] = cadeia[i,]
          
        }
        
        
      }
      
      return(cadeia)
      
    }
    
  # Estima para a o modelo criado
    
    valor_inicial = c(1,1,1)
    it = 10e4
    
    cadeia = metropolis_MCMC(valor_inicial,it)
    
  # Testa autocorrelação da amostra gerada
    
    acf(cadeia)
    
  # Para diminuir o viés do valor inicial e diminuir a autocorrelação da amostra, utilizamos
  # um período de burn-in e um step-size maior que 1
    
    bi = 10e3
    lg = 5
    nf = (it-bi)/lg
    
    cadeia_final = array(dim = c(nf,3))
    
    for(i in 1:nf){
      cadeia_final[i,] = cadeia[(bi+1)+(i-1)*de,]
      
    }
    
  # Testa novamente para autocorrelação
    
    acf(cadeia_final)
      
  # Resultados
    
    # Cria janela para os plots
      windows()
    
    # Altera estrutura para a janela com os plots
      par(mfrow = c(2,3))
    
    # Histograma de b0
      
      hist(cadeia_final[,1], nclass = 30, main="posteriori de b0", 
            xlab="Valor real = linha vermelha")
      
      abline(v = mean(cadeia_final[,1]))
      abline(v = b0, col="red" )
      
    # Histograma de b1
      
      hist(cadeia_final[,2], nclass = 30, main="posteriori de b1", 
            xlab="Valor real = linha vermelha")
      
      abline(v = mean(cadeia_final[,2]))
      abline(v = b1, col="red" )    
  
    # Histograma de de
      
      hist(cadeia_final[,3], nclass = 30, main="posteriori do desvio-padrão", 
            xlab="Valor real = linha vermelha")
      
      abline(v = mean(cadeia_final[,3]))
      abline(v = de, col="red" )  
    
    # Valores da amostragem (cadeia gerada)
        
      plot(cadeia_final[,1], type = "l", xlab="Valor real = linha vermelha", 
            main = "Valores da cadeia de b0", ylab="Cadeia")
      
      abline(h = b0, col="red" )
      
      plot(cadeia_final[,2], type = "l", xlab="Valor real = linha vermelha", 
            main = "Valores da cadeia de b1", ylab="Cadeia")
      
      abline(h = b1, col="red" )
      
      plot(cadeia_final[,3], type = "l", xlab="Valor real = linha vermelha", 
            main = "Valores da cadeia de de", ylab="Cadeia")
      
      abline(h = de, col="red" )
    