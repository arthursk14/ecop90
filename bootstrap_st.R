
# Bootstrap para modelagem de séries de tempo:
# Método "moving blocks", baseado na seção 8.6 do livro de Bradley Efron, R.J. Tibshirani (1994)
  
  # Seed
  
  seed = sample(seq.int(from=1, to=20000, by=1),1)
  set.seed(seed)
  
  
  # Utilizando série de preços das ações da IBM
  
    # Carrega pacotes
      require("tseries")
      require("quantmod")
      require("rugarch")
    
    # Pega série
      
      # Declara as datas
      startDate = as.Date("2007-01-03")
      endDate = as.Date("2019-11-27")
      
      # Pega dados
      getSymbols("IBM", from=startDate, to=endDate)
      
      head(IBM)
      IBM[length(IBM[,1]),]
      
      # Monta a série de preços de fechamento
      y = IBM[,4]
      
      # Série dos retornos
        n = length(y)
        
        r = diff(log(y))
        r = r[2:n]
      
      # Gráfico das séries
      chart_Series(y)
      chart_Series(r)
      
    # Testa para estacionariedade
    adf.test(y)
    adf.test(r)
    
    # Estima modelo AR(p) para as séries (y - precos; r - retornos)
    ar(y)
    ar(r)
    
    # Tamanho da série dos retornos original
    n = length(r)
    
    # Tamanho dos "moving blocks"
    k = 40
    
    # Número de repetições (cada repetição é uma reamostragem da amostra original)
    nr = 400
    
    # Cria matriz para salvar os coeficientes AR(p) estimados
    p_hat = array(dim = c(nr,8))
    
    # Loop para a reamostragem em blocos
    
    for(i in 1:nr){
      
      # Cria o vetor para as séries reamostradas
      r_bt = rep(NA,n)
      
      # Loop específico para lidar com os blocos
      for(j in 1:ceiling(n/k)){
        
        fim = sample(k:n, size=1)
        r_bt[(j-1)*k+(1:k)] = r[fim-(k:1)+1]
        
      }
      
      # Trunca para o caso em que k não divide n
      r_bt = r_bt[1:n]
      
      # Estima os coeficientes do AR(p)
      
        # Recria a série y_bt (preços), com base em r_bt (retornos)
        y_bt = y
        
        for (q in 2:n+1){
          y_bt[q] = exp(log(y_bt[q-1])+r_bt[q-1])
          
        }
        
        # Trunca a ordem do AR(p) para no máximo 8
        ordem_max = ar(y_bt)$order
        
        if (ordem_max > 8){
          ordem_max = 8
          
        }
        
        for (l in 1:ordem_max){
          
          p_hat[i,l] <- ar(y_bt)$p[l]
        
        }
      
    }
    
    # Histograma da distribuição do parâmetro AR(1)
    hist(p_hat[,1], prob = T)
    
    # Compara a média das estimativas bootstrap com a estimativa da série original
    mean(p_hat[,1])
    ar(y)$p[1]
