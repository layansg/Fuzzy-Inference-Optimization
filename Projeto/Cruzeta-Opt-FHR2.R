#SIMULACAO-OTIMIZACAO DA REGRA DE OPERACAO DE RESERVATORIO
#Regra de racionamento fuzzy-genetico


#Os arquivos que seguem devem estar salvos no Diretorio especificado
#Especifique a pasta onde estao os arquivos (Ctrl+Shift+H)

#1) ENTRADA DE DADOS

Ev <- read.table("EvCruzetaPRH.txt", header = T)    #Dados de evaporacao (mm/mes)
Q <- read.table("Vazoes.txt", header = T)           #Dados de Vazao (m3/s)
P <- read.table("Prec.txt", header = T)             #Precipitacao (mm/mes)
CAV <- read.table("CAV.txt", header = T)            #Relacao Cota-Area-Volume

options(scipen = 999)  #Desabilita notacao cientifica

#2) FUNCOES

#2.1) Funcao interpolacao da area correspondente ao volume
Area_cor <- function(S)
{
  for (linha in 1:nrow(CAV))
  {
    if (S > CAV[linha,3] & S < CAV[linha+1,3])
    {
      (V2 <- CAV[linha+1,3])
      (V1 <- CAV[linha,3])
      (A2 <- CAV[linha+1,2])
      (A1 <- CAV[linha,2])
      (A <- (((S - V1) * (A2 - A1)) / (V2 - V1)) + A1)
    }else{
      if (S == CAV[linha, 3])
      {
        A <- CAV[linha, 2]
      }
    }
  }
  A
}

#2.2) Funcao numero de dias no mes

n_dias <- function(Mes)
{
  if(Mes == 1 | Mes == 3 | Mes == 5 | Mes == 7 | Mes == 8 | Mes == 10 | Mes == 12)
  {
    31
  }else{
    if(Mes == 4 | Mes == 6 | Mes == 9 | Mes == 11)
    {
      30
    }else{
      28
    }
    
  }
}

#2.3) Funcao de determinacao do suprimento - REGRA DE RACIONAMENTO FUZZY

library(frbs)
Suprimento.fuzzy <- function(S, D, CR1, CR2, A1, A2, B1, B2)
{
  #2.3.1) Input Fuzzy Set Setting
  ## Define a matrix representing shape and parameters of membership functions of input variables.
  ## The matrix has 5 rows where the first row represent the type of the membership function whereas
  ## others are values of its parameters.                
  varinp.mf <- matrix(c(4, 0, 0, S_min, S_min,
                        4, S_min, S_min, CR1-B1, CR1+B1,
                        4, CR1-B1, CR1+B1, CR2-B2, CR2+B2,
                        4, CR2-B2, CR2+B2, S_max, S_max), 
                      nrow=5, byrow=F)
  ## Define number of fuzzy terms of input variables.
  num.varinput <- 1
  num.fvalinput <- matrix(c(4), nrow=1)           #1 input with 4 labels
  ## Give the names of the linguistic labels of each input variables.
  varinput.1 <- c("MUITO.BAIXO","BAIXO", "MEDIO", "ALTO")
  names.varinput <- c(varinput.1)
  
  #2.3.2) Inference setting
  ## Define inference parameters. 
  type.tnorm <- "MIN"
  type.snorm <- "MAX"
  type.model <- "MAMDANI"
  type.implication.func <- "ZADEH"
  type.defuz <- "WAM"  # Weighted average method
  
  #2.3.3) Rule base setting
  ##Define rulebase parameters.
  r1 <- c("MUITO.BAIXO", "and", "dont_care", "->", "ZERO")
  r2 <- c("BAIXO", "and", "dont_care", "->", "MUITO RACIONADO")
  r3 <- c("MEDIO", "and", "dont_care", "->", "POUCO RACIONADO")
  r4 <- c("ALTO", "and", "dont_care", "->", "MAXIMO")
  rule.a <- matrix(c(r1,r2,r3,r4), nrow=4, byrow = T)
  rule <- rulebase(type.model, rule.a, func.tsk = NULL)
  
  #2.3.4) Fuzzification
  ## Provide new data for predicting result
  newdata <- matrix(c(S), nrow=1)
  MF <- fuzzifier(newdata, num.varinput, num.fvalinput, varinp.mf)
  
  #2.3.5) Inference Implication
  miu.rule <- inference(MF, rule, names.varinput, type.tnorm, type.snorm)
  
  #2.3.6) Output Defuzzification Setting
  range.output <- matrix(c(0, D), nrow=2)
  ## Define number of fuzzy terms of output variable.
  num.fvaloutput <- matrix(c(4), nrow=1)              #1 output with 4 labels
  ## Give the names of the fuzzy terms of the output variable.
  varoutput.1 <- c("ZERO", "MUITO RACIONADO", "POUCO RACIONADO", "MAXIMO")
  names.varoutput <- c(varoutput.1)
  ## Define the shapes and parameters of the membership functions of the output variables.
  varout.mf <- matrix(c(1, 0, 0, A1*D, NA,
                        1, 0, A1*D, A2*D, NA,
                        1, A1*D, A2*D, D, NA,
                        1, A2*D, D, D, NA), 
                      nrow=5, byrow=F)
  ## the names of variables
  colnames.var <- c("input1", "output1")
  
  #2.3.7) Defuzzification
  result <- defuzzifier(newdata, rule, range.output, names.varoutput,
                        varout.mf, miu.rule, type.defuz, type.model, func.tsk = NULL)
  result
} #Fim da função de inferencia fuzzy

#2.4) Funcao objetivo da otimizacao
funcao_objetivo <- function(CR1.s, CR2.s, CR1.u, CR2.u, alfa1, alfa2, beta1, beta2)
{
  ### Simulação do balaço hídrico ###
  #_________________________________#
  # Atribuicoes iniciais
      #Demandas
  D_irriga_anual <- 27562464   #m3/ano
  Perc_irriga_mensal <- c(0.086433443,0.045772557,0.022052284,0.016802301,0.041870383,0.066710103,0.086716542,0.113809281,0.127722914,0.142399325,0.13529419,0.114416678)
  D_industrial <- 0.01        #m3/s
  D_prior <- 0.054            #m3/s -> abastecimento humano e dessedentacao animal
  D_prioritario <- 0
  
  meses_simulados <- 0
  falhas <- 0
  S_max <- 23545745.33
  S_min <- 1179186
  S_ti <- S_max / 2           #Reservatório inicial com metade da capacidade
  A_ti <- Area_cor(S_ti)
  S_tf <- 0
  A_tf <- 0
  A_med1 <- A_ti              #Estimativa inicial de area media
  A_med2 <- 0
  D_t <- 0
  D_total <- 0
  Y_t <- 0
  Y_total <- 0
  Sp_t <- 0
  Sp_total <- 0
  Y_prior <- 0
  Y_seco <- 0
  D_seco <- 0
  Rv_seco <- 0
  
  # Varredura do periodo de tempo em simulacao
  for(t_meses in 1:840)   #70 anos de dados historicos mensais para otimizacao
  {
    meses_simulados <- meses_simulados + 1      #Contagem de meses de iteracao
    mes <- Q[t_meses, 2]                        #Mes da iteração
    nDias <- n_dias(mes)                        #Numero de dias no mes
    D_t <- D_irriga_anual*Perc_irriga_mensal[mes] + ((D_prior+D_industrial)*nDias*24*60*60)       #Demanda no mes
    D_total <- D_total + D_t                    #Somatorio de demanda total
    D_prioritario <- D_prioritario + D_prior*nDias*24*60*60     #Somatorio de demanda prioritaria
    
    #Suprimento determinado pela regra de racionamento
    if (mes==1 | mes==2 | mes==7 | mes==8 | mes==9 | mes==10 | mes==11 | mes==12){
      Y_t <-as.double( Suprimento.fuzzy(S_ti, D_t, CR1.s, CR2.s, alfa1, alfa2, beta1, beta2) )
    }else{
      Y_t <-as.double( Suprimento.fuzzy(S_ti, D_t, CR1.u, CR2.u, alfa1, alfa2, beta1, beta2) )
    }
    
    #Condicao de disponibilidade de água abaixo de D_t*alfa1
    WA <- S_ti - S_min
    if (S_ti > S_min)
    {
      if (WA < (D_t * alfa1))
      {
        Y_t <- WA
      }
    }else{
      Y_t <- 0
    }
    
    Y_total <- Y_total + Y_t                    #Somatorio de suprimento total
    
    #Soma de suprimento e demanda no periodo seco do ano
    if (mes==1 | mes==2 | mes==7 | mes==8 | mes==9 | mes==10 | mes==11 | mes==12){
      Y_seco <- Y_seco + Y_t
      D_seco <- D_seco + D_t
    }
    
    #Soma de suprimento prioritario
    if (Y_t > 0){
      if (Y_t < (D_prior*nDias*24*60*60)){
        Y_prior <- Y_prior + Y_t
      }else{
        Y_prior <- Y_prior + (D_prior*nDias*24*60*60)
      }
    }
    
    #Minimizacao do erro da area de superficie media
    erro <- 1
    while (erro >= 0.05)
    {
      #Calculo do armazenamento final
      S_tf <- S_ti + (Q[t_meses, 1] * nDias * 24 * 60 * 60) +
        + (((P[mes, 1]/1000) - (Ev[mes, 1]/1000)) * A_med1) - Y_t 
      
      #Correcao se o armazenamento final for superior ao maximo
      if (S_tf > S_max)
      {
        #Sp_t <- S_tf - S_max
        S_tf <- S_max
      }
      
      #Correcao se o armazenamento final for inferior a zero
      if (S_tf <= 0)
      {
        S_tf = 0.000001 
      }
      A_tf <- Area_cor(S_tf)              #Interpolacao: Area final correspondente ao volume final
      A_med2 <- (A_ti + A_tf) / 2         #Calculo da area media
      erro <- abs(1 - (A_med2 / A_med1))  #Calculo do erro entre iteracoes da area media
      A_med1 <- A_med2                    #Atribuicao para proxima iteracao
    }#Fecha correcao do erro
    
    #Atribuicao para proximo mes, final torna-se inicial
    S_ti <- S_tf
    A_ti <- A_tf
    
  }#Fecha o tempo de simulacao
  #_____________________________________# 
  # Fim da simulacao do balanco hidrico #
  
  #Penalizacao para restricoes
  if(CR1.u>CR2.u){
    c1 <- 0
  }else{
    c1 <- 1
  }
  
  if(CR1.s>CR2.s){
    c2 <- 0
  }else{
    c2 <- 1
  }
  
  if(alfa1>alfa2){
    c3 <- 0
  }else{
    c3 <- 1
  }
  
  if(beta1<0 | beta1>2*(CR1.s-S_min) | beta1>2*(CR1.u-S_min) | beta1>(CR2.s-CR1.s) | beta1>(CR2.u-CR1.u)){
    c4 <- 0
  }else{
    c4 <- 1
  }
  if(beta2<0 | beta2>(CR2.s-CR1.s) | beta2>(CR2.u-CR1.u) | beta2>2*(S_max-CR2.s) | beta2>2*(S_max-CR2.u)){
    c5 <- 0
  }else{
    c5 <- 1
  }
  
  if((Y_prior / D_prioritario)*0.75){   #calibrar restricao
    c6 <- 0
  }else{
    c6 <- 1
  }
  
  #Valor resultante da funcao objetivo
  #Rt <- 1-(falhas/meses_simulados)         #Confiabilidade
  #Rv <- Y_total/D_total                    #Eficiencia geral
  Rv_seco <- Y_seco / D_seco                #Eficiencia no periodo seco
  #Rv_prior <- Y_prior / D_prioritario       #Eficiencia para demanda prioritaria
  
  Rv_seco * c1 * c2 * c3 * c4 * c5 * c6     #Resultado da funcao objetivo
} #Fecha a funcao objetivo

###3) OTIMIZACAO
library(GA)
library(parallel)
library(doParallel)

S_max <- 23545745.33
S_min <- 1179186
min<-c(S_min,S_min,S_min,S_min,0,0,0,0)
max<-c(S_max,S_max,S_max,S_max,1,1,2553737,2553737)

system.time({         #Inicia analise do tempo de processamento
  GA<-ga(type="real-valued", fitness=function(x){funcao_objetivo(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8])},
         lower=min, upper=max, popSize=100, pcrossover = 0.8, pmutation=0.1, maxiter=3000, run = 1000,
         names=c("CR1.s", "CR2.s", "CR1.u", "CR2.u", "alfa1", "alfa2", "beta1", "beta2"),
         #suggestions=c(S_max/2, S_max, S_max/2, S_max, 1, 1, 2553737, 2553737),
         parallel = T)
})  #Fecha analise do processamento

summary(GA)
#GA@solution

###4) PLOTAGENS

plot(GA)

#TESTES
CR1.s <- S_max/3
CR2.s <- S_max/2
CR1.u <- S_max/3
CR2.u <- S_max/2
alfa1 <- 1
alfa2 <- 1
beta1 <- 2553737
beta2 <- 2553737
funcao_objetivo(CR1.s, CR2.s, CR1.u, CR2.u, alfa1, alfa2, beta1, beta2)