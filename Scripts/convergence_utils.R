# ---------------------------------------------------------------------------
#           ANÁLISE DE CONVERGENCIA EM SAIDAS DO BLUPF90 
#         Biblioteca de funções para analise de convergencia
#
# Created at: 10-04-2015
# Vinicius S. Junqueira and Delvan Silva
#
# Modificated at: 08-05-2020
# Thales de Lima Silva 
# ---------------------------------------------------------------------------


# Obs: Caso tenha problemas com a acentuacao, consulte este link:
# https://support.rstudio.com/hc/en-us/articles/200532197-Character-Encoding

library(coda)
library(mcmcse)
library(boa)
library(miscTools)

# FUNCAO PARA GERAR ESTATISTICAS A POSTERIORI 
Get.Bayes.Stats <- function(dir, data_files, burnin, thin, outfile){
  nf = length(data_files)
  resul = paste(dir, '/', outfile, '/', sep='')
  dir.create(resul, showWarnings = F)
  direc = paste(resul, 'Estatisticas_completas.txt', sep = '')
  fn = file(direc, 'wt')
  grava = paste(resul, 'P_Values.txt', sep = '')
  gw = file(grava, 'wt')
  for(i in 1:nf){
    data = import_file(data = data_files, i = i, dir = dir)
    cat('Analise =', data_files[i], '\n', file = fn)
    cat('-------------------', '\n', file = fn)
    cat('Media a posteriori: ', colMeans(data[,4:11]), '\n', file = fn)
    cat('Mediana posteriori: ', colMedians(data[,4:11]), '\n', file = fn)
    cat('Moda a posteriori:  ', moda(data), '\n', file = fn)
    cat('Desvpad posteriori: ', PSD(data), '\n', file = fn)
    cat('Z-Geweke (p.valor): ', gewe(data), '\n', file = fn)
    cat('Z-Geweke_(p.valor):', gewe(data), '\n', file = gw)
    cat('Tamanho efet (ESS): ', ESS(data), '\n', file = fn)
    cat('Intervalo de credibilidade: ', '\n', file = fn)
    Intervalo_cred(data, fn = fn)
  }
  close(fn)
}

# FUNCAO PARA IMPUTAR DADOS DE FORMA RECURSIVA
import_file <- function(data, i, dir){
  # Import data files
  fold = paste(dir, '/', data[i], '/postgibbs_samples1', sep='')
  files = read.table(fold)
  # Burnin
  files2 = files[(burnin+1):nrow(files),]
  # Thin samples
  gethin = rep(1:thin, times = nrow(files2)/thin)
  idx = gethin == thin
  files3 = files2[idx,]
  return(files3)
}

# FUNCAO PARA CALCULAR A "MODA"
moda = function(data){
  val = NULL
  for(j in 1:(ncol(data)-3)){
    aaa = hist(data[[3+j]], plot=F)
    tmp = aaa$mids[aaa$counts == max(aaa$counts)]
    val = cbind(val, tmp)
  }
  return(val)
}

# FUNCAO PARA CALCULAR O "TAMANHO EFETIVO DAS AMOSTRAS"
ESS = function(data){
  val = NULL
  for(j in 1:(ncol(data)-3)){
    tmp = ess(data[[3+j]])
    val = cbind(val, tmp)
  }
  return(val)
}

# FUNCAO PARA CALCULAR O "DESVIO PADRAO A POSTERIORI"
PSD = function(data){
  val = NULL
  for(j in 1:(ncol(data)-3)){
    tmp = sd(data[[3+j]])
    val = cbind(val, tmp)
  }
  return(val)
}

# FUNCAO PARA CALCULAR O "P-VALUE" DO TESTE GEWEKE
gewe = function(data){
  tmp = boa.geweke(data[,4:11], p.first = .1, p.last = .5)
  return(tmp[,2])
}

# FUNCAO PARA CALCULAR O "INTERVALO DE CREDIBILIDADE"
Intervalo_cred = function(data, fn){
  val = NULL
  for(j in 1:(ncol(data)-3)){
    tmp1 = quantile(data[[3+j]], probs = .025)
    tmp2 = quantile(data[[3+j]], probs = .975)
    cat('Efeito ', j,' = ', tmp1, tmp2, '\n', file=fn)
  }
  return(val)
}