# ---------------------------------------------------------------------------
#           ANÁLISE DE CONVERGENCIA EM SAIDAS DO BLUPF90 
#
# Created by: Thales de Lima Silva 
# Colaborators: Vinicius S. Junqueira and Delvan Silva
#
# Updated at: 09-05-2020
# Thales de Lima Silva 
# ---------------------------------------------------------------------------

### Carregar pacotes
library(data.table)
library(dplyr)

### Definir diretorio de trabalho
rm(list = ls())
setwd('C:/Users/Thales.Thales-PC/Dropbox/UNESP GMA/Tese/Analises/Analise2/BiPEHT/')

### Carregar arquivo de dados (saida do postgibbs)
data <- fread('postgibbs_samples', 
              col.names = c("Ef1","Ef2","Ef3","Va1", "Cov1e2", "Va2", "Ve1", "Ve2"))

# ===========================================================================

### Criar colunas com herdabilidades e correlacao
data1 <- data %>% mutate(Ha1 = Va1/(Va1+Ve1),
                         Ha2 = Va2/(Va2+Ve2),
                         Cor1e2 = Cov1e2/(sqrt(Va1)+sqrt(Va2)))

### Analise Exploratoria 
summary(data1[,4:11])

plot(data1$Va1, type = "l", xlab = "Iteracoes", ylim = c(3.5,4.5))

#Salvar grafico: png("AP.png", width = 480, height = 480, units = "px", pointsize = 12,
#    bg = "white", res = NA, family = "", restoreConsole = TRUE,
#    type = c("windows", "cairo", "cairo-png"), antialias)


data2 <- data1[20000:200000,] #Da um chute no valor do burn-in
autocorr.plot(data2, lag.max= 50, auto.layout = T) #Serve para definir o thin


### Gravar arquivo editado na pasta
write.table(data1, 'postgibbs_samples1', sep = ' ', quote = F, row.names = F, 
            append = F, col.names = F)

# ===========================================================================

# Biblioteca de utilitários para teste de convergência
source("convergence_utils.R")

# Definir parâmetros para a função Bayes Stats
dir        = "C:/Users/Thales.Thales-PC/Dropbox/UNESP GMA/Tese/Analises/Analise2"
data_files = c('BiPEHT')
burnin     = 25000
thin       = 25
outfile    = "resultados_analises"

####Colocar um loop aqui
# Rodar a função 
Get.Bayes.Stats(dir, data_files, burnin, thin, outfile)
# Carregar p-valor dos parâmetros genéticos calculados
path <- paste(dir, outfile, "P_Values.txt", sep = "/")
p.values = fread(path, drop = 1, 
                 col.names = c("Va1", "Cov1e2", "Va2", "Ve1", "Ve2", "Ha1", "Ha2", "Cor1e2"))
ifelse(p.values[,] > 0.05, "convergiu!", "ainda não")


for(i in 1:length(p.values)){
  ifelse(p.values[i] > 0.05, "convergiu!", "ainda não")
}
if (p.values[,] > 0.05) {
  "convergiu!"
} else {"ainda não"}
  
#============================================================================

nro_iter <- nrow(data)

# Funcao para selecao de Burn-in e Thin
conv.parameter.selection(num.iters = 20, parameter, ) {
  set.seed(10)
  conv <- 0
  while(conv == 0) {
    for burnin 0:(nro.iter/2) {  #repete até 50% do número de iterações
      # Roda a funcao
      GetBayesStats(dir,data_files,burnin,thin,outfile)
      # resultado: se Geweke >= 0.05
      if # Se não convergiu
      thin <- thin + 1
      burnin <- burnin + 1000
      GetBayesStats(dir,data_files,burnin,thin,outfile)
      else
    }
  }
}  
  
  
}
