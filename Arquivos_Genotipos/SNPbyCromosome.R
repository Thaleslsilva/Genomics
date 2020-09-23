################################################################################ 
#
#  ANÁLISE DE DADOS GENÔMICOS 
#
#  Comparar a quantidade de SNP's por cromossomo em diferentes painéis
#
#  Atualização: 23/09/2020                          Dev.: Thales de Lima Silva
#
################################################################################


# 1. PREPARANDO O AMBIENTE #####################################################
# Limpando workspace
rm(list=ls()) 
# Não deixar o R converter automaticamente caracter em fator
options(stringsAsFactors=F) 

# Carregando pacotes
library(dplyr)
library(stringr)
library(gtools)

# Definindo o diretório onde estão os arquivos
setwd("C:/Users/Thales/Google Drive/UNESP_GMA/AnlDadosGen/aula1/ex3")


# 2. CARREGANDO OS DADOS #######################################################
# Criando objeto com o nome dos arquivos que serão lidos, iniciando com "SNP"
filenames <- dir(pattern = "^[SNP]") 
length(filenames)
head(filenames)


# 3. PROCESSAMENTO #############################################################
# Criando tabela com a quantidade de SNP's por cromossomo
table <- data.frame(Chr = character())
panel <- 0
cols  <- str_replace(filenames, ".txt", "") 
for (k in filenames){
  panel <- panel + 1
  mapa  <- read.table(k, header = T, sep = "\t")
  conta <- mapa %>% count(Chromosome)
  colnames(conta) <- c("Chr", cols[panel])
  table <- merge(table, conta, all = T)
}

str(table)

# Ordenando tabela pela identificação do cromossomo
ordena <- mixedorder(table$Chr)
table <- table[ordena,]
View(table)

# 4. SALVANDO A TABELA GERADA A PARTIR DOS DADOS ###############################
write.table(table, file="NroSNPporChr.txt", 
            quote=F, row.names=F, col.names=T, sep="\t")
