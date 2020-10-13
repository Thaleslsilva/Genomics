################################################################################ 
#
#  ANÁLISE DE DADOS GENÔMICOS 
#
#  Avalia a presença de estrutura (estratificação) populacional em um objeto em
#  R com amostras de dados genotípicos
#
#  Obs.: Este script necessita do "Utils_ParG_EstPop.R" na mesma pasta para ser 
#  executado corretamente
#
#  Versão 1.0                                       Dev.: Roberto Carvalheiro
#  Atualização: 07/10/2020                        Atual.: Thales de Lima Silva
#
################################################################################ 


# 1. PREPARANDO O AMBIENTE #####################################################
# Limpando workspace
rm(list=ls()) 
# Não deixar o R converter automaticamente caracter em fator
options(stringsAsFactors=F) 

# Definindo o diretório de trabalho
setwd("C:/Users/Thales/Google Drive/UNESP_GMA/AnlDadosGen/aula3/Ex3")
source("Utils_ParG_EstPop.R")


# 2. CARREGANDO OS DADOS #######################################################
# Objeto R salvo da leitura dos arquivos de genótipos
load("genostrat.Rdata")
dim(genostrat)


# 3. APLICANDO ANÁLISE NO CONUNTO DE DADOS ORIGINAL ############################
# Criando variáveis globais vazias
G.matrix <- matrix()     # Matriz G
G.list   <- data.frame() # Dataframe com elementos da triangular superior da G
rel.pc   <- matrix()     # Matriz de autovalores na base 100%
pc.ord   <- matrix()     # Matriz de componentes principais

# Calcula a matiz de parentesco genômico (G) pelo método 1 de VanRaden
# Converte os elementos da triangular superior da G em um dataframe (G.list)
VanRaden.1(genostrat)
G <- G.matrix

# Estratifica as amostras pelo método de componentes principais
PCA.geno(G.matrix)

# Imprime o gráfico PCplot
PC.plot(pc.ord)

# Definindo subgrupos manualmente a partir do resultado do PCA
g1 <- which(pc.ord[,1] < 0)
g2 <- which(pc.ord[,1] > 0)


# 4. APLICANDO ANÁLISE NO PRIMEIRO SUBGRUPO ####################################
# Calcula a matiz de parentesco genômico (G) pelo método 1 de VanRaden
VanRaden.1(genostrat[g1,])
G1 <- G.matrix

# Aplica o método de componentes principais no primeiro subgrupo
PCA.geno(G.matrix)

# Imprime o gráfico PCplot do primeiro subgrupo
PC.plot(pc.ord)


# 5. APLICANDO ANÁLISE NO SEGUNDO SUBGRUPO #####################################
# Calcula a matiz de parentesco genômico (G) pelo método 1 de VanRaden
VanRaden.1(genostrat[g2,])
G2 <- G.matrix

# Aplica o método de componentes principais no segundo subgrupo
PCA.geno(G.matrix)

# Imprime o gráfico PCplot do segundo subgrupo
PC.plot(pc.ord)


# 6. COMPARANDO GRUPOS #########################################################
# Verificar a relação de parentesco entre animais dos diferentes grupos (G e G1)
# na presença de estrutura populacional
summary(diag(G1))
summary(diag(G)[1:nrow(G1)])

plot(diag(G)[1:nrow(G1)],diag(G1), 
     main = "Distribuição da Relação de Parentesco",
     xlab = "Parentesco em G", ylab = "Parentesco em G1")

cor(diag(G)[1:nrow(G1)],diag(G1))

# Verificar a relação de parentesco entre animais dos diferentes grupos (G e G2)
# na presença de estrutura populacional
summary(diag(G2))
summary(diag(G)[(length(g1)+1):nrow(G)])

plot(diag(G)[(length(g1)+1):nrow(G)],diag(G2),
     main = "Distribuição da Relação de Parentesco",
     xlab = "Parentesco em G", ylab = "Parentesco em G2")

cor(diag(G)[(length(g1)+1):nrow(G)],diag(G2))