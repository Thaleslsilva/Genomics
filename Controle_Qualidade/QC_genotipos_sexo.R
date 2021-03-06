################################################################################ 
#
#  AN�LISE DE DADOS GEN�MICOS 
#
#  Inferindo o sexo de amostras (bovinos) com base no n�vel de heterozigosidade 
#  do cromossomo X
#
#  Vers�o 1.0                                       Dev.: Roberto Carvalheiro
#  Atualiza��o: 05/10/2020                        Atual.: Thales de Lima Silva
#
################################################################################ 


# 1. PREPARANDO O AMBIENTE #####################################################
# Limpando workspace
rm(list=ls()) 
# N�o deixar o R converter automaticamente caracter em fator
options(stringsAsFactors=F) 

# Carregando o pacote snpStats
library(snpStats)

# Definindo o diret�rio onde est�o os arquivos de gen�tipos gerados pelo 
# software GenomeStudio
setwd("C:/Users/Thales/Google Drive/UNESP_GMA/AnlDadosGen/aula2/ex3")


# 2. CARREGANDO OS DADOS #######################################################
# Objeto R salvo da leitura dos arquivos de gen�tipos
load("genoAN.Rdata")
dim(genoAN)

# Lendo "SNP Map"
snpmap <- read.table("SNP_Map_GGPHDv2.txt", sep = "\t", header = T)
str(snpmap)
table(snpmap$Chromosome)


# 3. LIMPEZA E PR�-PROCESSAMENTO ###############################################
# Filtrando apenas SNPs do cromossomo X
xchr <- which(snpmap$Chromosome == "X")
mapaX <- snpmap[xchr, 1:4]
head(mapaX)

# filtrar SNPs nos genotipos
dim(genoAN)
snp.x <- match(mapaX$Name, colnames(genoAN))
which(is.na(snp.x))
genoX <- genoAN[, snp.x]
dim(genoX)


# 4. CONFERINDO ESTAT�STICAS DAS AMOSTRAS PARA O CROMOSSOMO X ##################
sample.stat <- row.summary(genoX) 
summary(sample.stat)
plot(sample.stat$Call.rate,sample.stat$Heterozygosity)
hist(sample.stat$Heterozygosity)

# 5. INFERINDO O SEXO DOS ANIMAIS ##############################################
# F�meas
femea <- which(sample.stat$Heterozygosity >= 0.1)
length(femea)
rownames(genoX)[femea]
################################################################################
# A letra "F" no nome da amostra diz que ela foi originalmente identificada    #
# como f�mea                                                                   #
################################################################################

# Machos
macho <- which(sample.stat$Heterozygosity < 0.1)
length(macho)
rownames(genoX)[macho]
################################################################################
# A letra "M" no nome da amostra diz que ela foi originalmente identificada    #
# como macho                                                                   #
################################################################################