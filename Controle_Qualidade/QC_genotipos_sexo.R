################################################################################ 
#
#  ANÁLISE DE DADOS GENÔMICOS 
#
#  Inferindo o sexo de amostras (bovinos) com base no nível de heterozigosidade 
#  do cromossomo X
#
#  Versão 1.0                                       Dev.: Roberto Carvalheiro
#  Atualização: 05/10/2020                        Atual.: Thales de Lima Silva
#
################################################################################ 


# 1. PREPARANDO O AMBIENTE #####################################################
# Limpando workspace
rm(list=ls()) 
# Não deixar o R converter automaticamente caracter em fator
options(stringsAsFactors=F) 

# Carregando o pacote snpStats
library(snpStats)

# Definindo o diretório onde estão os arquivos de genótipos gerados pelo 
# software GenomeStudio
setwd("C:/Users/Thales/Google Drive/UNESP_GMA/AnlDadosGen/aula2/ex3")


# 2. CARREGANDO OS DADOS #######################################################
# Objeto R salvo da leitura dos arquivos de genótipos
load("genoAN.Rdata")
dim(genoAN)

# Lendo "SNP Map"
snpmap <- read.table("SNP_Map_GGPHDv2.txt", sep = "\t", header = T)
str(snpmap)
table(snpmap$Chromosome)


# 3. LIMPEZA E PRÉ-PROCESSAMENTO ###############################################
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


# 4. CONFERINDO ESTATÍSTICAS DAS AMOSTRAS PARA O CROMOSSOMO X ##################
sample.stat <- row.summary(genoX) 
summary(sample.stat)
plot(sample.stat$Call.rate,sample.stat$Heterozygosity)
hist(sample.stat$Heterozygosity)

# 5. INFERINDO O SEXO DOS ANIMAIS ##############################################
# Fêmeas
femea <- which(sample.stat$Heterozygosity >= 0.1)
length(femea)
rownames(genoX)[femea]
################################################################################
# A letra "F" no nome da amostra diz que ela foi originalmente identificada    #
# como fêmea                                                                   #
################################################################################

# Machos
macho <- which(sample.stat$Heterozygosity < 0.1)
length(macho)
rownames(genoX)[macho]
################################################################################
# A letra "M" no nome da amostra diz que ela foi originalmente identificada    #
# como macho                                                                   #
################################################################################