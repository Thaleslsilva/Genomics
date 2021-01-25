################################################################################ 
#
#  ANÁLISE DE DADOS GENÔMICOS 
#
#  GWAS utilizando o método BayesC e o software GS3
#  Gera arquivo de genótipos para o GS3, a partir de um objeto SnpMatrix
#
#  Versão 1.0                                       Dev.: Roberto Carvalheiro
#  Atualização: 25/01/2021                        Atual.: Thales de Lima Silva
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
setwd("D:/Unesp/disciplinas/genomica_2020/aula5/GS3/")


# 2. CARREGANDO OS DADOS #######################################################
# Upload do objeto SnpMatrix com genotipos (573 animais x 16072 SNPs)
load("genodata.RData")
dim(genodata)  


# 3. PRÉ-PROCESSAMENTO #########################################################
# Salvando nomes dos SNPs para gerar Manhattan plot (posGS3.R)
SNPnames <- colnames(genodata) 
write.table(SNPnames, file = "SNPnames.txt", quote = F, row.names = F, 
            col.names = F)


# 4. SALVANDO ARQUIVO DE GENOTIPOS ############################################# 
rownames(genodata)
# Sample IDs com formato fixo (10 caracteres)
sampleIDs <- sprintf('%-10s', rownames(genodata))
(rownames(genodata) <- sampleIDs)   
write.SnpMatrix(genodata, file = "genotipos.txt", quote = F, row.names = T,
                col.names = F, na = 5, sep = "")
