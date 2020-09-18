################################################################################ 
#
#  ANÁLISE DE DADOS GENÔMICOS
#
#  Leitura de genótipos com o programa illumina2preGS
#
#  Atualização: 18/09/2020                          Dev.: Thales de Lima Silva
#
################################################################################ 


# 1. PREPARANDO O AMBIENTE #####################################################
# Limpando workspace
rm(list=ls()) 
# Não deixar o R converter automaticamente caracter em fator
options(stringsAsFactors=F) 

# Definindo o diretório onde estão os arquivos de genótipos gerados pelo 
# software GenomeStudio
setwd("C:/Users/Thales/Google Drive/UNESP_GMA/AnlDadosGen/aula1/ex2")


# 2. CARREGANDO OS DADOS #######################################################
# Criando objeto com o nome dos arquivos que serão lidos. Arquivos iniciando com
# a letra "F" (FinalReport"i".txt)
filenames <- dir(pattern="^[F]") 
length(filenames)
head(filenames)

# Concatenando arquivos de genotipos
genodat <- NULL

for (k in filenames){
  report <- read.table(k, skip = 11, sep="\t")
  genodat<- rbind(genodat, report)
}

# Conferir a dimensão do objeto com os genótipos
dim(genodat)  

# Adicionando nomes às colunas do DataFrame
colnames(genodat) <- c("SNP Name", "Sample ID", "Allele1 - Forward",
                       "Allele2 - Forward",	"Allele1 - Top",	"Allele2 - Top",
                       "Allele1 - AB", "Allele2 - AB", "Log R Ratio",	
                       "B Allele Freq", "GC Score")
View(genodat)

# Salvando objeto (genótipos) que será utilizado nos próximos passos
write.table(genodat, file="FinalReport.txt", quote=F, row.names=F, col.names=T, 
            sep="\t")

################################################################################
#  DEVE-SE COPIAR O CABEÇALHO DE UMA DOS ARQUIVOS DO GenomeStudio E COLAR NO   #
#  ARQUIVO "FinalReport.txt" ANTES DE USÁ-LO COMO INPUT NO ILLUMINA2PREGS!!!   #
################################################################################
