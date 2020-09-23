################################################################################ 
#
#  ANÁLISE DE DADOS GENÔMICOS 
#
#  Leitura de genótipos com a função read.snps.long do pacote snpStats
#
#  Atualização: 17/09/2020                          Dev.: Thales de Lima Silva
#
################################################################################ 


# 1. PREPARANDO O AMBIENTE #####################################################
# Limpando workspace
rm(list=ls()) 
# Não deixar o R converter automaticamente caracter em fator
options(stringsAsFactors=F) 

# Instalando o pacote snpStats
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("snpStats")

# Carregando o pacote snpStats
library(snpStats)

# Definindo o diretório onde estão os arquivos de genótipos gerados pelo 
# software GenomeStudio
setwd("C:/Users/Thales/Google Drive/UNESP_GMA/AnlDadosGen/aula1/ex1")


# 2. CARREGANDO OS DADOS #######################################################
# Criando objeto com o nome dos arquivos que serão lidos. Arquivos iniciando com
# a letra "F" (FinalReport"i".txt)
filenames <- dir(pattern="^[F]") 
length(filenames)
head(filenames)

# Lendo "SNP map"
snpmap <- read.table("SNP_Map_GGPHDv2.txt",sep="\t",header=T)
str(snpmap)
head(snpmap)
snpids <- snpmap$Name
head(snpids)


# 3. LIMPEZA E PRÉ-PROCESSAMENTO ###############################################
# Lendo sample IDs dos proprios arquivos de genotipos (segunda coluna)
animids <- NULL
count   <- 0
for (k in filenames){count <- count+1
                     animids[count] <- scan(k, skip=11, nlines=1, sep="\t",
                                            what="character")[2]
}
animids

# Leitura de genótipos usando a função do snpStats (read.snps.long)
# vide referência snpStats_manual.pdf
genodat <- read.snps.long(file=filenames,
                          sample.id=animids,
                          snp.id=snpids,
                          diploid=NULL,
                          fields=c(sample=2,snp=1,allele1=3,allele2=4,confidence=5),
                          codes=c("A","B"),
                          threshold=0.5,#GC score
                          lower=TRUE,
                          skip=11,
                          sep="\t",
                          verbose=TRUE, 
                          in.order=FALSE,
                          every=10000)

# Conferir a dimensão do objeto com os genótipos
dim(genodat)  

# Conferir como são salvos os nomes das linhas e das colunas
head(rownames(genodat))
head(colnames(genodat))


# 4. SALVANDO OS OBJETOS GERADOS A PARTIR DOS DADOS ############################
# Salvando objeto (genótipos) que será utilizado nos próximos passos
save(genodat,file="genodat.Rdata")

# Salvando snpmatrix no formato texto
sampleIDs <- sprintf('%-10s',rownames(genodat))
rownames(genodat) <- sampleIDs   
write.SnpMatrix(genodat, file="genodat.txt", quote=F, row.names=T, col.names=F,
                na=5, sep="")
