################################################################################ 
#
#  ANÁLISE DE DADOS GENÔMICOS 
#
#  Preparação de genótipos para controle de qualidade utilizando o software
#  QCF90
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
setwd("C:/Users/Thales/Google Drive/UNESP_GMA/AnlDadosGen/aula2/ex2")


# 2. CARREGANDO OS DADOS #######################################################
# Objeto R salvo da leitura dos arquivos de genótipos     
load("genodat.Rdata")   
dim(genodat)

# Objeto adicional com outros genótipos
load("genodat2.Rdata")   
dim(genodat2)

# Concatenando genótipos em um único objeto - ***importante***: os objetos devem
# ter os mesmos SNPs e as colunas igualmente ordenadas
genotipo <- rbind(genodat, genodat2)
dim(genotipo)

# Removendo objetos anteriores para liberar memória
rm(genodat, genodat2)

# Lendo "SNP map"
snpmap <- read.table("SNP_Map_50K.txt", sep = "\t", header = T)
str(snpmap)

mapa <- subset(snpmap, select = c(Name, Chromosome, Position))
head(mapa)
table(mapa$Chromosome)


# 3. LIMPEZA E PRÉ-PROCESSAMENTO ###############################################
# Filtrando apenas autossomos
mapa <- subset(mapa, as.numeric(Chromosome) >= 1)
table(mapa$Chromosome)

# Ordenando por cromossomo e posição
ordena <- order(as.numeric(mapa$Chromosome), mapa$Position) 
mapa <- mapa[ordena, ]
head(mapa)
tail(mapa)

# Removendo SNPs com mesma posição
chrpos <- paste(mapa$Chromosome, mapa$Position, sep = "_")
## Listando os duplicados
mapa[duplicated(chrpos)|duplicated(chrpos, fromLast=TRUE), ] 
dupli <- which(duplicated(chrpos))
mapa  <- mapa[-dupli, ]

# Filtrando SNPs nos genótipos
dim(genotipo)
nrow(mapa)
head(colnames(genotipo))
snp.ok <- match(mapa$Name, colnames(genotipo))
which(is.na(snp.ok))
genotipo <- genotipo[ ,snp.ok]
dim(genotipo)


# 4. SALVANDO SNPMATRIX COMO TEXTO PARA RODAR QC NO SOFTWARE QCF90 #############
sampleIDs <- sprintf('%-10s', rownames(genotipo))
rownames(genotipo) <- sampleIDs   
write.SnpMatrix(genotipo, file = "geno_qcf90.txt", quote = F, row.names = T, 
                col.names = F, na = 5, sep = "")


# 5. FAZENDO O CONTROLE DE QUALIDADE NO SOFTWARE QCF90 #########################
# QCF90 --dry-run (Verifica inconsistências na execução do QC)
command1 <- "qcf90 --snpfile geno_qcf90.txt --crm 0.90 --cra 0.97 --maf 0.02 --dry-run"
system(command1)

# QCF90 --save-clean (Executa o QC)
command2<-"qcf90 --snpfile geno_qcf90.txt --crm 0.90 --cra 0.97 --maf 0.02 --save-clean"
system(command2)