################################################################################ 
#
#  ANÁLISE DE DADOS GENÔMICOS 
#
#  GWAS utilizando o método BayesC e o software GS3
#  Gera Trace plots e Manhattan plots
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
setwd("C:/Users/Thales/Dropbox/Repositories/Genomics/GWAS/GS3")


# 2. CARREGANDO OS DADOS #######################################################
# Arquivo com mapa dos SNPs (utilizar coordenadas ARS-UCD1.2: chrA e posA)
mapa <- read.table("HDmap_UMD_ARS.txt", header = T, sep = "\t")
head(mapa)
str(mapa)

# lista de SNPs do arquivo de genótipos
SNPnames <- read.table("SNPnames.txt")
head(SNPnames)


# 3. PRÉ-PROCESSAMENTO #########################################################
# conferir se localizou as coordenadas de todos os SNPs
linkSNPs <- match(SNPnames$V1,mapa$snpname)
length(which(is.na(linkSNPs)))

# mapa dos SNPs no arquivo de genótipos
map <- mapa[linkSNPs, c(1, 4, 5)]
colnames(map) <- c('Name', 'Chromosome', 'Position')
head(map)

# converter identificação do cromossomo para numérico
table(map$Chromosome)
map$Chromosome <- as.numeric(map$Chromosome)


# 4. ANÁLISE EXPLORATÓRIA - OUTPUTS GS3 ######################################## 
# Estimativas das variâncias e de pi
var <- read.table("GWAS_BayesC.var", header = T)
head(var)
summary(var)
# Trace plot variância aditiva dos SNPs
plot(var$vara)         
# Trace plot variância residual
plot(var$vare)         
# Trace plot variância genética aditiva
plot(var$X2varapqpi)   

plot(density(var$vara))
summary(var$vara)
plot(density(var$vare))
summary(var$vare)
plot(density(var$X2varapqpi))
summary(var$X2varapqpi)

# Estimativa de herdabilidade
h2 <- var$X2varapqpi / (var$X2varapqpi + var$vare) 
plot(density(h2))
summary(h2)

# Output: solucoes SNPs
sol <- read.table("GWAS_BayesC.sol", header = T)
head(sol)
table(sol$effect)
SNP <- subset(sol, effect == 2)
summary(SNP)
# Densidade a posteriori de pi
plot(density(SNP$p)) 

# Frequências alélicas
freq <- read.table("freq") 
freq <- freq$V1
summary(freq)

head(SNP)
# abs(solution)
betaSNP <- abs(SNP[,3])                     
summary(betaSNP)
# posterior inclusion probability (Berg, 2013)
pip <- SNP[,5]                              
summary(pip)
# Bayes Factor (Legarra, 2013)
BF <- 99 * SNP[,5] / (1 - SNP[,5])                
summary(BF)
var <- 100*(SNP[,3]**2) * (2*freq*(1-freq)) / sum((SNP[,3]**2) * (2*freq*(1-freq)))
summary(var)
sum(var)

# Listar top 10 SNPs com base nos efeitos
topSNP <- order(betaSNP, decreasing = T)
findSNP <- data.frame(Chr = map[topSNP,2], BPpos = map[topSNP,3], 
                      SNPeffect = SNP[topSNP,3], PIP = pip[topSNP], 
                      BF = BF[topSNP], var = var[topSNP], freq[topSNP])
head(findSNP, n = 10)


# 5. SALVANDO LISTA DE TOP SNPs PARA INDETIFICAR GENES E QTLs ##################
# salvar lista de top SNPs para identificar genes e QTLs
write.table(findSNP[1:10,], file = "topSNPs.txt", quote = F, sep = "\t", 
            row.names = F, col.names = T)


# 6. GERANDO MANHATTAN PLOTS COM FUNÇÕES AUXILIARES#############################
# Manhattan plot para abs(SNP effect)
jpeg(file = "Manhattan_BayesC_sol.jpeg", width = 550, height = 450)
source("manhattan_function_sol.R")
pro.man <- data.frame(CHR = map[,2], BP = map[,3], P = betaSNP)
manhattan(pro.man)     
dev.off()

# Manhattan plot para posterior inclusion probability
jpeg(file = "Manhattan_BayesC_pip.jpeg", width = 550, height = 450)
source("manhattan_function_pip.R")
pro.man <- data.frame(CHR = map[,2], BP = map[,3], P = pip)
manhattan(pro.man)     
dev.off()

# Manhattan plot para Bayes Factor
jpeg(file = "Manhattan_BayesC_BF.jpeg", width = 550, height = 450)
source("manhattan_function_BF.R")
pro.man <- data.frame(CHR = map[,2], BP = map[,3], P = BF)
manhattan(pro.man)     
dev.off()

# Manhattan plot para % variância
jpeg(file = "Manhattan_BayesC_var.jpeg", width = 550, height = 450)
source("manhattan_function_var.R")
pro.man <- data.frame(CHR = map[,2], BP = map[,3], P = var)
manhattan(pro.man)     
dev.off()
