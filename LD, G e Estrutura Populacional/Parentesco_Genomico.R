################################################################################ 
#
#  ANÁLISE DE DADOS GENÔMICOS 
#
#  Comparando diferentes métodos de cálculo da matriz de parentesco genômico (G)
#  Compara as matrizes G em relação aos elementos da diagonal e os elementos 
#  fora da diagonal. Faz a comparação com base em estatísticas descritivas 
#  (média, mínimo, máximo, ...) e correlação.
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

# Carregando o pacote snpStats
library(snpStats)

# Definindo o diretório onde estão os arquivos de genótipos gerados pelo 
# software GenomeStudio
setwd("C:/Users/Thales/Google Drive/UNESP_GMA/AnlDadosGen/aula3/Ex2")


# 2. CARREGANDO OS DADOS #######################################################
# Objeto R salvo da leitura dos arquivos de genótipos
load("geno1.Rdata")
dim(geno1)


# 3. PRÉ-PROCESSAMENTO E ANÁLISE EXPLORATÓRIA ##################################
# M: matriz de genotipos (n alelos B)
M <- as(geno1, 'numeric')
M[1:10, 1:5]

snpstat <- col.summary(geno1)
head(snpstat)
# RAF é a frequência do alelo B
fB <- snpstat$RAF 


# 4. CALCULADO A MATRIZ G POR DIFERENTES MÉTODOS ###############################

## 4.1. Matriz G pelo Método 1 VanRaden 
P <- 2 * rep(1, nrow(M)) %*% t(fB)
Z <- M - P
sum2pq <- sum(2 * fB * (1 - fB))
G <- (Z %*% t(Z)) / sum2pq
summary(diag(G))
plot(density(G))

# Convertendo os elementos da triangular superior da G em um dataframe (Glist)
Glist <- G
Glist[lower.tri(Glist,diag = FALSE)] <- NA
Glist <- na.omit(data.frame(as.table(Glist)))
colnames(Glist) <- c("ANi", "ANj", "Gij")
head(Glist)

# Estatísticas da diagonal
summary(Glist$Gij[Glist$ANi == Glist$ANj])
plot(density(Glist$Gij[Glist$ANi == Glist$ANj]))

# Estatísticas off-diagonal
summary(Glist$Gij[Glist$ANi != Glist$ANj])
plot(density(Glist$Gij[Glist$ANi != Glist$ANj]))


## 4.2. Matriz G pelo Método 2 VanRaden (G2)
D  <- 1 / (ncol(M) * (2 * fB * (1 - fB)))
G2 <- Z %*% (D * t(Z))

# Convertendo os elementos da triangular superior da G2 em um dataframe (G2list)
G2list <- G2
G2list[lower.tri(G2list,diag = FALSE)] <- NA
G2list <- na.omit(data.frame(as.table(G2list)))
colnames(G2list) <- c("ANi", "ANj", "G2ij")
head(G2list)

# Comparando diagonais G x G2
summary(diag(G))
summary(diag(G2))
cor(diag(G), diag(G2))
plot(diag(G), diag(G2))

# Comparando off-diagonais G x G2
summary(Glist$Gij[Glist$ANi != Glist$ANj])
summary(G2list$G2ij[G2list$ANi != G2list$ANj])
cor(Glist$Gij[Glist$ANi != Glist$ANj], G2list$G2ij[G2list$ANi != G2list$ANj])
plot(Glist$Gij[Glist$ANi != Glist$ANj], G2list$G2ij[G2list$ANi != G2list$ANj])


## 4.3. Matriz G Normalizada (GN)
trace <- (sum(diag(Z %*% t(Z))))
GN <- (Z %*% t(Z)) / (trace / nrow(M))

# Convertendo os elementos da triangular superior da GN em um dataframe (GNlist)
GNlist <- GN
GNlist[lower.tri(GNlist, diag = FALSE)] <- NA
GNlist <- na.omit(data.frame(as.table(GNlist)))
colnames(GNlist) <- c("ANi", "ANj", "GNij")
head(GNlist)

# Comparando diagonais G x GN
summary(diag(G))
summary(diag(GN))
cor(diag(G), diag(GN))
plot(diag(G), diag(GN))
trace / nrow(M)
sum2pq

# Comparando off-diagonais G x GN
summary(Glist$Gij[Glist$ANi != Glist$ANj])
summary(GNlist$GNij[GNlist$ANi != GNlist$ANj])
cor(Glist$Gij[Glist$ANi != Glist$ANj], GNlist$GNij[GNlist$ANi != GNlist$ANj])
plot(Glist$Gij[Glist$ANi != Glist$ANj], GNlist$GNij[GNlist$ANi != GNlist$ANj])


################################################################################
#                                                                              #
# REFERÊNCIAS BIBLIOGRÁFICAS:                                                  #
#                                                                              #
# VanRaden PM. Genomic measures of relationship and inbreeding. Interbull Bull # 
# 2007, 37:33-36.                                                              #
#                                                                              #
# VanRaden PM. Efficient methods to compute genomic predictions. J Dairy Sci   #
# 2008, 91:4414-4423.                                                          #
#                                                                              #
# Forni S, Aguilar I, Misztal I. Different genomic relationship matrices for   #
# single-step analysis using phenotypic, pedigree and genomic information.     #
# Genetics Selection Evolution 2011, 43:1.                                     #
#                                                                              #
# Strandén I, Christensen OF. Allele coding in genomic evaluation. Genetics    #
# Selection Evolution 2011, 43:25.                                             #
#                                                                              #
# Tier B, Meyer K, Ferdosi MH. Which genomic relationship matrix? Proc. Assoc. #
# Advmt. Breed. Genet. 2015, 21: 461-464.                                      #
#                                                                              #
# Toro MA, García-Cortés LA, Legarra A. A note on the rationale for estimating #
# genealogical coancestry from molecular markers. Genetics Selection Evolution #
# 2011, 43:27.                                                                 #
#                                                                              #
# Powell JE, Visscher PM, Goddard ME. Reconciling the analysis of IBD and IBS  #
# in complex trait studies. Nat Rev Genet 2010, 11:800-805.                    #
#                                                                              #
################################################################################
