################################################################################ 
#
#  AN�LISE DE DADOS GEN�MICOS 
#
#  M�dulo auxiliar de fun��es para c�lculo de parentesco gen�mico e 
#  estratifica��o populacional
#
#  Atualiza��o: 13/10/2020             Dev.: Roberto Carvalheiro & Thales Silva
#
################################################################################


# Carregando o pacote snpStats
library(snpStats)


# FUN��O QUE CALCULA A MATRIZ DE PARENTESCO GEN�MICO PELO M�TODO 1 DE VANRADEN
VanRaden.1 <- function(x) {
  # M: Matriz de gen�tipos (nro alelos B)
  M  <- as(x,'numeric')      
  snpstat <- col.summary(x)
  # fB: Frequ�ncia do alelo B
  fB <- snpstat$RAF                 
  P  <- 2 * rep(1, nrow(M)) %*% t(fB)
  Z  <- M - P
  sum2pq  <- sum(2 * fB * (1 - fB))
  G  <- (Z %*% t(Z)) / sum2pq
  G.matrix <<- G
  cat("\n", "Estat�sticas da diagonal G", "\n")
  print(summary(diag(G)))
  # Plota a densidade de G
  plot(density(G), main = "Densidade do Coeficiente de Endogamia")
  
  # Converte os elementos da triangular superior da G em um dataframe (Glist)
  Glist <- G
  Glist[lower.tri(Glist,diag = FALSE)] <- NA
  Glist <- na.omit(data.frame(as.table(Glist)))
  colnames(Glist) <- c("Animal i", "Animal j", "Gij")
  G.list <<- Glist
  cat("\n", "Elementos da matriz triangular superior de G", "\n")
  head(Glist)
}



# FUN��O QUE CALCULA A MATRIZ DE PARENTESCO GEN�MICO PELO M�TODO 2 DE VANRADEN
VanRaden.2 <- function() {
  D  <- 1 / (ncol(M) * (2 * fB * (1 - fB)))
  G2 <- Z %*% (D * t(Z))
  # Retorna estat�sticas da diagonal de G2
  summary(diag(G2))
  
  # Converte os elementos da triangular superior da G2 em um dataframe (G2list)
  G2list <- G2
  G2list[lower.tri(G2list,diag = FALSE)] <- NA
  G2list <- na.omit(data.frame(as.table(G2list)))
  colnames(G2list) <- c("ANi", "ANj", "G2ij")
  head(G2list)
}



# FUN��O DE ESTRATIFICA��O DAS AMOSTRAS PELO MET�DO DE AN�LISE DE COMPONENTES 
# PRINCIPAIS (PCA)
PCA.geno <- function(G) {
  # Gerando autovalores e autovetores a partir da matriz G
  pc.all <- eigen(G, symmetric=TRUE) 
  # Salvando os autovetores
  pc     <- pc.all$vectors           
  # Salvando autovalores na base 100%
  rel.pc <<- pc.all$values * 100/sum(pc.all$values) 
  cat("\n", "Variabiliade explicada pelos autovalores em %", "\n")
  print(head(rel.pc))
  # Salvando os 2 primeiros autovetores
  pc.ord <<- cbind(pc[,1], pc[,2])    
  cat("\n", "Autovetores 1 e 2", "\n")
  head(pc.ord)
}



# FUN��O PARA IMPRIMIR GR�FICO PCPLOT
PC.plot <- function(pc) {
  plot(pc, pch = 20, main = "An�lise de Componentes Principais",
       xlab = paste("PC1: ", round(rel.pc[1], digits=2), "%", sep=""),
       ylab = paste("PC2: ", round(rel.pc[2], digits=2), "%", sep=""),
       cex = 2) 
}