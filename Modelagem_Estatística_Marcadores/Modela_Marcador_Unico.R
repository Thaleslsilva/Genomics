################################################################################ 
#
#  ANÁLISE DE DADOS GENÔMICOS 
#
#  Como modelar estatisticamente os marcadores?
#  Marcador único - diferentes modelos
#
#  Versão 1.0                                       Dev.: Roberto Carvalheiro
#  Atualização: 19/01/2021                        Atual.: Thales de Lima Silva
#
################################################################################


# PREPARANDO O AMBIENTE ########################################################
rm(list=ls())

setwd("C:/Users/Thales/Google Drive/UNESP_GMA/ADGen_2020/aula4-Modelagem_Estatistica_Marcadores/lab1")

# Ler arquivo e conferir variáveis
dados<-read.table("aula4.dat",header=T)
str(dados)
View(dados)
table(dados$geno)
summary(dados$pesod)


# Modelo 1: Estimar Efeitos Genotípicos (g) ####################################
### Pressupondo mu=0
# Cria-se um vetor com os dados fenotípicos
y  <- dados$pesod
y
# Cria a matriz de incidência dos genótipos e concatena com o vetor de médias
Q  <- model.matrix(~as.factor(geno)+0,dados) 
Q
# Cria a mão esquerda das equações
QQ <- crossprod(Q)
QQ
# Cria a mão direita das equações
Qy <-crossprod(Q,y)
Qy
# Resolvendo o sistema de equações e calculando o vetor de solução dos efeitos
sol1a <- solve(QQ,Qy)
row.names(sol1a) <- c("mu+gAA","mu+gAB","mu+gBB")
sol1a

### Pressupondo gAA=0
# Cria a matriz de incidência dos genótipos e concatena com o vetor de médias
Q  <- (model.matrix(~as.factor(geno),dados)) 
Q
# Cria a mão esquerda das equações
QQ <- crossprod(Q)
QQ
# Cria a mão direita das equações
Qy <- crossprod(Q,y)
Qy
# Resolvendo o sistema de equações e calculando o vetor de solução dos efeitos
sol1b <- solve(QQ,Qy)
row.names(sol1b) <- c("mu+gAA","gAB-gAA","gBB-gAA")
sol1b


# Modelo 2: Estimar Efeitos Aditivo (a) e de Dominância (d) ####################
# Cria a matriz de incidência com efeitos aditivoa e de dominância
dados$geno
A <- dados$geno-1
A
D <- as.numeric(dados$geno==1)
D
Z <- cbind(1,A,D)
Z
# Cria a mão esquerda das equações
ZZ <- crossprod(Z)
ZZ
# Cria a mão direita das equações
Zy <- crossprod(Z,y)
Zy
# Resolvendo o sistema de equações e calculando o vetor de solução dos efeitos
sol2 <- solve(ZZ,Zy)
row.names(sol2)<-c("mu","a","d")
sol2


# Modelo 3: Regressão pelo número de alelos B ##################################
# Cria a matriz de incidência e concatena com o vetor de médias
dados$geno
W<-cbind(1,dados$geno)
W
# Cria a mão esquerda das equações
WW<-crossprod(W)
WW
# Cria a mão direita das equações
Wy<-crossprod(W,y)
Wy
# Resolvendo o sistema de equações e calculando o vetor de solução dos efeitos
sol3<-solve(WW,Wy)
row.names(sol3)<-c("alfa","beta")
sol3


# Modelo 4: Média e Efeito de Substituição Alélica #############################
# Cria a matriz de incidência e concatena com o vetor de médias
dados$geno-1
X<-cbind(1,dados$geno-1)
X
# Cria a mão esquerda das equações
XX<-crossprod(X)
XX
# Cria a mão direita das equações
Xy<-crossprod(X,y)
Xy
# Resolvendo o sistema de equações e calculando o vetor de solução dos efeitos
sol4<-solve(XX,Xy)
row.names(sol4)<-c("mu","tal")
sol4
class(row.names(sol4))
