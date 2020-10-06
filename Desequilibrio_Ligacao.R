################################################################################ 
#
#  ANÁLISE DE DADOS GENÔMICOS 
#
#  Comparar o Desequilíbrio de Ligação (LD) local de duas regiões do genoma que 
#  tiveram alteração expressiva de SNPs comparando os mapas (assemblies) de 
#  referência UMD3.1 e ARS-UCD1.2
#  As duas regiões são: Chr.21:59-61Mb e Chr.27:0-2Mb
#
#  Versão 1.0                                       Dev.: Roberto Carvalheiro
#  Atualização: 06/10/2020                        Atual.: Thales de Lima Silva
#
################################################################################


# 1. PREPARANDO O AMBIENTE #####################################################
# Limpando workspace
rm(list=ls()) 
# Não deixar o R converter automaticamente caracter em fator
options(stringsAsFactors=F) 

## Para instalar o pacote "LDHeatmap" baixar e instalar o Rtools 4.0: 
## https://cran.r-project.org/bin/windows/Rtools/

## Instalar o devtools
## install.packages("devtools")

## Instalar a última versão do LDheatmap a partir do GitHub 
## devtools::install_github("SFUStatgen/LDheatmap")

## Caso haja alguma incompatibilidade entre a versão do R instalada e a versão
## do LDHeatmap disponível no GitHub, versões anteriores deste pacote podem ser
## encontradas na página do archive do CRAN:
## https://cran.r-project.org/src/contrib/Archive/LDheatmap/

# Carregando os pacotes
library(LDheatmap)
library(snpStats)

# Definindo o diretório de trabalho
setwd("C:/Users/Thales/Google Drive/UNESP_GMA/AnlDadosGen/aula3/Ex1")


# 2. CARREGANDO OS DADOS #######################################################
# Leitura dos genótipos
load("genoHDaula3.Rdata")
dim(genoHDaula3)

# Leitura do mapa
mapa <- read.table("HDmap_UMD_ARS.txt", header = T, sep = "\t")
str(mapa)
head(mapa)


# 3. LIMPEZA E PRÉ-PROCESSAMENTO ###############################################
# Filtranto SNPs do arquivo de genótipos no mapa
aux    <- match(colnames(genoHDaula3), mapa$snpname)
which(is.na(aux))
mapaG  <- mapa[aux, 1:5]
chrU_A <- paste(mapaG$chrU, mapaG$chrA, sep="_")
table(chrU_A)

# Separando SNPs de interesse em listas
snplist1<-which(mapaG$chrU=="21" & mapaG$posU/1000000>=59 & mapaG$posU/1000000<=61)
snplist2<-which(mapaG$chrA=="21" & mapaG$posA/1000000>=59 & mapaG$posA/1000000<=61)
snplist3<-which(mapaG$chrU=="27" & mapaG$posU/1000000>=0 & mapaG$posU/1000000<=2)
snplist4<-which(mapaG$chrA=="27" & mapaG$posA/1000000>=0 & mapaG$posA/1000000<=2)

# Criando variável auxiliar no mapa para montar os heatmaps 
mapaG$w <- 0
mapaG$w[snplist1] <- 1
mapaG$w[snplist2] <- 2
mapaG$w[snplist3] <- 3
mapaG$w[snplist4] <- 4
table(mapaG$w)


# 4. GERANDO GRÁFICOS COM LD NAS ANELAS DE INTERESSE ###########################

################################################################################
# Comparando o desequilíbrio de ligação na janela do cromossomo 21 nas duas    #
# referências                                                                  #
################################################################################
# LD janela 1: Chr.21:59-61Mb (UMD3.1) 
LDheatmap(genoHDaula3[, mapaG$w==1], mapaG$posU[mapaG$w == 1], flip = TRUE, 
          color = heat.colors(20),
          title = paste("LD heatmap UMD3.1 Chr.", 
                        unique(mapaG$chrU[mapaG$w == 1]), ":", 
                        min(mapaG$posU[mapaG$w == 1]),"-", 
                        max(mapaG$posU[mapaG$w == 1]), 
                        "  (N_snps = ", length(which(mapaG$w == 1)), ")", sep = ""),
          add.map = T)


# LD janela 2: Chr.21:59-61Mb (ARS-UCD1.2) 
LDheatmap(genoHDaula3[, mapaG$w==2], mapaG$posA[mapaG$w == 2], flip = TRUE,
          color = heat.colors(20),
          title = paste("LD heatmap ARS-UCD1.2 Chr.", 
                        unique(mapaG$chrA[mapaG$w == 2]), ":", 
                        min(mapaG$posA[mapaG$w == 2]), "-",
                        max(mapaG$posA[mapaG$w == 2]), 
                        "  (N_snps = ", length(which(mapaG$w == 2)), ")", sep = ""),
          add.map=T)


################################################################################
# Comparando o desequilíbrio de ligação na janela do cromossomo 27 nas duas    #
# referências                                                                  #
################################################################################
# LD janela 3: Chr.27:0-2Mb (UMD3.1) 
LDheatmap(genoHDaula3[, mapaG$w == 3], mapaG$posU[mapaG$w == 3], flip = TRUE,
          color = heat.colors(20),
          title = paste("LD heatmap UMD3.1 Chr.",
                        unique(mapaG$chrU[mapaG$w == 3]), ":",
                        min(mapaG$posU[mapaG$w == 3]),"-",
                        max(mapaG$posU[mapaG$w == 3]),
                        "  (N_snps = ", length(which(mapaG$w == 3)), ")", sep = ""),
          add.map=T)


# LD janela 4: Chr.27:0-2Mb (ARS-UCD1.2) 
LDheatmap(genoHDaula3[, mapaG$w == 4], mapaG$posA[mapaG$w == 4], flip = TRUE,
          color = heat.colors(20),
          title = paste("LD heatmap ARS-UCD1.2 Chr.",
                        unique(mapaG$chrA[mapaG$w == 4]), ":",
                        min(mapaG$posA[mapaG$w == 4]), "-",
                        max(mapaG$posA[mapaG$w == 4]),
                        "  (N_snps = ", length(which(mapaG$w == 4)), ")", sep = ""),
          add.map=T)
