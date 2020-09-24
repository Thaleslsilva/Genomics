################################################################################ 
#
#  ANÁLISE DE DADOS GENÔMICOS 
#
#  Calcular o número de SNPs do painel HD que tiveram a posição alterada para 
#  outro cromossomo ao contrastar a posição das referências/assemblies UMD3.1 
#  e ARS-UCD1.2
#
#  Atualização: 24/09/2020                          Dev.: Thales de Lima Silva
#
################################################################################


# 1. PREPARANDO O AMBIENTE #####################################################
# Limpando workspace
rm(list=ls()) 
# Não deixar o R converter automaticamente caracter em fator
options(stringsAsFactors=F) 

# Carregando pacotes
library(gtools)
library(ggplot2)

# Definindo o diretório onde estão os arquivos
setwd("C:/Users/Thales/Google Drive/UNESP_GMA/AnlDadosGen/aula1/ex4")


# 2. CARREGANDO E VISUALIZANDO OS DADOS ########################################
# Criando objeto com o nome dos arquivos que serão lidos
filenames <- dir(pattern = "HD") # arquivos iniciando com "SNP"
head(filenames)

# SNP_Map_HD
map2 <- read.table(filenames[2], header = T, sep = "\t")
str(map2)
# 9913_ARS1.2_777962_HD_marker_name_180910
map1 <- read.table(filenames[1], header = F, sep = "\t",
                   col.names = c("Chromosome", "Name", "Y", "Position"))
str(map1)


# 3. EDIÇÃO E PRÉ-PROCESSAMENTO ################################################

## 9913_MAP_marker_names_180910.README
## The chromosomes are coded as integers according to the PLINK specification.
## Autosomes 1..29, X=30, Y=31, PAR=32, MT=33, ChrUn is coded as 0 (zero). 

## Nota: PAR = pseudoautosomal region of X

# Editando arquivos de mapas
map2.1 <- map2[,-c(1,5:9)]
names(map2.1)[1:3] <- c("Name", "Chr_UMD", "Pos_UMD")

map1.1 <- map1[,-c(3)]
names(map1.1)[1:3] <- c("Chr_ARS", "Name", "Pos_ARS")
## A coluna que indica o cromossomo no mapa 1.1 deve ser do tipo "caracter" para
## poder ser comparada com o mapa 2.1
map1.1$Chr_ARS = as.character(as.numeric(map1.1$Chr_ARS))
str(map1.1)

### Alterando a identificação dos cromossomos para parear como o mapa 2.1
table(map2.1$Chr_UMD)
table(map1.1$Chr_ARS)
map1.1$Chr_ARS[map1.1$Chr_ARS == "30"] <- "X"
map1.1$Chr_ARS[map1.1$Chr_ARS == "32"] <- "X"
map1.1$Chr_ARS[map1.1$Chr_ARS == "31"] <- "Y"
map1.1$Chr_ARS[map1.1$Chr_ARS == "33"] <- "MT"


# 4. COMPARAR MAPAS/ASSEMBLIES #################################################
# Unindo os dados em uma única tabela
table_maps <- merge(map1.1, map2.1, by = "Name", all = T)
View(table_maps)

length(which(table_maps$Pos_UMD<0))
# For map files, if there is not a REF/ALT identified then the position is set 
# to negative.

## 4.1 SNPs que mudaram de cromossomo ##########################################
dif <- which(table_maps$Chr_ARS != table_maps$Chr_UMD)
dif_Chr <- as.data.frame(table(paste(table_maps$Chr_ARS[dif],
                                     table_maps$Chr_UMD[dif], sep = "_")))
dif_Chr
dif_Chr[order(dif_Chr$Freq),]
sum(dif_Chr$Freq)

# RESPOSTA: #####################################################################
# 6.477 é o número de SNP's do painel HD que tiveram a posição alterada para    #
# outro cromossomo, ao contrastar a posição das referências UMD3.1              #
# (SNP_Map_HD.txt) e ARS- UCD1.2 (9913_ARS1.2_777962_HD_marker_name_180910.map) #
#################################################################################

## 4.2 SNPs que não mudaram de cromossomo ######################################
eq_Chr <- as.data.frame(table(paste(table_maps$Chr_ARS[-dif],
                                    table_maps$Chr_UMD[-dif], sep = "_")))
eq_Chr$Perc <- eq_Chr$Freq/table(map2.1$Chr_UMD)*100
eq_Chr
sum(eq_Chr$Freq)

# RESPOSTA: #####################################################################
# 771.485 SNP's do painel HD não tiveram a posição alterada para outro          #
# cromossomo, ao contrastar a posição das referências UMD3.1 (SNP_Map_HD.txt) e #
# ARS- UCD1.2 (9913_ARS1.2_777962_HD_marker_name_180910.map)                    #
#################################################################################

## 4.3 Mudanças de posição dentro de cromossomo ################################
difMb <- (table_maps$Pos_ARS[-dif] - abs(table_maps$Pos_UMD[-dif]))/1000000
summary(difMb)
boxplot(difMb)
plot(density(difMb))

maps2 <- table_maps[-dif,]
summary(maps2)

for (var in mixedsort(unique(maps2$Chr_UMD[maps2$Chr_UMD != "0"]))) {
  dev.new()
  print(ggplot(maps2[maps2$Chr_UMD == var,], 
               aes(x = Pos_UMD / 1000000, y = abs(Pos_ARS)/1000000)) 
        + geom_point() 
        + geom_smooth(method = lm, se = F) 
        + labs(title = paste("Chromosome", var), x = "UMD position (Mb)", 
               y = "ARS position (Mb)"))
}

# 5. SALVANDO TABELA COM SNPS QUE TIVERAM A POSIÇÃO ALTERADA  ##################
snp_pos_altr <- table_maps[table_maps$Chr_ARS != table_maps$Chr_UMD, ]

write.table(snp_pos_altr, file = "SNPs_alterados.txt", quote = F, row.names = F,
            col.names = T, sep = "\t")
