################################################################################ 
#
#  ANÁLISE DE DADOS GENÔMICOS 
#
#  Controle de qualidade de SNPs e amostras incluindo SNPs correlacionados e 
#  report do QC (documenta quantos SNPs foram excluídos por critérios de 
#  controle de qualidade)
#
#  Versão 2.0                                       Dev.: Roberto Carvalheiro
#  Atualização: 29/09/2020                        Atual.: Thales de Lima Silva
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
setwd("C:/Users/Thales/Google Drive/UNESP_GMA/AnlDadosGen/aula2/ex1")


# 2. CARREGANDO OS DADOS #######################################################
# Objeto R salvo da leitura dos arquivos de genótipos     
load("genodat.Rdata")   
dim(genodat)
head(genodat)

# Objeto adicional com outros genótipos
load("genodat2.Rdata")   
dim(genodat2)
head(genodat)

# Concatenando genótipos em um único objeto - ***importante***: os objetos devem
# ter os mesmos SNPs e as colunas igualmente ordenadas
genotipo<-rbind(genodat,genodat2)
dim(genotipo)

# Removendo objetos anteriores para liberar memória
rm(genodat, genodat2)

# Lendo "SNP map"
snpmap <- read.table("SNP_Map_50K.txt", sep = "\t", header = T)
str(snpmap)

mapa <- data.frame(snp = snpmap$Name, chr = snpmap$Chromosome, 
                   pos = snpmap$Position)
head(mapa)
table(mapa$chr)


# 3. LIMPEZA E PRÉ-PROCESSAMENTO ###############################################
# filtrar apenas autossomos
auto <- which(!is.na(match(mapa$chr, 1:29)))
mapa <- mapa[auto,]
table(mapa$chr)

# Documentando alterações nos SNPS com objeto Report
(report <- data.frame(situation = "all snps", n = nrow(snpmap)))
auxn <- data.frame(situation = "non-autosome", n = nrow(snpmap) - nrow(mapa))
(report <- rbind(report, auxn))

# Ordenando por cromossomo e posicao
ordena <- order(as.numeric(mapa$chr), mapa$pos) 
mapa <- mapa[ordena,]
head(mapa)
tail(mapa)

# Removendo SNPs com mesma posicao
dupli <- which(duplicated(data.frame(mapa$chr, mapa$pos)))
mapa[dupli,]
if (length(dupli) > 0) mapa <- mapa[-dupli,]

auxn <- data.frame(situation = "same position", n = length(dupli))
(report <- rbind(report, auxn))

# Filtrando SNPs no genotipo
dim(genotipo)
head(colnames(genotipo))

snp.ok <- match(mapa$snp, colnames(genotipo))
genotipo <- genotipo[,snp.ok]
dim(genotipo)

# Conferindo estatísticas por SNP
snpstat <- col.summary(genotipo)
summary(snpstat)

# Conferindo estatísticas por amostra
samplestat <- row.summary(genotipo)
summary(samplestat)


# 4. APLICANDO O CONTROLE DE QUALIDADE #########################################
# definindo parametros para QC dos SNPs
maf.thr     <-0.02 # MAF mínima
callsnp.thr <-0.90 # Call-rate mínimo (SNP)
hwe.thr     <-1e-5 # HWE mínimo (p valor para o z-score)
r2.thr      <-1.01 # r2 máximo com outro SNP dentro da janela (valor maior que 
                   ## 1 desativa este critério)
window.size <-100  # Tamanho da janela para calcular r2
N_high_r2   <-0    # Objeto para armazenar numero de SNPs que serao descartados
                   ## por alto r2

# definindo parametros para QC das amostras
callsamp.thr <-0.97 # Call-rate mínimo (amostra)

niter  <- 5         # numero de iteracoes do QC
genoQC <- genotipo  # objeto para armazenar genotipos pos QC

for (iter in 1:niter) {
  # SNP QC
  nsnps.old <- dim(genoQC)[2]
  snp.qc <- col.summary(genoQC)
  cols.ok <- which((snp.qc$MAF > maf.thr) & (snp.qc$Call.rate > callsnp.thr) &
                   ((1-pnorm(abs(snp.qc$z.HWE))) > hwe.thr))
  genoQC <- genoQC[,cols.ok]
  
  # QC de SNPs altamente correlacionados
  ref.map <- match(colnames(genoQC), mapa$snp)
  mapa <- mapa[ref.map,]
  mafes <- snp.qc$MAF[cols.ok]
  drop <- array(1, nrow(mapa))
  for (i in 1:29) {
    ind.cromo <- which(as.character(mapa$chr) == i)
    inew <- order(mapa$pos[ind.cromo])
    ind.cromo <- ind.cromo[inew]
    ####### ld é a funcao do snpStats para calcular r²
    ld1  <- ld(genoQC[,ind.cromo], depth = window.size, stats = c("R.squared"))
    ldmat <- summary(ld1)
    ind.cut <- which(ldmat[,3] > r2.thr)
    ir1 <- ldmat$i[ind.cut]
    ir2 <- ldmat$j[ind.cut]
    maf.teste1 <- mafes[ind.cromo[ir1]]
    maf.teste2 <- mafes[ind.cromo[ir2]]
    drop.r2 <- ir1
    ####### ficar com o SNP com a maior MAF, do par altamente correlacionado
    isub <- which(maf.teste2 > maf.teste1) 
    drop.r2[isub] <- ir2[isub]  
    drop.r2 <- unique(drop.r2)
    drop[ind.cromo[drop.r2]] <- 0
  }
  cols.ok <- which(drop == 1)
  N_high_r2 <- N_high_r2 + length(which(drop == 0))
  genoQC <- genoQC[,cols.ok]
  nsnps.new <- dim(genoQC)[2]
  ####### parar o processo se nao houver nenhum corte de SNPs
  if (nsnps.new == nsnps.old) {break}
  
  # Controle por filtro ####################################
  crit3 <- length(which((snp.qc$MAF <= maf.thr) & 
                          (snp.qc$Call.rate <= callsnp.thr) & 
                          ((1-pnorm(abs(snp.qc$z.HWE)))<= hwe.thr)))
  maf_cr <- length(which((snp.qc$MAF <= maf.thr) & 
                           (snp.qc$Call.rate <= callsnp.thr))) - crit3
  maf_hwe <- length(which((snp.qc$MAF <= maf.thr) & 
                            ((1-pnorm(abs(snp.qc$z.HWE))) <= hwe.thr))) - crit3
  cr_hwe <- length(which((snp.qc$Call.rate <= callsnp.thr) & 
                           ((1-pnorm(abs(snp.qc$z.HWE))) <= hwe.thr))) - crit3
  maf.crit <- length(which((snp.qc$MAF <= maf.thr))) - maf_cr - maf_hwe - crit3
  hwe.crit <- length(which((1-pnorm(abs(snp.qc$z.HWE))) <= hwe.thr)) - maf_hwe - cr_hwe - crit3
  cr.crit <- length(which((snp.qc$Call.rate <= callsnp.thr))) - cr_hwe - maf_cr - crit3
  auxn <- data.frame(situation = "niter", n = iter)
  report <- rbind(report, auxn)
  auxn <- data.frame(situation = "MAF only", n = maf.crit)
  report <- rbind(report, auxn)
  auxn <- data.frame(situation = "HWE only",n = hwe.crit)
  report <- rbind(report, auxn)
  auxn <- data.frame(situation = "SNP call rate only", n = cr.crit)
  report <- rbind(report, auxn)
  auxn <- data.frame(situation = "both MAF and HWE", n = maf_hwe)
  report <- rbind(report, auxn)
  auxn <- data.frame(situation = "both MAF and call rate", n = maf_cr)
  report <- rbind(report, auxn)
  auxn <- data.frame(situation = "both call rate and HWE", n = cr_hwe)
  report <- rbind(report, auxn)
  auxn <- data.frame(situation = "both MAF, call rate and HWE", n = crit3)
  report <- rbind(report, auxn)
  auxn <- data.frame(situation = "high r2", n = N_high_r2)
  report <- rbind(report, auxn)
  ####################################
  
  # Amostra do QC
  sample.qc <- row.summary(genoQC)
  rows.ok   <- which(sample.qc$Call.rate > callsamp.thr) 
  genoQC    <- genoQC[rows.ok,]
}

iter
report
dim(genotipo)
dim(genoQC)

# Conferindo estatísticas pós QC
snp.qc <- col.summary(genotipo)
sample.qc <- row.summary(genotipo)
summary(data.frame(snp.qc$Call.rate,snp.qc$MAF,snp.qc$z.HWE))
summary(sample.qc$Call.rate)

fsnp.qc <- col.summary(genoQC)
fsample.qc <- row.summary(genoQC)
summary(data.frame(fsnp.qc$Call.rate,fsnp.qc$MAF,fsnp.qc$z.HWE))
summary(fsample.qc$Call.rate)

# Plotando resultados comparados 
par(mfrow = c(1,2))
boxplot(snp.qc$Call.rate,main = "SNP Call Rate")
boxplot(fsnp.qc$Call.rate,main = "SNP Call Rate after QC")

boxplot(snp.qc$MAF,main="MAF")
boxplot(fsnp.qc$MAF,main="MAF after QC")

boxplot(snp.qc$z.HWE,main="zHWE")
boxplot(fsnp.qc$z.HWE,main="zHWE after QC")


# 5. SALVANDO O OBJETO COM GENÓTIPOS PÓS QC E GERANDO RELATÓRIO ################
save(genoQC, file = "genoQC.Rdata")

report
write.table(report, "reportQC.txt", quote = F, row.names = F, sep = "\t")

pro.venn <- by(report$n, report$situation, sum)
source("C:/Users/Thales/Google Drive/UNESP_GMA/AnlDadosGen/aula2/ex1/vennQC.R")
png(file = "C:/Users/Thales/Google Drive/UNESP_GMA/AnlDadosGen/aula2/ex1/SNP_QC_report.png",
    width = 1000, height = 800)
venndia(nA = pro.venn[which(names(pro.venn) == "MAF only")], 
        nB = pro.venn[which(names(pro.venn) == "HWE only")],
        nC = pro.venn[which(names(pro.venn) == "SNP call rate only")],
        nAB = pro.venn[which(names(pro.venn) == "both MAF and HWE")],
        nAC = pro.venn[which(names(pro.venn) == "both MAF and call rate")],
        nBC = pro.venn[which(names(pro.venn) == "both call rate and HWE")],
        nABC =pro.venn[which(names(pro.venn) == "both MAF, call rate and HWE")],
        nD = pro.venn[which(names(pro.venn) == "high r2")])
dev.off()