As ferramentas disponíveis nesta pasta foram desenvolvidas pelo Prof. Roberto Carvalheiro (roberto.carvalheiro@unesp.br)

# Controle de Qualidade de Genótipos 

Conteúdo prático da disciplina Análise de Dados Genômicos do Programa de Pós-Graduação em Genética e Melhoramento Animal da FCAV/UNESP de Jaboticabal - SP.

# Material disponível

* ## R Script - QC_genotipos.R
    * Faz o controle de qualidade de SNPs e amostras incluindo SNPs correlacionados e documenta a quantidade de SNPs excluídos por critérios de controle de qualidade
    
* ## R Script - vennQC.R
    * Script de funções auxiliares ao QC_genotipos.R
    
* ## R Script - QC_genotipos_qcf90.R
    * Faz o preparação dos genótipos e executa o controle de qualidade utilizando o software QCF90 (Yutaka Masuda)

* ## R Script - QC_genotipos_sexo.R
    * Faz a inferência do sexo de amostras (bovinos) com base no nível de heterozigosidade do cromossomo X

* ## Arquivos de Dados - genodat.Rdata e genodat2.Rdata
    * Objetos em R com exemplos de dados gentípicos para realização do controle de qualidade
    
* ## Arquivos de Dados - genoAN.Rdata
    * Objetos em R com exemplos de dados gentípicos para inferência do sexo pelo QC_genotipos_sexo.R
    
* ## Arquivo de mapa - SNP_Map_50K.txt
    * Mapa de referência dos SNPs.

* ## Arquivo de mapa - SNP_Map_GGPHDv2.txt
    * Mapa de referência dos SNPs para o QC_genotipos_sexo.R.
