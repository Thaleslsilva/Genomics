As ferramentas disponíveis nesta pasta foram desenvolvidas pelo Prof. Roberto Carvalheiro (roberto.carvalheiro@unesp.br)

# Desequilíbrio de Ligação, Parentesco Genômico e Estrutura Populacional

Conteúdo prático da disciplina Análise de Dados Genômicos do Programa de Pós-Graduação em Genética e Melhoramento Animal da FCAV/UNESP de Jaboticabal - SP

# Material disponível

* ## R Script - Desequilibrio_Ligacao.R
    * Compara o desequilíbrio de ligação (LD) local de duas regiões do genoma bovino que tiveram alteração expressiva de SNPs comparando os mapas (assembly) UMD3.1 e ARS-UCD1.2

* ## R Script - Parentesco_Genomico.R
    * Compara diferentes métodos de cálculo da matriz de parentesco genômico (G) utilizando dados reais
    
* ## R Script - Estrutura_Populacao.R
    * Avalia a presença de estrutura (estratificação) populacional entre genótipos de uma base de dados. Caso haja estratificação expressiva, pode-se recalcular a G e refazer o PC plot dentro de cada subgrupo.
    
* ## R Script - Utils_ParG_EstPop.R
    * Módulo de funções, auxiliar ao script Estrutura_Populacao.R

* ## Arquivo de Dados - genoHDaula3.Rdata
    * Objeto em R com exemplos de dados gentípicos para cálculo do LD via Desequilibrio_Ligacao.R
    
* ## Arquivo de Dados - geno1.Rdata
    * Objeto em R com dados gentípicos para cálculo da matriz de parentesco genômico (G) via Parentesco_Genomico.R

* ## Arquivo de Dados - genostrat.Rdata
    * Objeto em R com dados gentípicos para avaliação da presença de estrutura (estratificação) populacional via Estrutura_Populacao.R

* ## Arquivo de mapa - HDmap_UMD_ARS.rar
    * Mapa de referência dos SNPs (compactado) para uso com o Desequilibrio_Ligacao.R
