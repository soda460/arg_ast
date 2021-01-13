# 27 avril 2020
# Jean-Simon Brouard

install.packages('svglite')
install.packages('ggrepel')

library(tidyverse)
library(ggplot2)
library(svglite)
library(ggrepel)

#------------------------------------------------------------------------------
# Function str_to_quosure
#------------------------------------------------------------------------------
# This mini function convert string to quosure objects
# It is helpful in circumstances when using an external source to define
# dplyr arguments to be used in functions
# See https://github.com/r-lib/rlang/issues/116
# and https://dplyr.tidyverse.org/articles/programming.html
#------------------------------------------------------------------------------
str_to_quosure <- function(string) {
  out <- quo(!!sym(string))
  return (out)
}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Function pheno_geno
#------------------------------------------------------------------------------
# This function take a phenotype (AST) (e.g. SUL or FOT)
# and a list of associated genes gene_list <- c(quo(sul1), quo(sul2), quo(sul3))
# Note that the gene list is a list of quosure for ensuring that
# we can use dyplr's verbs inside functions
# see https://dplyr.tidyverse.org/articles/programming.html
# usage example: pheno_geno (df, SUL, sul_genes_list)
#------------------------------------------------------------------------------
pheno_geno <- function (df, AST, gene_list, antibio) {
  
  header_df <- c('freq_presence', 'freq_absence', 'var00', 'var11', 'var10', 'var01', 'pct_00', 'pct_11' )
  AST <- enquo(AST)
  AST_name <- paste0(quo_name(AST))  # Will be used to name the plot
  AST_name <- trimws(AST_name, "r")
  print (gene_list)
  
  # Apply our function to calculate var_11, var_00,  etc! to a list of genes implicated in the resistance
  out <- lapply(gene_list, calculate_concordance, df=df, AST_col=!!AST)
  
  df2 <- as.data.frame(do.call(rbind, out))  # bind the results in a dataframe
  
  names(df2) <- header_df  # name the columns dataframe
  
  # Get string objects from quosure?
  L <- lapply(gene_list, quo_name)
  
  df3 <- df2 %>% add_column(gene=unlist(L)) %>% select(gene, everything())
  
  filename <- paste("plots_2020/", AST_name, ".csv", sep = '')
  
  print(df3)
  
  write.csv(df3, filename, row.names = FALSE, quote = FALSE)

  # First save the plot in a variable
  p <- ggplot(df3, aes(x=pct_11, y=pct_00)) +
    xlim(c(0, 1)) +
    ylim(c(0, 1)) +
    geom_text_repel(aes(label=gene)) +
    geom_point(color = 'red') +
    theme_classic(base_size = 16) +
    ggtitle(paste("Resistance determinant to ", antibio, " (" , AST_name, ")", sep=""))

  # Then save the plot in a file
    ggsave(plot=p, paste("plots_2020/", AST_name, ".png", sep=""), width = 14, height = 7, dpi=600)

  # And finally, print the same plot in the RStudio console  
  # Alternative plot with the ggrepel library to avoid surimposed text
  ggplot(df3, aes(x=pct_11, y=pct_00)) +
    xlim(c(0, 1)) +
    ylim(c(0, 1)) +
    geom_text_repel(aes(label=gene)) +
    geom_point(color = 'red') +
    theme_classic(base_size = 16) +
    ggtitle(paste("Resistance determinant to ", antibio, " (" , AST_name, ")", sep=""))
}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Function calculate_concordance
#------------------------------------------------------------------------------
# This function take a dataframe containing genotypes phenotypes data
# and return a small dataframe with calculated PCT_00 and PCT_11 values
# Usage example: calculate_concordance(df, SUL, sul2)
#------------------------------------------------------------------------------
calculate_concordance <- function (df, AST_col, gene_col) {
  
  # Take the parameter AST_col and gene_col and just remember, R, that it is an expression.
  # Wn this case, the parameter of a function, denoting a column.
  # When using this NSE in function, prefix with !!
  AST <- enquo(AST_col)
  gene <- enquo(gene_col)
  
  # Give frequencies of the gene
  freq_presence <- df %>% filter (!!AST %in% c(0, 1) & !!gene != 0) %>% select (!!AST, !!gene) %>% count() %>% pull(n)
  freq_absence <- df %>% filter (!!AST %in% c(0, 1) & !!gene == 0) %>% select (!!AST, !!gene) %>% count() %>% pull(n)
  
  # Values corresponding to 00 (absence of resistance gene, absence of observed resistance)
  var00 <- df %>% filter (!!gene == 0 & !!AST == 0) %>% count() %>% pull(n)
  
  # Values corresponding to 11 (presence of resistance gene, presence of observed resistance)
  var11 <- df %>% filter (!!gene != 0 & !!AST == 1) %>% count() %>% pull(n)
  
  # Values corresponding to 10 (presence of resistance gene, absence of observed resistance)
  var10 <- df %>% filter (!!gene != 0 & !!AST == 0) %>% count() %>% pull(n)
  
  # Values corresponding to 01 (absence of resistance gene, presence of observed resistance)
  var01 <- df %>% filter (!!gene == 0 & !!AST == 1) %>% count() %>% pull(n)
  
  pct_00 = var00/freq_absence
  pct_11 = var11/freq_presence
  
  data <- c(freq_presence, freq_absence, var00, var11, var10, var01, pct_00, pct_11)
  return(data)
  print(data)
}
#------------------------------------------------------------------------------


# 1  - Loading the Phenotype dataset 
# Using the version where raw phenotypes (R, S, I) were interpreted as 1 or 0 by DPL

ast_carcass <- read.csv('data_140120/ast_carcass.csv', sep =',')
ast_feces_etc <- read.csv('data_140120/ast_feces_etc.csv', sep = ',')
ast <- ast_carcass %>% bind_rows(ast_feces_etc)
ast2 <- ast %>% rename(Strain = strain) # in the csv we have 'strain' not 'Strain' as in res dataset 


# 1.B - Adding the sensitive strains in the phenotype dataset

ast_3_meta <- read.csv('data_160420/metadata_sensitive_strains.csv')
ast_3_results <- read.csv('data_160420/AST_sensibles.csv')
ast_3_results <- ast_3_results %>% rename(Strain = strain)
ast_sensibles <- inner_join(ast_3_meta, ast_3_results, by = "Strain")
ast_final <- full_join(ast_sensibles, ast2)
ast_final <- ast_final %>% select(-Sample_ID)


# 2 - Loading the final genotype dataset
# Note that the merge of the first 131 strains and the new sensitive datasets
# has been done in the data_wrangling.R sript

res_final <- read.csv('data_160420/resistome_final_200420.csv')


# 3 - Merging the genotype and the phenotype datasets

df <- ast_final %>% inner_join(res_final, by = 'Strain')

# 3.1 Put in a variable the Antibiotic resistance genes actually found in the strains
ARG_found <- as.list(colnames(df))[-1:-47]   # cols 1 to 47 are other stuff or metadata


# 4 - Loading the resistance determinant dataset
res_det <- read.csv('data_210120/res_det.csv', sep =',')


# 5 --> Playing interactively with the dataset with a nice example, the SUL phenotype

# 5.1 Little sanity check : see the 4 Providencia strains that are resistant to MEM
# and that bear the IMP-7 resistance gene 
df %>% filter( MEM == 1) %>% select(Strain, genus, MEM, IMP_7)

df %>% filter (sul1 != 0) %>% select (Strain, SUL, sul1, sul2, sul3)

# See var_00, i.e. strains which are devoid of sul2 and sensitive to SUL (normal)
df %>% filter (sul2 == 0 & SUL == 0 ) %>% select (Strain, genus, SUL, sul1, sul2, sul3)

# See var_11, i.e the many strains carrying the sul2 gene and actually resistant to SUL (normal)
df %>% filter (sul2 != 0 & SUL == 1) %>% select (Strain, genus, SUL, sul1, sul2, sul3)

# See var_10, i.e. the strains carrying the sul2 gene and actually sensitive to SUL (weird)
df %>% filter (sul2 != 0 & SUL == 0) %>% select (Strain, genus, SUL, sul1, sul2, sul3)

# See var_01, i.e. the strains devoid of the sul2 gene, but actually resistant to SUL (other gene implicated)
df %>% filter (sul2 == 0 & SUL == 1) %>% select (Strain, genus, SUL, sul1, sul2, sul3)


# 6 Debuging and testing section, keep it for understanding it in the future :)

# 6.1 Performing analysis with home-made gene lists with quosures!
# SUL_list <- c(quo(sul2), quo(sul1), quo(sul3))

# Then calling the main function pheno_geno that notably produce graphs!
# pheno_geno(df, SUL, SUL_list)

# Alternatively, the quos (note the s) function take several arguments and will save time...
# MEM_list <- quos(IMP_7, ACT_2, ACT_17, ramA, Ecoli_acrR, Ecoli_marR, marA)
# pheno_geno(df, MEM, MEM_list)

# 6.2 Trying to construct these lists with R

# SUL <- res_det %>% filter (SUL == 1) %>% select (determinant)
# b <- as.list(as.character(SUL$determinant))
# Ce qui m'a permis de comprendre
#  b[[2]]
# quo(!!sym(b[[2]]))
# c <- map(b, str_to_quosure) # where str_to_quosure is a minimal functionn that use sym and !!
# pheno_geno(df, SUL, c)
# Got the same result!

# Another test with FOT
# FOT <- res_det %>% filter (FOT == 1) %>% select (determinant)
# b <- as.list(as.character(FOT$determinant))
# Problem can occur because some gene determinant are not present in the pheno-geno dataset
# Here we list the ARG found
# ARG_found <- as.list(colnames(df))[-1:-48]   # cols 1 to 47 are other stuff or metadata
# c <- keep(b, (b %in% ARG_found))
# d <- map(c, str_to_quosure) # map is a purrr fn
# pheno_geno(df, FOT, d, "fotantib")


# 7 - The true analysis with the 24 antibiotics!
#------------------------------------------------------------------------------
# 1 - AMP
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (AMP == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, AMP, d, "Ampicillin")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 2 - AUG
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (AUG == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, AUG, d, "Augmentin (amoxicillin-clavulanic acid)")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 3 - TZP
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (TZP == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, TZP, d, "Piperacillin-Tazobactam")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 4 - FAZ
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (FAZ == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, FAZ, d, "Cefazolin")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 5 - CEP
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (CEP == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, CEP, d, "Cephalothin")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 6 - CTX/FOT
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (FOT == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, FOT, d, "Cefotaxime")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 7 - POD
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (POD == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, POD, d, "Cefpodoxime")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 8 - TAZ
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (TAZ == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, TAZ, d, "Ceftazidime")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 9 - AXO
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (AXO == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, AXO, d, "Ceftriaxone")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 10 - FTC *
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (FTC == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, FTC, d, "Cefotaxime-Clavulanic acid")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 11 - CCV
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (CCV == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, CCV, d, "Ceftazidime-Clavulanic acid")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 12 - FEP
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (FEP == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, FEP, d, "Cefepime")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 13 - FOX
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (FOX == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, FOX, d, "Cefoxitin")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 14 - IMP
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (IMP == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, IMP, d, "Imipenem")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 15 - MER/MEM
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (MEM == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, MEM, d, "Meropenem")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 16 - SUL
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (SUL == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, SUL, d, "Sulfisoxazole")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 17 - SXT
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (SXT == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, SXT, d, "Trimethoprim-Sulfamethoxazole")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 18 - TET
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (TET == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, TET, d, "Tetracycline")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 19 - GEN
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (GEN == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, GEN, d, "Gentamicin")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 20 - STR
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (STR == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, STR, d, "Streptomycin")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 21 - AZI
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (AZI == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, AZI, d, "Azithromycin")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 22 - CIP
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (CIP == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, CIP, d, "Ciprofloxacin")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 23 - NAL
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (NAL == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, NAL, d, "Nalidixic acid")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 24 - CHL
#------------------------------------------------------------------------------
df2 <- res_det %>% filter (CHL == 1) %>% select (determinant)
b <- as.list(as.character(df2$determinant))
c <- keep(b, (b %in% ARG_found))
d <- map(c, str_to_quosure) # map is a purrr fn
pheno_geno(df, CHL, d, "Chloramphenicol")
#------------------------------------------------------------------------------


























