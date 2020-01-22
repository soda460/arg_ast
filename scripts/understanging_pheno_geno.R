
library(tidyverse)
library(ggplot2)


# --> The Phenotype dataset <-- 
# Using the version where raw phenotypes (R, S, I) were interprend as 1 or 0 by DPL

ast_carcass <- read.csv('data_140120/ast_carcass.csv', sep =',')

ast_feces_etc <- read.csv('data_140120/ast_feces_etc.csv', sep = ',')

ast <- ast_carcass %>% bind_rows(ast_feces_etc)


# Keep the 131 isolates (no NAs in steve_strain_name)
ast <- ast %>% filter(steve_strain_names != 'NAs')

# Rename the column steve_strain_names with strain
ast <- ast %>% rename(Strain = steve_strain_names)


# Genotypes dataset
res <- read.csv('data_271119/resistome_steve_271119.csv', sep =',')

# Merge Resistance and AST datasets
df <- ast %>% inner_join(res, by = 'Strain' )

# Global variables used by functions
header_df <- c('freq_presence', 'freq_absence', 'var00', 'var11', 'var10', 'var01', 'pct_00', 'pct_11' )


#Playing interactively with the dataset

# Freq i.e. le gene est present
df %>% filter (sul1 != 0) %>% select (SUL, sul1, sul2, sul3)

# Quelle est la proportion de sensible et de resistant
df %>% filter (sul1 != 0 & SUL == 'S') %>% select (SUL, sul1, sul2, sul3)
df %>% filter (sul1 != 0 & SUL == 'R') %>% select (SUL, sul1, sul2, sul3)


# On veut eliminer les NAs et selectionner les cas ou on aurait plus d'unne allele du gene

# Donne la frquence des genes sul2 (avec des resultats AST (pas de NAs))
df %>% filter (SUL %in% c('R', 'S') & sul2 != 0) %>% select (SUL, sul1, sul2, sul3) 
df %>% filter (SUL %in% c('R', 'S') & sul2 == 0) %>% select (SUL, sul1, sul2, sul3) 

df %>% filter (SUL %in% c('R', 'S') & sul2 == 0) %>% select (SUL, sul1, sul2, sul3)

# lignes corresondant aux 00 (absence du gene, absence de la resistance)
df %>% filter (SUL == 'S' & sul2 == 0) %>% select (SUL, sul1, sul2, sul3)

# lignes corresondant aux 11 (presence du gene, presence de la resistance)
df %>% filter (SUL == 'R' & sul2 != 0) %>% select (SUL, sul1, sul2, sul3)

# lignes corresondant aux 10 (presence du gene, ABSENCE de la resistance)
df %>% filter (SUL == 'S' & sul2 != 0) %>% select (SUL, sul1, sul2, sul3)

# lignes corresondant aux 01 (absence du gene, et resistance)
df %>% filter (SUL == 'R' & sul2 == 0) %>% select (SUL, sul1, sul2, sul3)



# Construction of tons of gene lists with quosures!
SUL_list <- c(quo(sul2), quo(sul1), quo(sul3))
print(SUL_list)
MEM_list <- c(quo(IMP_7), quo(ACT_2), quo(ACT_17), quo(RAMA))

# Call the main function pheno_geno that call another function and produce graphs
pheno_geno(df, SUL, SUL_list)
pheno_geno(df, MEM, MEM_list)





# now in a function

# apply our function to calculate 11 00 etc to a list of genes implicated in the resistance to 1 AST
out <- lapply(SUL_list, calculate_concordance, df=df, AST_col=SUL)

# bind the results in a dataframe
SUL_df <- as.data.frame(do.call(rbind, out))
names(SUL_df) <- header_df

gene_list <- lapply(SUL_list, quo_name)

df2 <- SUL_df %>% add_column(gene=L) %>% select(gene, everything())

# Call the plot
ggplot(df2, aes(x=pct_11, y=pct_00)) +
  geom_point(shape=21, fill='blue') +
  xlim(c(0, 1)) +
  ylim(c(0, 1)) +
  geom_text(aes(label=gene), size=3, hjust=1, vjust=-0.5)








##########################################################################0
# Function pheno_geno
##########################################################################0
# This function take a phenotype (AST) (e.g. SUL or FOT)
# and a list of associated genes gene_list <- c(quo(sul1), quo(sul2), quo(sul3))
# Note that the gene list is a list of quosure for ensuring that
# we can use dyplr's verbs inside functions
# see https://dplyr.tidyverse.org/articles/programming.html
# usage example: pheno_geno (df, SUL, sul_genes_list)
##########################################################################0
pheno_geno <- function (df, AST, gene_list) {
  
  AST <- enquo(AST)

  print('ok')
  print(AST)

  
  enquos(gene_list)
  print(gene_list)
  
  
  # apply our function to calculate 11 00 etc to a list of genes implicated in the resistance
  out <- lapply(gene_list, calculate_concordance, df=df, AST_col=!!AST)
  
  # bind the results in a dataframe
  SUL_df <- as.data.frame(do.call(rbind, out))
  names(SUL_df) <- header_df
  
  
  # 
  L <- lapply(gene_list, quo_name)
  
  df2 <- SUL_df %>% add_column(gene=L) %>% select(gene, everything())
  
  print(df2)
  # Call the plot
  ggplot(df2, aes(x=pct_11, y=pct_00)) +
    geom_point(shape=21, fill='blue') +
    xlim(c(0, 1)) +
    ylim(c(0, 1)) +
    geom_text(aes(label=gene), size=3, hjust=1, vjust=-0.5)
}
##########################################################################0









##########################################################################0
# Function calculate_concordance
##########################################################################0
# This function take a dataframe containing genotypes phenotypes data
# and return a small dataframe with calculated PCT_00 and PCT_11 values
# usage example: calculate_concordance(df, SUL, sul2)
##########################################################################0
calculate_concordance <- function (df, AST_col, gene_col) {
  
  # Take the parameter AST_col and gene_col and just remember, R, that it is an expression.
  # in this case, the parameter of a function, denoting a column.
  # When we want to use this NSE in function, prefix with !!
  AST <- enquo(AST_col)
  gene <- enquo(gene_col)
  
  # Give frequencies of the gene
  freq_presence <- df %>% filter (!!AST %in% c('R', 'S') & !!gene != 0) %>% select (!!AST, !!gene) %>% count() %>% pull(n)
  freq_absence <- df %>% filter (!!AST %in% c('R', 'S') & !!gene == 0) %>% select (!!AST, !!gene) %>% count() %>% pull(n)
  
  
  # Values corresponding to 00 (absence of resistance gene, absence of observed resistance)
  var00 <- df %>% filter (!!gene == 0 & !!AST == 'S') %>% count() %>% pull(n)
  
  # Values corresponding to 11 (presence of resistance gene, presence of observed resistance)
  var11 <- df %>% filter (!!gene != 0 & !!AST == 'R') %>% count() %>% pull(n)
  
  # Values corresponding to 10 (presence of resistance gene, absence of observed resistance)
  var10 <- df %>% filter (!!gene != 0 & !!AST == 'S') %>% count() %>% pull(n)
  
  # Values corresponding to 01 (absence of resistance gene, presence of observed resistance)
  var01 <- df %>% filter (!!gene == 0 & !!AST == 'R') %>% count() %>% pull(n)
  
  pct_00 = var00/freq_absence
  pct_11 = var11/freq_presence
  
  data <- c(freq_presence, freq_absence, var00, var11, var10, var01, pct_00, pct_11)
  
  return(data)
  print(data)
}
##########################################################################0
























