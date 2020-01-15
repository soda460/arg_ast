library(tidyverse)
library(ggplot2)
df1 <- read.csv("data_091219/ast.csv", sep ="\t")



# Chargement de deux listes equivalentes de strain_names 
steve_strain_names <- read.csv("data/steve_strain_names.list", sep = '\t')
DPL_strain_names <- read.csv("data/DPL_strain_names.list", sep = '\t')

# Les deux listes sont dans le meme ordre, on pourra donc les joindre cote a cote
steve_DPL_strain_names <- cbind(steve_strain_names, DPL_strain_names)

# Pour faire un join, il faut une colonne commune, renommons la colonne DPL_strain_names par strain
steve_DPL_strain_names <- rename(steve_DPL_strain_names,  'strain' = DPL_strain_names)



# dyplr full join et operation pour placer la nouvelle colonne en position 2
df2 <- full_join(df1, steve_DPL_strain_names, by='strain') %>% select (strain, steve_strain_names, everything())


# a ce moment-ci on a un dataset complet de resistome Res13, on produit un CSV!
write.csv(df2, 'data_091219/ast_edited.csv', row.names = FALSE, quote = FALSE)

