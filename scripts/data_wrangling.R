# 13 nov 2019 repris le 27 novembre car les genes de resistance exclusifs aux plasmides etaient absents!
# le but de ce script est de preparer un dataset pour Steve Methot
# le dataset final contiendra des noms de souches renommes pour faciliter le travail de Steve dans SAS
# e.g Res13-Abat-PEA15-P4-01 devient APEA15P401
# de meme le dataset final contiendra des noms de gene de resistance renommes
#install.packages("tidyverse")
library(tidyverse)
library(ggplot2)

plasmids <- read.csv("data_271119/SAS_plasmids.1.csv", sep ="\t")

chr <- read.csv("data_271119/SAS_chr.1.csv", sep ="\t")

# Cette etape permet de bien merger les datasets plasmides et chr
all <- bind_rows(chr, plasmids)

# je vais l'examiner pour etre certain
write.csv(all, 'data_271119/all.csv', row.names = FALSE, quote = FALSE)

# le probleme c'est que le bind_row a mis des NAs et on vourdrait des 0
# on va couper all en deux parties (colonnes 1,14 et colonnes 15,143 )
# dans la partie 2 on va remplacer les NAs par des 0 et recoller ensuite

A <- all[,1:14]
B <- all[,15:143]

# On remplace les NA dans B (effectif)
B[is.na(B)] <- 0

# df1 notre dataset avec les valeurs de resistome pour les plasmides et les chr 
df1 <- cbind(A,B)




# Chargement de deux listes equivalentes de strain_names 
steve_strain_names <- read.csv("data/steve_strain_names.list", sep = '\t')
DPL_strain_names <- read.csv("data/DPL_strain_names.list", sep = '\t')

# Les deux listes sont dans le meme ordre, on pourra donc les joindre cote a cote
steve_DPL_strain_names <- cbind(steve_strain_names, DPL_strain_names)

# Pour faire un join, il faut une colonne commune, renommons la colonne DPL_strain_names par strain
steve_DPL_strain_names <- rename(steve_DPL_strain_names,  'strain' = DPL_strain_names)



# dyplr full join
df2 <- full_join(df1, steve_DPL_strain_names, by='strain')

# a ce moment-ci on a un dataset complet de resistome Res13, on produit un CSV!
write.csv(df2, 'data_271119/resistome_271119.csv', row.names = FALSE, quote = FALSE)


# maintenant on va travailler en 'par souche'
# pour y arriver, j'ai utilise mon script amr_collapsed_by_strain.pl sur le csv precedent
# je charge ici ce csv pour le reste de l'analyse


df_collapsed <- read.csv("data_271119/resistome_271119_collapsed.csv", sep ="\t")


# Bien mais comme notre d1 avait plus de lignes (244 souches == datast complet)
# que le nombre de souches avec lesquelles Steve devrait faire l'analyse (130)
# il est normal que le full join ait mis des valeurs NA dans la colonne steve_gene_names

# Beaucoup d'operations sequentielles ici
# 1. On enleve les lignes ayant NA dans la colonne steve_strain_names (car pas de donnnees de sensititre)
# 2. On enelve la colonne Strain avec les anciens noms.
# 3. On renomme la colonne steve_strain_names par Strain et on la place dans la premiere colonne
#df_collapsed_2 <- df_collapsed %>% filter(!is.na(steve_strain_names)) %>% select(-Strain) %>% rename('Strain' = steve_strain_names) %>% select (Strain, everything())
# edited april 2020
df_collapsed_2 <- df_collapsed %>% filter(!is.na(steve_strain_names)) %>% select(-steve_strain_names)


# avant d'exporter enlever les genes qui ne figurent plus dans le dataset
C <- df_collapsed_2[,1:14]
D <- df_collapsed_2[,15:143]

D2 <- D[,colSums(D) > 0]
df_collapsed_3 <- cbind(C,D2)



#write.csv(df_collapsed_3 , 'data_271119/resistome_steve_271119.csv', row.names = FALSE, quote = FALSE)
write.csv(df_collapsed_3 , 'data_271119/resistome_steve_210420.csv', row.names = FALSE, quote = FALSE)


# Les etapes suivants ont ete ajoutees durant la pandemie de COVI-19
# en teletravail pour ajouter 60 souches sensibles a l'analyse


# vendredi trace ; load final resistome file with holoenzyme coding!
sensitive_strains_amr_genes <- read.csv("data_160420/sensitive_strains_collapsed_AMR_final_resistome.csv")


# Ne pas commettre l'erreur de dire by="Strain", on fait un natural join ici et ca marche!
all_amr <- full_join(df_collapsed_3, sensitive_strains_amr_genes)

# on veut remplacer les NAs par des 0 dans la section des genes de resistance
# on va remplacer les NAs par des 0 et recoller ensuite

C <- all_amr[,1:14]
D <- all_amr[,15:139]

# On remplace les NA dans D (effectif)
D[is.na(D)] <- 0

# df1 notre dataset avec les valeurs de resistome pour les plasmides et les chr 
all_amr_2 <- cbind(C,D)




write.csv(all_amr_2, 'data_160420/resistome_final_200420.csv', row.names = FALSE, quote = FALSE)









