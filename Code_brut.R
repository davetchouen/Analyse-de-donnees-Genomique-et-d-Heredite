
### QUESTION 1

install.packages('kinship2')

library(kinship2)

# Charger les données du fichier pedigree.tsv
pedigree <- read.table("pedigree.tsv", header = TRUE, sep = "\t")

# Calculer les coefficients de parenté
K <- kinship(id = pedigree$ID,dadid = pedigree$FATHER,momid = pedigree$MOTHER)

# Extraire les coefficients d'apparentement de la diagonale
phi <- diag(K)

# Calculer les coefficients de consanguinité pour chaque individu
f<-(phi*2)-1

# Compter le nombre d'individus avec phi > 0.25
nb_individus_f_025 <- sum(f >= 0.25)
nb_individus_f_025 ## 28 individus ont un coeficient d'apparentement supperieur ou egale a 0.25


##QUESTION 2
install.packages("gaston")
library(gaston)

# Modifier les options par défaut de Gaston
options(gaston.autosomes = 26, gaston.chr.x = 27, gaston.chr.y = NA)

# Lire les données à partir des fichiers .bed, .bim et .fam

data_bed <- read.bed.matrix("20180209_SoayPlates1-83")


length(data_bed@ped$famid)## nous avons 7268 individus 

length(data_bed@snps$id) ## on a 39176 snps
table(data_bed@snps$chr)

proportion_low_maf <- sum(data_bed@snps$maf < 0.05) / length(data_bed@snps$id)

proportion_low_maf ## proportion de la maf est de 0.09322034

callratemin<- min(data_bed@ped$callrate)
callratemin # le callrate minimal des individus est -1.731523 

callratemin_snp <- min(data_bed@snps$callrate)
callratemin_snp ## le callrate minimal des snp est 0.5939736


selected_data <- select.snps(data_bed, callrate >= 0.95)
selected_data <- select.snps(selected_data, maf >= 0.05)


#QUESTION 3

hwe_results <- set.hwe(selected_data)

#vérification de la structure de hwe_results
head(hwe_results@snps$hwe)
tail(hwe_results@snps$hwe)

#premiere methode 
n1 <- ncol(hwe_results)
plot( -log10(ppoints(n1)), -log10(sort(hwe_results@snps$hwe)), pch=".", cex=2)

#deuxiemme methode 
qqplot.pvalues(hwe_results@snps$hwe,  col.abline = "green")

#En conclusion, le QQ-plot des p-valeurs pour le test des proportions de Hardy-Weinberg montre que les 
#données sont approximativement en accord avec l'hypothèse de Hardy-Weinberg. Cependant, il y a quelques
#points qui s'écartent de la ligne rouge, ce qui pourrait indiquer que l'hypothèse n'est pas vraie pour tous les 




#QUESTION 4

# Appliquer le thinning LD

data_thinned <- LD.thin(selected_data, 0.1, which.snps = is.autosome(selected_data@snps$chr))
data_thinned
#A bed.matrix with 7268 individuals and 155 markers.
#snps stats are set
#ped stats are set

# Calculer la matrice d'apparentement génomique (GRM)
grm<- GRM(data_thinned)
ei_grm <- eigen(grm)

plot(ei_grm$vectors, xlab = 'PC1', ylab = 'PC2')

# Comparer les valeurs obtenues pour les individus communs aux deux jeux de données
# Supposons que vous avez déjà une matrice d'apparentement K calculée à partir d'une autre méthode et stockée dans un objet appelé "k_matrix"

?intersect
# Sélectionner les individus communs aux deux jeux de données
commun <- intersect(rownames(K), rownames(grm))

# Extraction des coefficients de la matrice K pour les individus communs
commun_K <- K[commun, commun]

# Extraction des coefficients de la matrice GRM pour les individus communs
commun_grm <- grm[commun, commun]

# Extraire les coefficients d'apparentement pour les individus communs dans la matrice K
coefficients_K <- diag(commun_K)
coefficients_K <- as.vector(coefficients_K)

# Extraire les coefficients d'apparentement pour les individus communs dans la matrice GRM
coefficients_grm <- diag(commun_grm)
coefficients_grm <- as.vector(coefficients_grm)
correlation <- cor(coefficients_K, coefficients_grm)

# Calculer la corrélation entre les coefficients d'apparentement des deux matrices
data_cor <- data.frame(coefficients_grm,coefficients_k)

 


#Question 5

beastx <- read.table("BEASTX.tsv", header = TRUE, sep = "\t")

beastx <- beastx[, c(1, 2, 3, 4, 12, 13)]

beastx$LambAdult<- as.factor(beastx$LambAdult)


hist(beastx$IgAmp[beastx$LambAdult == "Lambs"],
     main = "Histogramme des valeurs d'IgA chez les agneaux",
     xlab = "IgA")
hist(beastx$IgAmp[(beastx$LambAdult == "Adults")],
     main = "Histogramme des valeurs d'IgA chez les adultes", 
     xlab = "IgA")

beastIGA <- beastx$IgAmp[!is.na(beastx$IgAmp)]
length(beastIGA) # nous avons 6483 individus avec les IGAmp disponible

match_vector<- match(beastx$MumID, pedigree$MOTHER, nomatch = 0)
#match_vector<- match(pedigree$MOTHER, beastx$MumID nomatch = 0)
match_valeurs <- which(match_vector!=0)
valeur_dispo <- match_vector[match_valeurs]
length(valeur_dispo) ##on a 6496 valeurs concordantes 
#################################################
all(beastx$ID %in% data_bed@ped$id)## tous les individus génotypés n'ont pas phenotype 

###################################################
beast_sans_na <- subset(beastx, !is.na(beastx$IgAmp))
Bed_disponible <- select.inds(data_bed, data_bed@ped$id %in% beast_sans_na$ID)

##Question 6

moigneaux <- subset(beast_sans_na, beast_sans_na$LambAdult=="Lambs")
mouton <- subset(beast_sans_na, beast_sans_na$LambAdult=="Adults")

id_moigneaux <- moigneaux$ID ## liste des identifiant des agneaux sans Na
id_moigneaux<-as.character(id_moigneaux) #transformation en charactere
colonne_grm <- colnames(grm) #noms des colonnes de la matrice grm
length(id_moigneaux)
valeurs_grm  <- id_moigneaux[id_moigneaux %in% colonne_grm] # noms des colonne de la grm retrouvé dans les identifiant des agneaux 
#valeurs_grm <- intersect(colonne_grm, id_moigneaux)
length(valeurs_grm)
grm_moigneaux<- grm[valeurs_grm, valeurs_grm]

dim(grm_moigneaux)

valeurs_IgAmp <- subset(moigneaux, moigneaux$ID %in% valeurs_grm)$IgAmp

heritabilite_moigneaux<- lmm.aireml(valeurs_IgAmp, 
           X = matrix(1, nrow = length(valeurs_IgAmp)),
           K= grm_moigneaux)

##Une fois qu’on a obtenu des estimations de τ et σ on peut estimer h par τ/τ+σ2
sigma2 <- heritabilite_moigneaux$sigma2
tau <- heritabilite_moigneaux$tau

# Calculer l'héritabilité
heritabilite <- tau**2 / (tau**2 + sigma2)
heritabilite # l'heritabilité des agneaux est de 0.0007540357

#################################################################

mouton <- subset(beast_sans_na, beast_sans_na$LambAdult=="Adults")

id_mouton <- mouton$ID ## liste des identifiant des agneaux sans Na
id_mouton <-as.character(id_mouton) #transformation en charactere
colonne_grm <- colnames(grm) #noms des colonnes de la matrice grm

valeurs_retrouvees <- id_mouton[id_mouton %in% colonne_grm] ##identité des mouton qu'on retrouve dans les colonne de la matrice
length(valeurs_retrouvees)
##valeurs_grm_mouton <- intersect(colonne_grm, id_mouton)

grm_mouton<- grm[valeurs_retrouvees, valeurs_retrouvees]

dim(grm_mouton)

valeurs_IgAmp_mouton <- subset(mouton, mouton$ID %in% valeurs_retrouvees)$IgAmp

heritabilite_mouton<- lmm.aireml(valeurs_IgAmp_mouton, 
                                    X = matrix(1, nrow = length(valeurs_IgAmp_mouton)),
                                    K= grm_mouton)

sigma2_m <- heritabilite_mouton$sigma2
tau_m <- heritabilite_mouton$tau

heritabilite_mouton <- tau_m**2 / (tau_m**2 + sigma2_m)
heritabilite_mouton # l'heritabilité des mouton est de 0.1285155



##QUESTION 7
#beast_sans_na <- subset(beastx, !is.na(beastx$IgAmp))
#moigneaux <- subset(beast_sans_na, beast_sans_na$LambAdult=="Lambs")
#mouton <- subset(beast_sans_na, beast_sans_na$LambAdult=="Adults")

K2 <- 2 * K
nom_colonne_K <- colnames(K2)

id_moigneaux <- moigneaux$ID ## liste des identifiant des agneaux sans Na
id_moigneaux<-as.character(id_moigneaux) #transformation en charactere
nom_colonne_K <- colnames(K2)#noms des colonnes de la matrice K2

#valeurs_retrouvees_k1 <- id_moigneaux[id_moigneaux %in% nom_colonne_K] # noms des colonne de la grm retrouvé dans les identifiant des agneaux 
#length(valeurs_retrouvees_k1) autre methode

#valeurs_k <- match(nom_colonne_K, id_moigneaux, nomatch = 0)
#match_val <- which(valeurs_k!=0)
#valeurs_k <- valeurs_k[match_val]
#valeurs_k<- as.character(valeurs_k)


valeurs_K <- intersect(nom_colonne_K, id_moigneaux)
length(valeurs_K)

K_moigneaux<- K2[valeurs_K, valeurs_K]
dim(K_moigneaux)

valeurs_IgAmp_K <- subset(moigneaux, moigneaux$ID %in% valeurs_K)$IgAmp
length(valeurs_IgAmp_K)

heritabilite_moigneaux_K<- lmm.aireml(valeurs_IgAmp_K, 
                                    X = matrix(1, nrow = length(valeurs_IgAmp_K)),
                                    K= K_moigneaux)

sigma2_moi_K <- heritabilite_moigneaux_K$sigma2
tau_mmoi_K <- heritabilite_moigneaux_K$tau

heritabilite_agneaux_K <- tau_mmoi_K**2 / (tau_mmoi_K**2 + sigma2_moi_K)
heritabilite_agneaux_K
print(paste("Héritabilité chez les agneaux est :",heritabilite_agneaux_K) )
############################################################################

K2 <- 2 * K
nom_colonne_K <- colnames(K2)

id_mouton<- mouton$ID ## liste des identifiant des agneaux sans Na
id_mouton<-as.character(id_mouton) #transformation en charactere
nom_colonne_K <- colnames(K2)#noms des colonnes de la matrice K2

valeurs_mouton_K <- id_mouton[id_mouton %in% nom_colonne_K] # noms des colonne de la K2 retrouvé dans les identifiant des agneaux 
#length(valeurs_retrouvees_k1) autre methode
#valeurs_mouton_K <- intersect(nom_colonne_K, id_mouton) #ceci ne marche pas 
length(valeurs_mouton_K) #valeur 3935

K_mouton<- K2[valeurs_mouton_K, valeurs_mouton_K]

dim(K_mouton)

valeurs_IgAmp_mouton_K <- subset(mouton, mouton$ID %in% valeurs_mouton_K)$IgAmp
length(valeurs_IgAmp_mouton_K)

heritabilite_mouton_K<- lmm.aireml(valeurs_IgAmp_mouton_K, 
                                      X = matrix(1, nrow = length(valeurs_IgAmp_mouton_K)),
                                      K= K_mouton)

sigma2_mouton_K <- heritabilite_mouton_K$sigma2
tau_mou_K <- heritabilite_mouton_K$tau

heritabilite_mouton_K <- tau_mou_K**2 / (tau_mou_K**2 + sigma2_mouton_K)
heritabilite_mouton_K
print(paste("Héritabilité chez les adultes est :",heritabilite_mouton_K) )



###QUESTION 8
# Accédez aux identifiants des individus à partir de la bed matrix
individus_bed_matrix <- data_bed@ped$id

# Sélectionnez les individus dont les identifiants sont dans id_moigneaux
individus_selectionnes <- individus_bed_matrix[individus_bed_matrix %in% id_moigneaux]

# Utilisez la fonction select.inds pour extraire la sous-matrice correspondant aux individus sélectionnés
bed_matrix_moigneaux <- select.inds(data_bed, id %in% individus_selectionnes)

#obtenir le jeux de données donc les ID sont disponible dans la GRM 
moigneaux_disponible <- subset(moigneaux, id_moigneaux %in% colonne_grm)
##########################################"
#MODELE MIXTE
aa <- association.test(bed_matrix_moigneaux,
                       Y=moigneaux_disponible$IgAmp, 
                       K = grm_moigneaux, method='lmm')
qqplot.pvalues(aa, col.abline = "green")
manhattan(aa)
abline(-log10(5e-8), 0, col='red')

#MODELE LINEAIRE

aa_linear <- association.test(bed_matrix_moigneaux,
                       Y=moigneaux_disponible$IgAmp, 
                       K = grm_moigneaux, method='lm', test = "wald")
qqplot.pvalues(aa_linear, col.abline = "green")

manhattan(aa_linear)
abline(-log10(5e-8), 0, col='red')


### QUESTION 9

individus_bed_matrix <- data_bed@ped$id

# Sélectionnez les individus de la bed matrix dont les identifiants sont dans id_mouton
individus_selectionnes_mouton <- individus_bed_matrix[individus_bed_matrix %in% id_mouton]
length(individus_selectionnes_mouton)
# Utilisez la fonction select.inds pour extraire la sous-matrice correspondant aux individus sélectionnés
bed_matrix_mouton <- select.inds(data_bed, id %in% individus_selectionnes_mouton)

## liste des ID dont les valeurs ne sont pas disponible dans la GRM
valeurs_absentes <- id_moigneaux[!id_mouton %in% colonne_grm]
#obtenir le jeux de données donc les ID sont disponible dans la GRM 
mouton_disponible <- subset(mouton, mouton$ID%in% colonne_grm)
mouton_disponible <- subset(mouton_disponible, individus_selectionnes_mouton%in% colonne_grm)
mouton_aggregated <- aggregate(IgAmp ~ ID, data = mouton_disponible, FUN = mean)



#######################################################""


# Accédez aux identifiants des individus à partir de la bed matrix
individus_bed_matrix <- data_bed@ped$id
length(individus_bed_matrix)

# Sélectionnez les individus dont les identifiants sont dans id_mouton
individus_selectionnes_mouton <- individus_bed_matrix[individus_bed_matrix %in% id_mouton]
length(individus_selectionnes_mouton)
individus_selectionnes_mouton <- as.character(individus_selectionnes_mouton)

#nouvelle matrice grm_mouton avec les individus dont les id sont disponible dans la bed matrix
grm_mouton_lm <- grm_mouton[individus_selectionnes_mouton, individus_selectionnes_mouton]
dim(grm_mouton_lm)

aa_linear_mouton <- association.test(bed_matrix_mouton,
                              Y=mouton_aggregated$IgAmp, 
                              K = grm_mouton_lm, method='lm', test = "wald")


qqplot.pvalues(aa_linear_mouton, col.abline = "green")

manhattan(aa_linear_mouton)
abline(-log10(5e-8), 0, col='red')




aa_mixte_mouton <- association.test(bed_matrix_mouton,
                       Y=mouton_aggregated$IgAmp, 
                       K = grm_mouton_lm, method='lmm')

qqplot.pvalues(aa_mixte_mouton, col.abline = "green")
manhattan(aa)
abline(-log10(5e-8), 0, col='red')



# 1. Sélectionner les SNPs du chromosome 24
snps_chr24 <- select.snps(data_bed, data_bed@snps$chr == 24)

individuschr24_bed_matrix <- snps_chr24@ped$id

# Sélectionnez les individus dont les identifiants sont dans id_moigneaux
individuschr24_selectionnes <- individuschr24_bed_matrix[individuschr24_bed_matrix %in% id_moigneaux]

# Utilisez la fonction select.inds pour extraire la sous-matrice correspondant aux individus sélectionnés
bed_matrix_moigneauxchr24 <- select.inds(snps_chr24, id %in% individuschr24_selectionnes)

#obtenir le jeux de données donc les ID sont disponible dans la GRM 
subset_chr24 <- subset(moigneaux, moigneaux$ID %in% individuschr24_selectionnes )
individuschr24_selectionnes<- as.character(individuschr24_selectionnes)

grmchr24 <- grm_moigneaux[individuschr24_selectionnes, individuschr24_selectionnes]

dim(grmchr24)



aa_linear_chr24 <- association.test(bed_matrix_moigneauxchr24,
                                     Y=subset_chr24$ID, 
                                     K = grmchr24, method='lm')


qqplot.pvalues(aa_linear_chr24, col.abline = "green")
manhattan(aa_linear_chr24)
abline(-log10(5e-8), 0, col='red')





aa_mixte_chr24 <- association.test(bed_matrix_moigneauxchr24,
                                    Y=subset_chr24$ID, 
                                    K = grmchr24, method='lmm', test = "wald")

qqplot.pvalues(aa_mixte_chr24, col.abline = "green")
manhattan(aa)
abline(-log10(5e-8), 0, col='red')




















