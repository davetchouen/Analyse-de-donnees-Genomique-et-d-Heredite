# Analyse de données Génomique et d'Hérédité : Étude des Valeurs d'IgA et Héritabilité chez les Ovins 🐑
Ce projet a été réalisé dans le cadre de ma dernière année de master en sciences des données. Il vise à analyser des données génomiques et phénotypiques d'ovins pour étudier l'héritabilité des valeurs d'IgA, un biomarqueur d'immunité, et détecter des associations génétiques significatives. Les données utilisées incluent des informations généalogiques (pedigree) et des données génomiques (SNPs), permettant d'explorer la structure génétique de la population et d'identifier les facteurs génétiques influençant les niveaux d'IgA.
## Objectifs du projet
1. Étudier les coefficients de consanguinité et d'apparentement pour une population ovine.
2. Analyser les SNPs (polymorphismes nucléotidiques simples) pour identifier les régions génomiques d'intérêt.
3. Estimer l'héritabilité des valeurs d'IgA chez les agneaux et les moutons adultes.
4. Effectuer des tests d'association génétique en utilisant des modèles mixtes et linéaires.

## Méthodologie
### 1. Calcul des Coefficients de Parenté et Consanguinité

Données : Fichier pedigree.tsv contenant les informations généalogiques (ID, père, mère).

Méthode : Utilisation du package kinship2 pour calculer la matrice de parenté (K) et les coefficients de consanguinité (f).

Résultats : Identification de 28 individus avec un coefficient de consanguinité supérieur ou égal à 0.25.

### 2. Analyse des Données Génétiques
Données : Fichiers .bed, .bim, et .fam contenant les informations génomiques.

Méthode : Utilisation du package gaston pour lire et manipuler les données génétiques. Filtrage des SNPs basé sur le taux d'appel (callrate) et la fréquence allélique mineure (MAF).

Résultats : Sélection de 155 marqueurs après filtrage et calcul de la matrice d'apparentement génomique (GRM).

### 3. Test de Hardy-Weinberg
Méthode : Vérification de l'équilibre de Hardy-Weinberg pour les SNPs sélectionnés.

Résultats : Les données sont globalement en accord avec l'hypothèse de Hardy-Weinberg, bien que quelques écarts soient observés.


### 4. Estimation de l'Héritabilité
Méthode : Utilisation de modèles mixtes linéaires pour estimer l'héritabilité du trait IgA chez les agneaux et les adultes.

Résultats :

Héritabilité chez les agneaux : 0.0007540357

Héritabilité chez les adultes : 0.1285155

### 5. Analyse d'Association
Méthode : Tests d'association entre les SNPs et les niveaux d'IgA en utilisant des modèles linéaires et mixtes.

Résultats : Visualisation des résultats via des QQ-plots et des manhattan plots pour identifier les SNPs significatifs.

## Conclusion
Ce projet a permis de mieux comprendre la structure génétique de la population étudiée et d'identifier les facteurs génétiques influençant les niveaux d'IgA. Les résultats montrent une faible héritabilité chez les agneaux, mais une héritabilité modérée chez les adultes, suggérant que les facteurs génétiques jouent un rôle plus important dans la régulation des niveaux d'IgA chez les individus adultes.

## Utilisation du Code
Pour reproduire les analyses,utilisez le code_brut.R fourni
## Fichiers Inclus
pedigree.tsv : Fichier contenant les informations généalogiques.

20180209_SoayPlates1-83.bed, 20180209_SoayPlates1-83.bim, 20180209_SoayPlates1-83.fam : Fichiers contenant les données génomiques.

Le fichier Genetique1.RDataTmp etant tres volumineux n'a pas été fourni vous pouvez me contacter pour l'avoir (davetchouenkou9@gmail.com)
