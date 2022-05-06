setwd("/Users/MAEL/Documents/M1_BI/Stage/research/new_models")

data_cys_restreint = read.table("bdd_models_byCys.txt_restreint", header = T)

which(is.na(data_cys_restreint), arr.ind = TRUE)
#Toutes les na sont sur la colonne 6 : qmean

#Matrice dans laquelle on indique le nombre de cystéines modifiées par protéine
df_prot_nb_cys = matrix(ncol = 2, nrow = 0)

for (num_access in data_cys_restreint$access) {
  if( length(which(df_prot_nb_cys[,1] == num_access)) == 0 ){
    df_prot_nb_cys = rbind(df_prot_nb_cys, c(num_access, 1))
  }else{
    df_prot_nb_cys[which(df_prot_nb_cys == num_access, arr.ind = TRUE)[1],2] = 
      as.integer(df_prot_nb_cys[which(df_prot_nb_cys == num_access, arr.ind = TRUE)[1],2]) + 1
  }
}

#Conversion en data frame de la matrice et des scores en int
df_prot_nb_cys = as.data.frame(df_prot_nb_cys)
colnames(df_prot_nb_cys) = c("access", "nb_cys")
df_prot_nb_cys$nb_cys <- as.numeric(df_prot_nb_cys$nb_cys)

#Summary sur les datas concernant les plis de Rossman
summary(data_cys_restreint$Rossman.like)
summary(as.factor(data_cys_restreint$Rossman.like))
summary(as.factor(data_cys_restreint$modif))

summary(data_cys_restreint$modif[which(data_cys_restreint$Rossman.like == 0)])
summary(data_cys_restreint$modif[which(data_cys_restreint$Rossman.like == 1)])

#Test du chi2 pour savoir si les modif (avec plus de 5 itérations)
# sont influencées par le pli de Rossman
M = as.table(rbind(c(15, 35, 59, 210), c(13, 10, 31, 128)))
dimnames(M) = list(ross = c("0", "1"), modif = c("mini", "NRg", "Nrg", "nRg"))
Xsq = chisq.test(M)
Xsq
Xsq$observed
Xsq$expected
#test non significatif

#Vérifications sur la position des cystéines et le taux de conservation
summary(data_cys_restreint$cys_pos)
summary(data_cys_restreint$taux_conserv[which(data_cys_restreint$taux_conserv != -1000)])
boxplot(data_cys_restreint$cys_pos)
boxplot(data_cys_restreint$taux_conserv[which(data_cys_restreint$taux_conserv != -1000)])

#Numéro d'accession des cys dans un pli de rossman
#Sélection des protéines avec une cystéine modifiée et un pli de rossman
list1 = data_cys_restreint$access[which(data_cys_restreint$Rossman.like == 1)]
list2 = df_prot_nb_cys$access[which(df_prot_nb_cys$nb_cys < 2)]
list3 = data_cys_restreint$access[which(data_cys_restreint$taux_conserv > 0.5)]
cys_ross_0.5 = intersect(intersect(list1, list2), list3)

#Choix des 4 protéines à observer sur PyMol
data_cys_restreint[which(data_cys_restreint$access == "Q39580"),c(1:22)]
data_cys_restreint[which(data_cys_restreint$access == "A8JH72"),c(1:22)]
data_cys_restreint[which(data_cys_restreint$access == "Q4U0V8"),c(1:22)]
data_cys_restreint[which(data_cys_restreint$access == "P09205"),c(1:22)]

#Petite vérif
summary(data_cys_restreint$modif)
vec_pos = c()
i = 1
for (name in cys_ross_0.5) {
  vec_pos[i] = which(data_cys_restreint$access == name)
  i = i+1
}
summary(data_cys_restreint$modif[vec_pos])


#random forest

library(rpart)
library(rpart.plot)

#Visualisation via cmdscale

#Visualisation réelle

vec_col = ifelse(data_cys_restreint$modif == "NRG", 0, 
       ifelse(data_cys_restreint$modif == "NRg", 1,
              ifelse(data_cys_restreint$modif == "NrG", 2,
                     ifelse(data_cys_restreint$modif == "Nrg", 3,
                            ifelse(data_cys_restreint$modif == "nRG", 4,
                                   ifelse(data_cys_restreint$modif == "nRg", 5,6))))))


#En utilisant seulement les descripteurs quantitatifs

dim(data_cys_restreint)
desc_quant = data_cys_restreint[,c(2,12:14,18,19,23:202)]

dist_matrix = dist(desc_quant)
Xp = cmdscale(dist_matrix, k = 2)
plot(Xp)

cl = kmeans(desc_quant, centers = 5)
cl$cluster
Xp = cbind(Xp, cl$cluster)
plot(Xp[,1:2], col = Xp[,3]+1)

desc_quant = data.frame(desc_quant, cluster = as.factor(cl$cluster))


#Choix du nombre de cluster optimal
vec_expl = c()
for (i in c(1:14)) {
  cl = kmeans(Xp[,1:2], centers = i)
  vec_expl[i] = cl$tot.withinss
}
plot(c(1:length(vec_expl)),vec_expl, "l", col = "Blue")

#Anova sur tous les descripteurs quantitatifs
summary(desc_quant$cluster)
vec_p_val = c()

for (i in c(1:(ncol(desc_quant)-1))) {
  vec_p_val[i] = summary(aov(desc_quant[,i]~cluster, data = desc_quant))[[1]][5][1,1]
}

signif_p_val = which(vec_p_val<0.05)
colnames(desc_quant)[which(vec_p_val<0.05)]
new_desc_quant = desc_quant[,c(signif_p_val, ncol(desc_quant))]

#par(mfrow = c(3,3), mar = c(1,1,1,1))

total_test_post_hoc = c()

for (i in c(1:(ncol(new_desc_quant)-1))) {
  total_test_post_hoc[i] = TukeyHSD(aov(new_desc_quant[,i]~cluster, data = new_desc_quant))
}

df_pval = data.frame(rep(0,10))
for (list_pval in total_test_post_hoc) {
  df_pval = cbind(df_pval, list_pval[,4])
}
df_pval = df_pval[,-1]

ref_mat = which(df_pval <= 0.05, arr.ind = TRUE)

mat_moy_desc = matrix(nrow = ncol(new_desc_quant)-1, ncol = 5)
mat_moy_desc = as.data.frame(mat_moy_desc)

for (i in c(1:(ncol(new_desc_quant)-1) ) ) {
  for (j in c(1:5)) {
    mat_moy_desc[i,j] = mean(new_desc_quant[,i][which(new_desc_quant$cluster == j)])
    rownames(mat_moy_desc)[i] = colnames(new_desc_quant)[i]
  }
}

colnames(mat_moy_desc) = c(1,2,3,4,5)
mat_moy_desc_clone = mat_moy_desc

for (i in c(1:nrow(mat_moy_desc)) ) {
  print(rownames(mat_moy_desc[i,]))
  print(colnames(sort(mat_moy_desc[i,])))
  mat_moy_desc_clone[i,] = colnames(sort(mat_moy_desc[i,]))
}

setwd("/Users/MAEL/Documents/M1_BI/Stage/research/recherche_perso_cys/descripteurs_cys_mod")
write.table(mat_moy_desc_clone, "moyenne_triee.txt")
write.table(ref_mat, "result_post_hoc.txt")
write.table(new_desc_quant, "data_utilisees.txt")
setwd("/Users/MAEL/Documents/M1_BI/Stage/research/new_models")

total_test_post_hoc[[4]][,4]

#plot(TukeyHSD(aov(cys_pos~cluster, data = desc_quant)))
#pairwise.t.test(desc_quant$cys_pos, desc_quant$cluster, p.adjust.method = "bonferroni")

desc_quant_norm = desc_quant[,-187]
for (i in c(1:186)) {
  desc_quant_norm[,i] = scale(desc_quant_norm[,i])
}

dist_matrix = dist(desc_quant_norm)
Xp = cmdscale(dist_matrix, k = 2)
plot(Xp)

cl = kmeans(desc_quant_norm, centers = 3)
cl$cluster
Xp = cbind(Xp, cl$cluster)
plot(Xp[,1:2], col = Xp[,3]+1)

#a = aov(desc_quant[,7]~cluster, data = desc_quant) ANOVA
#b = summary(a) SUMMARY
#b[[1]][5][1,1] P-VALUE

#Suite de l'analyse sur les 5 groupes de cystéines
gr1 = which(cl$cluster == 1)
gr2 = which(cl$cluster == 2)
gr3 = which(cl$cluster == 3)
gr4 = which(cl$cluster == 4)
gr5 = which(cl$cluster == 5)
num_quant = c(2,12:14,18,19,23:202)

summary(as.factor(data_cys_restreint$Rossman.like[gr1]))
summary(as.factor(data_cys_restreint$Rossman.like[gr2]))
summary(as.factor(data_cys_restreint$Rossman.like[gr3]))
summary(as.factor(data_cys_restreint$Rossman.like[gr4]))
summary(as.factor(data_cys_restreint$Rossman.like[gr5]))

# test chi2

table_grp_ross = as.table(rbind(summary(as.factor(data_cys_restreint$Rossman.like[gr1])),
                                summary(as.factor(data_cys_restreint$Rossman.like[gr2])),
                                summary(as.factor(data_cys_restreint$Rossman.like[gr3])),
                                summary(as.factor(data_cys_restreint$Rossman.like[gr4])),
                                summary(as.factor(data_cys_restreint$Rossman.like[gr5]))))

Xsq = chisq.test(table_grp_ross)
Xsq
Xsq$observed
Xsq$expected

summary(as.factor(data_cys_restreint$modif[gr1]))
summary(as.factor(data_cys_restreint$modif[gr2]))
summary(as.factor(data_cys_restreint$modif[gr3]))
summary(as.factor(data_cys_restreint$modif[gr4]))
summary(as.factor(data_cys_restreint$modif[gr5]))

table_grp_modif = as.table(rbind(summary(as.factor(data_cys_restreint$modif[gr1])),
                                 summary(as.factor(data_cys_restreint$modif[gr2])),
                                 summary(as.factor(data_cys_restreint$modif[gr3])),
                                 summary(as.factor(data_cys_restreint$modif[gr4])),
                                 summary(as.factor(data_cys_restreint$modif[gr5]))))

table_grp_modif

Xsq = chisq.test(table_grp_ross)
Xsq
Xsq$observed
Xsq$expected

#En utilisant seulement les descripteurs acides aminés

desc_quant = data_cys_restreint[,c(23:202)]

dist_matrix = dist(desc_quant)
Xp = cmdscale(dist_matrix, k = 2)
plot(Xp)

cl = kmeans(desc_quant, centers = 3)
cl$cluster
Xp = cbind(Xp, cl$cluster)
plot(Xp[,1:2], col = Xp[,3]+1)

vec_expl = c()
for (i in c(1:14)) {
  cl = kmeans(Xp[,1:2], centers = i)
  vec_expl[i] = cl$tot.withinss
}
plot(c(1:length(vec_expl)),vec_expl, "l", col = "Blue")

#Hclust
#En utilisant seulement les descripteurs quantitatifs
desc_quant = data_cys_restreint[,c(2,12:14,18,19,23:202)]
dist_matrix = dist(desc_quant)

hc=hclust(dist_matrix,method="ward.D2")
plot(hc, main="gene")
cl2 = cutree(hc, k = 7)

Xp = cmdscale(dist_matrix, k = 2)
Xp_hclust = cbind(Xp, cl2)
plot(Xp_hclust[,1:2], col = Xp_hclust[,3]+1)

#En utilisant seulement les descripteurs acides aminés
desc_quant = data_cys_restreint[,c(23:202)]
dist_matrix = dist(desc_quant)
hc=hclust(dist_matrix,method="ward.D2")
plot(hc, main="gene")
cl2 = cutree(hc, k = 7)

Xp = cmdscale(dist_matrix, k = 2)
Xp_hclust = cbind(Xp, cl2)
plot(Xp_hclust[,1:2], col = Xp_hclust[,3]+1)



