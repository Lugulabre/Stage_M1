setwd("/Users/MAEL/Documents/M1_BI/Stage/research/new_models")
library(stringr)
library(corrplot)
library(nparcomp)

data_cys_restreint = read.table("bdd_models_byCys.txt_restreint", header = T)
data_cys_AA = data_cys_restreint[,23:202]

list_AA = c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V")
colnames(data_cys_AA)[1:20] = list_AA

data_AA = matrix(data = 0, nrow = 501, ncol = 20)
data_AA = as.data.frame(data_AA)
colnames(data_AA) = list_AA

for (i in c(1:20)) {
  col_aa = which(substr(colnames(data_cys_AA), 1,1) == list_AA[i])
  data_AA[,i] = apply(data_cys_AA[,col_aa], 1, sum)
}

for (i in c(1:nrow(data_AA))) {
  data_AA[i,] = data_AA[i,]/sum(data_AA[i,])
}

summary(apply(data_AA, 1, sum))

#aromatic residues
listAR = c('F', 'Y', 'H', 'W')
freqAR = apply(data_AA[,which(is.element(colnames(data_AA), listAR))], 1, sum)
data_AA$freqAR = freqAR

#polar residues
listPR = c('C','D','E','H','K','N','Q','R','S','T','W','Y')
freqPR = apply(data_AA[,which(is.element(colnames(data_AA), listPR))], 1, sum)
data_AA$freqPR = freqPR

#aliphatic residues
listAlR = c('I', 'L', 'V')
freqAlR = apply(data_AA[,which(is.element(colnames(data_AA), listAlR))], 1, sum)
data_AA$freqAlR = freqAlR

#charged residues
listCR = c('D', 'E', 'R', 'K', 'H')
freqCR = apply(data_AA[,which(is.element(colnames(data_AA), listCR))], 1, sum)
data_AA$freqCR = freqCR

#negative residues
listNR = c('D', 'E')
freqNR = apply(data_AA[,which(is.element(colnames(data_AA), listNR))], 1, sum)
data_AA$freqNR = freqNR

#positive residues
listPosR = c('H', 'K', 'R')
freqPosR = apply(data_AA[,which(is.element(colnames(data_AA), listPosR))], 1, sum)
data_AA$freqPosR = freqPosR

#hydrophobic residues
listHR = c('A','V','L','I','M','F','W', 'Y')
freqHR = apply(data_AA[,which(is.element(colnames(data_AA), listHR))], 1, sum)
data_AA$freqHR = freqHR

#small residues
listSR = c('C','V','T','G','A','S','D','N','P')
freqSR = apply(data_AA[,which(is.element(colnames(data_AA), listSR))], 1, sum)
data_AA$freqSR = freqSR

#tiny residues
listTR = c('A', 'C', 'G', 'S')
freqTR = apply(data_AA[,which(is.element(colnames(data_AA), listTR))], 1, sum)
data_AA$freqTR = freqTR

#matrice corrélation
MatCor = cor(data_AA)
#corrplot(MatCor)
a = which(MatCor > 0.8, arr.ind = TRUE)
a = a[-which(a[,1] == a[,2]),]
#0.81 entre freqNR et D -> conservation

#Kmeans
dist_matrix = dist(data_AA)
Xp = cmdscale(dist_matrix, k = 2)
#plot(Xp)

cl = kmeans(data_AA, centers = 4)
#Pour réutiliser les clusters du rapport
clust_util = read.table("../recherche_perso_cys/analyse_AA_env/bis/cluster_utilises.txt", header = TRUE)
cl$cluster = as.vector(clust_util[,1])
Xp = cbind(Xp, cl$cluster)
plot(Xp[,1:2], col = Xp[,3]+1, main = "cluster selon environnement AA",
     xlab = "x", ylab = "y")
legend("bottomright", c(as.character(c(1,2,3,4))),
       col = c(2,3,4,5),
       lty = 1, cex = 0.6, ncol = 2)

data_AA = cbind(data_AA, as.factor(cl$cluster))
colnames(data_AA)[30] = "cluster"

summary(as.factor(data_AA$cluster))
#représentation en boxplot des descripteurs
boxplot(data_AA[,-30])

vec_expl = c()
for (i in c(1:14)) {
  cl = kmeans(Xp[,1:2], centers = i)
  vec_expl[i] = cl$tot.withinss
}
plot(c(1:length(vec_expl)),vec_expl, "l", col = "Blue", xlab = 'nb k-means', ylab = "tot. withinss")

#conditions d'application
vec_bartl = c()
#vec_norm = c()
for (i in c(1:(ncol(data_AA)-1))) {
  vec_bartl[i] = bartlett.test(data_AA[,i]~cluster, data = data_AA)[[3]]
  #vec_norm[i] = by(data_AA[,i], data_AA$cluster, shapiro.test)[[2]]
}
colnames(data_AA)[which(vec_bartl > 0.05)]
#vec_norm < 0.05
#colnames(data_AA)

#ANOVA entre les 4 groupes
vec_p_val = c()
vec_p_val_kw = c()

for (i in c(1:(ncol(data_AA)-1))) {
  if(vec_bartl[i]>0.05){
    vec_p_val[i] = summary(aov(data_AA[,i]~cluster, data = data_AA))[[1]][5][1,1]
  }else{
    vec_p_val_kw[i] = kruskal.test(data_AA[,i]~cluster, data = data_AA)[[3]]
  }
}

#signif anova
signif_p_val = which(vec_p_val<0.05)
colnames(data_AA)[which(vec_p_val<0.05)]
vec_names_AOV = colnames(data_AA)[which(vec_p_val<0.05)]
#signif kruskal wallis
signif_p_val_kw = which(vec_p_val_kw<0.05)
colnames(data_AA)[which(vec_p_val_kw<0.05)]
vec_names_KW = colnames(data_AA)[which(vec_p_val_kw<0.05)]

new_desc_quant = data_AA[,c(sort(c(signif_p_val, signif_p_val_kw)), ncol(data_AA))]

new_desc_quant$cluster = as.factor(new_desc_quant$cluster)

#TEST POST HOC
total_test_post_hoc = c()
for (i in c(1:(ncol(data_AA)-1))) {
  total_test_post_hoc[i] = TukeyHSD(aov(data_AA[,i]~cluster, data = new_desc_quant))
}
total_test_post_hoc = total_test_post_hoc[signif_p_val]
#dim(total_test_post_hoc)

#TEST POST HOC NON PARAM
total_post_hoc_np = matrix(nrow = (ncol(data_AA)-1), ncol = 6)
for (i in c(1:(ncol(data_AA)-1))) {
  temp = nparcomp::nparcomp(data = data_AA, data_AA[,i] ~ cluster)
  temp = as.vector(temp[[3]][1][which(temp[[3]][6] < 0.05),])
  if(length(temp) != 0){
    total_post_hoc_np[i,c(1:length(temp))] = temp

  }
#  a = kw_comp[[3]][1][which(kw_comp[[3]][6] < 0.05),]
#  a = as.vector(a)
}

total_post_hoc_np = as.data.frame(total_post_hoc_np)
rownames(total_post_hoc_np) = colnames(data_AA[,-30])
total_post_hoc_np = total_post_hoc_np[signif_p_val_kw,]

len_df = length(total_test_post_hoc[[1]][,4])
df_pval = data.frame(rep(0,len_df))
for (list_pval in total_test_post_hoc) {
  df_pval = cbind(df_pval, list_pval[,4])
}
df_pval = df_pval[,-1]

ref_mat = which(df_pval <= 0.05, arr.ind = TRUE)

mat_moy_desc = matrix(nrow = ncol(new_desc_quant)-1, ncol = length(summary(new_desc_quant$cluster)))
mat_moy_desc = as.data.frame(mat_moy_desc)

for (i in c(1:(ncol(new_desc_quant)-1) ) ) {
  for (j in c(1:ncol(mat_moy_desc))) {
    mat_moy_desc[i,j] = mean(new_desc_quant[,i][which(new_desc_quant$cluster == j)])
    rownames(mat_moy_desc)[i] = colnames(new_desc_quant)[i]
  }
}

colnames(mat_moy_desc) = c(1,2,3,4)
mat_moy_desc_clone = mat_moy_desc

for (i in c(1:nrow(mat_moy_desc)) ) {
  print(rownames(mat_moy_desc[i,]))
  print(colnames(sort(mat_moy_desc[i,])))
  mat_moy_desc_clone[i,] = colnames(sort(mat_moy_desc[i,]))
}

namegrp1 = data.frame(name = data_cys_restreint$access[which(data_AA$cluster == 1)],
                      cys = data_cys_restreint$cys_pos[which(data_AA$cluster == 1)])
namegrp2 = data.frame(name = data_cys_restreint$access[which(data_AA$cluster == 2)],
                      cys = data_cys_restreint$cys_pos[which(data_AA$cluster == 2)])
namegrp3 = data.frame(name = data_cys_restreint$access[which(data_AA$cluster == 3)],
                      cys = data_cys_restreint$cys_pos[which(data_AA$cluster == 3)])
namegrp4 = data.frame(name = data_cys_restreint$access[which(data_AA$cluster == 4)],
                      cys = data_cys_restreint$cys_pos[which(data_AA$cluster == 4)])

setwd("/Users/MAEL/Documents/M1_BI/Stage/research/recherche_perso_cys/analyse_AA_env/bis")
write.table(mat_moy_desc_clone, "moyenne_triee.txt")
write.table(total_post_hoc_np, "result_post_hoc_NP.txt")
#write.table(total_post_hoc, "result_post_hoc_fin.txt")
write.table(ref_mat, "result_post_hoc.txt")
write.table(new_desc_quant, "data_utilisees.txt")
write.table(data_AA$cluster, "cluster_utilises.txt")
write.table(colnames(data_AA)[signif_p_val], "names_result_post_hoc.txt")
write.table(namegrp1, "names_grp1.txt")
write.table(namegrp2, "names_grp2.txt")
write.table(namegrp3, "names_grp3.txt")
write.table(namegrp4, "names_grp4.txt")
setwd("/Users/MAEL/Documents/M1_BI/Stage/research/new_models")


data_grp1 = data_cys_restreint[which(data_AA$cluster == 1),]
data_grp2 = data_cys_restreint[which(data_AA$cluster == 2),]
data_grp3 = data_cys_restreint[which(data_AA$cluster == 3),]
data_grp4 = data_cys_restreint[which(data_AA$cluster == 4),]

#chi2 sur modif selon les grp déterminés
M <- as.table(rbind(summary(as.factor(data_grp1$modif)),
                    summary(as.factor(data_grp2$modif)),
                    summary(as.factor(data_grp3$modif)),
                    summary(as.factor(data_grp4$modif))))
dimnames(M) <- list(group = c("1","2", "3", "4"),
                    modif = c("regroup.","NRg", "NrG", "Nrg", "nRG", "nRg", "nrG"))
M[,1] = M[,1] + M[,2] + M[,3] + M[,5] + M[,7]
M = M[,-c(2,3,5,7)]
print("chi2 : modif~cluster")
(Xsq <- chisq.test(M))  # Prints test summary
Xsq$observed   # observed counts (same as M)
Xsq$expected   # expected counts under the null
Xsq$residuals

M = as.table(rbind(c(111, 109, 66, 143),
                   c(17, 30, 10, 7)
                   ))

#chi2 sur rossman selon les grp déterminés
M <- as.table(rbind(summary(as.factor(data_grp1$Rossman.like)),
                    summary(as.factor(data_grp2$Rossman.like)),
                    summary(as.factor(data_grp3$Rossman.like)),
                    summary(as.factor(data_grp4$Rossman.like))))
dimnames(M) <- list(group = c("1","2", "3", "4"),
                    modif = c("0","1"))
print("chi2 : Rossman~cluster")
(Xsq <- chisq.test(M))  # Prints test summary
Xsq$observed   # observed counts (same as M)
Xsq$expected   # expected counts under the null
Xsq$residuals

#repartition position cys en fonction cluster
hist(data_grp1$cys_pos)
hist(data_grp2$cys_pos)
hist(data_grp3$cys_pos)
hist(data_grp4$cys_pos)

data_cys_restreint$cluster = as.factor(data_AA$cluster)
data_cys_restreint$AS
aov.cys_pos = aov(formula = ASA_cys ~ cluster, data = data_cys_restreint)
summary(aov.cys_pos)
plot(TukeyHSD(aov.cys_pos))

par(mfrow = c(2,2), mar = c(1,1,1,1))
data_grp1 = data_cys_restreint[which(data_AA$cluster == 1),]
data_grp2 = data_cys_restreint[which(data_AA$cluster == 2),]
data_grp3 = data_cys_restreint[which(data_AA$cluster == 3),]
data_grp4 = data_cys_restreint[which(data_AA$cluster == 4),]
#Analyse sur les groupes
#gr1
#data_grp1
vec_col = ifelse(data_grp1$modif == "NRG", 0, 
                 ifelse(data_grp1$modif == "NRg", 1,
                        ifelse(data_grp1$modif == "NrG", 2,
                               ifelse(data_grp1$modif == "Nrg", 3,
                                      ifelse(data_grp1$modif == "nRG", 4,
                                             ifelse(data_grp1$modif == "nRg", 5,6))))))

data_grp1 = data_grp1[,c(2,12:14,18,19,23:202)]
desc_quant_norm = apply(data_grp1, 2, scale)
desc_quant_norm = as.data.frame(desc_quant_norm)
dist_matrix = dist(desc_quant_norm)
Xp = cmdscale(dist_matrix, k = 2)
plot(Xp, col = vec_col+3)

#gr2
#data_grp2
vec_col = ifelse(data_grp2$modif == "NRG", 0, 
                 ifelse(data_grp2$modif == "NRg", 1,
                        ifelse(data_grp2$modif == "NrG", 2,
                               ifelse(data_grp2$modif == "Nrg", 3,
                                      ifelse(data_grp2$modif == "nRG", 4,
                                             ifelse(data_grp2$modif == "nRg", 5,6))))))

data_grp2 = data_grp2[,c(2,12:14,18,19,23:202)]
desc_quant_norm = apply(data_grp2, 2, scale)
desc_quant_norm = as.data.frame(desc_quant_norm)
dist_matrix = dist(desc_quant_norm)
Xp = cmdscale(dist_matrix, k = 2)
plot(Xp, col = vec_col+3)

#gr3
#data_grp3
vec_col = ifelse(data_grp3$modif == "NRG", 0, 
                 ifelse(data_grp3$modif == "NRg", 1,
                        ifelse(data_grp3$modif == "NrG", 2,
                               ifelse(data_grp3$modif == "Nrg", 3,
                                      ifelse(data_grp3$modif == "nRG", 4,
                                             ifelse(data_grp3$modif == "nRg", 5,6))))))

data_grp3 = data_grp3[,c(2,12:14,18,19,23:202)]
desc_quant_norm = apply(data_grp3, 2, scale)
desc_quant_norm = as.data.frame(desc_quant_norm)
dist_matrix = dist(desc_quant_norm)
Xp = cmdscale(dist_matrix, k = 2)
plot(Xp, col = vec_col+3)

summary(data_cys_restreint$modif)
#gr4
#data_grp4
vec_col = ifelse(data_grp4$modif == "NRG", 0, 
                 ifelse(data_grp4$modif == "NRg", 1,
                        ifelse(data_grp4$modif == "NrG", 2,
                               ifelse(data_grp4$modif == "Nrg", 3,
                                      ifelse(data_grp4$modif == "nRG", 4,
                                             ifelse(data_grp4$modif == "nRg", 5,6))))))

data_grp4 = data_grp4[,c(2,12:14,18,19,23:202)]
desc_quant_norm = apply(data_grp4, 2, scale)
desc_quant_norm = as.data.frame(desc_quant_norm)
dist_matrix = dist(desc_quant_norm)
Xp = cmdscale(dist_matrix, k = 2)
plot(Xp, col = vec_col+3)

#summary(data_grp4)
data_cys_restreint$cluster = data_AA$cluster
dev.off()
boxplot(cys_pos ~ cluster, data = data_cys_restreint)

#verification numéro d'accès - groupes

df_test = data.frame(cluster = data_AA$cluster, name = data_cys_restreint$access)
df_test$name = as.character(df_test$name)
df_test$cluster = as.character(df_test$cluster)

ref_name = df_test[1,2]
vec_grp = c(df_test[1,1])
rg_grp = 2
multi_modif_mm_grp = data.frame(name = NA, cluster = NA)
multi_modif_diff_grp = data.frame(name = NA, cluster = NA)
uni_modif = data.frame(name = NA, cluster = NA)
flag = TRUE

for (i in c(2:nrow(df_test))) {
  if(ref_name == df_test[i,2]){
    vec_grp[rg_grp] = df_test[i,1]
  
    if(length(vec_grp) > 1){
      if(vec_grp[rg_grp-1] != vec_grp[rg_grp]){
        flag = FALSE
      }
    }
    rg_grp = rg_grp+1
    
  }else{
    if(length(vec_grp) == 1){
      uni_modif = rbind(uni_modif,
                        c(df_test[i-1,2], df_test[i-1,1]))
    }else if(flag){
      multi_modif_mm_grp = rbind(multi_modif_mm_grp,
                                 c(df_test[i-1,2], 
                                   paste(vec_grp, collapse = " ")) )
    }else{
      multi_modif_diff_grp = rbind(multi_modif_diff_grp,
                                   c(df_test[i-1,2],
                                     paste(vec_grp, collapse = " ")) )
    }
    ref_name = df_test[i,2]
    vec_grp = c(df_test[i,1])
    rg_grp = 2
    flag = TRUE
  }
}

uni_modif
multi_modif_mm_grp
multi_modif_diff_grp

summary(as.factor(data_AA$cluster))

uni_modif$cluster = as.factor(uni_modif$cluster)
summary(uni_modif$cluster)

length(multi_modif_diff_grp[81,2])
