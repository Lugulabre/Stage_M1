setwd("/Users/MAEL/Documents/M1_BI/Stage/research/new_models")
#library(stringr)
#library(corrplot)
#library(nparcomp)

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

#Kmeans
#dist_matrix = dist(data_AA)
#Xp = cmdscale(dist_matrix, k = 2)
#plot(Xp)

cl = kmeans(data_AA, centers = 4)
#Pour réutiliser les clusters du rapport
clust_util = read.table("../recherche_perso_cys/analyse_AA_env/bis/cluster_utilises.txt", header = TRUE)
cl$cluster = as.vector(clust_util[,1])
#Xp = cbind(Xp, cl$cluster)

data_AA = cbind(data_AA, as.factor(cl$cluster))
colnames(data_AA)[30] = "cluster"

#verification numéro d'accès - groupes

df_test = data.frame(cluster = data_AA$cluster, name = data_cys_restreint$access)
df_test$name = as.character(df_test$name)
df_test$cluster = as.character(df_test$cluster)

ref_name = df_test[1,2]
vec_grp = c(df_test[1,1])
rg_grp = 2
multi_modif_mm_grp = matrix(data = NA, nrow = 0, ncol = 7)
multi_modif_diff_grp = matrix(data = NA, nrow = 0, ncol = 11)
uni_modif = matrix(data = NA, nrow = 0, ncol = 2)
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
      lign_sup = c(df_test[i-1,2], vec_grp)
      if(length(lign_sup) < 7){
        lign_sup[7] = NA
      }
      multi_modif_mm_grp = rbind(multi_modif_mm_grp, lign_sup)

    }else{
      lign_sup = c(df_test[i-1,2], vec_grp)
      if(length(lign_sup) < 11){
        lign_sup[11] = NA
      }
      multi_modif_diff_grp = rbind(multi_modif_diff_grp, lign_sup)

    }
    ref_name = df_test[i,2]
    vec_grp = c(df_test[i,1])
    rg_grp = 2
    flag = TRUE
  }
}

uni_modif = as.data.frame(uni_modif)
multi_modif_mm_grp = as.data.frame(multi_modif_mm_grp)
multi_modif_diff_grp = as.data.frame(multi_modif_diff_grp)

row.names(uni_modif) = c(1:nrow(uni_modif))
row.names(multi_modif_mm_grp) = c(1:nrow(multi_modif_mm_grp))
row.names(multi_modif_diff_grp) = c(1:nrow(multi_modif_diff_grp))


summary(as.factor(data_AA$cluster))

colnames(uni_modif) = c("prot", "cluster")
uni_modif$cluster = as.factor(uni_modif$cluster)
summary(uni_modif$cluster)

length(multi_modif_diff_grp[81,2])


#Nombre de modifications selon les groupes
summary(as.factor(multi_modif_mm_grp[3,c(2:7)]))

min(which(is.na(multi_modif_mm_grp[3,c(2:7)])))


mat_cys_grp = matrix(nrow = 2, ncol = 4, data = 0)

for (i in c(1:nrow(uni_modif))) {
  mat_cys_grp[1,as.integer(uni_modif[i,2])] = mat_cys_grp[1,as.integer(uni_modif[i,2])]+1
}

sauv = mat_cys_grp[1,]

for (j in c(1:nrow(multi_modif_mm_grp) )) {
  if( is.na(multi_modif_mm_grp[j,7]) ){
    nb_cys = which(is.na(multi_modif_mm_grp[j,c(2:7)]))[1] - 1
    mat_cys_grp[1,as.integer(multi_modif_mm_grp[j,2])] = mat_cys_grp[1,as.integer(multi_modif_mm_grp[j,2])]+nb_cys
  }else{
    mat_cys_grp[1,as.integer(multi_modif_mm_grp[j,2])] = mat_cys_grp[1,as.integer(multi_modif_mm_grp[j,2])]+6
  }
}

mat_cys_grp[1,] = mat_cys_grp[1,] - sauv

mat_cys_grp[2,] = summary(as.factor(cl$cluster)) - mat_cys_grp[1,] - sauv

mat_cys_grp = as.data.frame(mat_cys_grp)
colnames(mat_cys_grp) = c("1", "2", "3", "4")
rownames(mat_cys_grp) = c("mm_prot", "prot_diff")

(Xsq2 <- chisq.test(mat_cys_grp))
Xsq2$p.value
Xsq2$residuals


#matrice en fonction nombre de modif

mat_cys_nb = matrix(nrow = 3, ncol = 4, data = 0)

for (i in c(1:nrow(uni_modif))) {
  mat_cys_nb[1,as.integer(uni_modif[i,2])] = mat_cys_nb[1,as.integer(uni_modif[i,2])]+1
}

for (j in c(1:nrow(multi_modif_mm_grp) )) {
  if( is.na(multi_modif_mm_grp[j,7]) ){
    nb_cys = which(is.na(multi_modif_mm_grp[j,c(2:7)]))[1] - 1
    mat_cys_grp[1,as.integer(multi_modif_mm_grp[j,2])] = mat_cys_grp[1,as.integer(multi_modif_mm_grp[j,2])]+nb_cys
  }else{
    mat_cys_grp[1,as.integer(multi_modif_mm_grp[j,2])] = mat_cys_grp[1,as.integer(multi_modif_mm_grp[j,2])]+6
  }
}



