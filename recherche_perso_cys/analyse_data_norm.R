setwd("/Users/MAEL/Documents/M1_BI/Stage/research/new_models")
library(stringr)

data_cys_restreint = read.table("bdd_models_byCys.txt_restreint", header = T)

#En utilisant seulement les descripteurs quantitatifs
dim(data_cys_restreint)
desc_quant = data_cys_restreint[,c(2,12:14,18,19,23:202)]

#Centrer-réduire les données
desc_quant_norm = apply(desc_quant, 2, scale)
desc_quant_norm = as.data.frame(desc_quant_norm)

#Matrice de distance
dist_matrix = dist(desc_quant_norm)
Xp = cmdscale(dist_matrix, k = 2)
plot(Xp)

#Cluster par kmeans
cl = kmeans(desc_quant_norm, centers = 2)
cl$cluster
Xp = cbind(Xp, cl$cluster)
plot(Xp[,1:2], col = Xp[,3]+1)

desc_quant_norm = data.frame(desc_quant_norm, cluster = cl$cluster)

#Choix du nombre de cluster optimal
vec_expl = c()
for (i in c(1:14)) {
  cl = kmeans(Xp[,1:2], centers = i)
  vec_expl[i] = cl$tot.withinss
}
plot(c(1:length(vec_expl)),vec_expl, "l", col = "Blue")

#Recherche des différences entre les deux groupes définis par kmeans
vec_t_test = c()
mat_moy = matrix(nrow = ncol(desc_quant_norm)-1, ncol = 2)
mat_moy = as.data.frame(mat_moy)
for (i in c(1: (ncol(desc_quant_norm)-1) ) ) {
  ech1 = desc_quant_norm[,i][which(desc_quant_norm$cluster == 1)]
  ech2 = desc_quant_norm[,i][which(desc_quant_norm$cluster == 2)]
  vec_t_test[i] = t.test(ech1, ech2)[[3]]
  mat_moy[i,1] = t.test(ech1, ech2)[[5]][[1]]
  mat_moy[i,2] = t.test(ech1, ech2)[[5]][[2]]
}

#Descripteurs entre les 2 groupes kmeans avec moyenne différente à 5%
colnames(desc_quant_norm)[which(vec_t_test <= 0.05)]

mat_write = matrix(nrow = length(vec_t_test), ncol = 3)
mat_write = as.data.frame(mat_write)
for (i in which(vec_t_test <= 0.05)) {
  mat_write[i,1] = colnames(desc_quant_norm)[i]
  mat_write[i,2] = ifelse(mat_moy[i,1]<mat_moy[i,2], 1, 2)
  mat_write[i,3] = ifelse(mat_moy[i,1]<mat_moy[i,2], 2, 1)
}
mat_write = na.omit(mat_write)

aa_pol_chrg = c("D", "E", "H", "R", "Y", "C", "K")
aa_pol_chrg_full = c("ASP", "GLU", "HIS", "TYR", "CYS", "LYS")

aa_phobes = c("V", "L", "A", "M", "W", "I", "F")
aa_phobes_full = c("VAL", "LEU", "ALA", "MET","TRP", "ILE", "PHE")

mat_write = data.frame(mat_write, pol_chr = rep(NA, nrow(mat_write)), 
                       hydrophobe = rep(NA, nrow(mat_write)) )

for (i in c(1:nrow(mat_write))) {
  desc_name = mat_write[i,1]
  if( is.na(as.integer(substr(desc_name, 2,3))) ){
    mat_write[i,4] = ifelse( length(which(str_detect(aa_pol_chrg_full, desc_name))) == 0,
                             0, 1)
    mat_write[i,5] = ifelse( length(which(str_detect(aa_phobes_full, desc_name))) == 0,
                             0, 1)
  }else{
    desc_name =  substr(desc_name, 1,1)
    mat_write[i,4] = ifelse( length(which(str_detect(aa_pol_chrg, desc_name))) == 0,
                             0, 1)
    mat_write[i,5] = ifelse( length(which(str_detect(aa_phobes, desc_name))) == 0,
                             0, 1)
  }
  
}

mat_write_clone = mat_write[-c(1:5),]
summary(as.factor(mat_write$pol_chr[which(mat_write_clone$V2 == 1)]))
summary(as.factor(mat_write$pol_chr[which(mat_write_clone$V2 == 2)]))

summary(as.factor(mat_write$hydrophobe[which(mat_write_clone$V2 == 1)]))
summary(as.factor(mat_write$hydrophobe[which(mat_write_clone$V2 == 2)]))

setwd("/Users/MAEL/Documents/M1_BI/Stage/research/recherche_perso_cys/analyse_data_norm")
write.table(mat_write, "moy_comp.txt")
setwd("/Users/MAEL/Documents/M1_BI/Stage/research/new_models")

hist(data_cys_restreint$cys_pos)
summary(data_cys_restreint$cys_pos)
summary(data_cys_restreint$modif)

#Cystéines selon les différents types de résidus



#Dendrogramme
hc=hclust(dist_matrix,method="ward.D2")
plot(hc, main="gene")
cl2 = cutree(hc, k = 2)

Xp = cmdscale(dist_matrix, k = 2)
Xp_hclust = cbind(Xp, cl2)
plot(Xp_hclust[,1:2], col = Xp_hclust[,3]+1)
