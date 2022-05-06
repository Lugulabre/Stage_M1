setwd("/Users/MAEL/Documents/M1_BI/Stage/research/recherche_perso_cys/analyse_AA_env")

data_cys_restreint$cluster = data_AA$cluster
summary(data_cys_restreint$modif)

#data totale
nb_modif = ifelse(data_cys_restreint$modif == "NRG", 3,
                  ifelse(data_cys_restreint$modif == "NRg" |
                           data_cys_restreint$modif == "NrG" |
                           data_cys_restreint$modif == "nRG",2,1))

data_quant = data_cys_restreint[,c(2,12:14,18,19,23:202)]
desc_quant_norm = apply(data_quant, 2, scale)
desc_quant_norm = as.data.frame(desc_quant_norm)
dist_matrix = dist(desc_quant_norm)
Xp = cmdscale(dist_matrix, k = 2)
plot(Xp, col = nb_modif+1)

summary(as.factor(nb_modif))

par(mfrow = c(2,2), mar = c(1,1,1,1))
data_grp1 = data_cys_restreint[which(data_AA$cluster == 1),]
data_grp2 = data_cys_restreint[which(data_AA$cluster == 2),]
data_grp3 = data_cys_restreint[which(data_AA$cluster == 3),]
data_grp4 = data_cys_restreint[which(data_AA$cluster == 4),]
#Analyse sur les groupes
#gr1
#data_grp1
nb_modif = ifelse(data_grp1$modif == "NRG", 3,
                  ifelse(data_grp1$modif == "NRg" |
                           data_grp1$modif == "NrG" |
                           data_grp1$modif == "nRG",2,1))

data_grp1 = data_grp1[,c(2,12:14,18,19,23:202)]
desc_quant_norm = apply(data_grp1, 2, scale)
desc_quant_norm = as.data.frame(desc_quant_norm)
dist_matrix = dist(desc_quant_norm)
Xp = cmdscale(dist_matrix, k = 2)
plot(Xp, col = nb_modif+1)

#gr2
#data_grp2

nb_modif = ifelse(data_grp2$modif == "NRG", 3,
                  ifelse(data_grp2$modif == "NRg" |
                           data_grp2$modif == "NrG" |
                           data_grp2$modif == "nRG",2,1))

data_grp2 = data_grp2[,c(2,12:14,18,19,23:202)]
desc_quant_norm = apply(data_grp2, 2, scale)
desc_quant_norm = as.data.frame(desc_quant_norm)
dist_matrix = dist(desc_quant_norm)
Xp = cmdscale(dist_matrix, k = 2)
plot(Xp, col = nb_modif+1)

#gr3
#data_grp3

nb_modif = ifelse(data_grp3$modif == "NRG", 3,
                  ifelse(data_grp3$modif == "NRg" |
                           data_grp3$modif == "NrG" |
                           data_grp3$modif == "nRG",2,1))

data_grp3 = data_grp3[,c(2,12:14,18,19,23:202)]
desc_quant_norm = apply(data_grp3, 2, scale)
desc_quant_norm = as.data.frame(desc_quant_norm)
dist_matrix = dist(desc_quant_norm)
Xp = cmdscale(dist_matrix, k = 2)
plot(Xp, col = nb_modif+1)

#gr4
#data_grp4

nb_modif = ifelse(data_grp4$modif == "NRG", 3,
                  ifelse(data_grp4$modif == "NRg" |
                           data_grp4$modif == "NrG" |
                           data_grp4$modif == "nRG",2,1))

data_grp4 = data_grp4[,c(2,12:14,18,19,23:202)]
desc_quant_norm = apply(data_grp4, 2, scale)
desc_quant_norm = as.data.frame(desc_quant_norm)
dist_matrix = dist(desc_quant_norm)
Xp = cmdscale(dist_matrix, k = 2)
plot(Xp, col = nb_modif+1)

dev.off()
