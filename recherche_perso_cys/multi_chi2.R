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

chi2_p_val = c()
ross_p_val = c()

#Kmeans

for (i in c(1:1000)) {
  data_AA = data_AA[,-30]
  dist_matrix = dist(data_AA)
  Xp = cmdscale(dist_matrix, k = 2)
  #plot(Xp)
  
  cl = kmeans(data_AA, centers = 4)
  
  Xp = cbind(Xp, cl$cluster)
  
  data_AA$cluster = as.factor(cl$cluster)
  #colnames(data_AA)[30] = "cluster"
  
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
  M[,1] = M[,1] + M[,2] + M[,3] + M[,5]
  M[,2] = M[,4] + M[,6] + M[,7]
  M = M[,-c(3:7)]
  #print("chi2 : modif~cluster")
  (Xsq <- chisq.test(M))  # Prints test summary
  chi2_p_val[i] = Xsq$p.value
  
  #chi2 sur rossman selon les grp déterminés
  M2 <- as.table(rbind(summary(as.factor(data_grp1$Rossman.like)),
                      summary(as.factor(data_grp2$Rossman.like)),
                      summary(as.factor(data_grp3$Rossman.like)),
                      summary(as.factor(data_grp4$Rossman.like))))
  dimnames(M2) <- list(group = c("1","2", "3", "4"),
                      modif = c("0","1"))
  #print("chi2 : Rossman~cluster")
  (Xsq2 <- chisq.test(M2))
  ross_p_val[i] = Xsq2$p.value
  print(i)
}

length(chi2_p_val[which(chi2_p_val < 0.05)])
length(ross_p_val[which(ross_p_val < 0.05)])
hist(chi2_p_val, xlim = 0:1)
hist(ross_p_val, xlim = 0:1)
t.test()

