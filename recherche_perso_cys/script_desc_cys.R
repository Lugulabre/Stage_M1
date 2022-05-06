setwd("/Users/MAEL/Documents/M1_BI/Stage/research/new_models")

data_cys_corps = read.table("bdd_models_byCys_copie.txt", header = F)
data_cys_header = read.table("bdd_models_byCys_copie_header.txt", header = T)
data_cys_l1 = read.table("bdd_models_byCys_copie_ligne_un.txt", header = F)

data_cys_corps= na.omit(data_cys_corps)

library(FactoMineR)

head(data_cys_corps)

test_type = c()

for (i in 1:ncol(data_cys_corps)) {
  test_type[i] =  typeof(data_cys_corps[,i])
}

test_type
head(data_cys_corps)
var = apply(data_cys_corps, 2, var)
which(is.na(var))

data_cys_corps = data_cys_corps[,-c(which(is.na(var)))]

#PCA(data_cys_corps)

#------------RESTREINT------------
data_cys_corps_restreint = read.table("bdd_models_byCys.txt_restreint_copie")
data_cys_corps_restreint = na.omit(data_cys_corps_restreint)
head(data_cys_corps_restreint)

names_cys = colnames(data_cys_header)
cys_corps = data.frame(data_cys_corps_restreint)
names(cys_corps) = names_cys
names(data_cys_corps_restreint) = names_cys

summary(cys_corps$patch_score)
summary(cys_corps$patch_name)


var2 = apply(cys_corps, 2, var)
cys_corps = cys_corps[,-c(which(is.na(var2)))]
head(cys_corps)
summary(cys_corps$cys_modelled)

#Suppression des descripteurs à variance nulle
i = 1
while (i <= (ncol(cys_corps)-1)) {
  if (var(cys_corps[,i]) == 0){
    cys_corps = cys_corps[,-i]
    i = i-1;
  }
  i = i+1;
}

#PCA(cys_corps)
#boxplot(cys_corps)

summary(cys_corps$pKa)
summary(cys_corps$cys_pos)
summary(data_cys_corps_restreint$modif)

#cys_corps = scale(cys_corps)

cys_corps = data.frame(cys_corps, modif = data_cys_corps_restreint$modif)
cys_corps$modif = ifelse(cys_corps$modif == "nRg", 1,0)

boxplot(cys_corps)

#correlation
library(corrplot)
library(ggplot2)
library(lattice)
library(caret)
matCor = cor(cys_corps[,-192])
cys_corps = cys_corps[,-findCorrelation(matCor, cutoff = 0.9)]

cys_corps$modif = as.factor(cys_corps$modif)

#Echantillon d'apprentissage
vIndApp = sample(nrow(cys_corps), size = (2*nrow(cys_corps)/3),
                 replace = FALSE)
matApp = cys_corps[vIndApp,]
#Echantillon de test
vIndTest = setdiff(1:nrow(cys_corps), vIndApp)
matTest = cys_corps[vIndTest,]

#regression
fit = glm(modif ~., data = as.data.frame(matApp), family = "binomial")
attributes(fit)
summary(fit)
#Y11 W11 T11 S11 P11 F11 M11 K11 L11 I11 H11 G11 E11 Q11 C11 D11 N11
#R11 A11
summary(cys_corps$Y11)
