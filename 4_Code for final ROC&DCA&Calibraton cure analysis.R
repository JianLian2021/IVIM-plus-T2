library(tidyverse)
library(caret)
library(pROC)
library(glmnet)
library(DMwR2)
library(rmda)
library(ggpubr)
library(ModelGood)
library(rms)
library(mRMRe)
library(DescTools)


A<- read.csv("dt.csv")
a<- c("IVIM_")
b<- colnames(A)
c<- paste(a, b, sep ="")
c
colnames(A)<- c
write.csv(A, "IVIM_.csv")



B<- read.csv("T2.csv")
d<- c("T2_")
e<- colnames(B)
f<- paste(d, e, sep ="")
f
colnames(B)<- f
write.csv(B, "T2_.csv")




A<- read.csv("Radscore.csv")

roc1<- roc(A$Label,A$T2,ci= T)
roc2<- roc(A$Label,A$IVIM,ci= T)
roc3<- roc(A$Label,A$Combined,ci= T)
plot(roc1, col= "green", lty= 1, lwd= 3, print.auc= T, print.auc.pattern= "AUC: %.2f (%.2f - %.2f)",
     legacy.axes= T, print.auc.y= .5)
plot(roc2, col= "blue", lty= 1, lwd= 3, print.auc= T, print.auc.pattern= "AUC: %.2f (%.2f - %.2f)",
     legacy.axes= T, add= T, print.auc.y= .45)
plot(roc3, col= "red", lty= 1, lwd= 3, print.auc= T, print.auc.pattern= "AUC: %.2f (%.2f - %.2f)",
     legacy.axes= T, add= T, print.auc.y= .40)
legend(x= .5, y= .25, lty= 1, lwd= 3, col = c("green", "blue","red"), legend = c("T2","IVIM", "Combined"))
roc.test(roc1,roc2)
roc.test(roc1,roc3)
roc.test(roc2,roc3)

cal1<- glm(Label~T2,data = A, family = binomial)
cal2<- glm(Label~IVIM, data = A, family = binomial)
cal3<- glm(Label~Combined, data = A, family = binomial)
a<- calPlot2(cal1, col = "green", lty = 2, lwd = 3, legend = F)
b<- calPlot2(cal2, col = "blue", lty = 2, lwd = 3, legend = F,add = T)
c<- calPlot2(cal3, col = "red", lty = 2, lwd = 3, legend = F,add = T)

legend(x= .7, y= .2, lty= 3, lwd= 3, col = c("green", "blue", "red"), legend = c("T2","IVIM","Combined"))
HosmerLemeshowTest(a$Frame$glm,a$Frame$jack)
HosmerLemeshowTest(b$Frame$glm,b$Frame$jack)
HosmerLemeshowTest(c$Frame$glm,c$Frame$jack)



dca1<- decision_curve(Label~T2, data = A)
dca2<- decision_curve(Label~IVIM, data = A)
dca3<- decision_curve(Label~Combined, data = A)
plot_decision_curve(list(dca1,dca2,dca3), col = c("green","blue","red" ,"purple", "grey"), 
                    confidence.intervals = F, curve.names = c("T2", "IVIM", "Combined"),
                    lty= 1, lwd= 3)



plot_decision_curve(dca1, col = c("red","purple", "green"), 
                    confidence.intervals = F, curve.names = "V1", 
                    lty= 1, lwd= 2.5)
plot_decision_curve(dca2, col = c("blue","purple", "green"), 
                    confidence.intervals = F, curve.names = "V2", 
                    lty= 1, lwd= 2.5, add= T)















