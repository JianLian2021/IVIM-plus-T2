library(tidyverse)
library(haven)
library(caret)
library(rms)
library(Boruta)
library(pROC)
library(ggpubr)
library(gridExtra)
library(mRMRe)
library(DMwR2)
library(Publish)
library(qqman)


set.seed(1)

setwd("C:/Users/223024351/Desktop/Huashan2")


dt_sel <- read_csv('dt.csv')#%>% select(-c(1:3))


dt_sel$Label <- dt_sel$Label %>% 
  factor(ordered = T, labels = c('Low', 'High'))




if(!is_empty(nearZeroVar(dt_sel)))
{
  dt_sel <- dt_sel[, -nearZeroVar(dt_sel)]
}



s_pre <- preProcess(dt_sel, method = 'medianImpute')
dt_sel <- predict(s_pre, dt_sel)


##############???????????????:????????????


utest_str <- paste(colnames(dt_sel)[1], paste(colnames(dt_sel)[-1], collapse = '+'), sep = '~') %>% as.formula()



utest_res <- univariateTable(utest_str, data = dt_sel)%>% summary() 
utest_res1 <- univariateTable(utest_str, data = dt_sel)#%>% summary() 
p_name <- utest_res$Variable[which(utest_res$`p-value`< 0.05)]



#p_name <- utest_res$Variable[which(utest_res$p.values < 0.05)]


write_csv(utest_res, 'utest_res.csv')


#######################??????????????????????????????p???????????????



dt_man <- tibble(P = utest_res1$p.values, SNP = names(utest_res1$p.values), 
                 BP = c(1:42, 1:108, 1:226, 1:20), 
                 CHR = c(rep(1, 42), rep(2, 108), rep(3, 226), rep(4, 20)))
                                                        
#pdf('manhattan.pdf', height = 6, width = 10)
manhattan(dt_man, col = get_palette(palette = 'default', k = 4), suggestiveline = -log10(1e-2), 
          xlab = 'Feature', cex = 2)
legend(x = 10, y = 5, legend = c('Histogram', 'GLCM', 'GLRLM', 'Shape'), 
       col =  get_palette(palette = 'default', k = 4), pch = 16 )


#########################?????????????????????,??????????????????????????????????????????p<0.05?????????


dt_sel <- select(dt_sel, c('Label', p_name))



s_pre2 <- preProcess(dt_sel, method = c('center', 'scale'))
dt_sel <- predict(s_pre2, dt_sel)



res_ulogit <- glmSeries(Label~1, data = dt_sel, vars = colnames(dt_sel)[-1], 
                        family = 'binomial')



res_ulogit <- filter(res_ulogit, Pvalue < 0.05)
write.csv(res_ulogit, file = 'ulogit.csv')



dt_sel <- select(dt_sel, c('Label',res_ulogit$Variable))
dt_sel <- na.omit(dt_sel)



########????????????????????????,??????mrmr?????????????????????
# select the feature further.
dt_mrmr <- mRMR.data(dt_sel)
f_sel <- mRMR.classic(data = dt_mrmr, target_indices = c(1), feature_count = 5)



useful_name <- featureNames(f_sel)[unlist(solutions(f_sel))]



dt_sel <- select(dt_sel, c('Label', useful_name))



# Show the ROCs of the remaining features after ulogit test


for_str <- paste('Label', paste(colnames(dt_sel)[-1], collapse = '+'), sep = '~')


roc_list <- roc(as.formula(for_str), data = dt_sel, ci = T)

auc_vec <- lapply(roc_list, pROC::auc) %>% unlist

roc_list <- roc_list[order(auc_vec, decreasing = T)]
col_vec <- get_palette(palette = 'jco', k = 5)
y_vec <- seq(0.2, 0.02, length.out = 5)

res_tbl <- tibble()

varName_roc <- names(roc_list)
for(i_a in 1:length(roc_list))
{
  if(i_a == 1)
  {
    plot(roc_list[[i_a]], print.auc = T, print.auc.pattern = paste(varName_roc[i_a], 'AUC: %.2f (%.2f - %.2f)'), 
         print.auc.y = y_vec[i_a], col = col_vec[i_a], legacy.axes = T, 
         print.auc.x = 0.8)
  }
  else
  {
    plot(roc_list[[i_a]], print.auc = T, print.auc.pattern = paste(varName_roc[i_a], 'AUC: %.2f (%.2f - %.2f)'), 
         print.auc.y = y_vec[i_a], add = T, col = col_vec[i_a], 
         print.auc.x = 0.8)
  }
  print(varName_roc[i_a])
  cutoff <- coords(roc_list[[i_a]], x = 'best', transpose = T)
  if(!is_empty(dim(cutoff)))
  {
    cutoff <- cutoff[, 1]
  } 
  cutoff <- cutoff[1]
  print(cutoff)
  
  res_bin <- ifelse(dt_sel[[varName_roc[i_a]]] > cutoff, 'High', 'Low') %>% factor(ordered = T)
  cmat <- confusionMatrix(res_bin, dt_sel$Label, positive = 'High')
  
  if(cmat$overall[1] < 0.5)
  {
    res_bin <- ifelse(dt_sel[[varName_roc[i_a]]] < cutoff, 'High', 'Low') %>% factor(ordered = T)
    cmat <- confusionMatrix(res_bin, dt_sel$Label, positive = 'High')
  }
  cmat_res <- c(cutoff, cmat$overall[c(1, 3, 4)], cmat$byClass[c(1:4)])
  
  res_tbl <- bind_rows(res_tbl, cmat_res)
}



res_tbl$VarName <- varName_roc



write_csv(res_tbl, path = 'variable_res.csv')

################################????????????????????????,???????????????????????????10?????????????????????????????????
###############################3
mytwoclasssummary <- function(data, lev = NULL, model = NULL)
{
  
  cmat <- confusionMatrix(data$pred, data$obs, positive = lev[2])
  
  roc_ <- pROC::roc(data$obs, data[, lev[1]], quiet = T)
  
  out <- c(pROC::auc(roc_), cmat$overall[1], cmat$byClass[c(1:4)])
  
  names(out) <- c('ROC', 'ACC', 'Sensitivity', 'Specificity', 'PPV', 'NPV')
  out
}



trainC <- trainControl(method = 'repeatedcv', number = 10, repeats = 10,
                       summaryFunction = mytwoclasssummary, 
                       classProbs = T)
#make.names(dt_train_final$Label,unique = FALSE, allow_ = TRUE )
fit1 <- train(Label~., data = dt_sel, method = 'glm',
              trControl = trainC, metric = 'Accuracy')



# fit1 <- step(fit_glm)



res_resample <- write_csv(fit1$resample, "resample.csv")#####???10??????????????????????????????
#res_data <- read_csv('resample_ADC.csv')
#resample <- select(fit_glm$resample[,]
res_data <- gather(fit1$resample[-7])
res_data$Model = "T2"                    
write.csv(res_data, "res_data2.csv")

ggboxplot(data = res_data, x = 'Model', y = 'value', fill = 'Model', alpha = 0.5) %>% 
  facet(facet.by = 'key') + theme(axis.text.x = element_blank())
ggsave('metrics.pdf', width = 12, height = 10)



mlogit <- publish(fit1$finalModel)
write_csv(mlogit, "mlogit.txt")
#######################33



final_formula <- paste(fit1$finalModel$coefficients, names(fit1$finalModel$coefficients), sep = '*') %>% paste(collapse = '+')
final_formula <- paste('Radscore', final_formula, sep = '=')



print(final_formula)
write_file(final_formula, path = 'radscore.txt')



##
dt_final <- dt_sel
res <- predict(fit1, data = dt_final, type = "prob")
res <- res[,-1]



res_roc <- roc(dt_sel$Label, res, ci = T)
cutoff <- coords(res_roc, x = 'best')[1]
res_bin <- ifelse(res > cutoff, 'High', 'Low') %>% as.factor



#res1 <- cbind.data.frame(res, dt_final$Label)
write.csv(res, "res.csv")



cmat <- confusionMatrix(res_bin, dt_sel$Label, positive = 'High')



res_res <- c(cutoff, cmat$overall[c(1, 3, 4)], cmat$byClass[1:4])



res_tbl_final <- bind_rows(res_res)



write_csv(res_tbl_final, path = 'res_mlogit.csv')



plot(res_roc, print.auc = T, print.auc.pattern = paste('AUC: %.2f (%.2f - %.2f)'), 
     print.auc.y = 0.3, col = col_vec[1], legacy.axes = T)



rmarkdown::render('results_final_texture.Rmd')

