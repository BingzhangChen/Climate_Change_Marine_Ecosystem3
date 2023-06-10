library(tidymodels)
library(forcats)
tidymodels_prefer()
load('SCSPico.Rdata')
set.seed(1001)
rmse   <- function(x,y) sd(x-y)

#Split train and test data
pico_split <- initial_split(np, prop = 0.8, strata = logChl)
Train      <- training(pico_split)
Test       <- testing(pico_split)

#Create a validation set that uses 3/4 of the model data for model fitting
val <- validation_split(Train, prop = 3/4)

#Random forests
rf_model <- 
  rand_forest(trees = 1000) %>% 
  set_engine("ranger") %>% 
  set_mode("regression") %>%
  translate()

rf_wflow <- 
  workflow() %>% 
  add_formula(
    logChl ~ lon + lat + Depth + DOY + logChl0 + T0 + I0
      ) %>% 
  add_model(rf_model) 

rf_fit <- rf_wflow %>% 
  fit(data = Train)

## Evaluating model performance ----
test.results <- Test %>%
  select(logChl) %>%
  bind_cols(predict(rf_fit, Test)) 

ggplot(test.results, aes(logChl, .pred)) +
  geom_point(alpha = 1/2) +
  xlab('Observation') +
  ylab('Prediction') +
  geom_abline(col = 'blue') +
  coord_obs_pred() +
  theme_light()

#Correlation between observation and model prediction
cor(test.results$logChl, test.results$.pred)

#Calculate RMSE
rmse(test.results$logChl, test.results$.pred)

##10-fold cross-validation
c_folds <- vfold_cv(Train, v = 10)

#Save the prediction to visualize the model fit and residuals
keep_pred <- control_resamples(save_pred = T, save_workflow = T)

system.time(
  rf_res <- rf_wflow %>%
    fit_resamples(resamples = c_folds, control = keep_pred)
)

#Get the performance metrics; 
#use the option summarize = F to get the metrics for each resample
collect_metrics(rf_res)

#Get the assessment predictions
#To obtain summarized values (averages of the replicate predictions)
#Use the option 'summarize = T'
assess_res <- collect_predictions(rf_res)

#Compare the observed and held-out predicted values
assess_res %>%
  ggplot(aes(x = logChl, y = .pred)) +
  geom_point(alpha = .2) +
  geom_abline(color = 'red') +
  coord_obs_pred() +
  ylab('Predicted') +
  theme_light()

#Use validation set
val_res <- rf_wflow %>%
  fit_resamples(resamples = val)

collect_metrics(val_res)

##Parallel processing
library(doParallel)

#Create a cluster object and then register
cl <- makePSOCKcluster(2)
registerDoParallel(cl)

#Run cross-validation
system.time(
  rf_res <- rf_wflow %>%
    fit_resamples(resamples = c_folds, control = keep_pred)
)
#Not much gain
stopCluster(cl)

##Understanding random forest model -------
vip_features <- c("lon","lat","Depth","DOY","logChl0","T0","I0")
vip_train <- 
  Train %>% 
  select(all_of(vip_features))

library(DALEXtra)
explainer_rf <- 
  explain_tidymodels(
    rf_fit, 
    data = vip_train, 
    y = Train$logChl,
    label = "random forest",
    verbose = FALSE
)

#Global explanations
pico_rf <- model_parts(explainer_rf, loss_function = loss_root_mean_square)

ggplot_imp <- function(...) {
  obj <- list(...)
  metric_name <- attr(obj[[1]], "loss_name")
  metric_lab <- paste(metric_name, 
                      "after permutations\n(higher indicates more important)")
  
  full_vip <- bind_rows(obj) %>%
    filter(variable != "_baseline_")
  
  perm_vals <- full_vip %>% 
    filter(variable == "_full_model_") %>% 
    group_by(label) %>% 
    summarise(dropout_loss = mean(dropout_loss))
  
  p <- full_vip %>%
    filter(variable != "_full_model_") %>% 
    mutate(variable = fct_reorder(variable, dropout_loss)) %>%
    ggplot(aes(dropout_loss, variable)) 
  if(length(obj) > 1) {
    p <- p + 
      facet_wrap(vars(label)) +
      geom_vline(data = perm_vals, aes(xintercept = dropout_loss, color = label),
                 linewidth = 1.4, lty = 2, alpha = 0.7) +
      geom_boxplot(aes(color = label, fill = label), alpha = 0.2)
  } else {
    p <- p + 
      geom_vline(data = perm_vals, aes(xintercept = dropout_loss),
                 linewidth = 1.4, lty = 2, alpha = 0.7) +
      geom_boxplot(fill = "#91CBD765", alpha = 0.4)
    
  }
  p +
    theme(legend.position = "none") +
    labs(x = metric_lab, 
         y = NULL,  fill = NULL,  color = NULL)
}

#Show the importance of features
ggplot_imp(pico_rf)

## Building Global Explanations from Local Explanations
pdp_sst <- model_profile(explainer_rf, N = 500, variables = "T0")

#another function for plotting the underlying data in this object:
ggplot_pdp <- function(obj, x) {
  
  p <- 
    as_tibble(obj$agr_profiles) %>%
    mutate(`_label_` = stringr::str_remove(`_label_`, "^[^_]*_")) %>%
    ggplot(aes(`_x_`, `_yhat_`)) +
    geom_line(data = as_tibble(obj$cp_profiles),
              aes(x = {{ x }}, group = `_ids_`),
              linewidth = 0.5, alpha = 0.05, color = "gray50")
  
  num_colors <- n_distinct(obj$agr_profiles$`_label_`)
  
  if (num_colors > 1) {
    p <- p + geom_line(aes(color = `_label_`), linewidth = 1.2, alpha = 0.8)
  } else {
    p <- p + geom_line(color = "midnightblue", linewidth = 1.2, alpha = 0.8)
  }
  
  return(p)
}

ggplot_pdp(pdp_sst, T0)  +
  labs(x = "Temperature", 
       y = "Chl (log)", 
       color = NULL)

#parameters to be optimized
#ntree: number of trees to grow
#mtry: number of variables randomly selected
#nperm: number of permutations
NTs <- c(500, 1000, 2000)
MTs <- c(3,6)

MEANR2 <- array(NA, dim=c(12, length(NTs), 
                          length(MTs)))
SDR2   <- MEANR2

for (m in 1:length(NTs)){
  for (n in 1:length(MTs)){
      rf.result <- foreach(i=1:10,.combine='rbind') %do% {
        x      <- sample(rownames(np), 0.5*nrow(np))
        y      <- setdiff(rownames(np),x)         #Find the elements of rownames(NP1) which did not belong to x
        Train  <- np[x,]   #Data for training
        Test   <- np[y,]
        c_rf   <- randomForest(logChl ~ lon + lat + 
                                 Depth + DOY + 
                                 logChl0 + T0 + I0,
                               ntree=NTs[m],
                               mtry =MTs[n],
                               data =Train)
        p_rf   <- randomForest(logPro ~ lon + lat + 
                                 Depth + DOY + 
                                 logChl0 + T0 + I0,
                               ntree=NTs[m],
                               mtry =MTs[n],
                               data =Train)
        s_rf   <- randomForest(logSyn ~ lon + lat + 
                                 Depth + DOY + 
                                 logChl0 + T0 + I0,
                               ntree=NTs[m],
                               mtry =MTs[n],
                               data =Train)
        e_rf   <- randomForest(logPeuk ~ lon + lat + 
                                 Depth + DOY + 
                                 logChl0 + T0 + I0,
                               ntree=NTs[m],
                               mtry =MTs[n],
                               data =Train)
        # print(p_rf)
        # imp    <- varImpPlot(p_rf)
        # impvar <- rownames(imp)[order(imp[, 1], decreasing = TRUE)]
        c.p    <- predict(c_rf,Test)
        p.p    <- predict(p_rf,Test)
        s.p    <- predict(s_rf,Test)
        e.p    <- predict(e_rf,Test)
        c(cor(c.p,Test$logChl)**2,
          cor(p.p,Test$logPro)**2,
          cor(s.p,Test$logSyn)**2,
          cor(e.p,Test$logPeuk)**2,
          rmse(c.p,Test$logChl),
          rmse(p.p,Test$logPro),
          rmse(s.p,Test$logSyn),
          rmse(e.p,Test$logPeuk),
          mean(c.p-Test$logChl),
          mean(p.p-Test$logPro),
          mean(s.p-Test$logSyn),
          mean(e.p-Test$logPeuk))
      }
      MEANR2[,m,n]=apply(rf.result,2,mean)
      SDR2[,m,n]  =apply(rf.result,2,sd)	
      save(MEANR2, SDR2, file = 'RF_optimsetting3.Rdata')
  }
}

load('RF_optimsetting3.Rdata')
#plot out
pdf("FigS4_RF_optim.pdf", width=6, height=8)
op   <- par(font.lab=1,
            family ="serif",
            mar    = c(4,4,2,2),
            mgp    = c(2,1,0),
            mfcol  = c(4,2),
            cex.lab= 1.2,
            cex.axis=1.2)

txts <- c('Chl', 'Pro', 'Syn', 'Peuk')
txts <- paste0(letters[1:8],') ',rep(txts,2))
nt   <- 0
for (i in 1:length(MTs)){
  for (j in 1:4){
    
    nt <- nt + 1
    #Chl fitting
    plot(NTs, MEANR2[j,,i], 
         type = 'b',
         ylim = c(0.6,0.9),
         xlab = 'Number of trees',
         ylab = bquote(R^2))
    
    #Add CI
    segments(NTs, MEANR2[j,,i] - 2*SDR2[j,,i],
             NTs, MEANR2[j,,i] + 2*SDR2[j,,i])
    
    mtext(paste0(txts[nt], ', mtry=', MTs[i]), side=3, adj = 0)
  }
}
dev.off()

 
# print(rf.chl)
# varImpPlot(rf.chl)
# 
 
# print(rf.syn)
# imp    <- varImpPlot(rf.syn)
# impvar <- rownames(imp)[order(imp[, 1], decreasing = TRUE)]
# 
# print(rf.peu)
# imp    <- varImpPlot(rf.peuk)
# impvar <- rownames(imp)[order(imp[, 1], decreasing = TRUE)]
# 
# pdf("partialPlot_Chl_RF.pdf",
#     width = 12,
#     height = 12)
# op     <- par(
#   font.lab = 1,
#   family = "serif",
#   mar    = c(3, 3, 2, .2),
#   mgp    = c(2, 1, 0),
#   mfrow  = c(3, 3)
# )
# 
# imp    <- varImpPlot(rf.chl)
# for (i in 1:nrow(imp)) {
#   partialPlot(rf.chl,
#               np,
#               rownames(imp)[i],
#               xlab = rownames(imp)[i],
#               main = "")
# }
# dev.off()
# 
# pdf("partialPlot_Pro_RF.pdf",
#     width = 12,
#     height = 12)
# op     <- par(
#   font.lab = 1,
#   family = "serif",
#   mar    = c(3, 3, 2, .2),
#   mgp    = c(2, 1, 0),
#   mfrow  = c(3, 3)
# )
# 
# imp    <- varImpPlot(rf.pro)
# for (i in 1:nrow(imp)) {
#   partialPlot(rf.pro,
#               np,
#               rownames(imp)[i],
#               xlab = rownames(imp)[i],
#               main = "")
# }
# dev.off()
# 
# pdf("partialPlot_Syn_RF.pdf",
#     width = 12,
#     height = 12)
# op     <- par(
#   font.lab = 1,
#   family = "serif",
#   mar    = c(3, 3, 2, .2),
#   mgp    = c(2, 1, 0),
#   mfrow  = c(3, 3)
# )
# 
# imp    <- varImpPlot(rf.syn)
# for (i in 1:nrow(imp)) {
#   partialPlot(rf.syn,
#               np,
#               rownames(imp)[i],
#               xlab = rownames(imp)[i],
#               main = "")
# }
# dev.off()
# 
# pdf("partialPlot_Peuk_RF.pdf",
#     width = 12,
#     height = 12)
# op     <- par(
#   font.lab = 1,
#   family = "serif",
#   mar    = c(3, 3, 2, .2),
#   mgp    = c(2, 1, 0),
#   mfrow  = c(3, 3)
# )
# 
# imp    <- varImpPlot(rf.peu)
# for (i in 1:nrow(imp)) {
#   partialPlot(rf.peu,
#               np,
#               rownames(imp)[i],
#               xlab = rownames(imp)[i],
#               main = "")
# }
# dev.off()
