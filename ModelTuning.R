
srf_survival_prob_for_time = function(rfmodel, df_to_predict, fixed_time, oob = FALSE){
  #finds survival prediction from srf model 
  
    if (oob) { predicted_matrix  =  predict(rfmodel, newdata = df_to_predict, ensemble = "oob", outcome= "test")
      } else { predicted_matrix  =  predict(rfmodel, newdata = df_to_predict)}
    
  j_for_fixedtime = match(1, round(predicted_matrix$time.interest,1) == fixed_time, nomatch = -100);
    
  if (j_for_fixedtime == -100){#print("no fixed time match was found, using closest")
      j_for_fixedtime = which.min(abs(predicted_matrix$time.interest - fixed_time))}
  
  if (oob) { y_predicted = predicted_matrix$survival.oob[ ,j_for_fixedtime]
      }else{ y_predicted = predicted_matrix$survival[, j_for_fixedtime]}
  
  
  return(y_predicted)
}

model_stats = function(y_predict, y_actual){
  
  #creates a list of metrics as below for predict/actual vectors
  # c(AUCROC, k_youden, TN_pred_act_00, FN_pred_act_01, TP_pred_act_11, FP_pred_act_10, 
  #    Balanced Accuracy, Accuracy, Sensitivity, Specificity, Precision) 
  
  #AUC and 95%CI for AUC
  proc_auc = pROC::roc(y_actual, y_predict)
  ciauc = pROC::ci.auc(proc_auc)
  model_auc = proc_auc$auc[1]
  
  # confusion matrix at Youden point
  k_youden = proc_auc$thresholds[which.max(proc_auc$sensitivities + proc_auc$specificities)]
  confm = caret::confusionMatrix(as.factor(y_actual), as.factor(ifelse(y_predict>k_youden,1,0)))
  #alternative ModelMetrics::confusionMatrix(as.double(df_tune$event_f), y_predicted, cutoff= k_youden)
  #f1score = ModelMetrics::f1Score(as.double(df_tune$event_f), y_predicted, cutoff = k_youden)
  metrix_model =  c("AUCROC" = model_auc,
                    "AUCROC_95CIlow" = ciauc[1],
                    "AUCROC_95CIhigh" = ciauc[3],
                    "k_youden" = k_youden,
                    "Brier" = ModelMetrics::brier(y_actual, y_predict, cutoff = k_youden),
                    "TN_pred_act_00" = confm$table["0","0"],
                    "FN_pred_act_01" = confm$table["0","1"],
                    "TP_pred_act_11" = confm$table["1","1"],
                    "FP_pred_act_10" = confm$table["1","0"],
                    confm$byClass["Balanced Accuracy"] , 
                    confm$overall["Accuracy"], 
                    confm$byClass["Sensitivity"],
                    confm$byClass["Specificity"], 
                    confm$byClass["Precision"])
}

####################################### TUNE SRF ############################################

#FUNCTION
srf_tune = function(df_tune ,  cv_number =5, predict.factors, fixed_time = NaN, seed_to_fix = 100,
                     mtry= c(3,4,5), nodesize = c(10,20,50), nodedepth = c(100), verbose = FALSE, oob = TRUE){
  #fine tuning survival random forest
  # if oob = TRUE, there is no CV !!! as OOB does the job already 
  #set seed
  set.seed(seed_to_fix)
  # limit mtry with the number of predictors and nodesize by 1/4 of sample size
  if (sum(nodesize > dim(df_tune)[1]/6) > 0) {print ("Warning - some min nodesize is > 1/6 of the sample size (1/2 of CV fold)")}
  nodesize = nodesize[nodesize <= dim(df_tune)[1]/6]
  #if all numbers higher than number of factors, only check this factor
  if (sum(mtry > length(predict.factors)) == length(mtry)) { mtry = c(length(predict.factors))}
  mtry = mtry[mtry <= length(predict.factors)]

  #grid of values to tune
  grid_of_values = expand.grid("mtry" = mtry, "nodesize" = nodesize, "nodedepth" = nodedepth)
  print(paste("Grid size is", dim(grid_of_values)[1]))
  
  if  (dim(grid_of_values)[1]==0) {
        output = list()
        return (output)}
  
  # defining output for fixed_time
  if (is.nan(fixed_time)){fixed_time = round(quantile(df_tune$time, 0.8),1)}
  df_tune$time_f = ifelse(df_tune$time <= fixed_time ,df_tune$time, fixed_time)
  df_tune$event_f =ifelse(df_tune$event==1 & df_tune$time <= fixed_time ,1,0)

  #going through combinations
  modelstats = list( ) 
  
  if (oob==FALSE) { #we do CV instead of using OOB predictions to tune SRF 
          set.seed(seed_to_fix)
          cv_folds = caret::createFolds(df_tune$event, k=cv_number, list = FALSE) #use caret to split into k-folds = cv_number
          for (i in 1:dim(grid_of_values)[1]){
              print(grid_of_values[i,])  
              mtry_i = grid_of_values[i,"mtry"]
              nodesize_i = grid_of_values[i,"nodesize"]
              nodedepth_i = grid_of_values[i,"nodedepth"]

              #train grid combination for each cv_iteration
              modelstats_cv = list()
              for (cv_iteration in 1:cv_number) {
                 df_train_cv = df_tune[cv_folds != cv_iteration, ]
                 df_test_cv = df_tune[cv_folds == cv_iteration, ]
                 #train SRF
                 rf.dt = rfsrc(as.formula(paste("Surv(time, event) ~", 
                        paste(predict.factors, collapse="+"))),
                        data = df_train_cv,
                        nodesize = nodesize_i,  # this is AVERAGE size, so we want this to be quite high
                        ntree = 300,
                        mtry = mtry_i, 
                        nodedepth = nodedepth_i,  
                        nsplit = 50, 
                        splitrule = "logrank", statistics= FALSE, membership=TRUE,
                        importance = "none", #to speed up by switching off VIMP calculations
                        seed = seed_to_fix
                        )
                   #compute predicted event probability for all people 
                  y_predicted = 1-srf_survival_prob_for_time(rf.dt, df_test_cv, fixed_time, oob= FALSE)
                  validation_stats = method_any_validate(y_predicted, fixed_time, df_train_cv, df_test_cv, weighted = TRUE)
                  modelstats_cv[[cv_iteration]] =  c("mtry" = mtry_i,  "nodesize" = nodesize_i, "nodedepth" = nodedepth_i, 
                                                 "time" = validation_stats$times_to_predict,
                                                 "auc_score" = validation_stats$auc_score, 
                                                 "brier_score"= validation_stats$brier_score,
                                                 "brier_score_scaled" = validation_stats$brier_score_scaled, 
                                                 "c_score" = validation_stats$c_score, 
                                                 "calibration_alpha" = validation_stats$calibration_alpha, "calibration_slope" = validation_stats$calibration_slope)
                }#end k-fold CV for one grid combination
              
              #averaging over cv-steps, firs transform to data.frame to use mean()
              modelstats_cv_df = data.frame(t(modelstats_cv[[1]]))
              for (j in 2:cv_number) {modelstats_cv_df = rbind(modelstats_cv_df,t(modelstats_cv[[j]]))}
              modelstats[[i]] = c(modelstats_cv[[1]]["mtry"],  modelstats_cv[[1]]["nodesize"], 
                                   modelstats_cv[[1]]["nodedepth"], 
                                  "auc_score" = mean(modelstats_cv_df$auc_score, na.rm=1), 
                                  "brier_score" = mean(modelstats_cv_df$brier_score, na.rm=1),
                                  "brier_score_scaled" = mean(modelstats_cv_df$brier_score_scaled, na.rm=1),
                                  "c_score"= mean(modelstats_cv_df$c_score, na.rm=1),
                                  "calibration_alpha" =  mean(modelstats_cv_df$calibration_alpha, na.rm=1),
                                  "calibration_slope" =  mean(modelstats_cv_df$calibration_slope, na.rm=1),
                                  "time"= fixed_time)
              
            }#end for grid

  } else { # end if(oob==false) 
          for (i in 1:dim(grid_of_values)[1]){
            print(grid_of_values[i,])  
            mtry_i = grid_of_values[i,"mtry"]
            nodesize_i = grid_of_values[i,"nodesize"]
            nodedepth_i = grid_of_values[i,"nodedepth"]
            
            rf.dt = rfsrc(as.formula(paste("Surv(time, event) ~", 
                                           paste(predict.factors, collapse="+"))),
                          data = df_tune,
                          nodesize = nodesize_i,  # this is AVERAGE size, so we want this to be quite high
                          ntree = 300,
                          mtry = mtry_i, 
                          nodedepth = nodedepth_i,  
                          nsplit = 50, 
                          splitrule = "logrank", statistics= FALSE, membership=TRUE,
                          importance = "none", #to speed up by switching off VIMP calculations
                          seed = seed_to_fix)
            
            #compute predicted event probability for all people 
            y_predicted = 1-srf_survival_prob_for_time(rf.dt, df_tune, fixed_time, oob= TRUE)
            validation_stats = method_any_validate(y_predicted, fixed_time, df_tune, df_tune, weighted = TRUE)
            modelstats[[i]] =  c("mtry" = mtry_i,  "nodesize" = nodesize_i, "nodedepth" = nodedepth_i, 
                                 "time" = validation_stats$times_to_predict,
                                 "auc_score" = validation_stats$auc_score, 
                                 "brier_score"= validation_stats$brier_score,
                                 "brier_score_scaled" = validation_stats$brier_score_scaled, 
                                 "c_score" = validation_stats$c_score, 
                                 "calibration_alpha" = validation_stats$calibration_alpha, 
                                 "calibration_slope" = validation_stats$calibration_slope)
                                               
          } #end for (i in grid)
  }#end else 
    
    #reshaping into data frame 
    df_modelstats = data.frame("V1" = modelstats[[1]])
    #check if there was more than 1 grid search
    if (dim(grid_of_values)[1] >1) {for (i in 2:dim(grid_of_values)[1]){ df_modelstats[i]= modelstats[[i]]}}
    df_modelstats = data.frame(t(df_modelstats))
  
  if (verbose == TRUE) {
        print(paste("AUC varies from", round(min(df_modelstats$auc_score ),4), "to", round(max(df_modelstats$auc_score ),4)))
        print(paste("Brier score varies from", round(min(df_modelstats$brier_score ),4), "to", round(max(df_modelstats$brier_score ),4)))
        print("Combination with highest AUC")
        print(df_modelstats[which.max(df_modelstats$auc_score),  c("mtry", "nodesize", "nodedepth")])
        print("Combination with lowest Brier Score")
        print(df_modelstats[which.min(df_modelstats$brier_score),  c("mtry", "nodesize", "nodedepth")])
        print("Combination with lowest AUC")
        df_modelstats[which.min(df_modelstats$auc_score), c("mtry", "nodesize", "nodedepth")]
  }
  output = list() 
  output$modelstats = df_modelstats 
  output$bestbrier = df_modelstats[which.min(df_modelstats$brier_score), ]
  output$bestauc = df_modelstats[which.max(df_modelstats$auc_score), ]
  output$bestcindex = df_modelstats[which.max(df_modelstats$c_score), ]
  return(output)
}


method_3_tune = function(df_tune, params_for_tree= c("age", "bmi", "hyp"), predict.factors, fixed_time=fixed_time,
                  n_cv = 3, maxdepth = maxdepth, cp =0.001, minbucket = minbucket, seed_to_fix = 100){
  #tuning ensemble 3 - fudging the tree in maxdepth, minbucket and risk factors to use to CV,
  # this is to get an average test c-index for one combination, all combinations are checked in method_3_cv() function
  set.seed(seed_to_fix)
  cv_folds = caret::createFolds(df_tune$event, k=n_cv, list = FALSE) #use caret to split into k-folds = cv_steps

  cindex_train = vector(length = n_cv); cindex_test = vector(length = n_cv)
  for (cv_iteration in 1:n_cv){
    df_train_cv = df_tune[cv_folds != cv_iteration, ]
    df_test_cv  = df_tune[cv_folds == cv_iteration, ]
  
    #grow shallow tree
    tree_cv = rpart::rpart(as.formula(paste("Surv(time, event) ~", paste(params_for_tree, collapse="+"))),
                               data = df_train_cv, minbucket = minbucket, maxdepth = maxdepth, cp = 0.001)
    #add the factor- identificator of a tree leaf, taking a predicted value 
    df_train_cv$cluster_tree = as.factor(round(predict(tree_cv, df_train_cv),6))
    df_test_cv$cluster_tree = as.factor(round(predict(tree_cv, df_test_cv),6))
    
    if(unique((round(predict(tree_cv, df_train_cv),6)))==1) {return (c(NaN, NaN))} 
    #add them to the Cox model with all the risk factors
    modified_cox_cv = coxph(as.formula(paste("Surv(time, event) ~",
                               paste(c(predict.factors, "cluster_tree"), collapse="+"))),
                               data = df_train_cv, x=TRUE)
    #replace NA to 0 for the coefficients 
    modified_cox_cv$coefficients[is.na(modified_cox_cv$coefficients)] = 0

    df_train_cv$modcox_lp =  predict(modified_cox_cv, newdata = df_train_cv, type = "lp", se.fit = FALSE, reference = "zero") 
    df_test_cv$modcox_lp =   predict(modified_cox_cv, newdata = df_test_cv, type = "lp", se.fit = FALSE, reference = "zero") 
      
    #check c-index 
    cindex_test[cv_iteration] = survConcordance(Surv(df_test_cv$time, df_test_cv$event)~df_test_cv$modcox_lp)$concordance
    cindex_train[cv_iteration] = survConcordance(Surv(df_train_cv$time, df_train_cv$event)~df_train_cv$modcox_lp)$concordance
  }
  return(c(mean(cindex_train), mean(cindex_test)))
}


############################## TUNE Ensemble 3B ######################
method_3B_tune = function(df_tune, params_for_tree, predict.factors, fixed_time=fixed_time,
                          n_cv = 3, nodedepth = nodedepth, cp =0.001, nodesize = nodesize, seed_to_fix = 100){
  #tuning ensemble 2b - fudging the tree in maxdepth, minbucket and risk factors to use to CV,
  # this is to get an average test c-index for one combination, all combinations are checked in method_3_cv() function
  set.seed(seed_to_fix)
  cv_folds = caret::createFolds(df_tune$event, k=n_cv, list = FALSE) #use caret to split into k-folds = cv_steps
  
  cindex_train = vector(length = n_cv); cindex_test = vector(length = n_cv)
  for (cv_iteration in 1:n_cv){
    df_train_cv = df_tune[cv_folds != cv_iteration, ]
    df_test_cv  = df_tune[cv_folds == cv_iteration, ]
    
    #grow shallow tree
    rf.dt = rfsrc(as.formula(paste("Surv(time, event) ~", paste(params_for_tree, collapse="+"))),
                  data = df_train_cv,
                  nodesize = minbucket,  # this is AVERAGE size, so we want this to be quite high
                  ntree = 1, #only need one tree
                  mtry =  length(params_for_tree), #try all predictors
                  nodedepth = nodedepth,  # calibrated parameter, from 2 till 7
                  nsplit = 20, splitrule = "logrank", statistics=TRUE,membership=TRUE,
                  importance = "none", #to speed up by switching off VIMP calculations
                  seed = seed_to_fix
    )
    #tree_to_plot = get.tree(rf.dt, 1); plot(tree_to_plot)
    
    predictedscore_rfdt_train =predict(rf.dt, newdata = df_train_cv, node= TRUE)
    predictedscore_rfdt_test =predict(rf.dt, newdata = df_test_cv, node= TRUE)
    predictedscore_rfdt_train$leaf.count #number of clusters
    clusters = unique(predictedscore_rfdt_train$predicted)
    
    df_train_cv$cluster_tree = as.factor(predictedscore_rfdt_train$predicted)
    df_test_cv$cluster_tree = as.factor(predictedscore_rfdt_test$predicted)
    
    #add them to the Cox model with all the risk factors
    modified_cox_cv = coxph(as.formula(paste("Surv(time, event) ~",
                                             paste(c(predict.factors, "cluster_tree"), collapse="+"))),
                            data = df_train_cv, x=TRUE)
    
    df_train_cv$modcox_lp =  predict(modified_cox_cv, newdata = df_train_cv, type = "lp", se.fit = FALSE, reference = "zero") 
    df_test_cv$modcox_lp =   predict(modified_cox_cv, newdata = df_test_cv, type = "lp", se.fit = FALSE, reference = "zero") 
    
    #check c-index 
    cindex_test[cv_iteration] = survConcordance(Surv(df_test_cv$time, df_test_cv$event)~df_test_cv$modcox_lp)$concordance
    cindex_train[cv_iteration] = survConcordance(Surv(df_train_cv$time, df_train_cv$event)~df_train_cv$modcox_lp)$concordance
  }
  return(c(mean(cindex_train), mean(cindex_test)))
}

############################## TUNE Ensemble 2B RF tree ######################

method_2B_tune = function(df_tune, params_for_tree , predict.factors, fixed_time=fixed_time,
                         n_cv = 3, nodedepth = nodedepth, cp =0.001, nodesize = nodesize, seed_to_fix = 100){
  #tuning ensemble 2b - fudging the tree in maxdepth, minbucket and risk factors to use to CV,
  # this is to get an average test c-index for one combination, all combinations are checked in method_3_cv() function
  set.seed(seed_to_fix)
  cv_folds = caret::createFolds(df_tune$event, k=n_cv, list = FALSE) #use caret to split into k-folds = cv_steps
  
  cindex_train = vector(length = n_cv); cindex_test = vector(length = n_cv)
  for (cv_iteration in 1:n_cv){
    df_train_cv = df_tune[cv_folds != cv_iteration, ]
    df_test_cv  = df_tune[cv_folds == cv_iteration, ]
    
    #grow shallow tree
    rf.dt = rfsrc(as.formula(paste("Surv(time, event) ~", paste(params_for_tree, collapse="+"))),
                  data = df_train_cv,
                  nodesize = nodesize,  # this is AVERAGE size, so we want this to be quite high
                  ntree = 1, #only need one tree
                  mtry =  length(params_for_tree), #try all predictors
                  nodedepth = nodedepth,  # calibrated parameter, from 2 till 7
                  nsplit = 50, splitrule = "logrank", statistics=TRUE,membership=TRUE,
                  importance = "none", #to speed up by switching off VIMP calculations
                  seed = seed_to_fix
    )
    #tree_to_plot = get.tree(rf.dt, 1); plot(tree_to_plot)

    predictedscore_rfdt_train =predict(rf.dt, newdata = df_train_cv, node= TRUE)
    predictedscore_rfdt_test =predict(rf.dt, newdata = df_test_cv, node= TRUE)
    predictedscore_rfdt_train$leaf.count #number of clusters
    clusters = unique(round(predictedscore_rfdt_train$predicted,6))
    df_train_cv$cluster_tree = as.factor(round(predictedscore_rfdt_train$predicted,6))
    df_test_cv$cluster_tree = as.factor(round(predictedscore_rfdt_test$predicted,6))
    
    #calibrate Cox models in each final leaf and#predict event probability at fixed_time
    df_train_cv$eventprob = NaN
    df_test_cv$eventprob = NaN 
    for (i in seq(length(clusters))){
      #calibrate model in cluster[i]
      if (sum(df_train_cv$cluster_tree== clusters[i]) < 20){ #<20 in the cluster 
         if (sum(df_train_cv$cluster_tree== clusters[i]) <=5){ #if only one, then nothing works
                    return(c(NaN, NaN))
              }else{  #if 2-20 then just use a KM for survival prediction, no further model
              cox_i = coxph(as.formula(paste("Surv(time, event) ~ 1")),data = df_train_cv[df_train_cv$cluster_tree== clusters[i],], x= TRUE)
              }
      }else{
        cox_i = coxph(as.formula(paste("Surv(time, event) ~", paste(c(predict.factors), collapse="+"))),
                    data = df_train_cv[df_train_cv$cluster_tree== clusters[i],], x= TRUE)
        #if all parameter is of the same value in a cluster, its coeff is NA and it breaks predictSurvProb, so we replace with 0
        cox_i$coefficients[is.na(cox_i$coefficients)] = 0
      }
      #predict event prob in train and test from cox_i
      df_test_cv[df_test_cv$cluster_tree == clusters[i] ,"eventprob"] = 1- 
          pec::predictSurvProb(cox_i, df_test_cv[df_test_cv$cluster_tree == clusters[i] ,] , fixed_time)
      df_train_cv[df_train_cv$cluster_tree == clusters[i] ,"eventprob"] = 1- 
        pec::predictSurvProb(cox_i, df_train_cv[df_train_cv$cluster_tree == clusters[i] ,] , fixed_time)
    }
    #check c-index 
    cindex_test[cv_iteration] = survConcordance(Surv(df_test_cv$time, df_test_cv$event)~df_test_cv$eventprob)$concordance
    cindex_train[cv_iteration] = survConcordance(Surv(df_train_cv$time, df_train_cv$event)~df_train_cv$eventprob)$concordance
  }
  return(c(mean(cindex_train), mean(cindex_test)))
}

############################## TUNE Ensemble 2A rpart ######################

method_2A_tune = function(df_tune, predictors.rpart , predict.factors,
                    fixed_time=10, n_cv = 3, maxdepth = 10, minbucket = 0, cp =0.001, seed_to_fix = 100){

  set.seed(seed_to_fix)
  cv_folds = caret::createFolds(df_tune$event, k=n_cv, list = FALSE) #use caret to split into k-folds = cv_steps
  cindex_train = vector(length = n_cv); cindex_test = vector(length = n_cv)
  if (minbucket==0) {minbucket = max(50,dim(df_train_cox_rpart)[1]/10)}
  
  for (cv_iteration in 1:n_cv){
    df_train_cv = df_tune[cv_folds != cv_iteration, ]
    df_test_cv  = df_tune[cv_folds == cv_iteration, ]
    
    rpart_params = rpart.control(minbucket =minbucket ,  maxdepth = maxdepth, cp = cp)
    rpart.m1 = rpart::rpart(as.formula(paste("Surv(time, event) ~", paste(predictors.rpart, collapse="+"))),
                            data = df_train_cv, control = rpart_params)
    clusters = unique(round(predict(rpart.m1,df_train_cv),6)); length(clusters) #number of clusters
    df_train_cv$cluster_tree = as.factor(round(predict(rpart.m1,df_train_cv),6))
    df_test_cv$cluster_tree = as.factor(round(predict(rpart.m1,df_test_cv),6))

    #calibrate Cox models in each final leaf and#predict event probability at fixed_time
    df_train_cv$eventprob = NaN
    df_test_cv$eventprob = NaN 
    for (i in seq(length(clusters))){
      #calibrate model in cluster[i]
      if (sum(df_train_cv$cluster_tree== clusters[i]) < 20){
        cox_i = coxph(as.formula(paste("Surv(time, event) ~ 1")),
                      data = df_train_cv[df_train_cv$cluster_tree== clusters[i],], x= TRUE)
      }else{
        cox_i = coxph(as.formula(paste("Surv(time, event) ~", paste(c(predict.factors), collapse="+"))),
                      data = df_train_cv[df_train_cv$cluster_tree== clusters[i],], x= TRUE)
        #if all parameter is of the same value in a cluster, its coeff is NA and it breaks predictSurvProb, so we replace with 0
        cox_i$coefficients[is.na(cox_i$coefficients)] = 0
      }
      #predict event prob in train and test from cox_i
      df_test_cv[df_test_cv$cluster_tree == clusters[i] ,"eventprob"] = 1- 
        pec::predictSurvProb(cox_i, df_test_cv[df_test_cv$cluster_tree == clusters[i] ,] , fixed_time)
      df_train_cv[df_train_cv$cluster_tree == clusters[i] ,"eventprob"] = 1- 
        pec::predictSurvProb(cox_i, df_train_cv[df_train_cv$cluster_tree == clusters[i] ,] , fixed_time)
    }
    #check c-index 
    cindex_test[cv_iteration] = survConcordance(Surv(df_test_cv$time, df_test_cv$event)~df_test_cv$eventprob)$concordance
    cindex_train[cv_iteration] = survConcordance(Surv(df_train_cv$time, df_train_cv$event)~df_train_cv$eventprob)$concordance
  }
  return(c(mean(cindex_train), mean(cindex_test)))
}


