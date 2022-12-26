
################### Ensemble methods functions ###################
# 1) training function on df_train
# 2) predicting event probability for fixed_time from a method, either for df_train or for df_test
# 3) validating performance from the method on df_test

#' Methods:
#' Cox 
#' Cox MFP
#' SRF, cross validate to tune 
#' Ensemble 1a = Cox into SRF, cross-validate SRF to tune 
#' Ensemble 1b = SRF into Cox, cross-validate SRF to tune 
#' Ensemble 2a = VIMP from SRF, then single RPART TREE, then Cox in each leaf
#' Ensemble 2b = same with SRF tree 
#' Ensemble 3 = VIMP from SRF, then single RPART TREE, then modified Cox with leaves identificators

################## libraries ##################
library(psych)
library(glmnet)
library(survival)

library(survminer)
library(ggplot2)
library(randomForest)
library(timeROC) # for time-dep AUC
library(zoo) # for timeROC
library(MALDIquant) # for match.closest
library(randomForestSRC)
library(data.tree)
library(rpart)
library(treeClust)
library(rpart.plot) # to plot rpart trees
library(classifierplots)
library(mfp) #for fractional polynomials
library(pec) # for predictSurvProb
library(dplyr)

################# validation statistics #####################

bs_surv = function(y_predicted_newdata,
                   #these are probability of events! not survival prob!
                   df_brier_train,
                   df_newdata,
                   time_points,
                   weighted = TRUE){
  #'This function calculates time-dependent BrierScore for df_newdata,
  #' overall same as
  #' https://scikit-survival.readthedocs.io/en/stable/api/generated/sksurv.metrics.brier_score.html#sksurv.metrics.brier_score
  #' it uses IPCW (inverse probability of censoring weights), computed with K-M curve
  #' where events are censored events from train data. i.e. df_brier_train
  #' returns vector of BS for each time in time_points
  
  # K-M prob of censoring for each observation till its individual time
  df_newdata$p_km_obs = survival_prob_km(df_brier_train, df_newdata$time, estimate_censoring  = TRUE)
  df_newdata$p_km_obs = pmax(pmin(df_newdata$p_km_obs, 0.9999),0.0001)
  # !! impute with mean observations if can't estimate !!
  #!!!
  df_newdata[is.na(df_newdata$p_km_obs), "p_km_obs"] = mean(df_newdata$p_km_obs, na.rm=1)
  #one number P(t_i is censored) = P(T > t_i) = S(t_i):
  p_km_t = survival_prob_km(df_brier_train, time_points, estimate_censoring  = TRUE)
  p_km_t= pmax(pmin(p_km_t, 0.9999), 0.0001); p_km_t
  
  bs = c()
  for (j in 1:length(time_points)) {
    # assign t_point and predicted  probabilities
    if (length(time_points) ==1) { #only 1 time
      t_point = time_points
      ppp = pmax(pmin(y_predicted_newdata, 0.99999), 0.00001)
    } else { #many times
      t_point = time_points[j]
      ppp = pmax(pmin(y_predicted_newdata[,j], 0.99999), 0.00001)
    }
    #cases and controls by time t_point
    id_case = ((df_newdata$time <= t_point) & (df_newdata$event==1))
    id_control=  (df_newdata$time > t_point)
    
    # compute BS with weights which are 1/G(t) for controls and 1/G(obs_i) for cases
    # if weights == false, use w=1 for all
    if (weighted==TRUE){
      #brier score is average of weighted squared errors
      bs[j] = (sum(as.numeric(id_case)*(1-ppp)^2*1/df_newdata$p_km_obs, na.rm=1)+
                 sum(id_control*(0-ppp)^2*1/p_km_t[j], na.rm=1))/dim(df_newdata)[1]
    }else{
      #unweighted BS
      bs[j] = sum(id_case*(1-ppp)^2 + id_control*(0-ppp)^2, na.rm=1)/dim(df_newdata)[1]
    }
  }
  names(bs) = round(time_points,6)
  return(bs)
}

survival_prob_km = function(df_km_train, times, estimate_censoring = FALSE){
  #This function calculates survival probability by K-M estimate returns vector (times)
  # if estimate_censoring is TRUE, censoring distr is estimated
  # we also extrapolate it in survival function space using poly(n=3)
  if (estimate_censoring == FALSE) {
    km = survfit(Surv(time, event)~1, data = df_km_train)
    kmf <- approxfun(km$time, km$surv,  method = "constant")
    kmdf = data.frame(cbind("time"=km$time, "surv" = km$surv))
    extrap = lm(surv ~ poly(time,3, raw = TRUE), data = kmdf)
    km_extrap = function(x){(cbind(1, x, x**2, x**3) %*% extrap$coefficients[1:4])}
    
  }else{
    df_km_train$censor_as_event = 1 - df_km_train$event
    km = survfit(Surv(time, censor_as_event)~1, data = df_km_train)
    kmdf = data.frame(cbind("time"=km$time, "surv" = km$surv))
    kmf <- approxfun(km$time, km$surv, method = "constant")
    extrap = lm(surv ~ poly(time,3, raw = TRUE), data = kmdf)
    km_extrap = function(x){(cbind(1, x, x**2, x**3) %*% extrap$coefficients[1:4])}
  }
  return(km_extrap(times))
}

method_any_validate = function(y_predict, times_to_predict, df_train, 
                               df_test, weighted = TRUE, alpha = "mean"){
  
  #This function computes auc, brier score, c-index,
  # calibration slope and alpha for df_test
  #for apparent statistics use test  = train
  
  auc_score = c()
  brier_score = c()
  brier_score_scaled = c()
  c_score = c()
  calibration_slope = c()
  calibration_alpha = c()
  
  for (i in 1:length(times_to_predict)){
    t_i = times_to_predict[i]  
    if (length(times_to_predict)>1) { y_hat = y_predict[,i] } else {y_hat = unlist(y_predict)}
    #' time dependent AUC
    if(class(try(timeROC::timeROC(T = df_test$time,delta=df_test$event,
                                  marker= y_hat, times = t_i*0.9999999,cause=1,
                                  weighting = "marginal"), silent = TRUE))=="try-error"){
      auc_score[i] = NaN}else{ 
        auc_score[i] = timeROC::timeROC(T = df_test$time,delta=df_test$event,
                                        marker= y_hat, times = t_i*0.9999999,cause=1,
                                        weighting = "marginal")$AUC[2]}
    #' compute time-dependent Brier score:
    if (class(try(bs_surv(y_hat, df_train, df_test, t_i, weighted = weighted), silent=TRUE))=="try-error"){
      brier_score[i] = NaN}else{
        brier_score[i]= bs_surv(y_hat, df_train, df_test, t_i, weighted = weighted)}
    
    #' compute concordance - time-dependent in a sense that a predictor is event probability at t_i:
    #' for Cox model it is the same for each time, as event prob
    #' will be in the same order as LPs at any time point
    if (class(try(concordancefit(Surv(df_test$time, df_test$event), -1*y_hat), silent=TRUE))=="try-error"){
      c_score[i]= NaN }else{
        if (is.null(concordancefit(Surv(df_test$time, df_test$event), -1*y_hat)$concordance)) {c_score[i]= NaN
        }else{c_score[i]=concordancefit(Surv(df_test$time, df_test$event), -1*y_hat)$concordance}
      }
    
    # can add later confusion matrix, but also need to find Youden point instead of 0.5
    # confmatrix = timeROC::SeSpPPVNPV(0.5, T = df_test$time,delta=df_test$event,
    #              marker= y_hat, times = t_i*0.999,cause=1, weighting = "marginal")
    # c(as.double(confmatrix$TP[2]), as.double(confmatrix$FP[2]),
    # as.double(confmatrix$PPV[2]), as.double(confmatrix$NPV[2]))
    
    # compute scaled Brier score, and calibration slope and alpha:
    # 1/0 by t_i:
    df_test$event_ti = ifelse(df_test$time <= t_i & df_test$event ==1, 1, 0)
    
    # cut 0 and 1 predicted probabilities for the logit to work:
    df_test$predict_ti = pmax(pmin(y_hat, 0.99999), 0.00001)
    
    #scaled BS = 1 - BS / (BS @ all predictions == prevalence by t_i )
    brier_score_base_i = bs_surv(rep(mean(df_test$event_ti),dim(df_test)[1]),
                                 df_train, df_test, t_i, weighted = weighted)
    brier_score_scaled[i]= 1 - (brier_score[i]/brier_score_base_i)
    
    #Calibration slope and alpha. Take out censored observations before t_i,
    #ie leaving those which state we know:
    df_test_in_scope = df_test[(df_test$time >= t_i) | (df_test$time <t_i & df_test$event ==1), ]
    
    y_hat_hat = log(df_test_in_scope$predict_ti / (1-df_test_in_scope$predict_ti))
    y_actual_i = df_test_in_scope$event_ti
    
    if(class(
      try( glm(y_actual_i ~ y_hat_hat,family = binomial(link = "logit")), silent = TRUE)
    )[1] =="try-error" ){
      calibration_slope[i] = NaN
      calibration_alpha[i] = NaN
    }else{
      calibration_slope[i] = glm(y_actual_i ~ y_hat_hat,
                                 family = binomial(link = "logit"))$coefficients[2]
      if (alpha == "logit"){ #take alpha from alpha: logit(y)~ logit(y_hat) + alpha
        calibration_alpha[i] = glm(y_actual_i ~ offset(y_hat_hat),
                                   family = binomial(link = "logit"))$coefficients[1]
      }else{  #take alpha as alpha= mean(y) - mean(y_hat)
        calibration_alpha[i] =  mean(y_actual_i) - mean(df_test_in_scope$predict_ti)
      }
    }#end else
  } #end for
  
  output = data.frame("T" = times_to_predict,
                      "AUCROC" = auc_score,
                      "BS" = brier_score,
                      "BS_scaled" = brier_score_scaled,
                      "C_score" = c_score,
                      "Calib_slope" = calibration_slope,
                      "Calib_alpha" = calibration_alpha
  )
  return (output)
}



populationstats = function(df_stats, time_f,  namedf = "df"){
  #df_stats = results_object[[3]]
  statsdf = data.frame()
  df_stats$event_fixed_time = ifelse(((df_stats$event ==1)&(df_stats$time < time_f)), 1, 0)
  df_stats$time_fixed_time = ifelse(df_stats$time < time_f, df_stats$time, time_f )
  statsdf["data_set_name","value"] = namedf
  statsdf["sample_size_N","value"] = nrow(df_stats)
  statsdf["time_fixed","value"] =   time_f  # results_object[[8]][1]
  statsdf["time_max","value"] = round(max(df_stats$time),2)
  statsdf["time_mean","value"] = round(mean(df_stats$time),2)
  statsdf["time_mean_event","value"] = round(sum(df_stats[df_stats$event ==1, "time"])/sum(df_stats$event ==1),2)
  statsdf["time_mean_censored","value"] = round(sum(df_stats[df_stats$event ==0, "time"])/sum(df_stats$event ==0),2)
  statsdf["time_mean_T","value"] = round(mean(df_stats$time_fixed_time),2)
  statsdf["time_mean_event_T","value"] = round(sum(df_stats[df_stats$event_fixed_time ==1, "time"])/sum(df_stats$event_fixed_time ==1),2)
  statsdf["time_mean_censored_T","value"] = round(sum(df_stats[df_stats$event_fixed_time ==0, "time"])/sum(df_stats$event_fixed_time ==0),2)
  statsdf["event_%","value"] = round(sum(df_stats$event==1) / nrow(df_stats) *100,2)
  statsdf["event_before_T_%","value"] = round(sum(df_stats$event_fixed_time==1) / nrow(df_stats)*100,2)
  statsdf["censored_before_T_%","value"] = 
    round(sum(df_stats$event_fixed_time==0 & df_stats$time_fixed_time < statsdf["time_fixed","value"]) / nrow(df_stats)*100,2)
  return (statsdf)
}

############# Basic Cox Model functions (in the same format as other methods) ###############

method_cox_train = function(df_train, predict.factors){
  # wrapper for coxph() function returning a trained Cox model
  cox.m  = coxph(as.formula(paste("Surv(df_train$time, df_train$event) ~",
                                  paste(predict.factors, collapse="+"))), 
                 data =df_train, x = TRUE)
  #!!!!! We replace NA coefficients with 0 
  #!!!!! i.e. ignore predictors which the Cox model couldn't estimate)
  #!!!!! This way the model can still produce some predictions 
  cox.m$coefficients[is.na(cox.m$coefficients)] = 0
  return (cox.m)
}

method_cox_predict = function(model_cox, newdata, times){
    #' returns event probability from trained cox model model_cox
    #' for newdata at given times
    #' could use pec package, but it makes data table out of newdata
    #' predicted_event_prob = 1-pec::predictSurvProb(model_cox, newdata, times)
    
    #define bh - baseline hazard as dataframe with "time" and "hazard"
    # if baseline hazard can't be calibrated, # return mean(y) for all times
    # we take baseline hazard from K-M estimate and lp from Cox !!!! :((
    if(class(try(basehaz(model_cox, centered = 0), silent= TRUE)) =="try-error"){
      bh = summary(survfit(model_cox$y~1),times)$cumhaz  
      predicted_event_prob = matrix(nrow = dim(newdata)[1],ncol = length(times))
      for (i in seq(length(times))){
        predicted_event_prob[,i] = 1 -
          exp(-bh[i]*exp(predict(model_cox, newdata = newdata, type = "lp", reference = "zero")))  }
      colnames(predicted_event_prob) = round(times,6)
      return(predicted_event_prob)
    } else {bh = basehaz(model_cox, centered = 0)}
    
    #define bh as function to compute bh for any time
    bh_approx = approxfun(bh[,"time"], bh[,"hazard"], method = "constant")
    
    #define bh_extrap how to extrapolate outside of the times in the training data
    if (class(try( lm(hazard ~ poly(time,3, raw = TRUE), 
                      data = bh), silent = TRUE))!="try-error"){
      extrap = lm(hazard ~ poly(time,3, raw = TRUE), data = bh)  
      bh_extrap = function(x){sum(c(1, x, x**2, x**3) * extrap$coefficients[1:4])}
    } else {
      min_bh = min(bh[,"hazard"], na.rm=1) 
      max_bh = max(bh[,"hazard"], na.rm=1)
      l =dim(bh)[1]
      bh[1, c("hazard", "time")] = c(0.0000001, min_bh);
      bh[l+1, c("hazard", "time")] = c(bh[l,"time"]+100000, max_bh)
      bh_extrap = approxfun(bh[,"time"], bh[,"hazard"], method = "constant")
    }
    
    #compute event probability for times
    # create placeholder
    predicted_event_prob = matrix(nrow = dim(newdata)[1], ncol = length(times))
    #go over each time in times
    for (i in seq(length(times))){
      
      if (is.na(bh_approx(times[i]))){  #if interpolation doesn't work, take extrapolated value
        bh_time = bh_extrap(times[i]); if(is.na(bh_time)){bh_time = mean(bh[, "hazard"], na.rm=TRUE)}
      }else{  
        bh_time = bh_approx(times[i])  
      }
      
      #if baseline hazard is infinite, event probability is 1
      if (bh_time == Inf){
        predicted_event_prob[,i] = 1
        
        #if baseline hazard is ==0, event prob is 0 for all with survival==1 
        # (somehow "survival" calculates even with bh==0)
      }else if (bh_time == 0){ 
        predicted_event_prob[,i] = 0
        
        #for those with non-definite survival, use K-M estimator 
        #(as Cox doesn't work for conflicting bh=0 and lp = exp)
        try({ if (sum(predict(model_cox, newdata = newdata,
                              type = "survival", reference = "zero")<1)>0) {
          #define K-M estimator with extrapolation beyond times in training data
          km = survfit(model_cox$y~1)
          km = as.data.frame(cbind("time"= km$time,"surv"= km$surv))
          min_km = min(km$surv, na.rm=1) 
          max_km = max(km$surv, na.rm=1)
          l = dim(km)[1]
          km[1, c("time", "surv")] = c(0.0000001, max_km);
          km[l+1, c("time", "surv")] = c(km[l,"time"]+100000, min_km)
          km_fun = approxfun(km$time, km$surv, method = "constant")
          #replace 0 with K-M probabilities for those with survival <1
          predicted_event_prob[predict(model_cox, newdata = newdata,
                                       type = "survival", reference = "zero")<1,i] = km_fun(times[i])}
        },silent = TRUE)
        #if baseline hazard is a number, use the survival formula   
      }else{
        predicted_event_prob[,i] = 1 -exp(-bh_time*exp(predict(model_cox, 
                                                               newdata = newdata,type = "lp", reference = "zero")))
      }
    }
    # name columns by the time for which it predicts event prob
    colnames(predicted_event_prob) = round(times,6)
    return (predicted_event_prob)
}

method_cox_cv = function(df, predict.factors, fixed_time =10, cv_number = 5, seed_to_fix = 100){
  time_0 = Sys.time()
  set.seed(seed_to_fix)
  cv_folds = caret::createFolds(df$event, k=cv_number, list = FALSE) #use caret to split into k-folds = cv_steps
  modelstats_train = list(); modelstats_test = list()
  for (cv_iteration in 1:cv_number){
    df_train_cv = df[cv_folds != cv_iteration, ]
    df_test_cv  = df[cv_folds == cv_iteration, ]
    cox.model = method_cox_train(df_train_cv, predict.factors)
    y_predict_test = method_cox_predict(cox.model, df_test_cv, fixed_time)
    y_predict_train = method_cox_predict(cox.model, df_train_cv, fixed_time)
    modelstats_test[[cv_iteration]] = method_any_validate(y_predict_test,
              fixed_time, df_train_cv, df_test_cv, weighted = 1)
    modelstats_train[[cv_iteration]] = method_any_validate(y_predict_train, 
              fixed_time, df_train_cv, df_train_cv, weighted = 1)
  }
  df_modelstats_test = data.frame(modelstats_test[[1]])
  df_modelstats_train = data.frame(modelstats_train[[1]])
  for (i in 2:cv_number){df_modelstats_test[i,]= modelstats_test[[i]]; df_modelstats_train[i,]= modelstats_train[[i]]}
  df_modelstats_test$test = 1; df_modelstats_train$test = 0
  output = list()
  output$test = df_modelstats_test
  output$train = df_modelstats_train
  output$testaverage = sapply(df_modelstats_test,mean,na.rm=1)
  output$trainaverage = sapply(df_modelstats_train,mean,na.rm=1)
  time_1 = Sys.time()
  output$time = time_1 - time_0
  return(output)
}

############# Augmented CoxPH with Fractional Polynomials ###############

method_coxmfp_train = function(df_train, predict.factors, verbose = FALSE){
  #Cox with flarctional polynomials, returns final Cox model with the selected fp() risk factors
  if (length(predict.factors)>20) {
    print ("Too many factors (>20) to perform MFP, default to baseline Cox")
    cox.mfp  = coxph(as.formula(paste("Surv(df_train$time, df_train$event) ~",
                                      paste(predict.factors, collapse="+"))), data =df_train, x = TRUE)
  }else{
    # create predict.factors.mfp = paste("fp(", predict.factors, ")", sep="") for all to be fp'd
    predict.factors.mfp = predict.factors  
    for (i in 1:length(predict.factors)) {
      # only "fp" continuous factors, or those with more than 5 values i=1
      if (length(table(select(df_train, predict.factors[i]))) >= 5){   
        predict.factors.mfp[i]= paste("fp(", predict.factors[i], ")", sep="")}
    }
    #calculate mfp 
    mfp.compute <- mfp(as.formula(paste("Surv(time, event) ~",paste(predict.factors.mfp, collapse="+"))),
                       data = df_train, family = "cox", verbose = FALSE, maxits = 15,select = 0.1)
    if (verbose){  print ("Multiple polynomial formula for Cox model _"); print (mfp.compute$formula)}
    cox.mfp  = coxph(formula = mfp.compute$formula, data = df_train, x = TRUE)
  }
  return(cox.mfp)
}

method_coxmfp_predict = function(model_coxmfp, df_test, fixed_time){
  #predicting event probability from CoxPh for all observations in df_test by fixed_time
  predicted_event_prob = 1-pec::predictSurvProb(model_coxmfp, df_test, fixed_time)
  return (predicted_event_prob)
}

method_coxmfp_cv = function(df, predict.factors, fixed_time =10, cv_number = 5,seed_to_fix = 100){
  # cross-validating CoxPh model, returns performance measures for each CV and averaged metrics
  time_0 = Sys.time()
  set.seed(seed_to_fix)
  cv_folds = caret::createFolds(df$event, k=cv_number, list = FALSE) #use caret to split into k-folds = cv_steps
  modelstats_train = list(); modelstats_test = list()
  models_for_each_cv = list()
  
  for (cv_iteration in 1:cv_number){
    df_train_cv = df[cv_folds != cv_iteration, ]
    df_test_cv  = df[cv_folds == cv_iteration, ]
    cox.model = method_coxmfp_train(df_train_cv, predict.factors)
    
    models_for_each_cv[[cv_iteration]] = cox.model
    
    y_predict_test = method_coxmfp_predict(cox.model, df_test_cv, fixed_time)
    y_predict_train = method_coxmfp_predict(cox.model, df_train_cv, fixed_time)
    modelstats_test[[cv_iteration]] = method_any_validate(y_predict_test, fixed_time, df_train_cv, df_test_cv, weighted = 1)
    modelstats_train[[cv_iteration]] = method_any_validate(y_predict_train, fixed_time, df_train_cv, df_train_cv, weighted = 1)
  }
  df_modelstats_test = data.frame(modelstats_test[[1]])
  df_modelstats_train = data.frame(modelstats_train[[1]])
  for (i in 2:cv_number){df_modelstats_test[i,]= modelstats_test[[i]]; df_modelstats_train[i,]= modelstats_train[[i]]}
  df_modelstats_test$test = 1; df_modelstats_train$test = 0
  
  output = list()
  output$test = df_modelstats_test
  output$train = df_modelstats_train
  output$testaverage = sapply(df_modelstats_test,mean,na.rm=1)
  output$trainaverage = sapply(df_modelstats_train,mean,na.rm=1)
  output$tuned_cv_models = models_for_each_cv
  
  time_1 = Sys.time()
  output$time = time_1 - time_0
  return(output)
}

############# SRF ################

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

srf_tune = function(df_tune,  cv_number =3, 
                    predict.factors, fixed_time = NaN, 
                    seed_to_fix = 100,mtry= c(3,4,5), 
                    nodesize = c(10,20,50),nodedepth = c(100),
                    verbose = FALSE, oob = TRUE){
  #function to tune survival random forest by mtry, nodesize and nodedepth grid 
  # if oob = TRUE, there is no CV !!! as OOB does the job already 
  
  #set seed
  set.seed(seed_to_fix)
  
  # limit mtry with the number of predictors and nodesize by 1/6 of sample size
  if (sum(nodesize > dim(df_tune)[1]/6) > 0) {
    if(verbose) print ("Warning - some min nodesize is > 1/6 of the sample size (1/2 of CV fold)")}
  
  nodesize = nodesize[nodesize <= dim(df_tune)[1]/6]
  
  #if all numbers higher than number of factors, only check this factor
  if (sum(mtry > length(predict.factors)) == length(mtry)) { 
    mtry = c(length(predict.factors))}
  mtry = mtry[mtry <= length(predict.factors)]
  
  #grid of values to tune
  grid_of_values = expand.grid("mtry" = mtry, 
                               "nodesize" = nodesize, "nodedepth" = nodedepth)
  if(verbose) print(paste("Grid size is", dim(grid_of_values)[1]))
  
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
      if(verbose) print(grid_of_values[i,])  
      mtry_i = grid_of_values[i,"mtry"]
      nodesize_i = grid_of_values[i,"nodesize"]
      nodedepth_i = grid_of_values[i,"nodedepth"]
      
      #train grid combination for each cv_iteration
      modelstats_cv = list()
      for (cv_iteration in 1:cv_number) {
        print (paste("SRF tuning CV step", cv_iteration, "/out of", cv_number))
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
                                           "time" = validation_stats$T,
                                           "AUCROC" = validation_stats$AUCROC, 
                                           "BS"= validation_stats$BS,
                                           "BS_scaled" = validation_stats$BS_scaled, 
                                           "C_score" = validation_stats$C_score, 
                                           "Calib_alpha" = validation_stats$Calib_alpha,
                                           "Calib_slope" = validation_stats$Calib_slope)
      }#end k-fold CV for one grid combination
      
      #averaging over cv-steps, firs transform to data.frame to use mean()
      modelstats_cv_df = data.frame(t(modelstats_cv[[1]]))
      for (j in 2:cv_number) {modelstats_cv_df = rbind(modelstats_cv_df,t(modelstats_cv[[j]]))}
      modelstats[[i]] = c(modelstats_cv[[1]]["mtry"],  modelstats_cv[[1]]["nodesize"], 
                          modelstats_cv[[1]]["nodedepth"], 
                          "AUCROC" = mean(modelstats_cv_df$AUCROC, na.rm=1), 
                          "BS" = mean(modelstats_cv_df$BS, na.rm=1),
                          "BS_scaled" = mean(modelstats_cv_df$BS_scaled, na.rm=1),
                          "C_score"= mean(modelstats_cv_df$C_score, na.rm=1),
                          "Calib_alpha" =  mean(modelstats_cv_df$Calib_alpha, na.rm=1),
                          "Calib_slope" =  mean(modelstats_cv_df$Calib_slope, na.rm=1),
                          "time"= fixed_time)
    }#end for grid
    
  } else { # end if(oob==false) 
    if(verbose) {print ("No internal CV for training SRF, 
           instead out-of-bag predictions used to assess performance")}
    for (i in 1:dim(grid_of_values)[1]){
      if(verbose) print(grid_of_values[i,])  
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
      modelstats[[i]] =  c("mtry" = mtry_i,  "nodesize" = nodesize_i, 
                           "nodedepth" = nodedepth_i, 
                           "time" = validation_stats$T,
                           "AUCROC" = validation_stats$AUCROC, 
                           "BS"= validation_stats$BS,
                           "BS_scaled" = validation_stats$BS_scaled, 
                           "C_score" = validation_stats$C_score, 
                           "Calib_alpha" = validation_stats$Calib_alpha, 
                           "Calib_slope" = validation_stats$Calib_slope)
      
    } #end for (i in grid)
  }#end else 
  
  #reshaping into data frame 
  df_modelstats = data.frame("V1" = modelstats[[1]])
  #check if there was more than 1 grid search
  if (dim(grid_of_values)[1] >1) {for (i in 2:dim(grid_of_values)[1]){ df_modelstats[i]= modelstats[[i]]}}
  df_modelstats = data.frame(t(df_modelstats))
  
  if (verbose == TRUE) {
    print(paste("AUC varies from", round(min(df_modelstats$AUCROC ),4), "to", round(max(df_modelstats$AUCROC ),4)))
    print(paste("Brier score varies from", round(min(df_modelstats$BS ),4), "to", round(max(df_modelstats$BS ),4)))
    print("Combination with highest AUC")
    print(df_modelstats[which.max(df_modelstats$AUCROC),  c("mtry", "nodesize", "nodedepth")])
    print("Combination with lowest Brier Score")
    print(df_modelstats[which.min(df_modelstats$BS),  c("mtry", "nodesize", "nodedepth")])
    print("Combination with lowest AUC")
    df_modelstats[which.min(df_modelstats$AUCROC), c("mtry", "nodesize", "nodedepth")]
  }
  output = list() 
  output$modelstats = df_modelstats 
  output$bestbrier = df_modelstats[which.min(df_modelstats$BS), ]
  output$bestauc = df_modelstats[which.max(df_modelstats$AUCROC), ]
  output$bestcindex = df_modelstats[which.max(df_modelstats$C_score), ]
  return(output)
}


method_srf_train = function(df_train, predict.factors, 
                            fixed_time=4, cv_number = 3, 
                            seed_to_fix = 100, fast_version = TRUE, oob = TRUE, verbose = FALSE){
  #for now only for best AUC but can be amended for brier score or cindex
  #defining the tuning grid for SRF 
  p = length(predict.factors) #number of predictors
  n = dim(df_train)[1] 
  #mtry grid
  mtry_default = round(sqrt(p),0)
  
  if (p<=10) {mtry = c(2,3,4,5)}else{if(p<=25){mtry = c(3,5,7,10,15)}else{
    mtry = c(round(p/10,0),round(p/5,0), round(p/3,0), round(p/2,0),mtry_default)}}
  
  #minimum nodesize grid
  nodesize = seq(min(15, round(n/6-1,0)), max(min(n/10,50),30), 5)
  #nodedepth grid
  nodedepth = c(50) #we don't tune this so just a big number 
  if (verbose) {print (paste("mtry", mtry, "nodedepth", nodedepth, "nodesize", nodesize))}
  
  if (fast_version == TRUE) {
    #take recommended mtry and check the best depth and node size 
    
    tune1 = srf_tune(df_train, cv_number = cv_number, predict.factors, 
                     fixed_time = fixed_time, seed_to_fix = seed_to_fix, 
                     mtry = mtry_default,nodesize = nodesize, 
                     nodedepth = nodedepth, oob = oob)
    nodesize_best = as.integer(tune1$bestauc["nodesize"])
    nodedepth_best = as.integer(tune1$bestauc["nodedepth"])
    #using the depth and size check the best mtry  
    tune2 = srf_tune(df_train, cv_number = cv_number,predict.factors, 
                     fixed_time = fixed_time, seed_to_fix = seed_to_fix, 
                     mtry = mtry, nodesize = nodesize_best,  
                     nodedepth = nodedepth_best , oob = oob)
    mtry_best = tune2$bestauc["mtry"] 
    best_combo_stat = tune2$bestauc
    modelstatsall = rbind(tune1$modelstats, tune2$modelstats)
  }else{
    tuneall = srf_tune(df_train, cv_number = cv_number,predict.factors,
                       fixed_time = fixed_time, seed_to_fix = seed_to_fix, 
                       mtry = mtry,nodesize = nodesize, 
                       nodedepth = nodedepth, oob = oob)
    nodesize_best = tuneall$bestauc["nodesize"]
    nodedepth_best = tuneall$bestauc["nodedepth"]
    mtry_best = tuneall$bestauc["mtry"]
    best_combo_stat = tuneall$bestauc
    modelstatsall= tuneall$modelstats
  }
  
  final.rfs =  rfsrc(as.formula(paste("Surv(time, event) ~", 
                                      paste(predict.factors, collapse="+"))),
                     data = df_train,
                     nodesize = nodesize_best,  # this is AVERAGE size, so we want this to be quite high
                     ntree = 500,
                     mtry = as.integer(mtry_best), 
                     nodedepth = as.integer(nodedepth_best),  
                     nsplit = 50, 
                     splitrule = "logrank", statistics= FALSE, membership=TRUE,
                     importance = "none", #to speed up by switching off VIMP calculations
                     seed = seed_to_fix
  )
  output = list()
  output$beststats = best_combo_stat
  output$allstats = modelstatsall
  output$model = final.rfs
  #calibrate SRF with the best parameters
  return (output)
}


method_srf_predict = function(model_srf, df_test, fixed_time =8, seed_to_fix = 100, oob= FALSE){
  
  predicted_event_prob = 1-srf_survival_prob_for_time(model_srf, df_test, fixed_time, oob= oob)
  if (length(fixed_time)==1) {return(predicted_event_prob)
  } else {
    predicted_event_prob = data.frame(predicted_event_prob)
    names(predicted_event_prob) = round(fixed_time,3)
    return (predicted_event_prob)
  }
}

method_srf_cv = function(df, predict.factors, fixed_time = 10, 
                         cv_number = 3, 
                         internal_cv_k = 3, 
                         seed_to_fix = 100){
  
  time_0 = Sys.time()
  set.seed(seed_to_fix)
  #use caret to split into k-folds = cv_steps
  cv_folds = caret::createFolds(df$event, k=cv_number, list = FALSE) 
  modelstats_train = list(); modelstats_test = list()
  srf_models_for_each_cv = list() #saving trained best SRF to re-use in ensemble 1A
  
  print (paste("Cross-validating Survival Random Forest with", cv_number,
               "outer loops, and ",internal_cv_k,"inner loops for model tuning"))
  
  for (cv_iteration in 1:cv_number){
    print (paste('External loop CV, step number = ', cv_iteration, '/ out of', cv_number))
    
    df_train_cv = df[cv_folds != cv_iteration, ]
    df_test_cv  = df[cv_folds == cv_iteration, ]
    
    srf.model.tuned = method_srf_train(df_train_cv, predict.factors,fixed_time = fixed_time,
                                       cv_number =internal_cv_k, seed_to_fix=seed_to_fix,
                                       fast_version = TRUE, oob = TRUE)
    
    y_predict_test =  method_srf_predict(srf.model.tuned$model, df_test_cv,
                                         fixed_time, seed_to_fix = seed_to_fix, oob = FALSE)
    y_predict_train = method_srf_predict(srf.model.tuned$model, df_train_cv, 
                                         fixed_time, seed_to_fix = seed_to_fix, oob= FALSE)
    
    modelstats_test[[cv_iteration]] = method_any_validate(y_predict_test, 
                                                          fixed_time, df_train_cv, df_test_cv, weighted = 1)
    modelstats_train[[cv_iteration]] = method_any_validate(y_predict_train, 
                                                           fixed_time, df_train_cv, df_train_cv, weighted = 1)
    
    srf_models_for_each_cv[[cv_iteration]]= srf.model.tuned$model
  }
  df_modelstats_test = data.frame(modelstats_test[[1]])
  df_modelstats_train = data.frame(modelstats_train[[1]])
  
  for (i in 2:cv_number){df_modelstats_test[i,]= modelstats_test[[i]]; df_modelstats_train[i,]= modelstats_train[[i]]}
  df_modelstats_test$test = 1; df_modelstats_train$test = 0
  output = list()
  output$test = df_modelstats_test
  output$train = df_modelstats_train
  output$testaverage = sapply(df_modelstats_test,mean,na.rm=1)
  output$trainaverage = sapply(df_modelstats_train,mean,na.rm=1)
  output$pretrained_srf_models = srf_models_for_each_cv
  time_1 = Sys.time()
  print (time_1 - time_0)
  output$time = time_1 - time_0
  return(output)
}

######################### Ensemble 1A ########################
# Method uses Cox model predictions to pass to Survival Random Forest 

method_1A_train = function(df_train, predict.factors, fixed_time=10, cv_number = 3, 
                           seed_to_fix = 100, fast_version = TRUE, oob = TRUE){
  #the function trains Cox model, then adds its predictions into Survival Random Forest model
  # to mimic stacking procedure and reduce overfitting,
  # we train Cox model on 0.9 of the data and predict on the rest 0.1 for each 1/10s fold 
  # so we pass out-of-the-bag prediction to SRF 
  
  #creating folds
  cv_folds = caret::createFolds(df_train$event, k=10, list = FALSE) 
  cindex_train = vector(length = 10); cindex_test = vector(length = 10)
  for (cv_iteration in 1:10){
    
    cox_train = df_train[cv_folds != cv_iteration, ]
    cox_oob  =  df_train[cv_folds == cv_iteration, ]
    #train cox model on cox_train
    cox_m_cv = method_cox_train(cox_train, predict.factors)
    #predict for cox_oob
    cox_predict_oob = method_cox_predict(cox_m_cv, cox_oob, fixed_time)
    #adding Cox prediction to the df_train in the column "cox_predict"
    df_train[cv_folds == cv_iteration, "cox_predict"] = cox_predict_oob
  }
  
  ##alternatively - just use all the data and pass apparent predictions to SRF 
  cox_model_for1a = method_cox_train(df_train, predict.factors)
  ##df_train$cox_predict = method_cox_predict(cox_model_for1a, df_train, fixed_time )
  
  #adding new factor and tuning SRF model with this added factor using srf_train 
  predict.factors.1A = c(predict.factors, "cox_predict")
  srf_model_for1a = method_srf_train(df_train, predict.factors = predict.factors.1A, 
                                     fixed_time=fixed_time, cv_number = cv_number,
                                     seed_to_fix = seed_to_fix, fast_version = fast_version, 
                                     oob = oob)
  
  v = vimp(srf_model_for1a$model, importance = "permute", seed = seed_to_fix)
  var_importance = sort(v$importance, decreasing = TRUE)
  
  srf_model_for1a$vimp10  = var_importance[1:min(length(var_importance),15)]
  srf_model_for1a$model_base  = cox_model_for1a
  
  return(srf_model_for1a)
}


method_1A_predict = function(model_1a, df_test, fixed_time, seed_to_fix = 100, oob= FALSE){
  #use model_base with the base Cox model to find cox_predict 
  df_test$cox_predict = method_cox_predict(model_1a$model_base, df_test, fixed_time)
  # now use "model" which is SRF which needs additional risk facto "cox_predict" which was created in the previous row
  predicted_event_prob = 1-srf_survival_prob_for_time(model_1a$model, df_test, fixed_time, oob= oob)
  return (predicted_event_prob)
}


method_1A_cv = function(df, predict.factors, fixed_time = 10, cv_number = 3,
                        internal_cv_k = 3, seed_to_fix = 100){
  time_0 = Sys.time()
  set.seed(seed_to_fix)
  cv_folds = caret::createFolds(df$event, k=cv_number, list = FALSE) #use caret to split into k-folds = cv_steps
  modelstats_train = list(); modelstats_test = list()
  models_for_each_cv=list()
  for (cv_iteration in 1:cv_number){
    print (paste('External loop CV, step number = ', cv_iteration, '/ out of', cv_number))
    
    df_train_cv = df[cv_folds != cv_iteration, ]
    df_test_cv  = df[cv_folds == cv_iteration, ]
    model.tuned = method_1A_train(df_train_cv, predict.factors,fixed_time = fixed_time,
                                  cv_number = internal_cv_k, seed_to_fix=seed_to_fix, 
                                  fast_version = TRUE, oob = TRUE)
    #  calculating prediction of the final 1A model 
    y_predict_test =  method_1A_predict(model.tuned, df_test_cv, fixed_time, seed_to_fix = seed_to_fix, oob = FALSE)
    y_predict_train = method_1A_predict(model.tuned, df_train_cv, fixed_time, seed_to_fix = seed_to_fix, oob= FALSE)
    
    modelstats_test[[cv_iteration]] = method_any_validate(y_predict_test, fixed_time, df_train_cv, df_test_cv, weighted = 1)
    modelstats_train[[cv_iteration]] = method_any_validate(y_predict_train, fixed_time, df_train_cv, df_train_cv, weighted = 1)
    
    models_for_each_cv[[cv_iteration]] = model.tuned
  }
  
  
  df_modelstats_test = data.frame(modelstats_test[[1]])
  df_modelstats_train = data.frame(modelstats_train[[1]])
  for (i in 2:cv_number){df_modelstats_test[i,]= modelstats_test[[i]]; 
  df_modelstats_train[i,]= modelstats_train[[i]]}
  
  df_modelstats_test$test = 1; df_modelstats_train$test = 0
  
  output = list()
  output$test = df_modelstats_test
  output$train = df_modelstats_train
  output$testaverage = sapply(df_modelstats_test,mean,na.rm=1)
  output$trainaverage = sapply(df_modelstats_train,mean,na.rm=1)
  output$tuned_cv_models = models_for_each_cv
  
  time_1 = Sys.time()
  print (time_1 - time_0)
  output$time = time_1 - time_0
  return(output)
}

################## Ensemble 1B #########################
# trains SRF, then passes OOB predictions to Cox 

method_1B_train = function(df_train, predict.factors, fixed_time=10, cv_number = 3, 
                           seed_to_fix = 100, fast_version = TRUE, oob = TRUE, pretrained_srf_model = NULL){
  # if SRF has been trained as a baseline model, the function can receive a tuned SRF, 
  # from which we take OOB predictions straight away
  if (is.null(pretrained_srf_model)){
    pretrained_srf_model_output = method_srf_train(df_train, 
                                                   predict.factors, fixed_time, 
                                                   cv_number, seed_to_fix, fast_version=1, 
                                                   oob=TRUE)
    pretrained_srf_model= pretrained_srf_model_output$model
  }
  df_train$srf_predict = method_srf_predict(pretrained_srf_model, 
                                            df_train, fixed_time, oob = TRUE)
  predict.factors.1B = c(predict.factors, "srf_predict")
  ens_model_1B = method_cox_train(df_train, predict.factors = predict.factors.1B)
  
  output = list()
  output$model = ens_model_1B
  output$model_base = pretrained_srf_model
  return(output)
}

method_1B_predict = function(model_1b, df_test, fixed_time){
  #finding SRF predictions
  df_test$srf_predict = method_srf_predict(model_1b$model_base, df_test, fixed_time, oob = TRUE)
  #this is just a Cox model, but df_test should have SRF predictions in it already  as df_test$srf_predict
  predicted_event_prob = 1-pec::predictSurvProb(model_1b$model, df_test, fixed_time)
  #to return a vector not a matrix from this function if there is only one time 
  
  return (predicted_event_prob)
}

method_1B_cv = function(df, predict.factors, fixed_time = 10, cv_number = 3, 
                        seed_to_fix = 100,pretrained_srf_models = NULL){
  time_0 = Sys.time()
  set.seed(seed_to_fix)
  cv_folds = caret::createFolds(df$event, k=cv_number, list = FALSE) 
  #use caret to split into k-folds = cv_steps
  modelstats_train = list(); modelstats_test = list()
  for (cv_iteration in 1:cv_number){
    df_train_cv = df[cv_folds != cv_iteration, ]
    df_test_cv  = df[cv_folds == cv_iteration, ]
    #if there are pretrained models, we pass them to the train function
    if ((is.null(pretrained_srf_models)==FALSE)&
        (length(pretrained_srf_models)== cv_number)){ 
      model1b.tuned = method_1B_train(df_train_cv, 
                                      predict.factors,
                                      fixed_time = fixed_time, 
                                      cv_number =3, 
                                      seed_to_fix=seed_to_fix, 
                                      fast_version = TRUE, 
                                      oob = TRUE, 
                                      pretrained_srf_model = pretrained_srf_models[[cv_iteration]])
    }else{
      model1b.tuned = method_1B_train(df_train_cv, 
                                      predict.factors,fixed_time = fixed_time,
                                      cv_number =3, seed_to_fix=seed_to_fix, 
                                      fast_version = TRUE, oob = TRUE, 
                                      pretrained_srf_model = "none")
    }
    #  calculating prediction 
    y_predict_test =  method_1B_predict(model1b.tuned, df_test_cv, fixed_time)
    y_predict_train = method_1B_predict(model1b.tuned, df_train_cv, fixed_time)
    
    modelstats_test[[cv_iteration]] = method_any_validate(y_predict_test, 
                                                          fixed_time, 
                                                          df_train_cv, 
                                                          df_test_cv, 
                                                          weighted = 1)
    modelstats_train[[cv_iteration]] = method_any_validate(y_predict_train, 
                                                           fixed_time, 
                                                           df_train_cv, 
                                                           df_train_cv, 
                                                           weighted = 1)
  }
  df_modelstats_test = data.frame(modelstats_test[[1]])
  df_modelstats_train = data.frame(modelstats_train[[1]])
  
  for (i in 2:cv_number){
    df_modelstats_test[i,]= modelstats_test[[i]]
    df_modelstats_train[i,]= modelstats_train[[i]]}
  df_modelstats_test$test = 1; df_modelstats_train$test = 0
  
  output = list()
  output$test = df_modelstats_test
  output$train = df_modelstats_train
  output$testaverage = sapply(df_modelstats_test,mean,na.rm=1)
  output$trainaverage = sapply(df_modelstats_train,mean,na.rm=1)
  time_1 = Sys.time()
  print (time_1 - time_0)
  output$time = time_1 - time_0
  return(output)
}

########## Ensemble 2A ###########
# Shallow tree -> Cox model in each cluster

method_2A_tune = function(df_tune, predictors.rpart , predict.factors,
                          fixed_time=10, n_cv = 3, 
                          maxdepth = 10, minbucket = 0, 
                          cp =0.001, seed_to_fix = 100, verbose = FALSE){
  
  set.seed(seed_to_fix)
  #use caret to split into k-folds = cv_steps
  cv_folds = caret::createFolds(df_tune$event, k=n_cv, list = FALSE) 
  cindex_train = vector(length = n_cv); cindex_test = vector(length = n_cv)
  if (minbucket==0) {minbucket = max(50,dim(df_train_cox_rpart)[1]/10)}
  
  for (cv_iteration in 1:n_cv){
    if (verbose) print (paste('CV tuning step number = ', 
                              cv_iteration, '/out of ', n_cv))
    df_train_cv = df_tune[cv_folds != cv_iteration, ]
    df_test_cv  = df_tune[cv_folds == cv_iteration, ]
    
    rpart_params = rpart.control(minbucket =minbucket ,  maxdepth = maxdepth, cp = cp)
    rpart.m1 = rpart::rpart(as.formula(paste("Surv(time, event) ~", 
                                             paste(predictors.rpart, collapse="+"))),
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
        cox_i$coefficients[is.na(cox_i$coefficients)] = 0
        
      }else{
        cox_i = coxph(as.formula(paste("Surv(time, event) ~", paste(c(predict.factors), collapse="+"))),
                      data = df_train_cv[df_train_cv$cluster_tree== clusters[i],], x= TRUE)
        #if all parameter is of the same value in a cluster, its coeff 
        #is NA and it breaks predictSurvProb, so we replace with 0
        cox_i$coefficients[is.na(cox_i$coefficients)] = 0
        
      }
      #predict event prob in train and test from cox_i
      df_test_cv[df_test_cv$cluster_tree == clusters[i] ,"eventprob"] = 1- 
        pec::predictSurvProb(cox_i, df_test_cv[df_test_cv$cluster_tree == clusters[i] ,] , fixed_time)
      df_train_cv[df_train_cv$cluster_tree == clusters[i] ,"eventprob"] = 1- 
        pec::predictSurvProb(cox_i, df_train_cv[df_train_cv$cluster_tree == clusters[i] ,] , fixed_time)
    }
    #check c-index 
    #print (describe(df_test_cv$eventprob))
    c_test = concordancefit(Surv(df_test_cv$time, df_test_cv$event), -1*df_test_cv$eventprob)$concordance
    c_train = concordancefit(Surv(df_train_cv$time, df_train_cv$event), -1*df_train_cv$eventprob)$concordance
    cindex_test[cv_iteration] = ifelse(is.na(c_test), 0, c_test)
    cindex_train[cv_iteration] = ifelse(is.na(c_train), 0, c_train)
  }
  
  mean_train = ifelse( is.na(mean(cindex_train))==TRUE, 0, mean(cindex_train))
  mean_test = ifelse( is.na(mean(cindex_test))==TRUE, 0, mean(cindex_test))
  
  return(c(mean_train,mean_test))
  
}


method_2A_train = function(df_train, predict.factors, fixed_time=10, 
                           internal_cv_k =3, seed_to_fix = 100, verbose = FALSE){
  p = length(predict.factors)
  n = dim(df_train)[1]
  
  #1) choosing the VIMP variables
  vimp_rfs =  rfsrc(as.formula(paste("Surv(time, event) ~", paste(predict.factors, collapse="+"))),
                    data = df_train,
                    nodesize = 15, ntree = 500, mtry = sqrt(p), nodedepth = NULL,  
                    nsplit = 30,  #to save the time as we only need importance 
                    splitrule = "logrank", statistics= FALSE, membership=TRUE,
                    importance = "permute", block.size = 250 , #https://www.rdocumentation.org/packages/randomForestSRC/versions/2.10.0/topics/rfsrc
                    seed = seed_to_fix)
  #sorting by importance, take first 10 
  var_importance = sort(vimp_rfs$importance, decreasing = TRUE)
  #resulting risk factors are the names of var_importance object
  var_sorted = names(var_importance) 
  #alternative method could be holdout.vimp.rfsrc(as.formula(paste("Surv(time, event) ~", paste(predict.factors, collapse="+"))), df_train, splitrule = "logrank", importance = "permute", ntree= 1000,  seed = seed_to_fix)
  
  
  #2) build the shallow tree: cross-validate the method by the number of VIMP factors and depth of the single tree 
  # number of factors for the tree 3,4,...,10; max tree depth from 3 to 7
  if(p>=3) {p_cv = 3:min(5, p)} else{p_cv=3} #CV by 3,4,...,10 factors for a shallow tree
  maxdepthlist = 2:4
  
  if(n>=200){ minbucket_list = seq(25, min(round(n/4,0),150), by = 50)
  }else{minbucket_list = c(25,50)}
  
  grid_of_values = expand.grid("p_cv" = p_cv, "tree_depth" = maxdepthlist, 
                               "minbucket" = minbucket_list)
  if(verbose) print(paste("Grid size for single tree tuning is", dim(grid_of_values)[1]))
  
  bestcindex = vector(mode = "double", length = dim(grid_of_values)[1]); i=1
  for (i in 1:dim(grid_of_values)[1]){
    params = var_sorted[1:grid_of_values[i, "p_cv"]]
    maxdepth = grid_of_values[i, "tree_depth"]
    minbucket = grid_of_values[i, "minbucket"]
    bestcindex[i] = method_2A_tune(df_train, params, predict.factors, 
                                   fixed_time=fixed_time,
                                   n_cv = internal_cv_k, maxdepth = maxdepth, 
                                   cp =0.005, minbucket = minbucket, 
                                   seed_to_fix = seed_to_fix)[2]
    i=i+1
  }
  #print(grid_of_values[which.max(bestcindex),])
  maxdepth_use = grid_of_values[which.max(bestcindex),"tree_depth"]
  p_use = var_sorted[1:grid_of_values[which.max(bestcindex),"p_cv"]]
  minbucket_use = grid_of_values[which.max(bestcindex),"minbucket"]
  
  #3) GROWING A SINGLE RPART TREE, only use the top-vars for it
  #train tree
  set.seed(seed_to_fix)
  rpart.m = rpart::rpart(as.formula(paste("Surv(time, event) ~", paste(p_use, collapse="+"))),
                         data = df_train, minbucket = minbucket_use, maxdepth = maxdepth_use, cp = 0.001)
  
  clusters = unique(round(predict(rpart.m,df_train),6)); length(clusters) #number of clusters
  df_train$cluster_tree = as.factor(round(predict(rpart.m,df_train),6))
  
  #calculating Cox regressions for each cluster separately - it has different baseline function AND
  #regression parameters, so it allows for different baseline (~time-dependency) and non-linearity in the relationships
  cox_models_in_clusters = list()
  for (i in seq(length(clusters))){
    #calibrate model in cluster[i]
    #i=1
    if (sum(df_train$cluster_tree== clusters[i]) < 30){
      cox_i = coxph(as.formula(paste("Surv(time, event) ~ 1")),
                    data = df_train[df_train$cluster_tree== clusters[i],], x= TRUE)
      cox_i$coefficients[is.na(cox_i$coefficients)] = 0
      cox_models_in_clusters[[i]] = cox_i
    }else{
      cox_i = coxph(as.formula(paste("Surv(time, event) ~", paste(c(predict.factors), collapse="+"))),
                    data = df_train[df_train$cluster_tree== clusters[i],], x= TRUE)
      #if all parameter is of the same value in a cluster, its coeff is NA and it breaks predictSurvProb, so we replace with 0
      cox_i$coefficients[is.na(cox_i$coefficients)] = 0
      cox_models_in_clusters[[i]] = cox_i
    }
  }
  
  output = list()
  output$vimp10 = var_importance[1:min(length(var_importance),15)]
  output$treemodel = rpart.m
  output$coxmodels = cox_models_in_clusters
  output$clusters = clusters
  #calibrate SRF with the best parameters
  return (output)
}

method_2A_predict = function(model_2a , df_test, fixed_time){
  df_test$cluster_tree = as.factor(round(predict(model_2a$treemodel, newdata = df_test),6))
  df_test$eventprob2a = NaN 
  for (i in 1:length(model_2a$clusters)){
    #i=2
    cluster = model_2a$clusters[[i]]
    coxm = model_2a$coxmodels[[i]]
    df_test[df_test$cluster_tree == cluster ,"eventprob2a"] = 1- 
      pec::predictSurvProb(coxm, df_test[df_test$cluster_tree == cluster,] , fixed_time)
    # if all event times in cluster < fixed_time, this function doesnt work, so we calculate manually
    #assuming constant hazard beyond max(time)
    if (sum(is.na(df_test[df_test$cluster_tree == cluster ,"eventprob2a"]))>0) {
      if(max(basehaz(coxm, centered = FALSE)[, "time"]) < fixed_time){
        tmax = max(basehaz(coxm, centered = FALSE)[, "time"])
        b_h = basehaz(coxm, centered = FALSE)[which.max(basehaz(coxm, centered = FALSE)[, "time"]), "hazard"]
        df_test[df_test$cluster_tree == cluster ,"eventprob2a"] = 1- 
          exp(-b_h)^exp(predict(coxm, newdata = df_test[df_test$cluster_tree == cluster ,], 
                                type = "lp", reference = "zero"))
      }
    }
  }
  output = list();   output$predict = df_test$eventprob2a;  output$clusters = df_test$cluster_tree
  #return(output)
  return(df_test$eventprob2a)
}

method_2A_cv = function(df, predict.factors, 
                        fixed_time = 10, 
                        cv_number = 3, 
                        internal_cv_k = 3, 
                        seed_to_fix = 100, verbose = FALSE){
  time_0 = Sys.time()
  set.seed(seed_to_fix)
  cv_folds = caret::createFolds(df$event, k=cv_number, list = FALSE) #use caret to split into k-folds = cv_steps
  modelstats_train = list(); modelstats_test = list()
  models_for_each_cv = list() #saving trained cv models 
  for (cv_iteration in 1:cv_number){
    print (paste('External loop CV, step', cv_iteration, '/ ', cv_number))
    
    df_train_cv = df[cv_folds != cv_iteration, ]
    df_test_cv  = df[cv_folds == cv_iteration, ]
    
    model2a_tuned = method_2A_train(df_train_cv, predict.factors,fixed_time = fixed_time, 
                                    internal_cv_k = internal_cv_k, seed_to_fix=seed_to_fix)
    
    y_predict_test =  method_2A_predict(model2a_tuned, df_test_cv, fixed_time)
    y_predict_train = method_2A_predict(model2a_tuned, df_train_cv, fixed_time)
    
    modelstats_test[[cv_iteration]] = method_any_validate(y_predict_test, fixed_time, df_train_cv, df_test_cv, weighted = 1)
    modelstats_train[[cv_iteration]] = method_any_validate(y_predict_train, fixed_time, df_train_cv, df_train_cv, weighted = 1)
    
    models_for_each_cv[[cv_iteration]]= model2a_tuned
  }
  df_modelstats_test = data.frame(modelstats_test[[1]])
  df_modelstats_train = data.frame(modelstats_train[[1]])
  
  for (i in 2:cv_number){df_modelstats_test[i,]= modelstats_test[[i]]; df_modelstats_train[i,]= modelstats_train[[i]]}
  df_modelstats_test$test = 1; df_modelstats_train$test = 0
  output = list()
  output$test = df_modelstats_test
  output$train = df_modelstats_train
  output$testaverage = sapply(df_modelstats_test,mean,na.rm=1)
  output$trainaverage = sapply(df_modelstats_train,mean,na.rm=1)
  output$tuned_cv_models = models_for_each_cv
  time_1 = Sys.time()
  print (time_1 - time_0)
  output$time = time_1 - time_0
  
  return(output)
}

########## Ensemble 2B - NOT READY ###########
# Shallow tree -> as clusters to the Cox model

method_2B_tune = function(df_tune, params_for_tree , predict.factors, fixed_time=fixed_time,
                          n_cv = 3, nodedepth = nodedepth, cp =0.001, nodesize = nodesize, seed_to_fix = 100){
  #tuning ensemble 2b - fudging the tree in maxdepth, minbucket and risk factors to use to CV,
  # this is to get an average test c-index for one combination, all combinations are checked in method_3_cv() function
  set.seed(seed_to_fix)
  cv_folds = caret::createFolds(df_tune$event, k=n_cv, list = FALSE) #use caret to split into k-folds = cv_steps
  
  cindex_train = vector(length = n_cv); cindex_test = vector(length = n_cv)
  for (cv_iteration in 1:n_cv){
    print (paste('CV step number = ', cv_iteration))
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
          cox_i$coefficients[is.na(cox_i$coefficients)] = 0
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
    c_test = concordancefit(Surv(df_test_cv$time, df_test_cv$event), -1*df_test_cv$eventprob)$concordance
    c_train = concordancefit(Surv(df_train_cv$time, df_train_cv$event), -1*df_train_cv$eventprob)$concordance
    cindex_test[cv_iteration] = ifelse(is.na(c_test), 0, c_test)
    cindex_train[cv_iteration] = ifelse(is.na(c_train), 0, c_train)
  }
  return(c(mean(cindex_train), mean(cindex_test)))
}

method_2B_train = function(df_train, predict.factors, fixed_time=10, seed_to_fix = 100){
  p = length(predict.factors)
  n = dim(df_train)[1]
  #1) choosing the VIMP variables
  vimp_rfs =  rfsrc(as.formula(paste("Surv(time, event) ~", paste(predict.factors, collapse="+"))),
                    data = df_train,
                    nodesize = 15, ntree = 500, mtry = sqrt(p), nodedepth = NULL,  
                    nsplit = 30,  #to save the time as we only need importance 
                    splitrule = "logrank", statistics= FALSE, membership=TRUE,
                    importance = "permute", block.size = 250 , #https://www.rdocumentation.org/packages/randomForestSRC/versions/2.10.0/topics/rfsrc
                    seed = seed_to_fix)
  #sorting by importance, take first 10 
  var_importance = sort(vimp_rfs$importance, decreasing = TRUE)
  #resulting risk factors are the names of var_importance object
  var_sorted = names(var_importance) 
  #alternative method could be holdout.vimp.rfsrc(as.formula(paste("Surv(time, event) ~", paste(predict.factors, collapse="+"))), df_train, splitrule = "logrank", importance = "permute", ntree= 1000,  seed = seed_to_fix)
  
  #2) build the shallow tree: cross-validate the method by the number of VIMP factors and depth of the single tree 
  # number of factors for the tree 3,4,...,10; max tree depth from 3 to 7
  if(p>=3) {p_cv = 3:min(10, p)} else{p_cv=3} #CV by 3,4,...,10 factors for a shallow tree
  maxdepthlist = 2:4
  
  if(n>=200){ minbucket_list = seq(25, min(round(n/4,0),150), by = 50)
  }else{minbucket_list = c(25,50)}
  
  grid_of_values = expand.grid("p_cv" = p_cv, "tree_depth" = maxdepthlist, "nodesize" = nodesize_list)
  print(paste("Grid size for single tree tuning is", dim(grid_of_values)[1]))
  
  bestcindex = vector(mode = "double", length = dim(grid_of_values)[1]); i=1
  for (i in 1:dim(grid_of_values)[1]){
    params = var_sorted[1:grid_of_values[i, "p_cv"]]
    nodedepth = grid_of_values[i, "tree_depth"]
    nodesize = grid_of_values[i, "nodesize"]
    bestcindex[i] = method_2B_tune(df_train, params_for_tree = params, predict.factors = predict.factors, fixed_time=fixed_time,
                                   n_cv = 3, nodedepth = nodedepth, cp =0.001, nodesize = nodesize, seed_to_fix = 100)[2]
    i=i+1
  }
  
  print(grid_of_values[which.max(bestcindex),])
  nodedepth_use = grid_of_values[which.max(bestcindex),"tree_depth"]
  p_use = var_sorted[1:grid_of_values[which.max(bestcindex),"p_cv"]]
  nodesize_use = grid_of_values[which.max(bestcindex),"nodesize"]
  
  #3) GROWING A SINGLE RPART TREE, only use the top-vars for it
  #train tree
  set.seed(seed_to_fix)
  
  rf.dt = rfsrc(as.formula(paste("Surv(time, event) ~", paste(p_use, collapse="+"))),
                data = df_train,
                nodesize = nodesize_use,  # this is AVERAGE size, so we want this to be quite high
                ntree = 1, #only need one tree
                mtry =  length(p_use), #try all predictors
                nodedepth = nodedepth_use,  # calibrated parameter, from 2 till 7
                nsplit = 20, splitrule = "logrank", statistics=TRUE,membership=TRUE,
                importance = "none", #to speed up by switching off VIMP calculations
                seed = seed_to_fix)
  
  predictedscore_rfdt_train =predict(rf.dt, newdata = df_train, node= TRUE)
  predictedscore_rfdt_train$leaf.count #number of clusters
  clusters = unique(predictedscore_rfdt_train$predicted)
  df_train$cluster_tree = as.factor(predictedscore_rfdt_train$predicted)
  
  cox_models_in_clusters = list()
  
  for (i in seq(length(clusters))){
    #calibrate model in cluster[i]
    if (sum(df_train$cluster_tree== clusters[i]) < 30){
      cox_i = coxph(as.formula(paste("Surv(time, event) ~ 1")),
                    data = df_train[df_train$cluster_tree== clusters[i],], x= TRUE)
      cox_i$coefficients[is.na(cox_i$coefficients)] = 0
      cox_models_in_clusters[[i]] = cox_i
    }else{
      cox_i = coxph(as.formula(paste("Surv(time, event) ~", paste(c(predict.factors), collapse="+"))),
                    data = df_train[df_train$cluster_tree== clusters[i],], x= TRUE)
      #if all parameter is of the same value in a cluster, its coeff is NA and it breaks predictSurvProb, so we replace with 0
      cox_i$coefficients[is.na(cox_i$coefficients)] = 0
      cox_models_in_clusters[[i]] = cox_i
    }
  }
  
  output = list()
  output$treemodel = rf.dt
  output$coxmodels = cox_models_in_clusters
  output$clusters = clusters
  #calibrate SRF with the best parameters
  return (output)
}

method_2B_predict = function(model_2b, df_test, fixed_time){
  # get cluster number from RF tree
  predictedscore_rfdt_test = predict(model_2b$treemodel, newdata = df_test, node= TRUE)
  df_test$cluster_tree = as.factor(predictedscore_rfdt_test$predicted)
  #compute event probability from respective Cox model for the cluster
  df_test$eventprob = NaN 
  for (i in 1:length(model_2b$clusters)){
    cluster = model_2b$clusters[[i]]
    coxm = model_2b$coxmodels[[i]]
    df_test[df_test$cluster_tree == cluster ,"eventprob"] = 1- 
      pec::predictSurvProb(coxm, df_test[df_test$cluster_tree == cluster ,] , fixed_time)
  }
  return(df_test$eventprob)
}

method_2B_cv = function(df, predict.factors, fixed_time = 10, cv_number = 3, seed_to_fix = 100){
  time_0 = Sys.time()
  set.seed(seed_to_fix)
  cv_folds = caret::createFolds(df$event, k=cv_number, list = FALSE) #use caret to split into k-folds = cv_steps
  modelstats_train = list(); modelstats_test = list()
  models_for_each_cv = list() #saving trained best SRF to re-use in ensemble 1A
  for (cv_iteration in 1:cv_number){
    print (paste('External loop CV, step ', cv_iteration, "/", cv_number))
    df_train_cv = df[cv_folds != cv_iteration, ]
    df_test_cv  = df[cv_folds == cv_iteration, ]
    
    model2b_tuned = method_2B_train(df_train_cv, predict.factors,fixed_time = fixed_time,seed_to_fix=seed_to_fix)
    
    y_predict_test =  method_2B_predict(model2b_tuned, df_test_cv, fixed_time)
    y_predict_train = method_2B_predict(model2b_tuned, df_train_cv, fixed_time)
    
    modelstats_test[[cv_iteration]] = method_any_validate(y_predict_test, fixed_time, df_train_cv, df_test_cv, weighted = 1)
    modelstats_train[[cv_iteration]] = method_any_validate(y_predict_train, fixed_time, df_train_cv, df_train_cv, weighted = 1)
    
    models_for_each_cv[[cv_iteration]]= model2b_tuned
  }
  df_modelstats_test = data.frame(modelstats_test[[1]])
  df_modelstats_train = data.frame(modelstats_train[[1]])
  
  for (i in 2:cv_number){df_modelstats_test[i,]= modelstats_test[[i]]; df_modelstats_train[i,]= modelstats_train[[i]]}
  df_modelstats_test$test = 1; df_modelstats_train$test = 0
  output = list()
  output$test = df_modelstats_test
  output$train = df_modelstats_train
  output$testaverage = sapply(df_modelstats_test,mean,na.rm=1)
  output$trainaverage = sapply(df_modelstats_train,mean,na.rm=1)
  output$tuned_cv_models = models_for_each_cv
  time_1 = Sys.time()
  print (time_1 - time_0)
  output$time = time_1 - time_0
  return(output)
}

########## Ensemble 3 ###########
# Modified Cox model with the cluster IDs as factors 


method_3_tune = function(df_tune, params_for_tree= c("age", "bmi", "hyp"), predict.factors,
                         fixed_time=fixed_time,
                         n_cv = 3, maxdepth = maxdepth, cp =0.001, 
                         minbucket = minbucket, seed_to_fix = 100){
  #tuning ensemble 3 - fudging the tree in maxdepth, minbucket and risk factors to use to CV,
  # this is to get an average test c-index for one combination, all combinations are checked in method_3_cv() function
  set.seed(seed_to_fix)
  cv_folds = caret::createFolds(df_tune$event, k=n_cv, list = FALSE) #use caret to split into k-folds = cv_steps
  
  cindex_train = vector(length = n_cv); cindex_test = vector(length = n_cv)
  for (cv_iteration in 1:n_cv){
    print (paste('CV step ', cv_iteration, "/", n_cv))
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
    cindex_test[cv_iteration] = concordancefit(Surv(df_test_cv$time, df_test_cv$event), -1*df_test_cv$modcox_lp)$concordance
    cindex_train[cv_iteration] = concordancefit(Surv(df_train_cv$time, df_train_cv$event), -1*df_train_cv$modcox_lp)$concordance
  }
  return(c(mean(cindex_train), mean(cindex_test)))
}

method_3B_tune = function(df_tune, params_for_tree, predict.factors, fixed_time=fixed_time,
                          n_cv = 3, nodedepth = nodedepth, cp =0.001, nodesize = nodesize, seed_to_fix = 100){
  #tuning ensemble 2b - fudging the tree in maxdepth, minbucket and risk factors to use to CV,
  # this is to get an average test c-index for one combination, all combinations are checked in method_3_cv() function
  set.seed(seed_to_fix)
  cv_folds = caret::createFolds(df_tune$event, k=n_cv, list = FALSE) #use caret to split into k-folds = cv_steps
  
  cindex_train = vector(length = n_cv); cindex_test = vector(length = n_cv)
  for (cv_iteration in 1:n_cv){
    print (paste('CV step number = ', cv_iteration))
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
    modified_cox_cv$coefficients[is.na(modified_cox_cv$coefficients)] = 0
    
    df_train_cv$modcox_lp =  predict(modified_cox_cv, newdata = df_train_cv, type = "lp", se.fit = FALSE, reference = "zero") 
    df_test_cv$modcox_lp =   predict(modified_cox_cv, newdata = df_test_cv, type = "lp", se.fit = FALSE, reference = "zero") 
    
    #check c-index 
    cindex_test[cv_iteration] = concordancefit(Surv(df_test_cv$time, df_test_cv$event), -1*df_test_cv$modcox_lp)$concordance
    cindex_train[cv_iteration] = concordancefit(Surv(df_train_cv$time, df_train_cv$event), -1*df_train_cv$modcox_lp)$concordance
  }
  return(c(mean(cindex_train), mean(cindex_test)))
}

method_3_train = function(df_train, predict.factors, fixed_time=10, 
                          internal_cv_k = 3, 
                          seed_to_fix = 100){
  
  p = length(predict.factors)
  n = dim(df_train)[1]
  
  #1) choosing the VIMP variables
  vimp_rfs =  rfsrc(as.formula(paste("Surv(time, event) ~", paste(predict.factors, collapse="+"))),
                    data = df_train,
                    nodesize = 15, ntree = 500, mtry = sqrt(p), nodedepth = NULL,  
                    nsplit = 30,  #to save the time as we only need importance 
                    splitrule = "logrank", statistics= FALSE, membership=TRUE,
                    importance = "permute", block.size = 250 , #https://www.rdocumentation.org/packages/randomForestSRC/versions/2.10.0/topics/rfsrc
                    seed = seed_to_fix)
  #sorting by importance, take first 10 
  var_importance = sort(vimp_rfs$importance, decreasing = TRUE)
  #resulting risk factors are the names of var_importance object
  var_sorted = names(var_importance) 
  #alternative method could be holdout.vimp.rfsrc(as.formula(paste("Surv(time, event) ~", paste(predict.factors, collapse="+"))), df_train, splitrule = "logrank", importance = "permute", ntree= 1000,  seed = seed_to_fix)
  
  
  #2) build the shallow tree: cross-validate the method by the number of VIMP factors and depth of the single tree 
  # number of factors for the tree 3,4,...,10; max tree depth from 3 to 7
  if(p>=3) {p_cv = 3:min(5, p)} else{p_cv=3} #CV by 3,4,...,10 factors for a shallow tree
  maxdepthlist = 2:4
  
  if(n>=200){ minbucket_list = seq(25, min(round(n/4,0),150), by = 50)
  }else{minbucket_list = c(25,50)}
  
  grid_of_values = expand.grid("p_cv" = p_cv, "tree_depth" = maxdepthlist, 
                               "minbucket" = minbucket_list)
  print(paste("Grid size for single tree tuning is", dim(grid_of_values)[1]))
  
  bestcindex = vector(mode = "double", length = dim(grid_of_values)[1]); i=1
  for (i in 1:dim(grid_of_values)[1]){
    params = var_sorted[1:grid_of_values[i, "p_cv"]]
    maxdepth = grid_of_values[i, "tree_depth"]
    minbucket = grid_of_values[i, "minbucket"]
    bestcindex[i] = method_2A_tune(df_train, params, predict.factors, fixed_time=fixed_time,
                                   n_cv = internal_cv_k, maxdepth = maxdepth, cp =0.001, minbucket = minbucket, seed_to_fix = seed_to_fix)[2]
    i=i+1
  }
  #print(grid_of_values[which.max(bestcindex),])
  maxdepth_use = grid_of_values[which.max(bestcindex),"tree_depth"]
  p_use = var_sorted[1:grid_of_values[which.max(bestcindex),"p_cv"]]
  minbucket_use = grid_of_values[which.max(bestcindex),"minbucket"]
  
  #3) GROWING A SINGLE RPART TREE, only use the top-vars for it
  #train tree
  set.seed(seed_to_fix)
  rpart.m = rpart::rpart(as.formula(paste("Surv(time, event) ~", paste(p_use, collapse="+"))),
                         data = df_train, minbucket = minbucket_use, maxdepth = maxdepth_use, cp = 0.001)
  
  clusters = unique(round(predict(rpart.m,df_train),6)); length(clusters) #number of clusters
  df_train$cluster_tree = as.factor(round(predict(rpart.m,df_train),6))
  
  
  #calculating Cox regression with clusters as a new factor
  modified_cox_model = coxph(as.formula(paste("Surv(time, event) ~",
                                              paste(c(predict.factors, "cluster_tree"), collapse="+"))),
                             data = df_train, x=TRUE)
  modified_cox_model$coefficients[is.na(modified_cox_model$coefficients)] = 0
  
  output = list()
  output$treemodel = rpart.m
  output$modcoxmodel = modified_cox_model
  output$clusters = clusters
  #calibrate SRF with the best parameters
  return (output)
}

method_3_predict = function(model_3, df_test, fixed_time){
  #predicts probability of event by fixed_time for df_test from model_3 for ensemble method 3
  # calculate which cluster df_test falls into 
  df_test$cluster_tree = as.factor(round(predict(model_3$treemodel, newdata = df_test),6))
  #compute probabilities from the modified Cox model 
  predictedprob = 1- pec::predictSurvProb(model_3$modcoxmodel, df_test, fixed_time)
  if (sum(is.nan(predictedprob))>0){print ("some probs are missing")  }
  
  return(predictedprob)
}

method_3_cv = function(df, predict.factors, fixed_time = 10, cv_number = 3, internal_cv_k =3, seed_to_fix = 100){
  time_0 = Sys.time()
  set.seed(seed_to_fix)
  cv_folds = caret::createFolds(df$event, k=cv_number, list = FALSE) #use caret to split into k-folds = cv_steps
  modelstats_train = list(); modelstats_test = list()
  modcox_models_for_each_cv = list() 
  cv_iteration=1
  for (cv_iteration in 1:cv_number){
    print (paste('External loop CV, step number = ', cv_iteration))
    df_train_cv = df[cv_folds != cv_iteration, ]
    df_test_cv  = df[cv_folds == cv_iteration, ]
    model3_tuned = method_3_train(df_train_cv, predict.factors, fixed_time=fixed_time,
                                  internal_cv_k=internal_cv_k,seed_to_fix=seed_to_fix)
    
    y_predict_test =  method_3_predict(model3_tuned, df_test_cv,  fixed_time)
    y_predict_train = method_3_predict(model3_tuned, df_train_cv, fixed_time)
    
    modelstats_test[[cv_iteration]]  = method_any_validate(y_predict_test, fixed_time, df_train_cv, df_test_cv, weighted = 1)
    modelstats_train[[cv_iteration]] = method_any_validate(y_predict_train, fixed_time, df_train_cv, df_train_cv, weighted = 1)
    
    modcox_models_for_each_cv[[cv_iteration]] = model3_tuned
  }
  df_modelstats_test = data.frame(modelstats_test[[1]])
  df_modelstats_train = data.frame(modelstats_train[[1]])
  for (i in 2:cv_number){df_modelstats_test[i,]= modelstats_test[[i]]; 
  df_modelstats_train[i,]= modelstats_train[[i]]}
  df_modelstats_test$test = 1; df_modelstats_train$test = 0
  output = list()
  output$test = df_modelstats_test
  output$train = df_modelstats_train
  output$testaverage = sapply(df_modelstats_test,mean,na.rm=1)
  output$trainaverage = sapply(df_modelstats_train,mean,na.rm=1)
  output$tuned_cv_models = modcox_models_for_each_cv
  time_1 = Sys.time()
  print (time_1 - time_0)
  output$time = time_1 - time_0
  return(output)
}


