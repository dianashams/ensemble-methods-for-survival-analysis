
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
#library(pec) # for predictSurvProb replaced
library(dplyr)

library("doParallel")
library("doFuture") #for parallel calculations
num.cores <- detectCores()-1
cluztrr <- makeCluster(num.cores)
registerDoParallel(cl = cluztrr)
registerDoFuture()
plan(cluster, workers = cluztrr)

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

method_any_cv = function(df, predict.factors, train_function, predict_function, valuation_times, 
                         cv_number = 2, seed_for_cv = 2024, parallel=FALSE,
                         model_args = list(), predict_args = list()){  #all arguments apart from data
  
  pb <- txtProgressBar(0, cv_number+2, style=3)  
  
  time_0 = Sys.time()
  set.seed(seed_for_cv)
  #define the method to evaluate
  #mmm = call(train_function, model_args)
  predict.factors = eligible_params(predict.factors, df)
  if (length(predict.factors)==0){print ("No eligible predictors."); return (NULL)}
  cv_folds = caret::createFolds(df$event, k=cv_number, list = FALSE) #use caret to split into k-folds = cv_steps
  modelstats_train = list(); modelstats_test = list();modelstats=list()
  model_list = list()
  
  setTxtProgressBar(pb, 1)
  
  if (!parallel){
    for (cv_iteration in 1:cv_number){
      df_train_cv = df[cv_folds != cv_iteration, ]; dim(df_train_cv)
      df_test_cv  = df[cv_folds == cv_iteration, ]; dim(df_test_cv)
      model_i = do.call(train_function, append(list(df_train_cv, predict.factors), model_args))
      model_list[[cv_iteration]] = model_i
      y_predict_test = do.call(predict_function, append(list(model_i, df_test_cv, valuation_times), predict_args))
      y_predict_train = do.call(predict_function, append(list(model_i, df_train_cv, valuation_times), predict_args))
      modelstats[[cv_iteration]] =rbind(
        "test" = method_any_validate(y_predict_test,valuation_times, df_train_cv, df_test_cv, weighted = 1),
        "train"=method_any_validate(y_predict_train, valuation_times, df_train_cv, df_train_cv, weighted = 1))
      
      setTxtProgressBar(pb, cv_iteration+1)
      
    }
  }else{
    parallel_stats = foreach::foreach(cv_iteration= 1:cv_number, .packages = c("survival", "timeROC"), .errorhandling = "pass"
    )%dopar%{
      df_train_cv = df[cv_folds != cv_iteration, ]
      df_test_cv  = df[cv_folds == cv_iteration, ]
      model_i = do.call(train_function, append(list(df_train_cv, predict.factors), model_args))
      y_predict_test = do.call(predict_function, append(list(model_i, df_test_cv, valuation_times), predict_args))
      y_predict_train = do.call(predict_function, append(list(model_i, df_train_cv, valuation_times), predict_args))
      setTxtProgressBar(pb, cv_iteration+1)
      
      list("model" = model_i,
           "stat" = rbind("test" = method_any_validate(y_predict_test,valuation_times, df_train_cv, df_test_cv, weighted = 1), 
                          "train" = method_any_validate(y_predict_train, valuation_times, df_train_cv, df_train_cv, weighted = 1)))

    }
    model_list = foreach::foreach(cv_iteration = 1:cv_number)%do%{parallel_stats[[cv_iteration]]$model}
    modelstats = foreach::foreach(cv_iteration = 1:cv_number)%do%{parallel_stats[[cv_iteration]]$stat}
  }
  
  #print(modelstats)
  #create data frame with results:
  df_modelstats_test = t(matrix(unlist(lapply(X= modelstats, FUN = function(x) x[1,])),
                                nrow = 7, ncol = length(modelstats)))
  df_modelstats_train = t(matrix(unlist(lapply(X= modelstats, FUN = function(x) x[length(valuation_times)+1,])),
                                 nrow = 7, ncol = length(modelstats)))
  
  
  if(length(valuation_times)>1){ 
    for (t in 2:length(valuation_times)){
      df_modelstats_test = rbind(df_modelstats_test, t(matrix(unlist(lapply(X= modelstats, 
                                                                            FUN = function(x) x[t,])),nrow = 7, ncol = length(modelstats))))
      df_modelstats_train = 
        rbind(df_modelstats_train, t(matrix(unlist(lapply(X= modelstats, 
                                                          FUN = function(x) x[length(valuation_times)+t,])),nrow = 7, ncol = length(modelstats))))
    }
  }
  
  df_modelstats_test = as.data.frame(df_modelstats_test);    names(df_modelstats_test)= names(modelstats[[1]])
  df_modelstats_train = as.data.frame(df_modelstats_train);    names(df_modelstats_train)= names(modelstats[[1]])
  df_modelstats_test$test = 1; df_modelstats_train$test = 0
  df_modelstats_test$cv_n = c(1:cv_number); df_modelstats_train$cv_n = c(1:cv_number)
  
  setTxtProgressBar(pb, cv_number+2)
  close(pb)
  #comprise  output object
  output = list()
  output$test = df_modelstats_test
  output$train = df_modelstats_train
  output$testaverage = sapply(df_modelstats_test[,1:8],mean,na.rm=1)
  output$trainaverage = sapply(df_modelstats_train[,1:8],mean,na.rm=1)
  output$model_list = model_list
  time_1 = Sys.time()
  output$time = time_1 - time_0
  return(output)
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


eligible_params = function(params, df){
  #'This function checks eligible predictors from params list for split
  #' It deletes those which are
  #'1) not in df and
  #'2) taking only 1 value (constants)
  #'! Later may delete collinear factors
  #'
  if (length(params)==0) return (NULL)
  
  # take only columns which are in df
  z = params %in% names(df)
  if (sum(!z)==length(params)){
    return (NULL) #no eligible params
  }else{
    params = params[z]# there are some potentially eligible    
  }
  params_eligible = params
  for (i in 1:length(params)){
    if (length(unique(df[,params[i]])) < 2){
      params_eligible = params_eligible[params_eligible!=params[i]]}
  }
  return (params_eligible)
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


method_cox_train = function(df_train, predict.factors, useCoxLasso = FALSE, fixed_time = NaN){
  # wrapper for coxph() function returning a trained Cox model
  if (useCoxLasso==FALSE){
    cox.m = NULL
    try({
      cox.m  = coxph(as.formula(
        paste("Surv(df_train$time, df_train$event) ~",
              paste(predict.factors, collapse="+"))),
        data =df_train, x = TRUE)
      #!!! We replace NA coefficients with 0 i.e. ignore predictors which Cox couldn't estimate
      # such that  predictions don't brake
      cox.m$coefficients[is.na(cox.m$coefficients)] = 0}, 
      silent = TRUE)
    if(is.null(cox.m)){
      print(paste("Warning: cox.m == NULL, N/Events=", dim(df_train)[1], sum(df_train$event==1)))
      write.csv(df_train, "failing_cox_df_train.csv")}
    return (cox.m)
  }else{
    return(method_coxlasso_train(df_train, predict.factors))
  }
}

method_coxlasso_train <- function(df_train, predict.factors, fixed_time = NaN){
  #' This function trains cox-lasso and re-trains 
  #' cox on non-zero predictors
  #' 
  library(glmnet)
  cox.m = NULL
  try({
    cv10 =glmnet::cv.glmnet(as.matrix(df_train[predict.factors]), 
                            Surv(df_train$time, df_train$event), 
                            family="cox", nfold=5, alpha = 1)
    new.predictors = rownames(coef(cv10, s="lambda.min"))[as.matrix(coef(cv10, s="lambda.min"))!=0]
    if (length(new.predictors)==0) {
      print ("0 predictors in lasso!")
      cox.m  = survival::coxph(Surv(df_train$time, df_train$event) ~ 1,data = df_train, x = TRUE)
    }else{
      f = as.formula(paste("Surv(df_train$time, df_train$event) ~",paste(new.predictors, collapse="+")))
      cox.m  = survival::coxph(f,data = df_train, x = TRUE)
      #!!!!! We replace NA coefficients with 0 
      # i.e. ignore predictors which the Cox model couldn't estimate)
      cox.m$coefficients[is.na(cox.m$coefficients)] = 0
    }
  }, silent = TRUE)
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
  #compute event probability for times:
  #create placeholder
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

#system.time({method_cox_cv(d_xt, params, fixed_time = seq(1,10,length.out=50), parallel=FALSE)})
## user  system elapsed 
## 24.76    0.39   25.58 
#system.time({method_cox_cv(d_xt, params, fixed_time = seq(1,10,length.out=50), parallel=TRUE)})
## user  system elapsed 
## 0.47    0.12    7.05 

method_cox_cv = function(df, predict.factors, fixed_time = NaN, 
                         cv_number = 5, seed_to_fix = 2024, useCoxLasso= FALSE, parallel=FALSE){
  
  predict.factors = eligible_params(predict.factors, df)
  
  # defining output for fixed_time
  if (sum(is.nan(fixed_time))>0){fixed_time = round(quantile(df_train[df_train$event==1, "time"], 0.8),1)}
  
  if (length(predict.factors)==0){print ("No eligible predictors."); return (NULL)}
  time_0 = Sys.time()
  set.seed(seed_to_fix)
  cv_folds = caret::createFolds(df$event, k=cv_number, list = FALSE) #use caret to split into k-folds = cv_steps
  modelstats_train = list(); modelstats_test = list();modelstats=list()
  
  if (!parallel){
  for (cv_iteration in 1:cv_number){
    df_train_cv = df[cv_folds != cv_iteration, ]; dim(df_train_cv)
    df_test_cv  = df[cv_folds == cv_iteration, ]; dim(df_test_cv)
    cox.model = method_cox_train(df_train_cv, 
                                 eligible_params(predict.factors, df_train_cv), 
                                 useCoxLasso= useCoxLasso)
    y_predict_test = method_cox_predict(cox.model, df_test_cv, fixed_time)
    y_predict_train = method_cox_predict(cox.model, df_train_cv, fixed_time)
    modelstats[[cv_iteration]] =rbind(
      "test" = method_any_validate(y_predict_test,fixed_time, df_train_cv, df_test_cv, weighted = 1),
      "train"=method_any_validate(y_predict_train, fixed_time, df_train_cv, df_train_cv, weighted = 1))
    }
  }else{
    modelstats = foreach::foreach(cv_iteration= 1:cv_number, .packages = c("survival", "timeROC"), .errorhandling = "pass"
    )%dopar%{
      df_train_cv = df[cv_folds != cv_iteration, ]
      df_test_cv  = df[cv_folds == cv_iteration, ]
      cox.model = method_cox_train(df_train_cv, eligible_params(predict.factors, df_train_cv), 
                                   useCoxLasso= useCoxLasso)
      y_predict_test = method_cox_predict(cox.model, df_test_cv, fixed_time)
      y_predict_train = method_cox_predict(cox.model, df_train_cv, fixed_time)
      rbind("test" = method_any_validate(y_predict_test,fixed_time, df_train_cv, df_test_cv, weighted = 1), 
      "train" = method_any_validate(y_predict_train, fixed_time, df_train_cv, df_train_cv, weighted = 1))
    }
  }

  #create data frame with results:
  df_modelstats_test = t(matrix(unlist(lapply(X= modelstats, FUN = function(x) x[1,])),
                                nrow = 7, ncol = length(modelstats)))
  df_modelstats_train = t(matrix(unlist(lapply(X= modelstats, FUN = function(x) x[length(fixed_time)+1,])),
                                   nrow = 7, ncol = length(modelstats)))
  if(length(fixed_time)>1){ 
    for (t in 2:length(fixed_time)){
      df_modelstats_test = rbind(df_modelstats_test, t(matrix(unlist(lapply(X= modelstats, 
          FUN = function(x) x[t,])),nrow = 7, ncol = length(modelstats))))
      df_modelstats_train = 
        rbind(df_modelstats_train, t(matrix(unlist(lapply(X= modelstats, 
          FUN = function(x) x[length(fixed_time)+t,])),nrow = 7, ncol = length(modelstats))))
    }
  }
  df_modelstats_test = as.data.frame(df_modelstats_test);    names(df_modelstats_test)= names(modelstats[[1]])
  df_modelstats_train = as.data.frame(df_modelstats_train);    names(df_modelstats_train)= names(modelstats[[1]])
  df_modelstats_test$test = 1; df_modelstats_train$test = 0
  df_modelstats_test$cv_n = c(1:cv_number); df_modelstats_train$cv_n = c(1:cv_number)
  
  #comprise  output object
  output = list()
  output$test = df_modelstats_test
  output$train = df_modelstats_train
  output$testaverage = sapply(df_modelstats_test[,1:8],mean,na.rm=1)
  output$trainaverage = sapply(df_modelstats_train[,1:8],mean,na.rm=1)
  time_1 = Sys.time()
  output$time = time_1 - time_0
  return(output)
}

############# Augmented CoxPH with Fractional Polynomials ###############

method_coxmfp_train = function(df_train, predict.factors, verbose = FALSE){
  #Cox with fractional polynomials, returns final Cox model with the selected fp() risk factors
  if (length(predict.factors)>20) {
    print ("Too many factors (>20) to perform MFP, default to baseline Cox")
    cox.mfp  = coxph(as.formula(paste("Surv(df_train$time, df_train$event) ~",
                                      paste(predict.factors, collapse="+"))),
                     data =df_train, x = TRUE)
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
  predicted_event_prob = method_cox_predict(model_coxmfp, df_test, fixed_time)
  return (predicted_event_prob)
}

method_coxmfp_cv = function(df, predict.factors, fixed_time = NaN, cv_number = 5,
                            seed_to_fix = 100, parallel= FALSE){
  # cross-validating CoxPh model, returns performance measures for each CV and averaged metrics
  time_0 = Sys.time()
  set.seed(seed_to_fix)
  
  # defining output for fixed_time
  if (sum(is.nan(fixed_time))>0){fixed_time = round(quantile(df[df$event==1, "time"], 0.85),1)}
  
  cv_folds = caret::createFolds(df$event, k=cv_number, list = FALSE) #use caret to split into k-folds = cv_steps
  modelstats_train = list(); modelstats_test = list()
  models_for_each_cv = list()
  
  if (!parallel){
    for (cv_iteration in 1:cv_number){
      df_train_cv = df[cv_folds != cv_iteration, ]
      df_test_cv  = df[cv_folds == cv_iteration, ]
      cox.model = method_coxmfp_train(df_train_cv, eligible_params(predict.factors,df_train_cv))
    
      models_for_each_cv[[cv_iteration]] = cox.model
    
      y_predict_test = method_coxmfp_predict(cox.model, df_test_cv, fixed_time)
      y_predict_train = method_coxmfp_predict(cox.model, df_train_cv, fixed_time)
      modelstats_test[[cv_iteration]] = method_any_validate(y_predict_test, fixed_time, df_train_cv, df_test_cv, weighted = 1)
      modelstats_train[[cv_iteration]] = method_any_validate(y_predict_train, fixed_time, df_train_cv, df_train_cv, weighted = 1)
    }
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
  
  if (oob) { predicted_matrix  =  predict(rfmodel, newdata = df_to_predict, 
                                          ensemble = "oob", outcome= "test")
  } else { predicted_matrix  =  predict(rfmodel, newdata = df_to_predict)}
  
  j_for_fixedtime = match(1, round(predicted_matrix$time.interest,1) == fixed_time, nomatch = -100);
  
  if (j_for_fixedtime == -100){#print("no fixed time match was found, using closest")
    j_for_fixedtime = which.min(abs(predicted_matrix$time.interest - fixed_time))}
  
  if (oob) { y_predicted = predicted_matrix$survival.oob[ ,j_for_fixedtime]
  }else{ y_predicted = predicted_matrix$survival[, j_for_fixedtime]}
  
  
  return(y_predicted)
}

# s = method_srf_train(d_xt, params)
# p = method_srf_predict(s, d_xt, 5)

srf_tune = function(df_tune,  cv_number =3, 
                    predict.factors, fixed_time = NaN, 
                    seed_to_fix = 100,mtry= c(3,4,5), 
                    nodesize = c(10,20,50),nodedepth = c(100),
                    verbose = FALSE, oob = TRUE){
  #function to tune survival random forest by mtry, nodesize and nodedepth grid 
  # if oob = TRUE, there is no CV !!! as OOB does the job already 
  
  #take out the factors which are not in df_tune or the ones which take 1 value
  predict.factors = eligible_params(predict.factors,df_tune)
  
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
  if (sum(is.nan(fixed_time)>0)| length(fixed_time)>1){  #not implemented for multiple time tuning yet
    fixed_time = round(quantile(df_tune[df_tune$event==1, "time"], 0.85),1)}
  
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
                            fixed_time = NaN, inner_cv = 3, 
                            seed_to_fix = 100, fast_version = TRUE, oob = TRUE, verbose = FALSE){
  #take out predictors which are not in df_train or constant
  predict.factors = eligible_params(predict.factors, df_train)
  #for now only for best AUC but can be amended for brier score or cindex
  #defining the tuning grid for SRF 
  p = length(predict.factors) #number of predictors
  n = dim(df_train)[1] 
  #mtry grid
  mtry_default = round(sqrt(p),0)
  
  # defining output for fixed_time
  if ( sum(is.nan(fixed_time))>0| (length(fixed_time)>1)){
    fixed_time = round(quantile(df_train[df_train$event==1, "time"], 0.85),1)}
  
  if (p<=10) {mtry = c(2,3,4,5)}else{if(p<=25){mtry = c(3,5,7,10,15)}else{
    mtry = c(round(p/10,0),round(p/5,0), round(p/3,0), round(p/2,0),mtry_default)}}
  
  #minimum nodesize grid
  nodesize = seq(min(15, round(n/6-1,0)), max(min(n/10,50),30), 5)
  #nodedepth grid
  nodedepth = c(50) #we don't tune this so just a big number 
  if (verbose) {print (paste("mtry", mtry, "nodedepth", nodedepth, "nodesize", nodesize))}
  
  if (fast_version == TRUE) {
    #take recommended mtry and check the best depth and node size 
    
    tune1 = srf_tune(df_train, cv_number = inner_cv, eligible_params(predict.factors,df_train), 
                     fixed_time = fixed_time, seed_to_fix = seed_to_fix, 
                     mtry = mtry_default,nodesize = nodesize, 
                     nodedepth = nodedepth, oob = oob)
    nodesize_best = as.integer(tune1$bestauc["nodesize"])
    nodedepth_best = as.integer(tune1$bestauc["nodedepth"])
    #using the depth and size check the best mtry  
    tune2 = srf_tune(df_train, cv_number = inner_cv,eligible_params(predict.factors,df_train),
                     fixed_time = fixed_time, seed_to_fix = seed_to_fix, 
                     mtry = mtry, nodesize = nodesize_best,  
                     nodedepth = nodedepth_best , oob = oob)
    mtry_best = tune2$bestauc["mtry"] 
    best_combo_stat = tune2$bestauc
    modelstatsall = rbind(tune1$modelstats, tune2$modelstats)
  }else{
    tuneall = srf_tune(df_train, cv_number = inner_cv,
                       eligible_params(predict.factors,df_train),
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
                                      paste(eligible_params(predict.factors,df_train), collapse="+"))),
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

method_srf_predict = function(model_srf, df_test, fixed_time, oob= FALSE){
  if (class(model_srf)=="list") {srf = model_srf$model}else{srf = model_srf}
  
  if (length(fixed_time)==1) {return(1- srf_survival_prob_for_time(srf, df_test, fixed_time, oob= oob))}
  
  predicted_event_prob = matrix(nrow = dim(df_test)[1], ncol = length(fixed_time))
  for (t in 1:length(fixed_time)){
    predicted_event_prob[,t] = 1- srf_survival_prob_for_time(srf, df_test, fixed_time[t], oob= oob)
  }
  colnames(predicted_event_prob) = round(fixed_time,3)
  return (predicted_event_prob)
}

method_srf_cv = function(df, predict.factors, fixed_time = NaN, 
                         cv_number = 3, 
                         inner_cv = 3, 
                         seed_to_fix = 100, parallel=FALSE){
  
  time_0 = Sys.time()
  set.seed(seed_to_fix)
  
  # defining output for fixed_time
  if (sum(is.nan(fixed_time))>0){fixed_time = round(quantile(df[df$event==1, "time"], 0.85),1)}
  
  predict.factors = eligible_params(predict.factors,df)
  if(length(predict.factors)==0){print ("No eliible params"); return (NULL)}
  
  #use caret to split into k-folds = cv_steps
  cv_folds = caret::createFolds(df$event, k=cv_number, list = FALSE) 
  modelstats_train = list(); modelstats_test = list()
  srf_models_for_each_cv = list() #saving trained best SRF to re-use in ensemble 1A
  
  print (paste("Cross-validating Survival Random Forest with", cv_number,
               "outer loops, and ",inner_cv,"inner loops for model tuning"))
  
  if (!parallel){
  for (cv_iteration in 1:cv_number){
    print (paste('External loop CV, step number = ', cv_iteration, '/ out of', cv_number))
    
    df_train_cv = df[cv_folds != cv_iteration, ]
    df_test_cv  = df[cv_folds == cv_iteration, ]
    
    srf.model.tuned = method_srf_train(df_train_cv, predict.factors,fixed_time = fixed_time,
                                       inner_cv = inner_cv, seed_to_fix=seed_to_fix,
                                       fast_version = TRUE, oob = TRUE)
    
    y_predict_test =  method_srf_predict(srf.model.tuned, df_test_cv,
                                         fixed_time,  oob = FALSE)
    y_predict_train = method_srf_predict(srf.model.tuned, df_train_cv, 
                                         fixed_time,  oob= FALSE)
    
    modelstats_test[[cv_iteration]] = method_any_validate(y_predict_test, 
                                                          fixed_time, df_train_cv, df_test_cv, weighted = 1)
    modelstats_train[[cv_iteration]] = method_any_validate(y_predict_train, 
                                                           fixed_time, df_train_cv, df_train_cv, weighted = 1)
    
    srf_models_for_each_cv[[cv_iteration]]= srf.model.tuned$model
  }
  }else{ #parallel cv
    
    
  } 
  
  df_modelstats_test = data.frame(modelstats_test[[1]])
  df_modelstats_train = data.frame(modelstats_train[[1]])
  
  for (i in 2:cv_number){df_modelstats_test[i,]= modelstats_test[[i]]
                         df_modelstats_train[i,]= modelstats_train[[i]]}
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

method_1A_train = function(df_train, predict.factors, fixed_time=NaN, inner_cv = 3, 
                           seed_to_fix = 100, fast_version = TRUE, oob = TRUE,
                           useCoxLasso = FALSE, var_importance_calc=1){
  #the function trains Cox model, then adds its predictions into Survival Random Forest model
  # to mimic stacking procedure and reduce overfitting,
  # we train Cox model on 0.9 of the data and predict on the rest 0.1 for each 1/10s fold 
  # so we pass out-of-the-bag prediction to SRF 
  predict.factors = eligible_params(predict.factors,df_train)
  if(length(predict.factors)==0){print ("No eliible params"); return (NULL)}
  
  # defining output for fixed_time
  if (sum(is.nan(fixed_time))>0){fixed_time = round(quantile(df_train[df_train$event==1, "time"], 0.85),1)}
  
  #creating folds
  cv_folds = caret::createFolds(df_train$event, k=10, list = FALSE) 
  cindex_train = vector(length = 10); cindex_test = vector(length = 10)
  for (cv_iteration in 1:10){
    cox_train = df_train[cv_folds != cv_iteration, ]
    cox_oob  =  df_train[cv_folds == cv_iteration, ]
    #train cox model on cox_train
    cox_m_cv = method_cox_train(cox_train, eligible_params(predict.factors, cox_train),
                                useCoxLasso =useCoxLasso)
    #predict for cox_oob
    cox_predict_oob = method_cox_predict(cox_m_cv, cox_oob, fixed_time)
    #adding Cox prediction to the df_train in the column "cox_predict"
    df_train[cv_folds == cv_iteration, "cox_predict"] = cox_predict_oob
  }
  
  ##alternatively - just use all the data and pass apparent predictions to SRF 
  cox_model_for1a = method_cox_train(df_train, eligible_params(predict.factors, cox_train),
                                     useCoxLasso =useCoxLasso)
  ##df_train$cox_predict = method_cox_predict(cox_model_for1a, df_train, fixed_time )
  
  #adding new factor and tuning SRF model with this added factor using srf_train 
  predict.factors.1A = c(predict.factors, "cox_predict")
  srf_model_for1a = method_srf_train(df_train, predict.factors = predict.factors.1A, 
                                     fixed_time=fixed_time, inner_cv = inner_cv,
                                     seed_to_fix = seed_to_fix, fast_version = fast_version, 
                                     oob = oob)
  
  if(var_importance_calc){
    v = vimp(srf_model_for1a$model, importance = "permute", seed = seed_to_fix)
    var_importance = sort(v$importance, decreasing = TRUE)
    srf_model_for1a$vimp10  = var_importance[1:min(length(var_importance),15)]
  }else{srf_model_for1a$vimp10  = c(NaN)}
  srf_model_for1a$model_base  = cox_model_for1a
  
  return(srf_model_for1a)
}


method_1A_predict = function(model_1a, df_test, fixed_time, 
                             seed_to_fix = 100, oob= FALSE,useCoxLasso =FALSE){
  #use model_base with the base Cox model to find cox_predict 
  df_test$cox_predict = method_cox_predict(model_1a$model_base, 
                                           df_test, fixed_time)
  #' now use "model" which is SRF which needs additional risk factor 
  #' "cox_predict" which was created in the previous row
  predicted_event_prob = 
    1-srf_survival_prob_for_time(model_1a$model, df_test, fixed_time, oob= oob)
  return (predicted_event_prob)
}


method_1A_cv = function(df, predict.factors, fixed_time = NaN, cv_number = 3,
                        internal_cv_k = 3, seed_to_fix = 100,useCoxLasso = FALSE){
  time_0 = Sys.time()
  set.seed(seed_to_fix)
  
  if (sum(is.nan(fixed_time))>0){fixed_time = round(quantile(df[df$event==1, "time"], 0.85),1)}
  
  predict.factors = eligible_params(predict.factors,df)
  if(length(predict.factors)==0){print ("No eliible params"); return (NULL)}
  
  cv_folds = caret::createFolds(df$event, k=cv_number, list = FALSE) #use caret to split into k-folds = cv_steps
  modelstats_train = list(); modelstats_test = list()
  models_for_each_cv=list()
  for (cv_iteration in 1:cv_number){
    print (paste('External loop CV, step number = ', cv_iteration, '/ out of', cv_number))
    
    df_train_cv = df[cv_folds != cv_iteration, ]
    df_test_cv  = df[cv_folds == cv_iteration, ]
    model.tuned = method_1A_train(df_train_cv, predict.factors,fixed_time = fixed_time,
                                  cv_number = internal_cv_k, seed_to_fix=seed_to_fix, 
                                  fast_version = TRUE, oob = TRUE, useCoxLasso = useCoxLasso)
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

method_1B_train = function(df_train, predict.factors, fixed_time=NaN, inner_cv = 3, 
                           seed_to_fix = 100, fast_version = TRUE, 
                           oob = TRUE, pretrained_srf_model = NULL,
                           useCoxLasso = FALSE){
  
  if (sum(is.nan(fixed_time))>0){fixed_time = round(quantile(df_train[df_train$event==1, "time"], 0.85),1)}
  
  predict.factors = eligible_params(predict.factors,df_train)
  if(length(predict.factors)==0){print ("No eliible params"); return (NULL)}
  
  # if SRF has been trained as a baseline model, the function can receive a tuned SRF, 
  # from which we take OOB predictions straight away
  if (is.null(pretrained_srf_model)){
    pretrained_srf_model_output = method_srf_train(df_train, 
                                                   predict.factors, fixed_time, 
                                                   inner_cv, seed_to_fix, 
                                                   fast_version=1,oob=TRUE)
    pretrained_srf_model= pretrained_srf_model_output$model
  }
  df_train$srf_predict = method_srf_predict(pretrained_srf_model, 
                                            df_train, fixed_time, oob = TRUE)
  predict.factors.1B = c(predict.factors, "srf_predict")
  ens_model_1B = method_cox_train(df_train, predict.factors = predict.factors.1B,useCoxLasso =useCoxLasso)
  
  output = list()
  output$model = ens_model_1B
  output$model_base = pretrained_srf_model
  return(output)
}

method_1B_predict = function(model_1b, df_test, fixed_time){
  #finding SRF predictions
  df_test$srf_predict = method_srf_predict(model_1b$model_base, df_test, fixed_time, oob = TRUE)
  #this is just a Cox model, but df_test should have SRF predictions in it already  as df_test$srf_predict
  predicted_event_prob = method_cox_predict(model_1b$model, df_test, fixed_time)
  #to return a vector not a matrix from this function if there is only one time 
  
  return (predicted_event_prob)
}

method_1B_cv = function(df, predict.factors, fixed_time = NaN, cv_number = 3, inner_cv=3,
                        seed_to_fix = 100,pretrained_srf_models = NULL,
                        useCoxLasso = FALSE){
  time_0 = Sys.time()
  set.seed(seed_to_fix)
  
  if (sum(is.nan(fixed_time))>0){fixed_time = round(quantile(df[df$event==1, "time"], 0.85),1)}
  
  predict.factors = eligible_params(predict.factors,df)
  if(length(predict.factors)==0){print ("No eliible params"); return (NULL)}
  
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
                                      cv_number =inner_cv, 
                                      seed_to_fix=seed_to_fix, 
                                      fast_version = TRUE, 
                                      oob = TRUE, 
                                      pretrained_srf_model = pretrained_srf_models[[cv_iteration]],
                                      useCoxLasso = useCoxLasso)
    }else{
      model1b.tuned = method_1B_train(df_train_cv, 
                                      predict.factors,fixed_time = fixed_time,
                                      cv_number =inner_cv, seed_to_fix=seed_to_fix, 
                                      fast_version = TRUE, oob = TRUE, 
                                      pretrained_srf_model = "none",
                                      useCoxLasso = useCoxLasso)
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

method_2A_tune = function(df_tune, predictors.rpart, predict.factors,
                          fixed_time=NaN, inner_cv = 3, 
                          maxdepth = 3, minbucket = 25, 
                          cp =0.001, seed_to_fix = 100, verbose = FALSE,
                          useCoxLasso = FALSE){
  
  set.seed(seed_to_fix)
  if (sum(is.nan(fixed_time))>0){fixed_time = round(quantile(df_tune[df_tune$event==1, "time"], 0.85),1)}
  #use caret to split into k-folds = cv_steps
  cv_folds = caret::createFolds(df_tune$event, k=inner_cv, list = FALSE) 
  cindex_train = vector(length = inner_cv); cindex_test = vector(length = inner_cv)
  if (minbucket==0) {minbucket = max(50,dim(df_train_cox_rpart)[1]/10)}
  
  for (cv_iteration in 1:inner_cv){
    if (verbose) print (paste('CV tuning step number = ', 
                              cv_iteration, '/out of ', inner_cv))
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
      coxi=  method_cox_train(df_train_cv[df_train_cv$cluster_tree== clusters[i],],
                              predict.factors,useCoxLasso = useCoxLasso )
      if(!is.null(coxi)){
        #predict event prob in train and test from coxi
        df_test_cv[df_test_cv$cluster_tree == clusters[i] ,"eventprob"] = 
          method_cox_predict(coxi, df_test_cv[df_test_cv$cluster_tree == clusters[i] ,], fixed_time)
        df_train_cv[df_train_cv$cluster_tree == clusters[i] ,"eventprob"] = 
          method_cox_predict(coxi, df_train_cv[df_train_cv$cluster_tree == clusters[i] ,], fixed_time)
      }
    }
    #check c-index 
    #print (describe(df_test_cv$eventprob))
    c_test = concordancefit(Surv(df_test_cv$time, df_test_cv$event), -1*df_test_cv$eventprob)$concordance
    c_train = concordancefit(Surv(df_train_cv$time, df_train_cv$event), -1*df_train_cv$eventprob)$concordance
    if (is.null(c_test)|is.null(c_train)){
      cindex_test[cv_iteration] = NaN;   cindex_train[cv_iteration] = NaN
    }else{cindex_test[cv_iteration] = c_test;   cindex_train[cv_iteration] = c_train}
    
  }
  
  mean_train = ifelse( is.na(mean(cindex_train, na.rm=1)), 0, mean(cindex_train, na.rm=1))
  mean_test =  ifelse( is.na(mean(cindex_test, na.rm=1)), 0, mean(cindex_test, na.rm=1))
  
  return(c(mean_train,mean_test))
  
}


method_2A_train = function(df_train, predict.factors, fixed_time= NaN, 
                           inner_cv =3, seed_to_fix = 100, verbose = FALSE,
                           useCoxLasso = FALSE, var_importance = NULL){
  
  predict.factors = eligible_params(predict.factors,df_train)
  if(length(predict.factors)==0){print ("No eliible params"); return (NULL)}
  
  if (sum(is.nan(fixed_time))>0){fixed_time = round(quantile(df_train[df_train$event==1, "time"], 0.85),1)}
  
  p = length(predict.factors)
  n = dim(df_train)[1]
  
  #1) choosing the VIMP variables
  if (is.null(var_importance)){
    vimp_rfs =  rfsrc(as.formula(paste("Surv(time, event) ~", paste(predict.factors, collapse="+"))),
                      data = df_train,
                      nodesize = 15, ntree = 500, mtry = sqrt(p), nodedepth = NULL,  
                      nsplit = 30,  #to save the time as we only need importance 
                      splitrule = "logrank", statistics= FALSE, membership=TRUE,
                      importance = "permute", block.size = 250 , 
                      #https://www.rdocumentation.org/packages/randomForestSRC/versions/2.10.0/topics/rfsrc
                      seed = seed_to_fix)
    #sorting by importance, take first 10 
    var_importance = sort(vimp_rfs$importance, decreasing = TRUE)}
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
  
  bestcindex = vector(mode = "double", length = dim(grid_of_values)[1])
  for (i in 1:dim(grid_of_values)[1]){
    params = var_sorted[1:grid_of_values[i, "p_cv"]]
    maxdepth = grid_of_values[i, "tree_depth"]
    minbucket = grid_of_values[i, "minbucket"]
    bestcindex[i] = method_2A_tune(df_train, params, predict.factors, 
                                   fixed_time=fixed_time,
                                   inner_cv = inner_cv, maxdepth = maxdepth, 
                                   cp =0.005, minbucket = minbucket, 
                                   seed_to_fix = seed_to_fix, useCoxLasso = useCoxLasso)[2]
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
    coxi=  method_cox_train(df_train[df_train$cluster_tree== clusters[i],],
                            predict.factors,useCoxLasso = useCoxLasso )
    cox_models_in_clusters[[i]] = coxi
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
    df_test[df_test$cluster_tree == cluster ,"eventprob2a"] = 
      method_cox_predict(coxm, df_test[df_test$cluster_tree == cluster,] , fixed_time)
  }
  return(df_test$eventprob2a)
}

method_2A_cv = function(df, predict.factors, 
                        fixed_time = NaN, 
                        cv_number = 3, 
                        inner_cv = 3, 
                        seed_to_fix = 100, verbose = FALSE,useCoxLasso = FALSE){
  time_0 = Sys.time()
  set.seed(seed_to_fix)
  
  if (sum(is.nan(fixed_time))>0){fixed_time = round(quantile(df[df$event==1, "time"], 0.85),1)}
  
  predict.factors = eligible_params(predict.factors,df)
  if(length(predict.factors)==0){print ("No eligible params"); return (NULL)}
  
  cv_folds = caret::createFolds(df$event, k=cv_number, list = FALSE) #use caret to split into k-folds = cv_steps
  modelstats_train = list(); modelstats_test = list()
  models_for_each_cv = list() #saving trained cv models 
  for (cv_iteration in 1:cv_number){
    print (paste('External loop CV, step', cv_iteration, '/ ', cv_number))
    
    df_train_cv = df[cv_folds != cv_iteration, ]
    df_test_cv  = df[cv_folds == cv_iteration, ]
    
    model2a_tuned = method_2A_train(df_train_cv, predict.factors,fixed_time = fixed_time, 
                                    inner_cv = inner_cv, seed_to_fix=seed_to_fix,
                                    useCoxLasso = useCoxLasso)
    
    y_predict_test =  method_2A_predict(model2a_tuned, df_test_cv, fixed_time)
    y_predict_train = method_2A_predict(model2a_tuned, df_train_cv, fixed_time)
    
    modelstats_test[[cv_iteration]] = method_any_validate(y_predict_test, fixed_time, df_train_cv, df_test_cv, weighted = 1)
    modelstats_train[[cv_iteration]] = method_any_validate(y_predict_train, fixed_time, df_train_cv, df_train_cv, weighted = 1)
    
    models_for_each_cv[[cv_iteration]]= model2a_tuned
  }
  df_modelstats_test = data.frame(modelstats_test[[1]])
  df_modelstats_train = data.frame(modelstats_train[[1]])
  
  for (i in 2:cv_number){
    df_modelstats_test[i,]= modelstats_test[[i]]
    df_modelstats_train[i,]= modelstats_train[[i]]
  }
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

method_3_tune = function(df_tune, params_for_tree, predict.factors,
                         fixed_time=NaN,
                         inner_cv = 3, maxdepth, cp =0.001, 
                         minbucket, seed_to_fix = 100,
                         useCoxLasso = FALSE){
  #tuning ensemble 3 - fudging the tree in maxdepth, minbucket and risk factors to use to CV,
  # this is to get an average test c-index for one combination, all combinations are checked in method_3_cv() function
  set.seed(seed_to_fix)
  cv_folds = caret::createFolds(df_tune$event, k=inner_cv, list = FALSE) #use caret to split into k-folds = cv_steps
  
  if (sum(is.nan(fixed_time))>0){fixed_time = round(quantile(df_tune[df_tune$event==1, "time"], 0.85),1)}
  
  cindex_train = vector(length = inner_cv); cindex_test = vector(length = inner_cv)
  for (cv_iteration in 1:inner_cv){
    #print (paste('CV step ', cv_iteration, "/", inner_cv))
    df_train_cv = df_tune[cv_folds != cv_iteration, ]
    df_test_cv  = df_tune[cv_folds == cv_iteration, ]
    
    #grow shallow tree
    tree_cv = rpart::rpart(as.formula(paste("Surv(time, event) ~", paste(params_for_tree, collapse="+"))),
                           data = df_train_cv, minbucket = minbucket, maxdepth = maxdepth, cp = 0.001)
    
    #add the factor- id of a tree leaf, taking a predicted value 
    df_train_cv$cluster_tree = as.factor(round(predict(tree_cv, df_train_cv),6))
    df_test_cv$cluster_tree = as.factor(round(predict(tree_cv, df_test_cv),6))
    
    if(length(unique(df_train_cv$cluster_tree)==1)) {return (c(NaN, NaN))} 
    #add them to the Cox model with all the risk factors
    modified_cox_cv= method_cox_train(df_train_cv,c(predict.factors, "cluster_tree"), useCoxLasso = useCoxLasso )
    df_train_cv$modcox_lp =  predict(modified_cox_cv, newdata = df_train_cv, 
                                     type = "lp", se.fit = FALSE, reference = "zero") 
    df_test_cv$modcox_lp =   predict(modified_cox_cv, newdata = df_test_cv, 
                                     type = "lp", se.fit = FALSE, reference = "zero") 
    
    #check c-index 
    cindex_test[cv_iteration] = concordancefit(Surv(df_test_cv$time, df_test_cv$event), 
                                               -1*df_test_cv$modcox_lp)$concordance
    cindex_train[cv_iteration] = concordancefit(Surv(df_train_cv$time, df_train_cv$event), 
                                                -1*df_train_cv$modcox_lp)$concordance
  }
  return(c(mean(cindex_train, na.rm=1), mean(cindex_test, na.rm=1)))
}



method_3_train = function(df_train, predict.factors, fixed_time=NaN, 
                          inner_cv = 3, 
                          seed_to_fix = 100,
                          useCoxLasso = FALSE, var_importance = NULL, verbose = FALSE){
  p = length(predict.factors)
  n = dim(df_train)[1]
  
  if (sum(is.nan(fixed_time))>0){fixed_time = round(quantile(df_train[df_train$event==1, "time"], 0.85),1)}
  
  #1) choosing the VIMP variables
  
  if(is.null(var_importance)){
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
  }
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
    bestcindex[i] = method_3_tune(df_tune = df_train, predict.factors = params, params_for_tree = predict.factors,
                                  fixed_time=fixed_time,inner_cv = inner_cv, maxdepth = maxdepth, cp =0.001, 
                                  minbucket = minbucket, seed_to_fix = seed_to_fix, useCoxLasso = useCoxLasso)[2]
    
  }
  #print(grid_of_values[which.max(bestcindex),])
  if(sum(is.nan(bestcindex)) != dim(grid_of_values)[1]){
    
    maxdepth_use = grid_of_values[which.max(bestcindex),"tree_depth"]
    p_use = var_sorted[1:grid_of_values[which.max(bestcindex),"p_cv"]]
    minbucket_use = grid_of_values[which.max(bestcindex),"minbucket"]
  }else{#all failed, no split has been made, default to NuLL model
    
    maxdepth_use = 1
    p_use = var_sorted[1:grid_of_values[1,"p_cv"]]
    minbucket_use = grid_of_values[1,"minbucket"]
    
  }
  
  #3) GROWING A SINGLE RPART TREE, only use the top-vars for it
  #train tree
  set.seed(seed_to_fix)
  rpart.m = rpart::rpart(as.formula(paste("Surv(time, event) ~", paste(p_use, collapse="+"))),
                         data = df_train, minbucket = minbucket_use, maxdepth = maxdepth_use, cp = 0.001)
  clusters = unique(round(predict(rpart.m,df_train),6)); length(clusters) #number of clusters
  df_train$cluster_tree = as.factor(round(predict(rpart.m,df_train),6))
  #calculating Cox regression with clusters as a new factor
  modified_cox_model = method_cox_train(df_train,c(predict.factors, "cluster_tree"), useCoxLasso = useCoxLasso)
  
  output = list()
  output$treemodel = rpart.m
  output$modcoxmodel = modified_cox_model
  output$clusters = clusters
  return (output)
}

method_3_predict = function(model_3, df_test, fixed_time){
  
  #predicts probability of event by fixed_time for df_test from model_3 for ensemble method 3
  # calculate which cluster df_test falls into 
  df_test$cluster_tree = as.factor(round(predict(model_3$treemodel, newdata = df_test),6))
  #compute probabilities from the modified Cox model 
  predictedprob = method_cox_predict(model_3$modcoxmodel, df_test, fixed_time)
  if (sum(is.nan(predictedprob))>0){print ("some probs are missing")  }
  
  return(predictedprob)
}

method_3_cv = function(df, predict.factors, fixed_time = NaN, cv_number = 3, inner_cv =3, seed_to_fix = 100){
  time_0 = Sys.time()
  set.seed(seed_to_fix)
  
  # defining output for fixed_time
  if (sum(is.nan(fixed_time))>0){fixed_time = round(quantile(df[df$event==1, "time"], 0.85),1)}
  
  cv_folds = caret::createFolds(df$event, k=cv_number, list = FALSE) #use caret to split into k-folds = cv_steps
  modelstats_train = list(); modelstats_test = list()
  modcox_models_for_each_cv = list() 
  cv_iteration=1
  for (cv_iteration in 1:cv_number){
    print (paste('External loop CV, step number = ', cv_iteration))
    df_train_cv = df[cv_folds != cv_iteration, ]
    df_test_cv  = df[cv_folds == cv_iteration, ]
    model3_tuned = method_3_train(df_train_cv, predict.factors, fixed_time=fixed_time,
                                  inner_cv=inner_cv,seed_to_fix=seed_to_fix)
    
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

method_1_2_3_bs = function(trials = 2,n_train = 1000,n_test = 20000,
                           rho_w = 1, percentcensored = 0.75,distr = "Exp",
                           sim_type = 1, #1- linear, 2- nonl, 3-xt, 4- nonPH
                           predict.factors = c("age", "bmi","sex","hyp"), 
                           observe_time = 10, training_time = 8, valuation_time = c(4,7,9.9), 
                           rs=123){
  
  t_super_start = Sys.time()
  
  if (sim_type==1){
    simfun = function(x){
      simulatedata_linear(N=x[1], randomseed = x[2], rho_w = rho_w,
                          observe_time = observe_time, distr = distr,
                          percentcensored = percentcensored)}
  }else if (sim_type==2){
    simfun = function(x){
      simulatedata_nonlinear(N=x[1], randomseed = x[2], rho_w = rho_w,
                             observe_time = observe_time, distr = distr,
                             percentcensored = percentcensored)}
  }else if (sim_type==3){
    simfun = function(x){
      simulatedata_crossterms(N=x[1], randomseed = x[2], rho_w = rho_w,
                              observe_time = observe_time, distr = distr,
                              percentcensored = percentcensored)}
  }else if (sim_type ==4) {
    simfun = function(x){
      simulatedata_nonPH(N=x[1], randomseed = x[2],
                         observe_time = observe_time, 
                         percentcensored = percentcensored)}
  }else {print ("Unknown simulation type"); return (NULL)}
  
  nn = c("Cox.T","Cox.AUCROC","Cox.BS","Cox.BS_scaled","Cox.C_score", "Cox.Calib_slope","Cox.Calib_alpha",
         "CoxMFP.T","CoxMFP.AUCROC","CoxMFP.BS","CoxMFP.BS_scaled","CoxMFP.C_score", "CoxMFP.Calib_slope","CoxMFP.Calib_alpha",
         "Ens1.T","Ens1.AUCROC","Ens1.BS","Ens1.BS_scaled","Ens1.C_score", "Ens1.Calib_slope","Ens1.Calib_alpha",
         "Ens2.T","Ens2.AUCROC","Ens2.BS","Ens2.BS_scaled","Ens2.C_score", "Ens2.Calib_slope","Ens2.Calib_alpha",
         "Ens3.T","Ens3.AUCROC","Ens3.BS","Ens3.BS_scaled","Ens3.C_score","Ens3.Calib_slope","Ens3.Calib_alpha")
  #place holders for 3 times 
  results_val_t1 = data.frame(matrix(NA, nrow = trials, ncol = 35))
  names(results_val_t1) = nn
  results_val_t2 = results_val_t1
  results_val_t3 = results_val_t1
  
  results_app_t1 =  results_val_t1
  results_app_t2 =  results_val_t1
  results_app_t3 =  results_val_t1
  names(results_val_t1)
  times_iter = c()
  n_leaves= c() #vector with n of final leaves
  def_leaves = list() #list with conditions for final leaves
  
  #simulate 1 external data
  dtest = lapply(list(c(n_test, rs)), simfun)[[1]]
  for (i in 1:trials){
    t0 = Sys.time()
    
    print (paste("Trial", i, "/", trials))
    dtrain = lapply(list(c(n_train, rs+i)), simfun)[[1]]
    #round(summary(method_cox_train(dtrain, params))$coeff,4)
    
    try({
      m1 = method_cox_train(dtrain, predict.factors)
      #m2 = method_srf_train(df_train = dtrain,predict.factors = predict.factors,fixed_time = training_time,
      #                        cv_number = 3, seed_to_fix = rs+i+2, fast_version = TRUE, verbose = FALSE)
      m2 = method_coxmfp_train(df_train = dtrain,predict.factors = predict.factors,verbose = FALSE)
      
      m3 = method_1A_train(df_train = dtrain, predict.factors = predict.factors, fixed_time = training_time,
                           cv_number = 3,seed_to_fix = rs+3,fast_version = TRUE, var_importance_calc =0)
      
      m4 = method_2A_train(df_train = dtrain,predict.factors = predict.factors, fixed_time = training_time, 
                           inner_cv = 3,seed_to_fix = rs+4,verbose = FALSE,useCoxLasso = FALSE)
      
      m5 = method_3_train(dtrain,predict.factors,fixed_time = training_time,inner_cv = 3,
                          seed_to_fix = rs+5,useCoxLasso =  FALSE,var_importance = m4$vimp10)
      
      #test performance t1
      p1t1 = method_cox_predict(m1, dtest,times = valuation_time[1])
      p2t1 = method_coxmfp_predict(m2, dtest, fixed_time = valuation_time[1])
      p3t1 = method_1A_predict(m3, dtest,valuation_time[1] )
      p4t1 = method_2A_predict(m4, dtest, valuation_time[1])
      p5t1 = method_3_predict(m5, dtest, valuation_time[1])
      
      p1t2 = method_cox_predict(m1, dtest,times = valuation_time[2])
      p2t2 = method_coxmfp_predict(m2, dtest, fixed_time = valuation_time[2])
      p3t2 = method_1A_predict(m3, dtest,valuation_time[2] )
      p4t2 = method_2A_predict(m4, dtest, valuation_time[2])
      p5t2 = method_3_predict(m5, dtest, valuation_time[2])
      
      p1t3 = method_cox_predict(m1, dtest,times = valuation_time[3])
      p2t3 = method_coxmfp_predict(m2, dtest, fixed_time = valuation_time[3])
      p3t3 = method_1A_predict(m3, dtest,valuation_time[3] )
      p4t3 = method_2A_predict(m4, dtest, valuation_time[3])
      p5t3 = method_3_predict(m5, dtest, valuation_time[3])
      
      p1t3_train = method_cox_predict(m1, dtrain,times = valuation_time[3])
      p2t3_train = method_coxmfp_predict(m2, dtrain, fixed_time = valuation_time[3])
      p3t3_train = method_1A_predict(m3, dtrain,valuation_time[3] )
      p4t3_train = method_2A_predict(m4, dtrain, valuation_time[3])
      p5t3_train = method_3_predict(m5, dtrain, valuation_time[3])
      
      stats1t1 = method_any_validate(p1t1, valuation_time[1], dtrain, dtest)
      stats2t1 = method_any_validate(p2t1, valuation_time[1], dtrain, dtest)
      stats3t1 = method_any_validate(p3t1, valuation_time[1], dtrain, dtest)
      stats4t1 = method_any_validate(p4t1, valuation_time[1], dtrain, dtest)
      stats5t1 = method_any_validate(p5t1, valuation_time[1], dtrain, dtest)
      
      stats1t2 = method_any_validate(p1t2, valuation_time[2], dtrain, dtest)
      stats2t2 = method_any_validate(p2t2, valuation_time[2], dtrain, dtest)
      stats3t2 = method_any_validate(p3t2, valuation_time[2], dtrain, dtest)
      stats4t2 = method_any_validate(p4t2, valuation_time[2], dtrain, dtest)
      stats5t2 = method_any_validate(p5t2, valuation_time[2], dtrain, dtest)
      
      stats1t3 = method_any_validate(p1t3, valuation_time[3], dtrain, dtest)
      stats2t3 = method_any_validate(p2t3, valuation_time[3], dtrain, dtest)
      stats3t3 = method_any_validate(p3t3, valuation_time[3], dtrain, dtest)
      stats4t3 = method_any_validate(p4t3, valuation_time[3], dtrain, dtest)
      stats5t3 = method_any_validate(p5t3, valuation_time[3], dtrain, dtest)
      
      stats1t3_train = method_any_validate(p1t3_train, valuation_time[3], dtrain, dtrain)
      stats2t3_train = method_any_validate(p2t3_train, valuation_time[3], dtrain, dtrain)
      stats3t3_train = method_any_validate(p3t3_train, valuation_time[3], dtrain, dtrain)
      stats4t3_train = method_any_validate(p4t3_train, valuation_time[3], dtrain, dtrain)
      stats5t3_train = method_any_validate(p5t3_train, valuation_time[3], dtrain, dtrain) 
      
      results_val_t1[i,]= cbind("Cox" = stats1t1, "CoxMFP" = stats2t1, "Ens1" = stats3t1,
                                "Ens2"=stats4t1, "Ens3"=stats5t1)
      results_val_t2[i,]= cbind("Cox" = stats1t2, "CoxMFP" = stats2t2, "Ens1" = stats3t2,
                                "Ens2"=stats4t2, "Ens3"=stats5t2)
      results_val_t3[i,]= cbind("Cox" = stats1t3, "CoxMFP" = stats2t3, "Ens1" = stats3t3,
                                "Ens2"=stats4t3, "Ens3"=stats5t3)
      results_app_t3[i,]= cbind("Cox" = stats1t3_train, "CoxMFP" = stats2t3_train, "Ens1" = stats3t3_train,
                                "Ens2"=stats4t3_train, "Ens3"=stats5t3_train)
      
    }, silent=TRUE)
    t1 = Sys.time()
    times_iter[i]= t1- t0
  }
  
  mean_na = function(x){mean(x, na.rm=1)}
  sd_na = function(x){SD(x, na.rm=1)}
  validated_stats_t1 = as.data.frame(cbind("mean" = apply(results_val_t1, 2, mean_na),
                                           "SD"= apply(results_val_t1, 2, sd_na)))
  validated_stats_t2 = as.data.frame(cbind("mean" = apply(results_val_t2, 2, mean_na),
                                           "SD"= apply(results_val_t2, 2, sd_na)))
  validated_stats_t3 = as.data.frame(cbind("mean" = apply(results_val_t3, 2, mean_na),
                                           "SD"= apply(results_val_t3, 2, sd_na)))
  
  apparent_stats = as.data.frame(cbind("mean" = apply(results_app_t3, 2, mean_na), 
                                       "SD"= apply(results_app_t3, 2, sd_na)))
  
  together = as.data.frame(cbind("train" = apparent_stats , "test_t1"= validated_stats_t1,
                                 "test_t2"= validated_stats_t2,"test_t3"= validated_stats_t3))
  
  together["compute_per_i_sec", ] = c(mean(times_iter, na.rm=1),sd(times_iter, na.rm=1),
                                      mean(times_iter, na.rm=1), sd(times_iter, na.rm=1),
                                      mean(times_iter, na.rm=1), sd(times_iter, na.rm=1),
                                      mean(times_iter, na.rm=1), sd(times_iter, na.rm=1))
  
  output = list()
  output$mean_stats = together
  output$train_stats = results_app_t3
  output$test_stats =  as.data.frame(rbind(results_val_t1, results_val_t2, results_val_t3))
  output$n_final_leaves =  n_leaves
  output$def_final_leaves =  def_leaves
  output$total_time =  Sys.time() - t_super_start
  return(output)
}

#######################################################################

############################### Bootstrap validation for one data ,method_1_2_3_bs_real ####################################
method_1_2_3_bs_real = function(trials = 2, mydata,  predict.factors, 
                                training_time, valuation_time=c(6,7,8),rs=123){
  
  t_super_start = Sys.time()
  
  nn = c("Cox.T","Cox.AUCROC","Cox.BS","Cox.BS_scaled","Cox.C_score", "Cox.Calib_slope","Cox.Calib_alpha",
         "CoxMFP.T","CoxMFP.AUCROC","CoxMFP.BS","CoxMFP.BS_scaled","CoxMFP.C_score", "CoxMFP.Calib_slope","CoxMFP.Calib_alpha",
         "Ens1.T","Ens1.AUCROC","Ens1.BS","Ens1.BS_scaled","Ens1.C_score", "Ens1.Calib_slope","Ens1.Calib_alpha",
         "Ens2.T","Ens2.AUCROC","Ens2.BS","Ens2.BS_scaled","Ens2.C_score", "Ens2.Calib_slope","Ens2.Calib_alpha",
         "Ens3.T","Ens3.AUCROC","Ens3.BS","Ens3.BS_scaled","Ens3.C_000000score","Ens3.Calib_slope","Ens3.Calib_alpha")
  #place holders for 3 times 
  
  results_val_t1 = data.frame(matrix(NA, nrow = trials, ncol = 35))
  names(results_val_t1) = nn
  results_val_t2 = results_val_t1
  results_val_t3 = results_val_t1
  
  results_app_t1 =  results_val_t1
  results_app_t2 =  results_val_t1
  results_app_t3 =  results_val_t1
  names(results_val_t1)
  times_iter = c()
  n_leaves= c() #vector with n of final leaves
  def_leaves = list() #list with conditions for final leaves
  
  #simulate 1 external data
  dtest = mydata
  for (i in 1:trials){
    t0 = Sys.time()
    print (paste("Bootstrap i", i, "/", trials))
    
    #sampled version of the data
    {set.seed(rs+i*37); s=sample(seq(dim(mydata)[1]),replace=TRUE)}
    dtrain = mydata[s, ]
    
    try({
      m1 = method_cox_train(dtrain, predict.factors)
      #m2 = method_srf_train(df_train = dtrain,predict.factors = predict.factors,fixed_time = training_time,
      #                        cv_number = 3, seed_to_fix = rs+i+2, fast_version = TRUE, verbose = FALSE)
      m2 = method_coxmfp_train(df_train = dtrain,predict.factors = predict.factors,verbose = FALSE)
      
      m3 = method_1A_train(df_train = dtrain, predict.factors = predict.factors, fixed_time = training_time,
                           cv_number = 3,seed_to_fix = rs+3,fast_version = TRUE, var_importance_calc =0)
      
      m4 = method_2A_train(df_train = dtrain,predict.factors = predict.factors, fixed_time = training_time, 
                           inner_cv = 3,seed_to_fix = rs+4,verbose = FALSE,useCoxLasso = FALSE)
      
      m5 = method_3_train(dtrain,predict.factors,fixed_time = training_time,inner_cv = 3,
                          seed_to_fix = rs+5,useCoxLasso =  FALSE,var_importance = m4$vimp10)
      
      #test performance t1
      p1t1 = method_cox_predict(m1, dtest,times = valuation_time[1])
      p2t1 = method_coxmfp_predict(m2, dtest, fixed_time = valuation_time[1])
      p3t1 = method_1A_predict(m3, dtest,valuation_time[1] )
      p4t1 = method_2A_predict(m4, dtest, valuation_time[1])
      p5t1 = method_3_predict(m5, dtest, valuation_time[1])
      
      p1t2 = method_cox_predict(m1, dtest,times = valuation_time[2])
      p2t2 = method_coxmfp_predict(m2, dtest, fixed_time = valuation_time[2])
      p3t2 = method_1A_predict(m3, dtest,valuation_time[2] )
      p4t2 = method_2A_predict(m4, dtest, valuation_time[2])
      p5t2 = method_3_predict(m5, dtest, valuation_time[2])
      
      p1t3 = method_cox_predict(m1, dtest,times = valuation_time[3])
      p2t3 = method_coxmfp_predict(m2, dtest, fixed_time = valuation_time[3])
      p3t3 = method_1A_predict(m3, dtest,valuation_time[3] )
      p4t3 = method_2A_predict(m4, dtest, valuation_time[3])
      p5t3 = method_3_predict(m5, dtest, valuation_time[3])
      
      p1t3_train = method_cox_predict(m1, dtrain,times = valuation_time[3])
      p2t3_train = method_coxmfp_predict(m2, dtrain, fixed_time = valuation_time[3])
      p3t3_train = method_1A_predict(m3, dtrain,valuation_time[3] )
      p4t3_train = method_2A_predict(m4, dtrain, valuation_time[3])
      p5t3_train = method_3_predict(m5, dtrain, valuation_time[3])
      
      stats1t1 = method_any_validate(p1t1, valuation_time[1], dtrain, dtest)
      stats2t1 = method_any_validate(p2t1, valuation_time[1], dtrain, dtest)
      stats3t1 = method_any_validate(p3t1, valuation_time[1], dtrain, dtest)
      stats4t1 = method_any_validate(p4t1, valuation_time[1], dtrain, dtest)
      stats5t1 = method_any_validate(p5t1, valuation_time[1], dtrain, dtest)
      
      stats1t2 = method_any_validate(p1t2, valuation_time[2], dtrain, dtest)
      stats2t2 = method_any_validate(p2t2, valuation_time[2], dtrain, dtest)
      stats3t2 = method_any_validate(p3t2, valuation_time[2], dtrain, dtest)
      stats4t2 = method_any_validate(p4t2, valuation_time[2], dtrain, dtest)
      stats5t2 = method_any_validate(p5t2, valuation_time[2], dtrain, dtest)
      
      stats1t3 = method_any_validate(p1t3, valuation_time[3], dtrain, dtest)
      stats2t3 = method_any_validate(p2t3, valuation_time[3], dtrain, dtest)
      stats3t3 = method_any_validate(p3t3, valuation_time[3], dtrain, dtest)
      stats4t3 = method_any_validate(p4t3, valuation_time[3], dtrain, dtest)
      stats5t3 = method_any_validate(p5t3, valuation_time[3], dtrain, dtest)
      
      stats1t1_train = method_any_validate(p1t3_train, valuation_time[1], dtrain, dtrain)
      stats2t1_train = method_any_validate(p2t3_train, valuation_time[1], dtrain, dtrain)
      stats3t1_train = method_any_validate(p3t3_train, valuation_time[1], dtrain, dtrain)
      stats4t1_train = method_any_validate(p4t3_train, valuation_time[1], dtrain, dtrain)
      stats5t1_train = method_any_validate(p5t3_train, valuation_time[1], dtrain, dtrain) 
      
      stats1t2_train = method_any_validate(p1t3_train, valuation_time[2], dtrain, dtrain)
      stats2t2_train = method_any_validate(p2t3_train, valuation_time[2], dtrain, dtrain)
      stats3t2_train = method_any_validate(p3t3_train, valuation_time[2], dtrain, dtrain)
      stats4t2_train = method_any_validate(p4t3_train, valuation_time[2], dtrain, dtrain)
      stats5t2_train = method_any_validate(p5t3_train, valuation_time[2], dtrain, dtrain) 
      
      stats1t3_train = method_any_validate(p1t3_train, valuation_time[3], dtrain, dtrain)
      stats2t3_train = method_any_validate(p2t3_train, valuation_time[3], dtrain, dtrain)
      stats3t3_train = method_any_validate(p3t3_train, valuation_time[3], dtrain, dtrain)
      stats4t3_train = method_any_validate(p4t3_train, valuation_time[3], dtrain, dtrain)
      stats5t3_train = method_any_validate(p5t3_train, valuation_time[3], dtrain, dtrain) 
      
      results_val_t1[i,]= cbind("Cox" = stats1t1, "CoxMFP" = stats2t1, "Ens1" = stats3t1,
                                "Ens2"=stats4t1, "Ens3"=stats5t1)
      results_val_t2[i,]= cbind("Cox" = stats1t2, "CoxMFP" = stats2t2, "Ens1" = stats3t2,
                                "Ens2"=stats4t2, "Ens3"=stats5t2)
      results_val_t3[i,]= cbind("Cox" = stats1t3, "CoxMFP" = stats2t3, "Ens1" = stats3t3,
                                "Ens2"=stats4t3, "Ens3"=stats5t3)
      results_app_t1[i,]= cbind("Cox" = stats1t1_train, "CoxMFP" = stats2t1_train, "Ens1" = stats3t1_train,
                                "Ens2"=stats4t1_train, "Ens3"=stats5t1_train)
      results_app_t2[i,]= cbind("Cox" = stats1t2_train, "CoxMFP" = stats2t2_train, "Ens1" = stats3t2_train,
                                "Ens2"=stats4t2_train, "Ens3"=stats5t2_train)
      results_app_t3[i,]= cbind("Cox" = stats1t3_train, "CoxMFP" = stats2t3_train, "Ens1" = stats3t3_train,
                                "Ens2"=stats4t3_train, "Ens3"=stats5t3_train)
      
    }, silent=TRUE)
    t1 = Sys.time()
    times_iter[i]= t1- t0
  }
  
  mean_na = function(x){mean(x, na.rm=1)}
  sd_na = function(x){SD(x, na.rm=1)}
  validated_stats_t1 = as.data.frame(cbind("mean" = apply(results_val_t1, 2, mean_na),
                                           "SD"= apply(results_val_t1, 2, sd_na)))
  validated_stats_t2 = as.data.frame(cbind("mean" = apply(results_val_t2, 2, mean_na),
                                           "SD"= apply(results_val_t2, 2, sd_na)))
  validated_stats_t3 = as.data.frame(cbind("mean" = apply(results_val_t3, 2, mean_na),
                                           "SD"= apply(results_val_t3, 2, sd_na)))
  
  apparent_stats_t1 = as.data.frame(cbind("mean" = apply(results_app_t1, 2, mean_na), 
                                       "SD"= apply(results_app_t1, 2, sd_na)))
  apparent_stats_t2 = as.data.frame(cbind("mean" = apply(results_app_t2, 2, mean_na), 
                                          "SD"= apply(results_app_t2, 2, sd_na)))
  apparent_stats_t3 = as.data.frame(cbind("mean" = apply(results_app_t3, 2, mean_na), 
                                          "SD"= apply(results_app_t3, 2, sd_na)))
  
  together = as.data.frame(cbind("train_t1" = apparent_stats_t1,
                                 "train_t2" = apparent_stats_t2,
                                 "train_t3" = apparent_stats_t3,
                                 "test_t1"= validated_stats_t1,
                                 "test_t2"= validated_stats_t2,
                                 "test_t3"= validated_stats_t3))
  
  together["compute_per_i_sec", ] = c(mean(times_iter, na.rm=1),sd(times_iter, na.rm=1),
                                      mean(times_iter, na.rm=1), sd(times_iter, na.rm=1),
                                      mean(times_iter, na.rm=1), sd(times_iter, na.rm=1),
                                      mean(times_iter, na.rm=1), sd(times_iter, na.rm=1),
                                      mean(times_iter, na.rm=1), sd(times_iter, na.rm=1),
                                      mean(times_iter, na.rm=1), sd(times_iter, na.rm=1))
  
  output = list()
  output$mean_stats = together
  output$train_stats = as.data.frame(rbind(results_app_t1, results_app_t2, results_app_t3))
  output$test_stats =  as.data.frame(rbind(results_val_t1, results_val_t2, results_val_t3))
  output$n_final_leaves =  n_leaves
  output$def_final_leaves =  def_leaves
  output$total_time =  Sys.time() - t_super_start
  return(output)
}


method_1_2_3_apparent = function(mydata,  predict.factors, 
                                training_time, valuation_time=c(6,7,8),rs=123){
  
  t_super_start = Sys.time()
  
  nn = c("Cox.T","Cox.AUCROC","Cox.BS","Cox.BS_scaled","Cox.C_score", "Cox.Calib_slope","Cox.Calib_alpha",
         "CoxMFP.T","CoxMFP.AUCROC","CoxMFP.BS","CoxMFP.BS_scaled","CoxMFP.C_score", "CoxMFP.Calib_slope","CoxMFP.Calib_alpha",
         "Ens1.T","Ens1.AUCROC","Ens1.BS","Ens1.BS_scaled","Ens1.C_score", "Ens1.Calib_slope","Ens1.Calib_alpha",
         "Ens2.T","Ens2.AUCROC","Ens2.BS","Ens2.BS_scaled","Ens2.C_score", "Ens2.Calib_slope","Ens2.Calib_alpha",
         "Ens3.T","Ens3.AUCROC","Ens3.BS","Ens3.BS_scaled","Ens3.C_score","Ens3.Calib_slope","Ens3.Calib_alpha")
  #place holders for 3 times 
  
  results_val_t1 = data.frame(matrix(NA, nrow = 1, ncol = 35))
  names(results_val_t1) = nn
  results_val_t2 = results_val_t1
  results_val_t3 = results_val_t1
  
  results_app_t1 =  results_val_t1
  results_app_t2 =  results_val_t1
  results_app_t3 =  results_val_t1
  names(results_val_t1)
  times_iter = c()
  n_leaves= c() #vector with n of final leaves
  def_leaves = list() #list with conditions for final leaves
  
  #simulate 1 external data
  dtest = mydata
  dtrain = mydata
  
  try({
      m1 = method_cox_train(dtrain, predict.factors)
      #m2 = method_srf_train(df_train = dtrain,predict.factors = predict.factors,fixed_time = training_time,
      #                        cv_number = 3, seed_to_fix = rs+i+2, fast_version = TRUE, verbose = FALSE)
      m2 = method_coxmfp_train(df_train = dtrain,predict.factors = predict.factors,verbose = FALSE)
      
      m3 = method_1A_train(df_train = dtrain, predict.factors = predict.factors, fixed_time = training_time,
                           cv_number = 3,seed_to_fix = rs+3,fast_version = TRUE, var_importance_calc =0)
      
      m4 = method_2A_train(df_train = dtrain,predict.factors = predict.factors, fixed_time = training_time, 
                           inner_cv = 3,seed_to_fix = rs+4,verbose = FALSE,useCoxLasso = FALSE)
      
      m5 = method_3_train(dtrain,predict.factors,fixed_time = training_time,inner_cv = 3,
                          seed_to_fix = rs+5,useCoxLasso =  FALSE,var_importance = m4$vimp10)
      
      #test performance t1
      p1t1 = method_cox_predict(m1, dtest,times = valuation_time[1])
      p2t1 = method_coxmfp_predict(m2, dtest, fixed_time = valuation_time[1])
      p3t1 = method_1A_predict(m3, dtest,valuation_time[1] )
      p4t1 = method_2A_predict(m4, dtest, valuation_time[1])
      p5t1 = method_3_predict(m5, dtest, valuation_time[1])
      
      p1t2 = method_cox_predict(m1, dtest,times = valuation_time[2])
      p2t2 = method_coxmfp_predict(m2, dtest, fixed_time = valuation_time[2])
      p3t2 = method_1A_predict(m3, dtest,valuation_time[2] )
      p4t2 = method_2A_predict(m4, dtest, valuation_time[2])
      p5t2 = method_3_predict(m5, dtest, valuation_time[2])
      
      p1t3 = method_cox_predict(m1, dtest,times = valuation_time[3])
      p2t3 = method_coxmfp_predict(m2, dtest, fixed_time = valuation_time[3])
      p3t3 = method_1A_predict(m3, dtest,valuation_time[3] )
      p4t3 = method_2A_predict(m4, dtest, valuation_time[3])
      p5t3 = method_3_predict(m5, dtest, valuation_time[3])
      
      p1t3_train = method_cox_predict(m1, dtrain,times = valuation_time[3])
      p2t3_train = method_coxmfp_predict(m2, dtrain, fixed_time = valuation_time[3])
      p3t3_train = method_1A_predict(m3, dtrain,valuation_time[3] )
      p4t3_train = method_2A_predict(m4, dtrain, valuation_time[3])
      p5t3_train = method_3_predict(m5, dtrain, valuation_time[3])
      
      stats1t1 = method_any_validate(p1t1, valuation_time[1], dtrain, dtest)
      stats2t1 = method_any_validate(p2t1, valuation_time[1], dtrain, dtest)
      stats3t1 = method_any_validate(p3t1, valuation_time[1], dtrain, dtest)
      stats4t1 = method_any_validate(p4t1, valuation_time[1], dtrain, dtest)
      stats5t1 = method_any_validate(p5t1, valuation_time[1], dtrain, dtest)
      
      stats1t2 = method_any_validate(p1t2, valuation_time[2], dtrain, dtest)
      stats2t2 = method_any_validate(p2t2, valuation_time[2], dtrain, dtest)
      stats3t2 = method_any_validate(p3t2, valuation_time[2], dtrain, dtest)
      stats4t2 = method_any_validate(p4t2, valuation_time[2], dtrain, dtest)
      stats5t2 = method_any_validate(p5t2, valuation_time[2], dtrain, dtest)
      
      stats1t3 = method_any_validate(p1t3, valuation_time[3], dtrain, dtest)
      stats2t3 = method_any_validate(p2t3, valuation_time[3], dtrain, dtest)
      stats3t3 = method_any_validate(p3t3, valuation_time[3], dtrain, dtest)
      stats4t3 = method_any_validate(p4t3, valuation_time[3], dtrain, dtest)
      stats5t3 = method_any_validate(p5t3, valuation_time[3], dtrain, dtest)
      

      results_val_t1= cbind("Cox" = stats1t1, "CoxMFP" = stats2t1, "Ens1" = stats3t1,
                                "Ens2"=stats4t1, "Ens3"=stats5t1)
      results_val_t2= cbind("Cox" = stats1t2, "CoxMFP" = stats2t2, "Ens1" = stats3t2,
                                "Ens2"=stats4t2, "Ens3"=stats5t2)
      results_val_t3= cbind("Cox" = stats1t3, "CoxMFP" = stats2t3, "Ens1" = stats3t3,
                                "Ens2"=stats4t3, "Ens3"=stats5t3)
      
    }, silent=TRUE)

  together = as.data.frame(rbind("Train_t1"= results_val_t1,
                                 "Train_t2"= results_val_t2,
                                 "Train_t3"= results_val_t3))
  output = list()
  output$mean_stats = t(together)
  output$total_time =  Sys.time() - t_super_start
  return(output)
}

