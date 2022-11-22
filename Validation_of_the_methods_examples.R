#'
#'  Examples of how to apply the functions from "EnsembleMethods_SeparateCodesByMethods.R"
#'  
#'  using simulated data sets from the "Simulated_data.R" file in the same folder at github
#'  
#' To apply these functions to another sample, make sure that survival outcome is 
#'             defined using "event" and "time" variables 
#'             

#'#############################################################################
##### Example of validating the methods for non-linear simulated data set ##### 
#'#############################################################################

# generate data
sim_data_nonlinear = simulatedata_nonlinear(500, 10, 0, 100)

# we will use the predictors in the simulated data, for another sample please change accordingly
predict.factors = c("age", "bmi", "hyp", "gender")

# use populationstats() to display statistics for the sample 
populationstats(sim_data_nonlinear, time_f = 5)

###### Cross-validation of the methods for prediction performance:######

# Cox model
cox = method_cox_cv(sim_data_nonlinear, predict.factors, fixed_time = 5, cv_number = 5, seed_to_fix = 10)

# Cox + polynomials
coxmfp = method_coxmfp_cv(sim_data_nonlinear, predict.factors, fixed_time = 5, cv_number = 5, seed_to_fix = 100)

# Survival Random Forest 
srf = method_srf_cv(sim_data_nonlinear, predict.factors, fixed_time = 5, cv_number = 3, internal_cv_k = 3,seed_to_fix = 100)

#Ensemble method 1a "Cox predictors to -> SRF"
ens1a = method_1A_cv(sim_data_nonlinear, predict.factors, fixed_time = 5, cv_number = 3, internal_cv_k = 3,seed_to_fix = 100)

#Ensemble method 2a RPart Tree -> Cox in Clusters 
ens2a = method_2A_cv(sim_data_nonlinear, predict.factors, fixed_time = 5, cv_number = 3, internal_cv_k = 3,seed_to_fix = 100)

#Ensemble method 3 RPart TREE -> Modified Cox with clusters IDs 
ens3 = method_3_cv(sim_data_nonlinear, predict.factors, fixed_time = 5, cv_number = 3, seed_to_fix = 100)

#test results - mean and std across CV steps
results_on_test_sets = cbind(  
                    "Cox" = cox$testaverage,   
                     "SRF" = srf$testaverage,   
                     "CoxMFP" = coxmfp$testaverage,
                     "1A" = ens1a$testaverage,  
                     "2A" = ens2a$testaverage,   
                      "3" = ens3$testaverage)
results_on_test_sets

std_on_test_sets =   cbind(  
                     "Cox" = apply(cox$test,  2, sd),
                     "SRF" = apply(srf$test,  2, sd),  
                     "CoxMFP" = apply(coxmfp$test,  2, sd),
                     "1A" = apply(ens1a$test,  2, sd),  
                     "2A" =apply(ens2a$test,  2, sd),   
                     "3" = apply(ens3$test,  2, sd))
std_on_test_sets

###### Train -> predict [-> externally validate] ######

#simulate train data
sim_data_nonlinear = simulatedata_nonlinear(N=1000, observe_time =10, randomseed = 100)

#simulate validation set (to demonstrate external validation)
sim_data_nonlinear_exttest = simulatedata_nonlinear(N=20000, observe_time =10, randomseed = 100)

#samples stats
populationstats(sim_data_nonlinear, 5)
populationstats(sim_data_nonlinear_exttest, 5)

# Train-predict-validate for Cox  
mcox = method_cox_train(sim_data_nonlinear, predict.factors)
  #predict on train set 
predictcoxext = method_cox_predict(mcox, sim_data_nonlinear_exttest, 5)
  # predict survival probability for the validation set 
predictcox = method_cox_predict(mcox, sim_data_nonlinear, 5)
  # compute predictive performance for the train set 
valcox = method_any_validate(predictcox, 5, sim_data_nonlinear, sim_data_nonlinear)
  # compute predictive performance for the validation set 
valcoxext = method_any_validate(predictcoxext, 5, sim_data_nonlinear, sim_data_nonlinear_exttest)

# Train-predict-validate for Cox mfp
mcoxmfp = method_coxmfp_train(sim_data_nonlinear, predict.factors)
predictcoxmfpext = method_cox_predict(mcoxmfp, sim_data_nonlinear_exttest, 5)
predictcoxmfp = method_cox_predict(mcoxmfp, sim_data_nonlinear, 5)
valmfp = method_any_validate(predictcoxmfp, 5, sim_data_nonlinear, sim_data_nonlinear)
valmfpext = method_any_validate(predictcoxmfpext, 5, sim_data_nonlinear, sim_data_nonlinear_exttest)

# Train-predict-validate for Survival Random Forest
msrf = method_srf_train(sim_data_nonlinear, predict.factors, fixed_time = 5,  cv_number = 3, seed_to_fix = 100)
predictsrfext = method_srf_predict(msrf$model, sim_data_nonlinear_exttest, fixed_time=8)
predictsrf = method_srf_predict(msrf$model, sim_data_nonlinear, fixed_time=8)
valsrf = method_any_validate(predictsrf, 5, sim_data_nonlinear, sim_data_nonlinear)
valsrfext = method_any_validate(predictsrfext, 5, sim_data_nonlinear, sim_data_nonlinear_exttest)

# Train-predict-validate for Ensemble 1A
m1a = method_1A_train(sim_data_nonlinear, predict.factors, cv_number = 3, fixed_time = 5, seed_to_fix = 100)

predict1bext = method_1A_predict(m1a, sim_data_nonlinear_exttest, 5)
predict1b = method_1A_predict(m1a, sim_data_nonlinear, 5)
val1b =method_any_validate(predict1b, 5, sim_data_nonlinear, sim_data_nonlinear)
val1bext =method_any_validate(predict1bext, 5, sim_data_nonlinear, sim_data_nonlinear_exttest)

# Train-predict-validate for Ensemble 2A
m2a = method_2A_train(sim_data_nonlinear, predict.factors, fixed_time = 5, 
                      internal_cv_k =3, seed_to_fix = 100)

survConcordance(Surv(sim_data_nonlinear$time, sim_data_nonlinear$event)~predict2a)$concordance
concordancefit(Surv(sim_data_nonlinear$time, sim_data_nonlinear$event), -1*predict2a)$concordance

predict2aext = method_2A_predict(m2a, sim_data_nonlinear_exttest, fixed_time=5)
predict2a = method_2A_predict(m2a, sim_data_nonlinear, fixed_time=5)
val2a = method_any_validate(predict2a, 5, sim_data_nonlinear, sim_data_nonlinear)
val2aextf =method_any_validate(predict2aext, 5, sim_data_nonlinear, sim_data_nonlinear_exttest)
  
  #Additionally: display single tree
rpart.plot(m2a$treemodel, nn=TRUE, roundint = FALSE, digits = 4)
  
  #display Cox models in each cluster 
  m2a$clusters #clusters correspond to the top values in cluster plot 
  m2a$coxmodels


# Train-predict-validate for Ensemble 3
m3 = method_3_train(sim_data_nonlinear, predict.factors, fixed_time = 5,  
                    internal_cv_k = 3, seed_to_fix = 100)
predict3ext = method_3_predict(m3, sim_data_nonlinear_exttest, fixed_time=8)
predict3 = method_3_predict(m3, sim_data_nonlinear, fixed_time=8)
val3 =method_any_validate(predict3, 5, sim_data_nonlinear, sim_data_nonlinear)
val3ext  =method_any_validate(predict3ext, 5, sim_data_nonlinear, sim_data_nonlinear_exttest)
  #additionally:
  #display Cox coefficients for the modified cox model 
round(summary(m3$modcoxmodel)$coefficients,4)
  # display the single tree with clusters 
rpart.plot(m3$treemodel, nn=TRUE, roundint = FALSE, digits = 4)

# display  results for all models 
rbind(valcoxext,valmfpext,valsrfext, val1bext, val2aextf, val3ext)
rbind(valcox,valmfp,valsrf, val1b, val2a, val3)
#library(clipr) Copy results to clip 
clipr::write_clip(rbind(valcoxext,valmfpext,valsrfext, val1bext, val2aextf, val3ext))
clipr::write_clip(rbind(valcox,valmfp,valsrf, val1b, val2a, val3))

#########################  hnscc_file  ##################################

hnscc_file = "~/Desktop/Study_KCL/PhD Projects/Ensemblemethods/LassoNet_Exercise/hnscc_merged.csv"
ddd = read.csv(hnscc_file)
predict.factors.hnscc = names(ddd)[4:26]
fixed_time_hns = 4
m2a = method_2A_train(ddd, predict.factors.hnscc, 
                      fixed_time = fixed_time_hns, 
                      internal_cv_k = 4, seed_to_fix = 50)

mcox = method_cox_train(ddd, predict.factors.hnscc)
mcox$coefficients

predictcox = method_cox_predict(mcox, ddd, fixed_time_hns)
concordancefit(Surv(df_test_cv$time, df_test_cv$event), -1*df_test_cv$eventprob)$concordance
valcox = method_any_validate(predictcox, fixed_time_hns, sim_data_nonlinear, sim_data_nonlinear)


