
#'#####################################################################
# Here we run and compute performance metrics from different methods # 
#'#####################################################################

#'############################################################################
########################### Non-Linear ################################
#'############################################################################

predict.factors = c("age", "bmi", "hyp", "gender")
sim_data_nonlinear = simulatedata_nonlinear(5000, 10, 0, 100)
#sim_data_nonlinear_exttest = simulatedata_nonlinear(10000, 10, 0, 101)

sim1  = simulatedata_crossterms(2000, 10, 0, 100)
populationstats(sim1, time_f = 5)
populationstats(sim_data_nonlinear_exttest, time_f = 5)
m1 = coxmfp$tuned_cv_models[[1]]
m1$means

#Cox
cox = method_cox_cv(sim_data_nonlinear, predict.factors, fixed_time = 5, cv_number = 5)
#Coxmfp
coxmfp = method_coxmfp_cv(sim_data_nonlinear, predict.factors, fixed_time = 5, cv_number = 5, seed_to_fix = 100)
#SRF
srf = method_srf_cv(sim_data_nonlinear, predict.factors, fixed_time = 5, cv_number = 3, seed_to_fix = 100)
#1a Cox -> SRF
ens1a = method_1A_cv(sim_data_nonlinear, predict.factors, fixed_time = 5, cv_number = 3, seed_to_fix = 100)
#1b SRF -> Cox 
ens1b = method_1B_cv(sim_data_nonlinear, predict.factors, fixed_time = 5, cv_number = 3, seed_to_fix = 100, 
                     pretrained_srf_models = srf$pretrained_srf_models)
#2a RPart Tree -> Cox in Clusters 
ens2a = method_2A_cv(sim_data_nonlinear, predict.factors, fixed_time = 5, cv_number = 3, seed_to_fix = 100)
#2b SRF Tree -> Cox in Clusters 
ens2b = method_2B_cv(sim_data_nonlinear, predict.factors, fixed_time = 5, cv_number = 3, seed_to_fix = 100)
#3 RPart TREE -> Modified Cox with clusters IDs 
ens3 = method_3_cv(sim_data_nonlinear, predict.factors, fixed_time = 5, cv_number = 3, seed_to_fix = 100)
#test results 

cbind(cox$testaverage, ens1a$testaverage)

res_test =   cbind(  "Cox" = cox$testaverage,   "SRF" = srf$testaverage,   "CoxMFP" = coxmfp$testaverage,
                     "1A" = ens1a$testaverage,  "1B" = ens1b$testaverage,  
                     "2A" = ens2a$testaverage,   "2B" = ens2b$testaverage,  "3" = ens3$testaverage)

std_test =   cbind(  "Cox" = apply(cox$test,  2, sd),   "SRF" = apply(srf$test,  2, sd),  
                     "CoxMFP" = apply(coxmfp$test,  2, sd),
                     "1A" = apply(ens1a$test,  2, sd),  "1B" = apply(ens1b$test,  2, sd),  
                     "2A" =apply(ens2a$test,  2, sd),   "2B" = apply(ens2b$test,  2, sd),  "3" = apply(ens3$test,  2, sd))

clipr::write_clip(res_test)
clipr::write_clip(std_test)
apply(cox$test,  2, sd)

#train results 
res_train = t(cbind(  "Cox" = cox$trainaverage,  "SRF" = srf$trainaverage,  "1A" = ens1a$trainaverage,
                      "1B" = ens1b$trainaverage,  "2A" = ens2a$trainaverage,  "2B" = ens2b$trainaverage,  "3" = ens3$trainaverage))
res_train


m2a = method_2a_train(sim_data_nonlinear, predict.factors, fixed_time = 8, seed_to_fix = 100)
#plot Rpart tree
rpart.plot(m2a$tree, digits =2)
#plot SRF tree 

############ Theoretical best linear #############

#theoretical best linear
cc = coxph(formula = as.formula(paste("Surv(sim_data$time, sim_data$event) ~",paste(c("age_", "bmi_", "hyp", "gender"), collapse = "+"))), x=TRUE, data = sim_data)
cc$coefficients = c(0.4,1,0.7,0) # change coefficients to the true ones
sim_data$eventprob = 1- pec::predictSurvProb(cc, sim_data, 8)
timeROC(T=sim_data$time,delta=sim_data$event,  marker=sim_data$eventprob,times=8*0.999,cause=1, weighting = "marginal")$AUC[2] 
#0.8579

#theoretical best non-linear
cc = coxph(formula = as.formula(paste("Surv(sim_data$time, sim_data$event) ~",
                                      paste(c("age_", "bmi_", "hyp", "gender"), collapse = "+"))), x=TRUE, data = sim_data)
cc$coefficients = c(0.4,1,0.7,0) # change coefficients to the true ones
sim_data$eventprob = 1- pec::predictSurvProb(cc, sim_data, 8)
timeROC(T=sim_data$time,delta=sim_data$event,  marker=sim_data$eventprob,times=8*0.999,cause=1, weighting = "marginal")$AUC[2] 
#0.8579

#theoretical best x-term
cc = coxph(formula = as.formula(paste("Surv(sim_data$time, sim_data$event) ~",
                                      paste(c("age_", "bmi_", "hyp", "gender"), collapse = "+"))), x=TRUE, data = sim_data)
cc$coefficients = c(0.4,1,0.7,0) # change coefficients to the true ones
sim_data$eventprob = 1- pec::predictSurvProb(cc, sim_data, 8)
timeROC(T=sim_data$time,delta=sim_data$event,  marker=sim_data$eventprob,times=8*0.999,cause=1, weighting = "marginal")$AUC[2] 
#0.8579


#'############################################################################
################################ CROSS TERMS #################################
#'############################################################################

predict.factors = c("age", "bmi", "hyp", "gender")
sim_data_xt = simulatedata_crossterms(5000, 10, 0, 100) #5000
sim_data_xt_exttest = simulatedata_crossterms(20000, 10, 0, 100)

populationstats(sim_data_xt, time_f = 8)
populationstats(sim_data_xt, time_f = 5)

#Cox
coxxt = method_cox_cv(sim_data_xt, predict.factors, fixed_time = 5, cv_number = 5, seed_to_fix = 100)
#CoxMFP
coxxtmfp = method_coxmfp_cv(sim_data_xt, predict.factors, fixed_time = 5, cv_number = 5)
#SRF
srfxt = method_srf_cv(sim_data_xt, predict.factors, fixed_time = 5, cv_number = 5, seed_to_fix = 100)
#1a Cox -> SRF
ens1axt = method_1A_cv(sim_data_xt, predict.factors, fixed_time = 5, cv_number = 5, seed_to_fix = 100)
#1b SRF -> Cox 
ens1bxt = method_1B_cv(sim_data_xt, predict.factors, fixed_time = 5, cv_number = 5, seed_to_fix = 100,
                       pretrained_srf_models = srfxt$pretrained_srf_models)
#2a RPart Tree -> Cox in Clusters 
ens2axt = method_2A_cv(sim_data_xt, predict.factors, fixed_time = 5, cv_number = 5, seed_to_fix = 100)
#2b SRF Tree -> Cox in Clusters 
ens2bxt = method_2B_cv(sim_data_xt, predict.factors, fixed_time = 5, cv_number = 3, seed_to_fix = 100)
#3 RPart TREE -> Modified Cox with clusters IDs 
ens3xt = method_3_cv(sim_data_xt, predict.factors, fixed_time = 5, cv_number = 5, seed_to_fix = 100)
#test results 

m2a = method_2A_train(sim_data_nonlinear, predict.factors, fixed_time = 8,  seed_to_fix = 100)
predict2a = method_2A_predict(m2a, sim_data_nonlinear_exttest, fixed_time=8)
method_validate

method_3_train(sim_data_nonlinear, predict.factors, fixed_time = 8,  seed_to_fix = 100)

res_xt = cbind("Cox" = coxxt$testaverage, 
               "SRF" = srfxt$testaverage,
               #"CoxMFP"= coxxtmfp$testaverage, 
               "1A" = ens1axt$testaverage,
               "1B" = ens1bxt$testaverage,   "2A" = ens2axt$testaverage, 
               "2B" = ens2bxt$testaverage, 
               "3" = ens3xt$testaverage)
#train results 
t(cbind(  "Cox" = coxxt$trainaverage,   "SRF" = srfxt$trainaverage,  
          "1A" = ens1axt$trainaverage,
          "1B" = ens1bxt$trainaverage,   "2A" = ens2axt$trainaverage,  
          "2B" = ens2bxt$trainaverage,
          "3" = ens3xt$trainaverage))

clipr::write_clip(res_xt)
clipr::write_clip(std_test_xt)
apply(cox$test,  2, sd)

std_test_xt =   cbind(  "Cox" = apply(coxxt$test,  2, sd),   "SRF" = apply(srfxt$test,  2, sd),  
                        "CoxMFP" = apply(coxxtmfp$test,  2, sd),
                        "1A" = apply(ens1axt$test,  2, sd),  "1B" = apply(ens1bxt$test,  2, sd),  
                        "2A" =apply(ens2axt$test,  2, sd),   "2B" = apply(ens2bxt$test,  2, sd),  "3" = apply(ens3xt$test,  2, sd))




#res_xt
#                       Cox     SRF      1A       1B       2A       2B
#times_to_predict    5.00000 5.00000 5.00000  5.00000 5.000000 5.000000
#auc_score           0.64771 0.71593 0.71545  0.70596 0.720765 0.715451
#brier_score         0.43422 0.42208 0.42269  0.41911 0.418767 0.417919
#brier_score_scaled  0.03496 0.06219 0.06069  0.06913 0.069543 0.073542
#c_score             0.61007 0.65693 0.65687  0.65092 0.661672 0.656760
#calibration_slope   0.95760 0.97948 0.98638  0.93901 0.986649 0.964383
#calibration_alpha  -0.01235 0.01015 0.01332 -0.01599 0.003473 0.003972
#test                1.00000 1.00000 1.00000  1.00000 1.000000 1.000000
#                        3
#times_to_predict    5.000000
#auc_score           0.721165
#brier_score         0.416939
#brier_score_scaled  0.073554
#c_score             0.662289
#calibration_slope   0.994691
#calibration_alpha  -0.001042
#test                1.000000

############################end Cross-terms##################################################

############# For paper - Cross-terms trees and Cox clusters display ####################

predict.factors = c("age", "bmi", "hyp", "gender")
sim_data_xt = simulatedata_crossterms(1000, 10, 0, 101) #5000

m2a_one = method_2A_train(sim_data_xt, predict.factors, fixed_time = 5,  seed_to_fix = 101)
#m2b_one = method_2B_train(sim_data_xt, predict.factors, fixed_time = 5,  seed_to_fix = 100)

m3_one = method_3_train(sim_data_xt, predict.factors, fixed_time = 5,  seed_to_fix = 101)
clipr::write_clip(summary(m3_one$modcoxmodel)$coefficients);summary(m3_one$modcoxmodel)$coefficients

a = summary(m3_one$modcoxmodel)$coefficients
rownames(a) = c("age", "bmi", "hyp", "gender", "cluster_2", "cluster_3", "cluster_4", "cluster_5") # "cluster_6", "cluster_7")
round(a,4)
summary(m3_one$modcoxmodel)$coefficients

m3_one$treemodel
m2a_one


rpart.plot(m3_one$treemodel, nn=TRUE, roundint = FALSE, digits = 2)
#plot RPart tree
rpart.plot(m2a_one$treemodel, nn=TRUE, roundint = FALSE, digits = 2)



#'############################################################################
################################ Linear ######################################
#'############################################################################
predict.factors = c("age", "bmi", "hyp", "gender")

sim_data_l = simulatedata_linear(5000, 10, 0, 100)
sim_data_l_exttest = simulatedata_linear(20000, 10, 0, 100)

populationstats(sim_data_xt, time_f = 8)
populationstats(sim_data_xt, time_f = 5)

#Cox
coxl = method_cox_cv(sim_data_l, predict.factors, fixed_time = 5, cv_number = 5)
#CoxMFP
coxlmfp = method_coxmfp_cv(sim_data_l, predict.factors, fixed_time = 5, cv_number = 5)
#SRF
srfl = method_srf_cv(sim_data_l, predict.factors, fixed_time = 5, cv_number = 5, seed_to_fix = 100)
#1a Cox -> SRF
ens1al = method_1A_cv(sim_data_l, predict.factors, fixed_time = 5, cv_number = 5, seed_to_fix = 100)
#1b SRF -> Cox 
ens1bl = method_1B_cv(sim_data_l, predict.factors, fixed_time = 5, cv_number = 5, seed_to_fix = 100,
                      pretrained_srf_models = srfl$pretrained_srf_models)
#2a RPart Tree -> Cox in Clusters 
ens2al = method_2A_cv(sim_data_l, predict.factors, fixed_time = 5, cv_number = 5, seed_to_fix = 100)
#2b SRF Tree -> Cox in Clusters 
ens2bl = method_2B_cv(sim_data_l, predict.factors, fixed_time = 5, cv_number = 3, seed_to_fix = 100)
#3 RPart TREE -> Modified Cox with clusters IDs 
ens3l = method_3_cv(sim_data_l, predict.factors, fixed_time = 5, cv_number = 5, seed_to_fix = 100)
#test results 

m2a = method_2A_train(sim_data_nonlinear, predict.factors, fixed_time = 8,  seed_to_fix = 100)
predict2a = method_2A_predict(m2a, sim_data_nonlinear_exttest, fixed_time=8)
method_validate

method_3_train(sim_data_nonlinear, predict.factors, fixed_time = 8,  seed_to_fix = 100)

res_l = cbind(  "Cox" = coxl$testaverage, "CoxMFP" = coxlmfp$testaverage,
                "SRF" = srfl$testaverage,   "1A" = ens1al$testaverage,
                "1B" = ens1bl$testaverage,   "2A" = ens2al$testaverage,   "2B" = ens2bl$testaverage,  "3" = ens3l$testaverage)
#train results 
cbind(  "Cox" = coxl$trainaverage,   "SRF" = srfl$trainaverage,   "1A" = ens1al$trainaverage,
        "1B" = ens1bl$trainaverage,   "2A" = ens2al$trainaverage,   "2B" = ens2bl$trainaverage,  "3" = ens3l$trainaverage)

clipr::write_clip(res_l)
clipr::write_clip(std_test_l)
apply(cox$test,  2, sd)

std_test_l =  cbind(  "Cox" = apply(coxl$test,  2, sd),   "SRF" = apply(srfl$test,  2, sd),  
                      "CoxMFP" = apply(coxlmfp$test,  2, sd),
                      "1A" = apply(ens1al$test,  2, sd),  "1B" = apply(ens1bl$test,  2, sd),  
                      "2A" =apply(ens2al$test,  2, sd),   "2B" = apply(ens2bl$test,  2, sd),  "3" = apply(ens3l$test,  2, sd))


########################### end non-linear######################################


#'#############################################################################
################################# ELSA ########################################
#' ############################################################################

data_diabetes = read.csv("C:/Users/dinab/Desktop/PhD Projects/Elsa Olesya/stata/stata13_se/diabetes_data_for_method.csv")
dim(data_diabetes)
predictors_diabetes = c("age_0_", "genderdum", "bmi_0_", "cvd_0", "hyp_0",
                        "baseline_exercise",  "B_wealth","baseline_B_smokstatus",
                        "t2dm_", "pc1_", "pc2_", "pc3_")

set.seed(101)
data_diabetes$time = data_diabetes$time_1

data_diabetes[predictors_diabetes]

#Cox
coxd = method_cox_cv(data_diabetes, predictors_diabetes, fixed_time = 7, cv_number = 5)
#CoxMFP
coxmfpd = method_coxmfp_cv(data_diabetes, predictors_diabetes, fixed_time = 7, cv_number = 5)
#SRF
srfd = method_srf_cv(data_diabetes, predictors_diabetes, fixed_time = 7, cv_number = 3, seed_to_fix = 100)
#1a Cox -> SRF
ens1ad = method_1A_cv(data_diabetes, predictors_diabetes, fixed_time = 7, cv_number = 3, seed_to_fix = 100)
#1b SRF -> Cox 
ens1bd = method_1B_cv(data_diabetes, predictors_diabetes, fixed_time = 7, cv_number = 3, seed_to_fix = 100,
                      pretrained_srf_models = srfd$pretrained_srf_models)
#2a RPart Tree -> Cox in Clusters 
ens2ad = method_2A_cv(data_diabetes, predictors_diabetes, fixed_time = 7, cv_number = 3, seed_to_fix = 101)
#2b SRF Tree -> Cox in data_diabetes 
ens2bd = method_2B_cv(data_diabetes, predictors_diabetes, fixed_time = 7, cv_number = 3, seed_to_fix = 100)
#3 RPart TREE -> Modified Cox with clusters IDs 
ens3d = method_3_cv(data_diabetes, predictors_diabetes, fixed_time = 7, cv_number = 3, seed_to_fix = 100)
#test results 

res_d = cbind(  "Cox" = coxd$testaverage, "CoxMFP" = coxmfpd$testaverage,
                "SRF" = srfd$testaverage,   "1A" = ens1ad$testaverage,
                "1B" = ens1bd$testaverage,      "3" = ens3d$testaverage)
#train results 
cbind(  "Cox" = coxl$trainaverage,   "SRF" = srfl$trainaverage,   "1A" = ens1al$trainaverage,
        "1B" = ens1bl$trainaverage,   "2A" = ens2al$trainaverage,   "2B" = ens2bl$trainaverage,  "3" = ens3l$trainaverage)

clipr::write_clip(res_d)
clipr::write_clip(std_test_d)
apply(cox$test,  2, sd)

std_test_d =  cbind(  "Cox" = apply(coxl$test,  2, sd),   "SRF" = apply(srfl$test,  2, sd),  
                      "CoxMFP" = apply(coxlmfp$test,  2, sd),
                      "1A" = apply(ens1al$test,  2, sd),  "1B" = apply(ens1bl$test,  2, sd),  
                      "3" = apply(ens3l$test,  2, sd))



#'#############################################################################
#'
##################### Validating performance on large 20 000 samples ###########
#'
#'#############################################################################

sim_data_nonlinear = simulatedata_linear(N=1000, observe_time = 10, percentcensored = 0.0, randomseed = 100)
sim_data_nonlinear_exttest = simulatedata_nonlinear(N=20000, observe_time =10, percentcensored = 0.1, randomseed = 100)
sim_data_nonlinear1 = simulatedata_nonlinear(N=1000, observe_time =10, percentcensored = 0.0, randomseed = 100)

mean(sim_data_nonlinear$event_time)
mean(sim_data_nonlinear1$event_time)
populationstats(sim_data_nonlinear, 5)

mcox = method_cox_train(sim_data_nonlinear, predict.factors)
predictcoxext = method_cox_predict(mcox, sim_data_nonlinear_exttest, 5)
predictcox = method_cox_predict(mcox, sim_data_nonlinear, 5)
valcox = method_any_validate(predictcox, 5, sim_data_nonlinear, sim_data_nonlinear)
valcoxext = method_any_validate(predictcoxext, 5, sim_data_nonlinear, sim_data_nonlinear_exttest)

mcoxmfp = method_coxmfp_train(sim_data_nonlinear, predict.factors)
predictcoxmfpext = method_cox_predict(mcoxmfp, sim_data_nonlinear_exttest, 5)
predictcoxmfp = method_cox_predict(mcoxmfp, sim_data_nonlinear, 5)
valmfp = method_any_validate(predictcoxmfp, 5, sim_data_nonlinear, sim_data_nonlinear)
valmfpext = method_any_validate(predictcoxmfpext, 5, sim_data_nonlinear, sim_data_nonlinear_exttest)

msrf = method_srf_train(sim_data_nonlinear, predict.factors, fixed_time = 8,  seed_to_fix = 100)
predictsrfext = method_srf_predict(msrf$model, sim_data_nonlinear_exttest, fixed_time=8)
predictsrf = method_srf_predict(msrf$model, sim_data_nonlinear, fixed_time=8)
valsrf = method_any_validate(predictsrf, 5, sim_data_nonlinear, sim_data_nonlinear)
valsrfext = method_any_validate(predictsrfext, 5, sim_data_nonlinear, sim_data_nonlinear_exttest)

m1b = method_1B_train(sim_data_nonlinear, predict.factors, fixed_time = 8,  seed_to_fix = 100, pretrained_srf_model = msrf$model)
predict1bext = method_1B_predict(m1b, sim_data_nonlinear_exttest, 8)
predict1b = method_1B_predict(m1b, sim_data_nonlinear, 8)
val1b =method_any_validate(predict1b, 5, sim_data_nonlinear, sim_data_nonlinear)
val1bext =method_any_validate(predict1bext, 5, sim_data_nonlinear, sim_data_nonlinear_exttest)

m2a = method_2A_train(sim_data_nonlinear, predict.factors, fixed_time = 8,  seed_to_fix = 100)
predict2aext = method_2A_predict(m2a, sim_data_nonlinear_exttest, fixed_time=8)
predict2a = method_2A_predict(m2a, sim_data_nonlinear, fixed_time=8)
val2a =method_any_validate(predict2a, 5, sim_data_nonlinear, sim_data_nonlinear)
val2aextf =method_any_validate(predict2aext, 5, sim_data_nonlinear, sim_data_nonlinear_exttest)

m3 = method_3_train(sim_data_nonlinear, predict.factors, fixed_time = 8,  seed_to_fix = 100)
predict3ext = method_3_predict(m3, sim_data_nonlinear_exttest, fixed_time=8)
predict3 = method_3_predict(m3, sim_data_nonlinear, fixed_time=8)
val3 =method_any_validate(predict3, 5, sim_data_nonlinear, sim_data_nonlinear)
val3ext  =method_any_validate(predict3ext, 5, sim_data_nonlinear, sim_data_nonlinear_exttest)

library(clipr)
rbind(valcoxext,valmfpext,valsrfext, val1bext, val2aextf, val3ext)
rbind(valcox,valmfp,valsrf, val1b, val2a, val3)
clipr::write_clip(rbind(valcoxext,valmfpext,valsrfext, val1bext, val2aextf, val3ext))
clipr::write_clip(rbind(valcox,valmfp,valsrf, val1b, val2a, val3))
