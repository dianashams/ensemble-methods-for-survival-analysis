# ################################################################################################################
#' Illutrative example 
#' Training and internally cross-validating projects' ensemble methods on
#' GBSG2 data 
#' https://scikit-survival.readthedocs.io/en/stable/api/generated/sksurv.datasets.load_gbsg2.html
# ################################################################################################################
# Jan 2023 version
# ################################################################################################################


# load GBSG2 data 
library(survival)
library(pec) #for GBSG2 data
data("GBSG2")

# re-format the data - hot-coding binary variables horTh and menostat, 
# also assuming that tgrade is ordinary variable 1<2<3
gbsg_data = GBSG2
gbsg_data$horTh = ifelse(gbsg_data$horTh== "no", 0, 1)
gbsg_data$menostat = ifelse(gbsg_data$menostat == "Post", 1, 0)
gbsg_data$tgrade = ifelse(gbsg_data$tgrade == "I", 1, ifelse(gbsg_data$tgrade=="II", 2, 3))
for (i in 1:dim(gbsg_data)[2]) print(c(names(gbsg_data)[i], class(gbsg_data[,i])))
params = c("age", "horTh", "menostat", "tsize", "tgrade", "pnodes", "progrec", "estrec")
gbsg_data$event = gbsg_data$cens
gbsg_data$time = gbsg_data$time/365 #convert into years 

#final analytical data
gbsg_data = gbsg_data[c("time", "event", params)]
quantile(gbsg_data[gbsg_data$event==1, "time"],0.95) #4.96
quantile(gbsg_data[gbsg_data$event==0, "time"],0.95) #6.34
final_time = 5

#----------- 1) COXPH ----------------------------------------------

#apparent validation: 
# train 
m3 = method_cox_train(gbsg_data, params)
#predict for t= 10
p10 = method_cox_predict(m3, gbsg_data, final_time)
#c-index
concordancefit(Surv(gbsg_data$time, gbsg_data$event),1-p10)$concordance #0.6879
#other stats
method_any_validate(p10, final_time, gbsg_data, gbsg_data)
#  T      AUCROC       BS  BS_scaled   C_score  Calib_slope  Calib_alpha
#  5   0.7440482 0.210829 0.1909064 0.6879283    1.409447   0.7512217

#internal 3-fold CV
coxcv = method_any_cv(gbsg_data, params, method_cox_train, method_cox_predict, 
                      valuation_times = final_time, #c(final_time*0.7,final_time),
                      cv_number = 10, seed_for_cv = 2024, parallel=TRUE)
round(coxcv$testaverage,3)
#     T      AUCROC          BS   BS_scaled     C_score Calib_slope Calib_alpha        test 
# 5.000       0.734       0.216       0.168       0.677       1.464       0.773       1.000   

# round(coxcv$test,3)
#    T AUCROC    BS BS_scaled C_score Calib_slope Calib_alpha test cv_n
# 1  5  0.757 0.192     0.238   0.656       1.689       1.055    1    1
# 2  5  0.702 0.201     0.078   0.691       0.412       0.417    1    2
# 3  5  0.806 0.160     0.314   0.677       2.004       1.103    1    3
# 4  5  0.622 0.220     0.061   0.629       0.811       1.442    1    4
# 5  5  0.675 0.272     0.095   0.705       1.195       0.650    1    5
# 6  5  0.790 0.204     0.228   0.721       2.065       0.706    1    6
# 7  5  0.794 0.212     0.208   0.666       2.374       0.583    1    7
# 8  5  0.753 0.241     0.131   0.698       0.860       0.422    1    8
# 9  5  0.680 0.231     0.119   0.665       1.797       0.448    1    9
# 10 5  0.758 0.224     0.208   0.666       1.434       0.906    1   10

#---------- 2) SRF--------------------------------------------------

# train 
srf_m = method_srf_train(gbsg_data, params, final_time, inner_cv = 3, seed_to_fix = 42)
#predict for t= 10
p10srf = method_srf_predict(srf_m, gbsg_data, final_time)
#c-index
concordancefit(Surv(gbsg_data$time, gbsg_data$event),1-p10srf)$concordance #0.8239
#other stats
round(method_any_validate(p10srf, final_time, gbsg_data, gbsg_data),3)
#   T  AUCROC    BS BS_scaled C_score Calib_slope Calib_alpha
#   5  0.916  0.145     0.444   0.824       2.855       0.661

#internal 3-fold CV
srf_cv = method_any_cv(gbsg_data, params, method_srf_train, method_srf_predict, 
                       valuation_times = final_time, #c(final_time*0.7,final_time),
                       cv_number = 5, seed_for_cv = 2023, parallel=TRUE, 
                       model_args = list("inner_cv"=3, "verbose"=TRUE),
                       predict_args = list())
#results by cv step
srf_cv$test
#   T    AUCROC        BS  BS_scaled   C_score Calib_slope Calib_alpha test cv_n
# 1 5 0.7593788 0.2204363 0.20179222 0.7021088    1.469467  0.8412348    1    1
# 2 5 0.7081672 0.2353548 0.08759203 0.7164307    1.323320  0.2847982    1    2
# 3 5 0.7405126 0.2097591 0.16009892 0.7091582    1.409517  0.6843569    1    3
# 4 5 0.6933945 0.2379000 0.12700981 0.6440329    1.232670  1.0690240    1    4
# 5 5 0.7271822 0.1833009 0.23616684 0.7235089    1.627641  0.9654919    1    5

#average CV results
round(srf_cv$testaverage,3)
#    T      AUCROC          BS   BS_scaled     C_score Calib_slope Calib_alpha        test 
#5.000       0.726       0.217       0.163       0.699       1.413       0.769       1.000


#---------- 3) Cox + fractional polynomials ----------------------------------------------

m_coxf = method_coxmfp_train(gbsg_data, params)
m_coxf$coefficients
#       hyp      I((age + 1.8)^1) I(((bmi + 3.5)/10)^1) I(((bmi + 3.5)/10)^2) 
# 1.0270031             0.2656635           -22.6543214            32.1645691 

p10_coxf = method_coxmfp_predict(m_coxf,gbsg_data, final_time)
round(method_any_validate(p10_coxf, final_time, gbsg_data, gbsg_data),3)
#    T AUCROC    BS BS_scaled C_score Calib_slope Calib_alpha
# 1  5  0.731 0.21     0.192   0.702       1.013       0.765

coxf_cv2 = method_any_cv(df = gbsg_data, params, method_coxmfp_train, method_coxmfp_predict, 
                         final_time , cv_number = 5, seed_for_cv = 42)
round(coxf_cv2$testaverage,3)
#     T      AUCROC          BS   BS_scaled     C_score Calib_slope Calib_alpha        test 
# 5.000       0.714       0.226       0.123       0.687       0.893       0.811       1.000 

#---------- 4) Ensemble 1 (COX -> SRF) ----------------------------------------------

m1a = method_1A_train(gbsg_data, params, final_time, inner_cv=3, seed_to_fix = 42)
p10_1a = method_1A_predict(model_1a = m1a,df_test = gbsg_data,fixed_time = final_time)
round(method_any_validate(p10_1a, final_time, gbsg_data, gbsg_data),3)
#     T   AUCROC     BS BS_scaled C_score  Calib_slope Calib_alpha
# 1   5  0.859    0.171     0.343   0.779       2.078       0.707
ens1a_cv = method_any_cv(df = gbsg_data, params, method_1A_train, method_1A_predict, 
                         valuation_times = final_time , cv_number = 5, seed_for_cv = 42)
round(ens1a_cv$testaverage,3)
#       T      AUCROC          BS   BS_scaled     C_score Calib_slope Calib_alpha        test 
#   5.000       0.726       0.222       0.138       0.696       1.452       0.641       1.000 

#var importance from SRF
round(m1a$vimp10,3)
# pnodes cox_predict     progrec         age      estrec       horTh      tgrade    menostat  tsize
#  0.037       0.035       0.021       0.016       0.006       0.002       0.001       0.000  0.000

#---------- 4) Ensemble 2 (Tree -> Cox in leaves) ----------------------------------------------

m2a = method_2A_train(gbsg_data, params, final_time, inner_cv=3, seed_to_fix = 42)

#plot the tree
rpart.plot(m2a$treemodel)
# n= 686 
# node), split, n, deviance, yval
# * denotes terminal node
# 1) root 686 956.4365 1.0000000  
# 2) pnodes< 3.5 376 410.6060 0.6362878  
# 4) progrec>=89.5 137 117.3362 0.3838470 *
#   5) progrec< 89.5 239 279.2723 0.8067135 *
#   3) pnodes>=3.5 310 481.4419 1.6102330 *

# cox models in th final leaves
m2a$coxmodels

p10_2a = method_2A_predict(model_2a = m2a,df_test = gbsg_data, fixed_time = final_time)
round(method_any_validate(p10_2a, final_time, gbsg_data, gbsg_data),3)
# T   AUCROC    BS     BS_scaled   C_score  Calib_slope Calib_alpha
# 5  0.749  0.204     0.216      0.71       1.059       0.144
ens2a_cv = method_any_cv(df = gbsg_data, params, method_2A_train, method_2A_predict, 
                         valuation_times = final_time , cv_number = 5, seed_for_cv = 42, 
                         model_args = list("inner_cv"=3, verbose= FALSE))
round(ens2a_cv$testaverage,3)
#       T      AUCROC          BS   BS_scaled     C_score Calib_slope Calib_alpha        test 
#   5.000       0.687       0.236       0.085       0.672       0.759       0.894       1.000 

#---------- 4) Ensemble 2 (Tree -> Cox in leaves) ----------------------------------------------

m3 = method_3_train(gbsg_data, params, final_time,inner_cv = 3, seed_to_fix = 42,useCoxLasso = FALSE)
#plot the tree 
rpart.plot(m3$treemodel, roundint = FALSE)

# n= 686 
# node), split, n, deviance, yval
# * denotes terminal node
# 1) root 686 956.4365 1.0000000  
# 2) pnodes< 3.5 376 410.6060 0.6362878 *
#   3) pnodes>=3.5 310 481.4419 1.6102330 *
  
ens3_cv = method_any_cv(gbsg_data, params, method_3_train, method_3_predict, cv_number = 5,
                        seed_for_cv = 42, valuation_times = final_time, parallel=FALSE,
                        model_args = list("inner_cv"=3, verbose=FALSE,"seed_to_fix" = 601))
round(ens3_cv$testaverage,3)
#      T      AUCROC          BS   BS_scaled     C_score Calib_slope Calib_alpha        test 
#  5.000       0.719       0.222       0.136       0.680       1.015       0.809       1.000 

#------------------------ Altogerther ----------------------------------------------
final_results = round(
  cbind("cox" = coxcv$testaverage,
        "cox_mfp" = coxf_cv2$testaverage,
        "srf" = srf_cv$testaverage,
        "ens1" = ens1a_cv$testaverage,
        "ens2" = ens2a_cv$testaverage,
        "ens3" = ens3_cv$testaverage), 4)

final_results

#                cox  cox_mfp   srf   ens1   ens2   ens3
# T           5.0000  5.0000 5.0000 5.0000 5.0000 5.0000
# AUCROC      0.7338  0.7143 0.7257 0.7263 0.6867 0.7187
# BS          0.2157  0.2259 0.2174 0.2223 0.2357 0.2224
# BS_scaled   0.1681  0.1225 0.1625 0.1376 0.0847 0.1364
# C_score     0.6773  0.6871 0.6990 0.6959 0.6718 0.6797
# Calib_slope 1.4641  0.8934 1.4125 1.4520 0.7587 1.0155
# Calib_alpha 0.7732  0.8114 0.7690 0.6408 0.8943 0.8085
# test        1.0000  1.0000 1.0000 1.0000 1.0000 1.0000

#' =>  WHAT DOES IT TELL ABOUT THE DATA? 
#' 1) Cox's calibration is quite bad (look at  coxcv$test for individual values for alpha and slope)
#' 2) some transformation of the coefficients as in cox-mfp seem to make calibration better,
#' however, the terms are not easy to interpret (check coxf_cv2$model_list with (pnodes/10)^-2)
#' 3) Ens 3 looks best on balance of discrimination and calibration
#' 4) Let's check how the tree divides the sample and which non-linear terms are added to the model in 5 CV trials:

for (i in 1:5) print (ens3_cv$model_list[[i]]$treemodel)
# => it always splits by pnodes at 3.5 or 4.5 threschold 

#5) check what is the impact size of the added non-linear I(pnodes>threschold ) term
for (i in 1:5) print (summary(ens3_cv$model_list[[i]]$modcoxmodel))
#=> it suggests  x2 the risk for pnodes>3.5, while pnodes coefficient becomes negligible. 
# i.e. pnodes acts as a step-function rather than a monotonous risk increase per each additional pnode

#Finally, we can check the trees in ensemlbe 2 to see if there are any other non-linearities being picked up
for (i in 1:5) print (ens2a_cv$model_list[[i]]$treemodel)
# similarly, it first splits by pnodes at 3.5, but then it picks some interaction with pnodes x estrec and pnodes x progrec
# overall, pnodes split seems stable across the ensemble methods, other non-linearity may still be sporadic 


