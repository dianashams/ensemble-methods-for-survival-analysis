
d1 = simulatedata_linear(1000, percentcensored = 0.5, randomseed=42, distr= "Exp")
d2 = simulatedata_nonlinear(1000, percentcensored = 0.5, randomseed=42, distr= "Weibull", rho_w=1.5)
d3 = simulatedata_crossterms(1000, percentcensored = 0.5, randomseed=42, distr= "Weibull", rho_w=0.8)
d4 = simulatedata_lin_nonPH(1000, percentcensored = 0.5,  randomseed=42)
params = names(d_xt)[1:4] # "age", "bmi", "hyp", "sex

populationstats(d1, time_f = 10, "linear")
populationstats(d2, time_f = 10, "linear")
populationstats(d3, time_f = 10, "linear")
populationstats(d4, time_f = 10, "linear")

# Train and Cross-validate all methods for cross-term data d3 (change it to d1-d4 or your data)

#----------- 1) COXPH ----------------------------------------------

#apparent validation: 
# train 
m3 = method_cox_train(d3, params)
#predict for t= 10
p10 = method_cox_predict(m3, d3, 10)
#c-index
concordancefit(Surv(d3$time, d3$event),1-p10)$concordance #0.6138
#other stats
method_any_validate(p10, 10, d3, d3)
#  T      AUCROC       BS  BS_scaled   C_score  Calib_slope  Calib_alpha
#1 10 0.64125592 0.1193591 0.2093106 0.6138496    0.9008747  0.05405041

#internal 3-fold CV
coxcv = method_any_cv(d3, params, method_cox_train, method_cox_predict, 
                      valuation_times = c(8,10),
                      cv_number = 3, seed_for_cv = 2024, parallel=TRUE)
round(coxcv$testaverage,3)
#     T      AUCROC          BS   BS_scaled     C_score Calib_slope Calib_alpha        test 
# 9.000       0.640       0.174       0.147       0.614       0.924       0.047       1.000 

round(coxcv$test,3)
#    T  AUCROC    BS BS_scaled C_score Calib_slope Calib_alpha test cv_n
# 1  8  0.666 0.219     0.112   0.632       1.098       0.030    1    1
# 2  8  0.613 0.235     0.057   0.602       0.724       0.027    1    2
# 3  8  0.639 0.229     0.091   0.607       1.007       0.066    1    3
# 4 10  0.672 0.110     0.278   0.632       1.068       0.040    1    1
# 5 10  0.614 0.123     0.194   0.602       0.692       0.051    1    2
# 6 10  0.634 0.126     0.147   0.607       0.953       0.070    1    3


#---------- 2) SRF--------------------------------------------------

# train 
srf_m = method_srf_train(d3, params, 10, inner_cv = 3, seed_to_fix = 42)
#predict for t= 10
p10srf = method_srf_predict(srf_m, d3, 10)
#c-index
concordancefit(Surv(d3$time, d3$event),1-p10srf)$concordance #0.7765
#other stats
round(method_any_validate(p10srf, 10, d3, d3),3)
#    T AUCROC    BS BS_scaled C_score Calib_slope Calib_alpha
# 1 10  0.862 0.086     0.431   0.777        1.55       0.044

#internal 3-fold CV
srf_cv = method_any_cv(d3, params, method_srf_train, method_srf_predict, 
                      valuation_times = c(8,10),
                      cv_number = 3, seed_for_cv = 2023, parallel=TRUE, 
                      model_args = list("inner_cv"=2, "verbose"=TRUE),
                      predict_args = list())
srf_cv$test
#    T    AUCROC         BS   BS_scaled   C_score Calib_slope Calib_alpha test cv_n
# 1  8 0.7506296 0.21026186  0.19343012 0.6920613   0.8020854  0.08842607    1    1
# 2  8 0.7851353 0.18311813  0.23915556 0.7436045   1.0051205 -0.03037066    1    2
# 3  8 0.7929125 0.17768885  0.27774247 0.7230778   1.0311116  0.05119165    1    3
# 4 10 0.7365965 0.15274265 -0.06153357 0.6924767   0.7331603  0.12930842    1    1
# 5 10 0.8163315 0.07055782  0.55167129 0.7487876   1.2199739 -0.02459155    1    2
# 6 10 0.7845998 0.09838314  0.34080398 0.7273556   0.9334690  0.04282694    1    3
round(srf_cv$testaverage,3)
#    T      AUCROC          BS   BS_scaled     C_score Calib_slope Calib_alpha        test 
#9.000       0.778       0.149       0.257       0.721       0.954       0.043       1.000 


#---------- 3) Cox + fractional polynomials ----------------------------------------------

m_coxf = method_coxmfp_train(d3, params)
m_coxf$coefficients
  #       hyp      I((age + 1.8)^1) I(((bmi + 3.5)/10)^1) I(((bmi + 3.5)/10)^2) 
  # 1.0270031             0.2656635           -22.6543214            32.1645691 

p10_coxf = method_coxmfp_predict(m_coxf,d3, 10)
round(method_any_validate(p10_coxf, 10, d3, d3),3)
  #    T AUCROC    BS BS_scaled C_score Calib_slope Calib_alpha
  # 1 10  0.755 0.105     0.302   0.711       1.255       0.045

coxf_cv2 = method_any_cv(df = d3, params, method_coxmfp_train, method_coxmfp_predict, 
                        10 , cv_number =3, seed_for_cv = 42)
round(coxf_cv2$testaverage,3)
  #     T      AUCROC          BS   BS_scaled     C_score Calib_slope Calib_alpha        test 
  #10.000       0.748       0.108       0.282       0.707       1.197       0.046       1.000

#---------- 4) Ensemble 1 (COX -> SRF) ----------------------------------------------

m1a = method_1A_train(d3, params, 10, inner_cv=3, seed_to_fix = 42)
p10_1a = method_1A_predict(model_1a = m1a,df_test = d3,fixed_time = 10)
round(method_any_validate(p10_1a, 10, d3, d3),3)
#     T   AUCROC    BS BS_scaled C_score  Calib_slope Calib_alpha
# 1   10  0.854  0.089     0.413    0.770       1.386       0.046
ens1a_cv = method_any_cv(df = d3, params, method_1A_train, method_1A_predict, 
                         valuation_times = 10 , cv_number =3, seed_for_cv = 42)
round(ens1a_cv$testaverage,3)
  #       T      AUCROC          BS   BS_scaled     C_score Calib_slope Calib_alpha        test 
  #  10.000       0.770       0.085       0.431       0.717       0.888      -0.008       1.000 

#var importance from SRF
m2a$vimp10
#bmi          hyp          age          sex 
#0.103548439  0.042983940  0.022202066 -0.002699149  

#---------- 4) Ensemble 2 (Tree -> Cox in leaves) ----------------------------------------------

m2a = method_2A_train(d3, params, 10, inner_cv=3, seed_to_fix = 42)

#plot the tree
rpart.plot(m2a$treemodel)
# cox models in th final leaves
m2a$coxmodels

p10_2a = method_2A_predict(model_2a = m2a,df_test = d3, fixed_time = 10)
round(method_any_validate(p10_2a, 10, d3, d3),3)
#         T   AUCROC    BS BS_scaled C_score  Calib_slope Calib_alpha
# 10  0.767    0.109           0.275   0.722       0.966       0.047
ens2a_cv = method_any_cv(df = d3, params, method_2A_train, method_2A_predict, 
                       valuation_times = 10 , cv_number =3, seed_for_cv = 42, 
                       model_args = list("inner_cv"=3, verbose= FALSE))
round(ens2a_cv$testaverage,3)
#       T      AUCROC          BS   BS_scaled     C_score Calib_slope Calib_alpha        test 
#  10.000       0.750       0.111       0.258       0.705       0.890       0.056       1.000

#---------- 4) Ensemble 2 (Tree -> Cox in leaves) ----------------------------------------------

m3 = method_3_train(d3, params, 10,inner_cv = 5, seed_to_fix = 42,useCoxLasso = FALSE)
#plot the tree 
rpart.plot(m3$treemodel, roundint = FALSE)

ens3_cv = method_any_cv(d3, params, method_3_train, method_3_predict, cv_number = 3,
                         seed_for_cv = 14, valuation_times = 10, parallel=FALSE,
                         model_args = list("inner_cv"=3, verbose=FALSE,"seed_to_fix" = 601))
round(ens3_cv$testaverage,3)
#      T      AUCROC          BS   BS_scaled     C_score Calib_slope Calib_alpha        test 
# 10.000       0.694       0.112       0.257       0.657       0.862       0.045       1.000 

#---------- Altogerther ----------------------------------------------
final_results = round(
  cbind("cox" = coxcv$testaverage,
  "cox_mfp" = coxf_cv2$testaverage,
  "srf" = srf_cv$testaverage,
  "ens1" = ens1a_cv$testaverage,
  "ens2" = ens2a_cv$testaverage,
  "ens3" = ens3_cv$testaverage), 4)

final_results

# T           9.0000 10.0000 9.0000 10.0000 10.0000 10.0000
# AUCROC      0.6398  0.7478 0.7777  0.7696  0.7502  0.6937
# BS          0.1738  0.1076 0.1488  0.0852  0.1110  0.1117
# BS_scaled   0.1466  0.2822 0.2569  0.4309  0.2583  0.2573
# C_score     0.6135  0.7070 0.7212  0.7169  0.7046  0.6573
# Calib_slope 0.9238  1.1969 0.9542  0.8875  0.8897  0.8616
# Calib_alpha 0.0474  0.0458 0.0428 -0.0081  0.0559  0.0452
# test        1.0000  1.0000 1.0000  1.0000  1.0000  1.0000
