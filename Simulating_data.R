
###########################################################################
populationstats = function(df_stats, time_f,  namedf = "df"){
  # function to show event statistics of a simulated data
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

###########################################################################

# These functions return simulated data with exponential or weibull distributions
# 0.5 / 0.75 / 0.9 - by time = 10 (scaled) expected number of events
# distr = "Exp" or "Weibull"  
# (rho_w =0.5, lambda = 0.447) hazard slopes down
# (rho_w=1.5, lambda 0.027) hazard up and down;
# drop-out - additional drop out before the end of study (expected, independent from event)


simulatedata_linear = function (N=1000, observe_time =10, 
                      percentcensored = 0.75,  
                      randomseed = 100,lambda = 0.1,
                      distr = "Exp",  rho_w = 1, 
                      drop_out = 0.3){ 
  #simulate the data
  df = simulate_population(N, randomseed)
  # calculate betas
  exp_beta = linear_beta(df)
  df["exp_beta"] = exp_beta
  # simulate censored and event times 
  df = df_event_times(exp_beta=exp_beta, df=df, N=N, observe_time=observe_time, 
                      percentcensored=percentcensored,randomseed= randomseed+1,
                      lambda = lambda, distr = distr, rho_w = rho_w, drop_out=drop_out)
  return(df)
}
  
simulatedata_nonlinear = function (N=1000, observe_time =10, 
                                percentcensored = 0.75,  
                                randomseed = 100,lambda = 0.1,
                                distr = "Exp",  rho_w = 1, 
                                drop_out = 0.3){ 
  #simulate the data
  df = simulate_population(N, randomseed)
  # calculate betas
  exp_beta = nonlinear_beta(df)
  # simulate censored and event times 
  df = df_event_times(exp_beta=exp_beta, df=df, N=N, observe_time=observe_time, 
                      percentcensored=percentcensored,randomseed= randomseed+1,
                      lambda = lambda, distr = distr, rho_w = rho_w, drop_out=drop_out)
  return(df)
}


simulatedata_crossterms = function (N=1000, observe_time =10, 
                                percentcensored = 0.75,  
                                randomseed = 100,lambda = 0.1,
                                distr="Exp",  rho_w = 1, 
                                drop_out = 0.3){ 
  #simulate the data
  df = simulate_population(N, randomseed)
  # calculate betas
  exp_beta = xt_beta(df)
  # simulate censored and event times 
  df = df_event_times(exp_beta=exp_beta, df=df, N=N, observe_time=observe_time, 
                      percentcensored=percentcensored,randomseed= randomseed+1,
                      lambda = lambda, distr = distr, rho_w = rho_w, drop_out=drop_out)
  return(df)
}

simulatedata_lin_nonPH = function (N=1000, observe_time =10, 
                                percentcensored = 0.75,  
                                randomseed = 100,lambda = 0.1,
                                drop_out = 0.3){ 
  #simulate the data
  
  df = simulate_population(N, randomseed)
  # calculate betas
  exp_beta = exp(0.4*df$age + 1.0*df$bmi+1*df$hyp+0.5*df$sex)
  df["exp_beta"] = exp_beta
  df["shape_rho"] = ifelse(df$sex==1, 2.5, 0.5)
  #simulate event times
  {set.seed(randomseed); v <- runif(n=N)}
  event_time = (- log(v) / (lambda * df$exp_beta))^(1/df$shape_rho) #Weibull density

  # re-scale the time to have 1-percentcensored of events=1
  # by the observe_time
  final_time = quantile(event_time, 1-percentcensored)
  df$event_time = 0.001 + pmin(round(event_time/final_time*observe_time, 3),observe_time) #scale to observe_time
  # generate drop-out times for random drop_out % observations
  if (drop_out >0 ){
    set.seed(randomseed+1)
    randcentime = runif(round(N*drop_out,0), 0, observe_time)
    cens_obs = sample.int(nrow(df), round(N*drop_out,0))
    df[cens_obs, "cens_time"] = randcentime
    df[-cens_obs, "cens_time"] = observe_time
  } else {df[, "cens_time"] = observe_time}
  
  #final time and event definition 
  #event =1 if event time < cens_time and observe_time
  df$time = pmin(df$event_time, df$cens_time, observe_time)
  df$event = ifelse(df$event_time == df$time, 1, 0)
  
  return (df)
}

mean_weibull=function(lambda, rho_w){ 
  # h0(t)=??*rho_w*t^(rho_w???1)    S(t) = exp(-lambda*t^rho_w)
  # ET = (1/lambda)^(1/rho_w) * gamma(1+1/rho_w)
  return(lambda^(-1/rho_w)* gamma(1+1/rho_w))}

plot_weibull_h= function(lambda =0.027, rho_w=1.5, title = "Baseline hazard"){
  # h0(t)=lambda*rho_w*t^(rho_w???1)    S(t) = exp(-lambda*t^rho_w)
  t = seq(0,1,0.02)
  h=function(x){return(lambda*rho_w*x^(rho_w-1))}
  plot(t, h(t), main =title)
}

linear_beta = function(df){
  return(exp(0.4*df$age + 1.0*df$bmi + 0.7*df$hyp))}

nonlinear_beta = function(df){
  #BMI impact is 2 for very low and high levels, 1 for high/ low level, 0 for normal range
  bmi_beta = ifelse((df$bmi < -1.5)|(df$bmi>1.5), 2, ifelse( (df$bmi < -1)|(df$bmi>1), 1, 0))
  #Age impact is 1 for age>=55; linear age impact is also present, but is smaller than in linear simulation
  age_beta = ifelse( (df$age>=1), 1, 0)
  return(exp(bmi_beta + (df$hyp*0.7)+ df$age*0.2 + age_beta))
  }

xt_beta = function(df){
  #BMI impact is 2 for very low and high levels, 1 for high/ low level, 0 for normal range
  bmi_beta = ifelse( (df$bmi < -1.5)|(df$bmi>1.5), 2, 
                     ifelse( (df$bmi < -1)|(df$bmi>1), 1, 0))
  
  # hypertension x age interaction
  hyp_beta = ifelse((df$age>=1 & df$hyp == 1), 2, 
                    ifelse((df$age<1 & df$hyp == 1),1,0))
  #simulating event time
  return(exp(bmi_beta + hyp_beta+ df$age*0.2))
}

simulate_population = function(N=1000, randomseed=42){
 set.seed(randomseed)
 df = data.frame(age = round(runif(N, -1.73, 1.73),1),
                bmi = round(rnorm(N, 0, 1),1),  
                hyp = rbinom(N,1,0.20),
                sex = rbinom(N,1,0.5))
 return(df)}

df_event_times = function(exp_beta, df, N, observe_time, 
   percentcensored,randomseed,lambda,distr, rho_w, drop_out){
  #simulate event times
  if (distr =="Exp") {
    # Exponential time distribution, h0(t)=??=0.1  with shape ??>0 and scale ??>0
    {set.seed(randomseed);  event_time = rexp(N, lambda * exp_beta)}
  }else{
    # simulate Weibull distribution 
    # h0(t)=lambda*rho*t^(rho???1), shape rho_w>0,  scale lambda>0. 
    # If rho=1, it is exponential
    # rho>1 => upward sloping, rho<1 => downward, rho=1 constant (exponential)
    {set.seed(randomseed); v <- runif(n=N)}
    event_time = (- log(v) / (lambda * exp_beta))^(1/rho_w) #Weibull density
  }
  
  # re-scale the time to have 1-percentcensored of events=1
  # by the observe_time
  final_time = quantile(event_time, 1-percentcensored)
  df$event_time = 0.001 + pmin(round(event_time/final_time*observe_time, 3),observe_time) #scale to observe_time
  # generate drop-out times for random drop_out % observations
  if (drop_out >0 ){
    set.seed(randomseed+1)
    randcentime = runif(round(N*drop_out,0), 0, observe_time)
    cens_obs = sample.int(nrow(df), round(N*drop_out,0))
    df[cens_obs, "cens_time"] = randcentime
    df[-cens_obs, "cens_time"] = observe_time
  } else {df[, "cens_time"] = observe_time}
  
  #final time and event definition 
  #event =1 if event time < cens_time and observe_time
  df$time = pmin(df$event_time, df$cens_time, observe_time)
  df$event = ifelse(df$event_time == df$time, 1, 0)
  
  return (df)
}

