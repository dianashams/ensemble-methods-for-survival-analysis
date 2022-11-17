
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

simulatedata_linear = function (N=1000, observe_time =10, 
                                percentcensored = 0.0, randomseed = 100){
  set.seed(randomseed)
  df = data.frame(age = round(runif(N, -1.73, 1.73),1),
                  bmi = round(rnorm(N, 0, 1),1),  
                  hyp = rbinom(N,1,0.20),
                  gender = rbinom(N,1,0.5))
  
  df$event_time = 0.01+ round((rexp(N, 0.1*exp(0.4*df$age + 1.0*df$bmi + 0.7*df$hyp))),2)

  #add censored / drop-out observations
  if (percentcensored >0 ){
    randcentime = runif(round(N*percentcensored,0),0,0 + observe_time)
    cens_obs = sample.int(nrow(df), round(N*percentcensored,0))
    df[cens_obs, "cens_time"] = randcentime
    df[-cens_obs, "cens_time"] = observe_time
  } else {df[, "cens_time"] = observe_time}

  df$time = pmin(df$event_time, df$cens_time, observe_time)
  df$event = ifelse(df$event_time <= df$time, 1, 0); sum(df$event ==1)/N

  return (df)
}

simulatedata_nonlinear = function (N=1000, observe_time =10, percentcensored = 0.5, randomseed = 100){
  set.seed(randomseed)
  df = data.frame(age = round(runif(N, -1.73,1.73),1),
                  bmi = round(rnorm(N, 0, 1),1),  
                  hyp = rbinom(N,1,0.20),
                  gender = rbinom(N,1,0.5))
  
  #BMI impact is 2 for very low and high levels, 1 for high/ low level, 0 for normal range
  bmi_beta = ifelse((df$bmi < -1.5)|(df$bmi>1.5), 2, 
                     ifelse( (df$bmi < -1)|(df$bmi>1), 1, 0))
  #Age impact is 1 for age>=55; linear age impact is also present, but is smaller than in linear simulation
  age_beta = ifelse( (df$age>=1), 1, 0)
  #simulating event time
  df$event_time = 0.01+ round((rexp(N, 0.08*exp(bmi_beta + 
                        (df$hyp*0.7)+ df$age*0.2 + age_beta))),2)
  #add censored / drop-out observations
  if (percentcensored >0 ){
    randcentime = runif(round(N*percentcensored,0),0,0 + observe_time)
    cens_obs = sample.int(nrow(df), round(N*percentcensored,0))
    df[cens_obs, "cens_time"] = randcentime
    df[-cens_obs, "cens_time"] = observe_time
  } else {df[, "cens_time"] = observe_time}
  
  df$time = pmin(df$event_time, df$cens_time, observe_time)
  df$event = ifelse(df$event_time <= df$time, 1, 0); sum(df$event ==1)/N
  return (df)
}

simulatedata_crossterms = function (N=1000, 
                                    observe_time =10, 
                                    percentcensored = 0.0, 
                                    randomseed=100){
  set.seed(randomseed)
  df = data.frame(age = round(runif(N, -1.73,1.73),1),
                  bmi = round(rnorm(N, 0, 1),1),  
                  hyp = rbinom(N,1,0.2),
                  gender = rbinom(N,1,0.5))
  
  #BMI impact is 2 for very low and high levels, 1 for high/ low level, 0 for normal range
  bmi_beta = ifelse( (df$bmi < -1.5)|(df$bmi>1.5), 2, 
                     ifelse( (df$bmi < -1)|(df$bmi>1), 1, 0))
  
  # hypertension x age interaction
  hyp_beta = ifelse((df$age>=1 & df$hyp == 1), 2, 
                    ifelse((df$age<1 & df$hyp == 1),1,0))
  #simulating event time
  df$event_time = 0.01 + round((rexp(N, 0.01+ 0.07*exp(bmi_beta + 
                                    (hyp_beta)+ df$age*0.2))),2)
  
  #add censored / drop-out observations
  if (percentcensored >0 ){
    randcentime = runif(round(N*percentcensored,0),0,0 + observe_time)
    cens_obs = sample.int(nrow(df), round(N*percentcensored,0))
    df[cens_obs, "cens_time"] = randcentime
    df[-cens_obs, "cens_time"] = observe_time
  } else {df[, "cens_time"] = observe_time}
  
  df$time = pmin(df$event_time, df$cens_time, observe_time)
  df$event = ifelse(df$event_time <= df$time, 1, 0); sum(df$event ==1)/N
  
  return (df)
}

#################################################################

simulatedata_linear_0 = function (N=1000,  observe_time =10, 
                                  percentcensored = 0.0, randomseed = 100){
  #old versions
  set.seed(randomseed)
  df = data.frame(age = round(runif(N, 20,75),1),
                  bmi = round(rnorm(N, 26, 3),1),  
                  hyp = rbinom(N,1,0.10),
                  gender = rbinom(N,1,0.5))
  
  df$bmi_ = (df$bmi-26)/3
  df$age_ = (df$age - 50)/15
  #BMI impact is 1 for low and high levels, 0 in-between
  bmi_beta = 1
  #Age impact is 1 for age>=55; linear age impact is also present, but is smaller than in linear simulation
  age_beta = ifelse( (df$age>=55), 1, 0)
  
  df$event_time = 0.01+ round((rexp(N, 0.15*exp(1.0*df$bmi_ + 0.7*df$hyp + 0.4*df$age_))),2)
  #make sure event-time is not too small, fudge a bit if so
  #df$event_time = ifelse(df$event_time<0.5, round(runif(1,0.2, 0.5),1), df$event_time)  
  #add censored / drop-out observations
  if (percentcensored >0 ){
    randcentime = runif(round(N*percentcensored,0),0,0 + observe_time)
    cens_obs = sample.int(nrow(df), round(N*percentcensored,0))
    df[cens_obs, "cens_time"] = randcentime
    df[-cens_obs, "cens_time"] = observe_time
  } else {df[, "cens_time"] = observe_time}
  
  df$time = pmin(df$event_time, df$cens_time, observe_time)
  df$event = ifelse(df$event_time <= df$time, 1, 0)
  
  return (df)
}

simulatedata_nonlinear_0 = function(N=1000, observe_time =10, percentcensored = 0.1, randomseed = 100){
  #old versions
  set.seed(randomseed)
  df = data.frame(age = round(runif(N, 20,75),1),
                  bmi = round(rnorm(N, 26, 3),1),  
                  hyp = rbinom(N,1,0.10),
                  gender = rbinom(N,1,0.5))
  
  df$bmi_ = (df$bmi-26)/3
  #BMI impact is 1 for low and high levels, 0 in-between
  bmi_beta = ifelse( (df$bmi_ < -1)|(df$bmi_>1), 1, ifelse( (df$bmi_ < -1.75)|(df$bmi_>1.75), 2, 0))
  #Age impact is 1 for age>=55; linear age impact is also present, but is smaller than in linear simulation
  age_beta = ifelse( (df$age>=55), 1, 0)
  
  df$event_time = 0.01+ round((rexp(N,0.15*exp(bmi_beta + (df$hyp*0.7)+ (df$age-50)/15*0.4 + age_beta))),2)
  #make sure event-time is not too small, fudge a bit if so
  #df$event_time = ifelse(df$event_time<0.5, round(runif(1,0.2, 0.5),1), df$event_time)  
  #add censored / drop-out observations
  if (percentcensored >0 ){
    randcentime = runif(round(N*percentcensored,0),0,0 + observe_time)
    cens_obs = sample.int(nrow(df), round(N*percentcensored,0))
    df[cens_obs, "cens_time"] = randcentime
    df[-cens_obs, "cens_time"] = observe_time
  } else {df[, "cens_time"] = observe_time}
  
  
  df$time = pmin(df$event_time, df$cens_time, observe_time)
  df$event = ifelse(df$event_time <= df$time, 1, 0); sum(df$event ==1)/N
  
  df$age_ = (df$age - 50)/15
  return (df)
}

simulatedata_crossterms_0 = function (N=1000, leftcenstime=0, observe_time =10, percentcensored = 0.1, randomseed=100){
  #old versions
  set.seed(randomseed)
  df = data.frame(age = round(runif(N, 25,75),1),
                  bmi = round(rnorm(N, 26, 3),1),  
                  hyp = rbinom(N,1,0.25),
                  gender = rbinom(N,1,0.5))
  df$bmi_ = (df$bmi-26)/3
  #BMI impact is 1 for low and high levels, 0 in-between
  bmi_beta = ifelse( (df$bmi_ < -1)|(df$bmi_>1), 1, ifelse( (df$bmi_ < -1.75)|(df$bmi_>1.75), 2, 0))
  
  # hypertension x age interaction
  hyp_beta = ifelse((df$age>=50 & df$hyp == 1), 2, ifelse((df$age<50 & df$hyp == 1),1,0))
  
  df$event_time = 0.01 + round((rexp(N, 0.01+ 0.15*exp(bmi_beta + (hyp_beta)+ (df$age - 50)/15*0.2))),2)
  #make sure event-time is not too small, fudge a bit if so
  #df$event_time = ifelse(df$event_time<0.5, round(runif(1,0.2, 0.5),1), df$event_time)  
  
  #add censored / drop-out observations
  if (percentcensored >0 ){
    randcentime = runif(round(N*percentcensored,0),0,0 + observe_time)
    cens_obs = sample.int(nrow(df), round(N*percentcensored,0))
    df[cens_obs, "cens_time"] = randcentime
    df[-cens_obs, "cens_time"] = observe_time
  } else {df[, "cens_time"] = observe_time}
  
  df$time = pmin(df$event_time, df$cens_time, observe_time)
  df$event = ifelse(df$event_time <= df$time, 1, 0); sum(df$event ==1)/N
  
  df$age_ = (df$age - 50)/15
  return (df)
}

