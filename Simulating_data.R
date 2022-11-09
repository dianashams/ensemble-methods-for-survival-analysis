
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

simulatedata_nonlinear = function (N=1000, observe_time =10, percentcensored = 0.1, randomseed = 100){
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
                                                  (df$hyp*0.7)+ df$age*0.2 
                                                + age_beta))),2)
  df$observe_time = observe_time
  
  #3_ add censored observations with a shorter observation time (drop-outs)
  # the time is randomly drawn uniformly from observe_time/20 to 
  # observe_time ('/20' to exclude very short observations)
  df$early_censored = 0 #marker if an observation dropped out early (1) 
  if (percentcensored >0 ){
    #assume that nobody drops out in the first 1/20th of the observation time 
    randcentime = runif(round(N*percentcensored,0),0,0 + observe_time)
    cens_obs = sample.int(nrow(df), round(N*percentcensored,0))
    # censored time is the end of observation in this simulation
    df[cens_obs, "cens_time"] = randcentime
    df[-cens_obs, "cens_time"] = observe_time
  } else {df[, "cens_time"] = observe_time}
  
  #4_ defining the outcome and time 
  # time is the first from event, censoring, or end of observation
  df$time = pmin(df$event_time, df$cens_time, df$observe_time)
  
  df$event = ifelse(df$event_time <= df$time, 1, 0)
  
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
  df$event = ifelse(df$event_time <= df$time, 1, 0)
  
  return (df)
}
