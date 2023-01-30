
library(stringr)

qnan = function(column){
  z = is.na(column) | (column==NaN)
  return (sum(z, na.rm=1)) }

notin = function(datafr, colstocheck){
  z = colstocheck %in% names(datafr)
  print(colstocheck[!z])
}

normalize_v = function (dfc){ return((dfc - mean(dfc, na.rm=1))/sd(dfc, na.rm=1))}

missing_values = function(data, colstocheck){
  df_missing= data.frame()
  for (n in colstocheck) { df_missing[n, 1]= qnan(data[n])}  #missing 418 B_alcohol,18 B_dep_bi, 27 B_smokstatus
  names(df_missing) = c("N missing")
  df_missing["N present"] = dim(data)[1] - df_missing["N missing"]
  df_missing["% missing"] = round(df_missing["N missing"]/dim(data)[1]*100,2)
  df_missing["N total"] = dim(data)[1]
  return (df_missing)
}

Clean_String <- function(string){
  
  # Remove everything that is not a number or letter (may want to keep more 
  # stuff in your actual analyses). 
  temp <- stringr::str_replace_all(temp,",", " ")
  temp <- stringr::str_replace_all(temp,"'", " ")
  temp <- stringr::str_replace_all(temp,'["]', '')
  # Shrink down to just one white space
  temp <- stringr::str_replace_all(temp,"[\\s]+", " ")
  # Split it
  temp <- stringr::str_split(temp, " ")[[1]]
  # Get rid of trailing "" if necessary
  indexes <- which(temp == "")
  if(length(indexes) > 0){
    temp <- temp[-indexes]
  } 
  return(temp)
}
