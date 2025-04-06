## Add recruits to current population

isEmpty <- function(x) {
  return(length(x)==0)
}


new.dp.recruits = function(rec) {
  recheight=runif(rec, 1, 2)
  
  newRpop = data.table(
    ColonyNo = 1:rec,
    Age = 1,
    Height = 0,
    stringsAsFactors = FALSE
  )
  
  newRpop$Height<-recheight
  newRpop[isEmpty(Height), Height := 0]
  
  newRpop=newRpop[ColonyNo!=0]
  
  return(newRpop)
}
