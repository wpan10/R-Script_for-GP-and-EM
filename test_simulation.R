library('dplyr')
source('simulation_4.R')

test_simulation<- function(tbl, 
                           testmatrix,
                           bb){
  num_id <- n_distinct(tbl$ID)
  df <- data.frame()
  
  for (i in 1:num_id){
    type <- tbl %>% filter(ID == i) %>% select(subtype)
    group <- type[["subtype"]][1]
    
    years <- tbl %>% filter(ID == i) %>% select( years_seen)
    years_seen <- years[["years_seen"]]
    
    pfvc_tmp <- tbl %>% filter (ID == i) %>% select(pfvc_gp)
    pfvc_gp <- pfvc_tmp[["pfvc_gp"]]
    
    predict_pfvc <- design(years_seen,bb) %*% testmatrix[,group]
   
    Group <- rep(group, length(years_seen))
    Index <- rep(i,length(years_seen))
    tmp_df <- data.frame(cbind(Index,years_seen,pfvc_gp,predict_pfvc,Group))
    df <- rbind(df,tmp_df)
  }
  
  colnames(df) = c('ID','years_seen','pfvc_gp','predict_pfvc_gp','subtype')
  return (df)   
}

