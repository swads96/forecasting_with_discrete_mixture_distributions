MakeDist <- function(distsdf){
  
  distdf <- distsdf[distsdf[,1] != 'Lst',]
  tdist <- distsdf[distsdf[,1] == 'Lst',]
  
  fun_dist <- 
    apply(distdf, FUN=function(x) {
      paste('distr::',x[1], '(', 
            ifelse(!is.na(x[2]) & (!is.na(x[3]) | !is.na(x[4])),
                   paste(x[2],',',sep=''),
                   ifelse(!is.na(x[2]) & is.na(x[3]) & is.na(x[4]),
                          x[2], '')), 
            ifelse(!is.na(x[3]) & !is.na(x[4]),
                   paste(x[3],',',sep=''),
                   ifelse(!is.na(x[3]) & is.na(x[4]), x[3],'')), 
            ifelse(!is.na(x[4]),x[4],''), ')',sep='')
    }, MARGIN = 1
    )
  
  fun_tdist <- apply(tdist, FUN=function(x) {
    paste0('distr::Td(',x[4],')*',x[3], '+', x[2])
  }, MARGIN = 1
  )
  
  dist_args <- paste(fun_dist, collapse=',',sep='')
  tdist_args <- paste0(fun_tdist,collapse=',')
  args <- ifelse(tdist_args!='',paste(dist_args,tdist_args,sep=','),dist_args)
  
  weights <- c(distdf[,5],tdist[,5])
  mixString <- paste('distr::UnivarMixingDistribution(',
                     args,',mixCoeff=weights)',sep='')
  mixDist <- eval(parse(text=mixString))
  
  return(mixDist)
}