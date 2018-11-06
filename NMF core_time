NMF_coretest <- function(){
  time1 <- Sys.time()
  for(i in 1:10){
    print(i)
    res <- nmf(rank_result, 2:10, nrun=3)
  }
  time2 <- Sys.time()
  cat("the P1 time =", difftime(time2,time1,units = "mins"), "\n")
  
  time3 <- Sys.time()
  for(i in 1:10){
    print(i)
    res <- nmf(rank_result, 2:10, nrun=3, .options="P4")
  }
  time4 <- Sys.time()
  cat("the P4 time =", difftime(time4,time3, units = "mins"), "\n")
}

NMF_coretest()
