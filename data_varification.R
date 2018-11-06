x <- "test_matirx.txt"
data_varification <- function(x){
  options(warn = 2)
  test_matrix <- try(read.table(x, header = T, stringsAsFactors = F, check.names = F, row.names = 1))
  if("try-error" %in% class(test_matrix))  #通过上述操作，就可以把信息两种错误可能情况全部爆出来（NULL 和 重复的行名）
  {
    cat("There is a Error in reading data!\n",test_matrix[1] )
    messages1 <- c("There is a Error in reading data.", test_matrix[1])
    write.table(messages1, file = "error.txt", quote = F, sep = "\t", row.names = F, col.names = F)
  }else{
    class_all <- c(apply(test_matrix, 1, class), apply(test_matrix, 2, class))
    if("character" %in% class_all)
    {
      cat("Error:The number in the matrix is in character format. Please modify the number format and upload again.\n")
      messages2 <- "The number in the matrix is in character format. Please modify the number format and upload again."
      write.table(messages2, file = "error.txt", quote = F, sep = "\t", row.names = F, col.names = F)
    }else 
    {
      samples_original <- colnames(test_matrix)
      if(0 %in% apply(test_matrix, 1, sum))  #判断是否有全0 unit行
      {
        zero_judge <- function(x){
          return(length(which(x==0)))
        }
        zero_position <- as.vector(which(apply(test_matrix, 1, zero_judge) == ncol(test_matrix)))
        test_matrix <- test_matrix[-zero_position, ]
        cat("Warning: Several full zero rows in the matrix have been deleted. \n")
        messages3 <- "Warning: Several full zero rows in the matrix have been deleted"
      }else{
        messages3 <- NA
      }
      
      if(0 %in% apply(test_matrix, 2, sum)) #判断是否有全0 unit列
      {
        zero_judge <- function(x){
          return(length(which(x==0)))
        }
        zero_position <- as.vector(which(apply(test_matrix, 2, zero_judge) == nrow(test_matrix)))
        test_matrix <- test_matrix[-zero_position, ]
        cat("Warning: Several full zero cols in the matrix have been deleted. \n")
        messages4 <- "Warning: Several full zero cols in the matrix have been deleted"
      }else{
        messages4 <- NA 
      }
      
      if(sum(is.na(test_matrix)))  #判断是否有NA数据
      {
        #RPPA_processing####
        NA_func <- function(x) 
        {
          #判断数据中的NA值缺失的比例，如果其中基因的缺失比例多于20%，则直接删除
          row_NANum<-which(rowMeans(is.na(x))>0.2)
          if(length(row_NANum)!=0){   
            x<-x[-row_NANum,]
          }
          col_NANum<-which(colMeans(is.na(x))>0.2)
          if(length(col_NANum)!=0){   
            x<-x[,-col_NANum]
          }
          #将处理之后的矩阵用均值填充之后计算
          x_STAD_meanfill <- function(a)
          {
            a[which(is.na(a))] <- mean(as.numeric(a), na.rm = T)
            return(a)
          }
          x<- apply(x, 1, x_STAD_meanfill)
          x<- as.data.frame(t(x))
          return(x)
        }
        test_matrix <- NA_func(test_matrix)
        samples_final <- colnames(test_matrix)
        Removed_samples <- setdiff(samples_original, samples_final)
        cat("Warning: There is some NA in the matrix, We delete or fill in the mean data for the corresponding NA data. \n")
        messages5 <- "Warning: There is some NA in the matrix, We delete or fill in the mean data for the corresponding NA data"
      }else{
        messages5 <- NA
      }
      
      if(sum(is.na(c(messages3,messages4,messages5)) == 3))
      {
        cat("Congratulations, there's no problem with the data")
        messages_success <- "Congratulations, there's no problem with the data"
        write.table(messages_success, file = "success.txt", quote = F, sep = "\t", row.names = F, col.names = F)
      }else{
        output_messages_position <- which(!is.na(c(messages3, messages4, messages5)))
        messages_all <- c(messages3, messages4, messages5)[output_messages_position]
        write.table(messages_all, file = "Warning.txt", quote = F, sep = "\t", row.names = F, col.names = F)
      }
    }
  }
}
