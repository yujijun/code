#输入数据的必备要求，样本的名称必须一样，按照各平台示例数据的要求上传数据，没有进行特征选择的数据
#接受参数传递####
argv <- commandArgs(TRUE)

#read the message of FileInfo####
Dirpath <- as.character(argv[1]) #which type of cancer?
setwd(Dirpath)
FileInfo <- read.table("COMSUC_FileInfo.txt", header = T, check.names = F)

#read the message of FileInfo####(temporary)
Filename <- paste0(c(1,2,3,4,5), ".txt")
Omicstype <- c("mRNA","miRNA","methylation","copy_number", "RPPA")
FileInfo <- data.frame(Filename=Filename, Omicstype=Omicstype)

#create the output matrix####
filename <- FileInfo$filename
error <- rep(NA, nrow(FileInfo))
warning <- rep(NA, nrow(FileInfo))
success <- rep(NA, nrow(FileInfo))
Output <- data.frame(Filename=Filename, Error=Error, Warning=Warning, success=success)

#读入数据校验函数####
DataCheck <- function(x){
  options(warn = 2)
  test_matrix <- try(read.table(x, header = T, stringsAsFactors = F, check.names = F, row.names = 1))
  if("try-error" %in% class(test_matrix))  #通过上述操作，就可以把信息两种错误可能情况全部爆出来（NULL 和 重复的行名）
  {
    cat("There is a Error in reading data!\n",test_matrix[1] )
    messages1 <- c("There is a Error in reading data.", test_matrix[1])
    #write.table(messages1, file = "error.txt", quote = F, sep = "\t", row.names = F, col.names = F)
  }else{
    messages1 <- NA
    class_all <- c(apply(test_matrix, 1, class), apply(test_matrix, 2, class))
    if("character" %in% class_all)
    {
      cat("Error:The number in the matrix is in character format. Please modify the number format and upload again.\n")
      messages2 <- "The number in the matrix is in character format. Please modify the number format and upload again."
      #write.table(messages2, file = "error.txt", quote = F, sep = "\t", row.names = F, col.names = F)
    }else{ 
        messages2 <- NA
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
          #create the function of NA judgement#
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
        
        #判断输出三种信息####
        #error#
        if(sum(is.na(c(messages1, messages2))) == 2){
          messages_error <- NA
        }else{
          output_messages_position <- which(!is.na(c(messages1, messages2)))
          messages_error <- c(messages1, messages2)[output_messages_position]
        }
        
        #warning#
        if(sum(is.na(c(messages3,messages4,messages5))) == 3){
          messages_warning <- NA
        }else{
          output_messages_position <- which(!is.na(c(messages3, messages4, messages5)))
          messages_warning <- c(messages3, messages4, messages5)[output_messages_position]
        }
        
        #success#
        if(sum(is.na(messages_error)) == 2){
          messages_success <- "success"
        }else{
          messages_success <- NA
        }
        
        #将结果导出####
        results <- list(messages_error=messages_error, messages_warning=messages_warning, messages_success=messages_success, test_matrix=test_matrix)
        return(results)
     }
  }
}

#read all of uploaded data####
UserdataList <- NULL
for(i in 1:nrow(FileInfo)){
  UserdataList[[i]] <- DataCheck(FileInfo[i,1])
  Output[i,2] <- UserdataList[[i]]$messages_error
  Output[i,3] <- UserdataList[[i]]$messages_warning
  Output[i,4] <- UserdataList[[i]]$messages_success
}
write.table(Output, file = "Output.txt", quote = F, sep = "\t", col.names=T, row.names = F)

#对数据进行样本对齐####
success_length <- length(which(Output[,4] == "success"))
if(success_length == nrow(FileInfo)){
  if(nrow(FileInfo) > 1){
    samples_align <- rownames(UserdataList[[1]]$test_matrix)
    for(i in 2:nrow(FileInfo)){
       samples_align <- interset(rownames(UserdataList[[i]]$test_matrix))
    }
  }
}

#输出align信息####
if(length(samples_align) == 0){
  cat("The sample cannot be aligned. Please check if the sample name of each group is consistent")
  messages_align <- "The sample cannot be aligned. Please check if the sample name of each group is consistent"
}else{
  cat("After alignment, there are ", length(samples_align), "samples left. Do you want to continue the following analysis?")
  messages_align <- paste0("After alignment, there are only ", length(samples_align), " samples left. Do you want to continue the following analysis?")
}
write.table(messages_align, file = "messages_align.txt", quote = F, sep = "\t", row.names = F, col.names = F)

#确认之后进行样本拼接####
