#----------------------Get extension of the input files-----------------
getExtension <- function(file){ 
  ex <- strsplit(basename(file), split="\\.")[[1]]
  return(ex[-1])
} 

#--------------------Train the mean-variance relationship from the user provided replicates-----------
VSS_train <- function(args)
{
  start_time <- Sys.time()
  path <- paste(getwd(),"/",args[4],sep="")
  dir.create(path)
  #------------ Getting the interval intersections of two input replicates--------------
  if(getExtension(args[2])=="bed" | getExtension(args[2])=="bedGraph"){
    system(paste("bedops --partition " ,args[2] ," ", args[3]," >", path, "/intersect_intervals.bed",sep=""))
    system(paste("bedtools closest -a ", path,"/intersect_intervals.bed -b ", args[2]," >", path, "/rep1_scores_closest.bed",sep=""))
    system(paste("bedtools closest -a ", path,"/intersect_intervals.bed -b ", args[3]," >", path, "/rep2_scores_closest.bed",sep=""))
  }else if(getExtension(args[2])=="bigWig"){
    system(paste("./bigWigToBedGraph " ,args[2] ," ",path,"/rep1_train_model.bedGraph",sep=""))
    system(paste("./bigWigToBedGraph " ,args[3] ," ",path,"/rep2_train_model.bedGraph",sep=""))
    system(paste("bedops --partition ",path, "/rep1_train_model.bedGraph " ,path, "/rep2_train_model.bedGraph >" , path, "/intersect_intervals.bed",sep = ""))
    system(paste("bedtools closest -a ", path,"/intersect_intervals.bed -b ", path,"/rep1_train_model.bedGraph >", path, "/rep1_scores_closest.bed",sep=""))
    system(paste("bedtools closest -a ", path,"/intersect_intervals.bed -b ", path,"/rep2_train_model.bedGraph >", path, "/rep2_scores_closest.bed",sep=""))
    system(paste("rm ",path,"/rep1_train_model.bedGraph",sep=""))
    system(paste("rm ",path,"/rep2_train_model.bedGraph",sep=""))
  }else if(getExtension(args[2])=="bam"){
    system(paste("bedtools genomecov -ibam " ,args[2], " -bga >",path, "/rep1_train_model.bedGraph",sep=""))
    system(paste("bedtools genomecov -ibam " ,args[3], " -bga >",path, "/rep2_train_model.bedGraph",sep=""))
    system(paste("bedops --partition ",path,"/rep1_train_model.bedGraph ",path,"/rep2_train_model.bedGraph >", path, "/intersect_intervals.bed",sep=""))
    system(paste("bedtools closest -a ",path,"/intersect_intervals.bed -b ",path,"/rep1_train_model.bedGraph >",path,"/rep1_scores_closest.bed",sep=""))
    system(paste("bedtools closest -a ",path,"/intersect_intervals.bed -b ",path,"/rep2_train_model.bedGraph >",path,"/rep2_scores_closest.bed",sep=""))
    system(paste("rm ",path,"/rep1_train_model.bedGraph",sep=""))
    system(paste("rm ",path,"/rep2_train_model.bedGraph",sep=""))
  }
  
  
  
  intersect_intervals <- fread(paste(path,"/intersect_intervals.bed",sep=""))
  rep1_scores_closest<- fread(paste(path,"/rep1_scores_closest.bed",sep=""))
  rep2_scores_closest<- fread(paste(path,"/rep2_scores_closest.bed",sep=""))
  #system(paste("rm ",path,"/intersect_intervals.bed",sep = ""))
  #system(paste("rm ",path,"/rep1_scores_closest.bed",sep = ""))
  #system(paste("rm ",path,"/rep2_scores_closest.bed",sep = ""))
  
  rep1_scores_closest[,c(1,4:6)] <- NULL
  rep1_scores_closest <- data.frame(rep1_scores_closest)
  rep1_scores_closest[,4] <- rep1_scores_closest[,2]-rep1_scores_closest[,1]
  rep1_scores_closest[,1:2] <- NULL
  rep2_scores_closest[,c(1,4:6)] <- NULL
  rep2_scores_closest <- data.frame(rep2_scores_closest)
  rep2_scores_closest[,4] <- rep2_scores_closest[,2]-rep2_scores_closest[,1]
  rep2_scores_closest[,1:2] <- NULL
  a <- rbind(rep1_scores_closest,rep2_scores_closest)
  b <- rbind(rep2_scores_closest,rep1_scores_closest)
  scores <- cbind(a,b)
  names(scores)[1:4] <- c("Base_score","Base_interval","Aux_score","Aux_interval")
  ordered_scores <- scores[order(scores$Base_score), ]
  L <- sum(as.numeric(ordered_scores$Base_interval))
  zero_signals=subset(ordered_scores,ordered_scores$Base_score==0)
  #non_zero_signals=subset(ordered_scores,ordered_scores$Base_score!=0)
  threshold <- 0.25*L
  #--------------Getting mean-variance--------------------------------
  if(sum(as.numeric(zero_signals$Base_interval))>=threshold){
    #---------mean-var curve smoothing parameters ---------
    beta_values=1000
    bin_sizes=100000
    alpha_values=2^(1/as.integer(beta_values))
    width_values=ceiling(log(1/0.01)/log(alpha_values))
    #------------------------------------------------------
    mean_ordered <- sum(as.numeric(ordered_scores$Base_score*ordered_scores$Base_interval))/sum(as.numeric(ordered_scores$Base_interval))
    less_than_mean <- subset(ordered_scores,ordered_scores$Base_score<=mean_ordered)
    higher_than_mean <- subset(ordered_scores,ordered_scores$Base_score>mean_ordered)
    Mean_1_over_sigma <- Weighted_Mean_1_over_sigma(width_values,as.integer(bin_sizes),alpha_values,higher_than_mean)
    Mean_1_over_sigma <- data.frame(Mean_1_over_sigma)
    mean_less_than_mean <- sum(as.numeric(less_than_mean$Aux_score*less_than_mean$Aux_interval))/sum(as.numeric(less_than_mean$Aux_interval))
    Mean_1_over_sigma[(nrow(Mean_1_over_sigma)+1),1] <- mean_less_than_mean
    sd_less_than_mean <- (sqrt(sum(as.numeric(((less_than_mean$Aux_score-mean_less_than_mean)^2)*less_than_mean$Aux_interval))/sum(as.numeric(less_than_mean$Aux_interval))))
    Mean_1_over_sigma[(nrow(Mean_1_over_sigma)),2] <- 1/sd_less_than_mean
    Mean_1_over_sigma[(nrow(Mean_1_over_sigma)+1),1] <- 0
    Mean_1_over_sigma[(nrow(Mean_1_over_sigma)),2] <- 1/sd_less_than_mean
    Mean_1_over_sigma <- data.frame(Mean_1_over_sigma)
    
  }else if(sum(as.numeric(zero_signals$Base_interval))<threshold){
    #---------mean-var curve smoothing parameters ---------
    beta_values=1e+07
    bin_sizes=1000
    alpha_values=2^(1/as.integer(beta_values))
    width_values=ceiling(log(1/0.01)/log(alpha_values))
    Mean_1_over_sigma <- Weighted_Mean_1_over_sigma(width_values,as.integer(bin_sizes),alpha_values,ordered_scores)
    Mean_1_over_sigma <- data.frame(Mean_1_over_sigma)
  }
  #path <- paste(getwd(),"/Trained_model",sep="")
  #dir.create(path)
  save(Mean_1_over_sigma,file=paste(path,"/Trained_mean_variance_model.Rdata",sep=""))
  write.table(Mean_1_over_sigma, file=paste(path,"/Trained_mean_variance_model.bedGraph",sep=""), quote=F, sep="\t", row.names=F, col.names=F)
  write.table(Mean_1_over_sigma, file=paste(path,"/Trained_mean_variance_model.bed",sep=""), quote=F, sep="\t", row.names=F, col.names=F)
  end_time <- Sys.time()
  print(end_time - start_time)
}
#------------------Identifying the mean-variance relationship-------------
Weighted_Mean_1_over_sigma<-function(distance,bin,alpha,ordered_scores)
{
  ordered_scores[,1:2] <- NULL
  ordered_scores[,3] <- ordered_scores$Aux_score*ordered_scores$Aux_interval
  ordered_scores[,4] <- (ordered_scores$Aux_score^2)*ordered_scores$Aux_interval
  names(ordered_scores)[3:4] <- c("sum","squared")
  L=ceiling(sum(ordered_scores$Aux_interval)/bin)
  l=bin
  temp <- ordered_scores
  temp[,5] <- floor(temp$Aux_interval/l)
  names(temp)[5] <- c("portion")
  Y=matrix(0,nrow=L,ncol = 1)
  a=matrix(0,nrow=L,ncol = 1)
  counter_bin <- 1
  z=matrix(nrow=((2*distance)+1),ncol = 1)
  weighted_z=matrix(nrow=((2*distance)+1),ncol=1)
  n=(2*distance)+1
  mean_1_over_sigma=matrix(nrow=L,ncol=2)
  y_value <- rep(0,nrow(Y))
  y_bins_counter <- rep(l,nrow(Y))
  a_value <- rep(0,nrow(Y))
  temp_score <- c()
  temp_score <- temp$Aux_score
  temp_interval <- c()
  temp_interval <- temp$Aux_interval
  temp_sum <- c()
  temp_sum <- temp$sum
  temp_squared <- c()
  temp_squared <- temp$squared
  temp_portion <- c()
  temp_portion <- temp$portion
  
for(i in 1:length(ordered_scores$Aux_score))
{
    main_interval <- temp_interval[i]
    if(temp_portion[i]==0 & temp_interval[i]<=y_bins_counter[counter_bin]){
      y_value[counter_bin] <- y_value[counter_bin] + temp_sum[i]
      a_value[counter_bin] <- a_value[counter_bin] + temp_squared[i]
      y_bins_counter[counter_bin] <- y_bins_counter[counter_bin] - temp_interval[i]
      temp_interval[i] <- 0
      if(y_bins_counter[counter_bin]==0){
        counter_bin <- counter_bin+1
      }
    }else if(temp_portion[i]==0 & temp_interval[i]>y_bins_counter[counter_bin]){
      #--------Fill the previous bin----------
      y_value[counter_bin] <- y_value[counter_bin] + (y_bins_counter[counter_bin]/main_interval)*temp_sum[i]
      a_value[counter_bin] <- a_value[counter_bin] + (y_bins_counter[counter_bin]/main_interval)*temp_squared[i]
      temp_interval[i] <- temp_interval[i] - y_bins_counter[counter_bin]
      y_bins_counter[counter_bin] <- 0
      counter_bin <- counter_bin+1
      
      #--------Fill the current bin---------
      y_value[counter_bin] <- y_value[counter_bin] + (temp_interval[i]/main_interval)*temp_sum[i]
      a_value[counter_bin] <- a_value[counter_bin] + (temp_interval[i]/main_interval)*temp_squared[i]
      y_bins_counter[counter_bin] <- y_bins_counter[counter_bin] - temp_interval[i]
      temp_interval[i] <- 0
      if(y_bins_counter[counter_bin]==0){
        counter_bin <- counter_bin+1
        }
    }else if(temp_portion[i]!=0){
      #--------Fill the previous bin---------
      left <- y_bins_counter[counter_bin]
      y_value[counter_bin] <- y_value[counter_bin] + (left/main_interval)*temp_sum[i]
      a_value[counter_bin] <- a_value[counter_bin] + (left/main_interval)*temp_squared[i]
      y_bins_counter[counter_bin] <- 0
      temp_interval[i] <- temp_interval[i] - left
      counter_bin <- counter_bin + 1
      
      #--------Fill the current bin----------
      num_score <- floor(temp_interval[i]/l)
      counter <- counter_bin
      first_bin <- counter_bin
      for(cc in counter:(num_score+first_bin-1))
        {
        y_value[cc] <- (l/main_interval)*temp_sum[i]
        a_value[cc] <- (l/main_interval)*temp_squared[i]
        y_bins_counter[counter_bin] <- 0
        counter_bin <- counter_bin+1
        temp_interval[i] <- temp_interval[i] - l
        }
      
      #--------Fill the next bin --------------
      num_score <- floor(temp_interval[i]/l)
      y_value[counter_bin] <- y_value[counter_bin] + (temp_interval[i]/main_interval)*temp_sum[i]
      a_value[counter_bin] <- a_value[counter_bin] + (temp_interval[i]/main_interval)*temp_squared[i]
      y_bins_counter[counter_bin] <- y_bins_counter[counter_bin] - temp_interval[i]
      temp_interval[i] <- 0
    }
}
  Y[,1] <- y_value
  a[,1] <- a_value
  ##########################
  for (j in 1:n)
  {
    z[j,1] <- (1/(alpha^(abs(distance+1-j)*l)))
    weighted_z[j,1] <- l*(1/(alpha^(abs(distance+1-j)*l)))
  }
  ##########################
  for(i in 1:L)
  {
    output=Preprocess_Weighted_Mean_1_over_sigma(Y,z,weighted_z,a,i,distance)
    mean_1_over_sigma[i,1]=output[1]
    mean_1_over_sigma[i,2]=output[2]
  }
  mean_1_over_sigma
}
#--------------------- Preprocessing for indetifying the mean-varaince relationship-------------------
Preprocess_Weighted_Mean_1_over_sigma<- function(weighted_sum,weights,sum_of_weights,sum_square,center_bin,width)
{
  epsilon=0.0001
  y=weighted_sum
  W=weighted_sum
  L=length(y)
  z=weights
  weighted_z=sum_of_weights
  k=center_bin
  w=width
  a=sum_square
  A=sum_square
  
  lower=max(1,(k-w))
  upper=min(L,(k+w))
  xx=upper-lower+1
  center_z=ceiling(length(z)/2)
  lower_z=(center_z-(k-lower))
  counter=lower_z
  upper_z=(center_z+(upper-k))
  for(i in lower:upper)
  {
    y[i,1]=y[i,1]*z[counter,1]
    counter=counter+1
  }
  weighted_mean=sum(y[lower:upper])/sum(weighted_z[lower_z:upper_z])
  counter=lower_z
  for(i in lower:upper)
  {
    a[i,1]=a[i,1]*z[counter,1]
    counter=counter+1
  }
  weighted_var=sum(a[lower:upper])/sum(weighted_z[lower_z:upper_z])
  weighted_var=weighted_var -(weighted_mean^2)
  weighted_var=weighted_var+epsilon
  mean_1_over_sigma_values <- c(weighted_mean,1/sqrt(weighted_var))
  mean_1_over_sigma_values
}

#-------------------- Transform the varaince-instabilized signals to variance-stabilized signals--------
VSS_transform <- function(args)
{
  path <- paste(getwd(),"/",args[4],sep="")
  dir.create(path)
  #------------ Getting the interval intersections of two input replicates--------------
  if(getExtension(args[2])=="bed" | getExtension(args[2])=="bedGraph"){
    #instabilized_replicate_signals <- read.table(args[2])
    instabilized_replicate_signals <- fread(args[2])
    instabilized_replicate_signals <- data.frame(instabilized_replicate_signals)
  }else if(getExtension(args[2])=="bigWig"){
    system(paste("./bigWigToBedGraph " ,args[2] ," ",path,"/instabilized_replicate_signals.bedGraph",sep=""))
    #instabilized_replicate_signals <- read.table(paste(path,"/instabilized_replicate_signals.bedGraph",sep=""))
    instabilized_replicate_signals <- fread(paste(path,"/instabilized_replicate_signals.bedGraph",sep=""))
    instabilized_replicate_signals <- data.frame(instabilized_replicate_signals)
    system(paste("rm ",path,"/instabilized_replicate_signals.bedGraph",sep=""))
  }else if(getExtension(args[2])=="bam"){
    system(paste("bedtools genomecov -ibam " ,args[2], " -bga >",path, "/instabilized_replicate_signals.bedGraph",sep=""))
    #instabilized_replicate_signals <- read.table(paste(path,"/instabilized_replicate_signals.bedGraph",sep=""))
    instabilized_replicate_signals <- fread(paste(path,"/instabilized_replicate_signals.bedGraph",sep=""))
    instabilized_replicate_signals <- data.frame(instabilized_replicate_signals)
    system(paste("rm ",path,"/instabilized_replicate_signals.bedGraph",sep=""))
  }
  
  
  chromosomes <- unique(instabilized_replicate_signals$V1)
  load(paste(getwd(),"/",args[3],"/Trained_mean_variance_model.Rdata",sep=""))
  Replicate1_scores=c()
  for(chr in chromosomes)
  {
    Mean_1_over_sigma <- data.frame(Mean_1_over_sigma)
    test_rep1<- subset(instabilized_replicate_signals, instabilized_replicate_signals[,1]==chr)
    replicate1_scores <- calculating_vss_signals(Mean_1_over_sigma,test_rep1)
    Replicate1_scores <- rbind(Replicate1_scores,replicate1_scores)
  }
  
  write.table(Replicate1_scores, file=paste(path,"/Variance_stabilized_signals.bedGraph",sep=""), quote=F, sep="\t", row.names=F, col.names=F)
  write.table(Replicate1_scores, file=paste(path,"/Variance_stabilized_signals.bed",sep=""), quote=F, sep="\t", row.names=F, col.names=F)
  
}

#------------------------Calculating variance-stabilized signals----------
calculating_vss_signals <- function(Mean_1_over_sigma,rep1)
{
  
  replicate1_scores <- rep1
  replicate1_scores <- replicate1_scores[order(replicate1_scores[,4]),]
  unstab_signals <- c()
  unstab_signals <- replicate1_scores[,4]
  unique_scores <- unique(replicate1_scores[,4])
  
  Mean_vs_sigma <- Mean_1_over_sigma
  Mean_vs_sigma <- data.frame(Mean_vs_sigma)
  Mean_vs_sigma[,1] <- Mean_1_over_sigma[,1]
  Mean_vs_sigma[,2] <- 1/Mean_1_over_sigma[,2]
  
  #-----------spline -----------
  fit_1 <- smooth.spline(Mean_vs_sigma[,1], (Mean_vs_sigma[,2]), cv=TRUE)
  predicted_sigma <- stats:::predict.smooth.spline(fit_1,replicate1_scores[,4])$y
  
  #-----------------------
  predicted_score=matrix(nrow=length(predicted_sigma),ncol=2)
  predicted_score[,1]=replicate1_scores[,4]
  predicted_score[,2]=predicted_sigma
  
  if(min(predicted_score[,2])<0){
    less <- subset(predicted_score,predicted_score[,1]<min(Mean_vs_sigma[,1]))
    if(min(less[,2])<0){
      less[,2] <- min(Mean_vs_sigma[,2])
    }
    
    middle <- subset(predicted_score,predicted_score[,1]>=min(Mean_vs_sigma[,1]) & predicted_score[,1]<=max(Mean_vs_sigma[,1]))
    high <- subset(predicted_score,predicted_score[,1]>max(Mean_vs_sigma[,1]))
    if(min(high[,2])<0){
      high[,2] <- max(Mean_vs_sigma[,2])
    }
    
    predicted_score <- rbind(less,middle,high)
  }
  predicted_score[predicted_score[,2]<0,2] <- mean(Mean_vs_sigma[,2])
  predicted_score[,2]=1/predicted_score[,2]
  
  x=predicted_score[,1]
  y=predicted_score[,2]
  ordered_data <- predicted_score[order(x,y),]
  ordered_data=ordered_data[!duplicated(ordered_data),]
  ordered_data=data.frame(ordered_data)
  stab_scores <- c()
  for(i in 1:length(unique_scores))
  {
    scores_list <- which(unstab_signals==unique_scores[i])
    start <- scores_list[1]
    end <- scores_list[length(scores_list)]
    sub <- subset(ordered_data, ordered_data[,1]<=unique_scores[i]) #recently
    score <-trapz(sub[,1],sub[,2])
    stab_scores[start:end] <- score
  }
  replicate1_scores[,5] <- stab_scores
  #replicate1_scores <- replicate1_scores[order(as.numeric(row.names(replicate1_scores))),]
  replicate1_scores<- replicate1_scores[order(replicate1_scores[,2]),]
  replicate1_scores[,4] <- NULL
  replicate1_scores
}

#-------------------Getting inputs ------------------
args = commandArgs(trailingOnly=TRUE)
if(args[1]=="train"){
  library(data.table)
  VSS_train(args)
}else if(args[1]=="transform"){
  library(pracma)
  library(data.table)
  VSS_transform(args)
}
#args=c("train","rep1.bed","rep2.bed","traindir")
#args=c("train","whole_rep1.bedGraph","whole_rep2.bedGraph","traindir")
#args=c("transform","ENCFF231JPM.bam","gm","gm2")
#args=c("train","rep1.bed","traindir","transformdir")
#Rscript Optimized_VSS.R train rep1.bed rep2.bed traindir
#Rscript Optimized_VSS.R transform rep1.bed traindir transformdir

#Rscript Optimized_VSS.R train whole_rep1.bedGraph whole_rep2.bedGraph traindir
#Rscript Optimized_VSS.R transform whole_rep1.bedGraph traindir transformdir


# Rscript Optimized_VSS.R train raw_rep1.bam raw_rep1.bam rawtraindir
# Rscript Optimized_VSS.R transform raw_rep1.bam rawtraindir rawtransformdir
# 
# Rscript Optimized_VSS.R train fold_rep1.bigWig fold_rep1.bigWig foldtraindir
# Rscript Optimized_VSS.R transform fold_rep1.bam foldtraindir foldtransformdir
# 
# Rscript Optimized_VSS.R train pval_rep1.bigWig pval_rep1.bigWig pvaltraindir
# Rscript Optimized_VSS.R transform pval_rep1.bam pvaltraindir pvaltransformdir
# 
