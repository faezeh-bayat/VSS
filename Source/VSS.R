##  Load these libraries 
library(data.table)
library(argparse)
library(bigmemory)
library(pracma)


## Gets extension of the input replicates (.bam, .bed, .bedGraph, .bigWig)
getExtension <- function(file){ 
  
  ex <- strsplit(basename(file), split="\\.")[[1]]
  return(ex[-1])
} 


## Gets the name of the replicate
getBasename <- function(file){
  
  ex <- strsplit(basename(file), split="\\.")[[1]]
  return(ex[1])
} 


## Parses input arguments for "load_inputs" and "train" queries (when using VSS tool)
parse_arguments <- function(){

  parser <- ArgumentParser()
  parser$add_argument("vss_function", type="character", nargs="+", help="Functions for load/train/transform data")
  parser$add_argument("tas1", type="character", nargs="+", help="Path for first TAGALIGN/Replicate1 file")
  parser$add_argument("--fraglen1", type="integer", help="Replicate 1 Fragment length.")
  parser$add_argument("tas2",type="character", nargs="+", help="Path for second TAGALIGN/Replicate2 file")
  parser$add_argument("--fraglen2", type="integer",help="Replicate 2 Fragment length.")
  parser$add_argument("--chrsz", type="character", help="2-col chromosome sizes file.")
  parser$add_argument("--gensz", type="character", help="Genome size (sum of entries in 2nd column of chr. sizes file, or hs for human, ms for mouse).")
  parser$add_argument("--signal", default="", type="character", help="Output signal (raw/fc/pval)")
  parser$add_argument("--inputdir", default="", type="character", help="Directory for loading data.")
  parser$add_argument("--traindir", default="", type="character", help="Directory for training the model.")
  parser$add_argument("--transformdir", default="", type="character", help="Directory for transformed signals.")
  args <- parser$parse_args()
  return(args)
}


## Parses input arguments for "transform" query (when using VSS tool)
parse_arguments_transform <- function(){
  
  parser <- ArgumentParser()
  parser$add_argument("vss_function", type="character", nargs="+", help="Functions for load/train/transform data")
  parser$add_argument("tas1", type="character", nargs="+", help="Path for first TAGALIGN/Replicate1 file")
  parser$add_argument("--signal", default="", type="character", help="Output signal (raw/fc/pval)")
  parser$add_argument("--inputdir", default="", type="character", help="Directory for loading data.")
  parser$add_argument("--traindir", default="", type="character", help="Directory for training the model.")
  parser$add_argument("--transformdir", default="", type="character", help="Directory for transformed signals.")
  args <- parser$parse_args()
  return(args)
}


###################################
## load_inputs
## 
## Input replicates can be in any of bed, bedGraph, bigWig or bam format. 
## In case that data are in the bam format (tag alignment data), you have multiple options:
## 1: you can either convert the bam file to raw signals or 
## 2: you can convert them to any of "Fold enrichment (fc)" or "p-value (pval)" signals. 
##
## We seperate these two conditions as you need to provide more arguments to pipeline to convert
## the bam file to either of fc or pval signals. We use ENCODE's default parameters for calculating the fc/pval signals.
## See "https://github.com/faezeh-bayat/VSS" for more examples on this.
##
## See "parse_arguments" for the arguments passed to this function.
## 
###################################

load_inputs <- function(args){
  
  dir.create(args$inputdir)
  path <- args$inputdir
  rep1 <- getBasename(args$tas1)
  rep2 <- getBasename(args$tas2)
  
  if(getExtension(args$tas1)=="bam" & args$signal!="raw"){
    
    
    
    system(paste("bedtools bamtobed -i ","'",args$tas1,"'",">","'",path,"'","/",rep1,sep=""))
    system(paste("bedtools bamtobed -i ","'",args$tas2,"'",">","'",path,"'","/",rep2,sep=""))
    
    ## Converts the bam files to fc/pval signals (outputs are in bigWig format)
    
    system(paste("python bam_to_fc.py ","'",path,"'","/",rep1 , " --fraglen ", args$fraglen1, " --shift 0", " --chrsz ", args$chrsz," --gensz ", args$gensz," --pval-thresh 0.01 --ctl-subsample 1 --out-dir ","'",path,"'"," --log-level INFO",sep=""))
    system(paste("python bam_to_fc.py ","'",path,"'","/",rep2 , " --fraglen ", args$fraglen2, " --shift 0", " --chrsz ", args$chrsz," --gensz ", args$gensz," --pval-thresh 0.01 --ctl-subsample 1 --out-dir ","'",path,"'"," --log-level INFO",sep=""))
    
    
    ## Converts the fc/pval signals (bigWig format) produced from previous "bam_to_fc.py" step to bed files
  
    system(paste("bigWigToBedGraph ","'",path,"'","/", rep1, ".fc.signal.bigwig ","'",path,"'","/",rep1,"_fc.bed",sep=""))
    system(paste("bigWigToBedGraph ","'",path,"'","/", rep2, ".fc.signal.bigwig ","'",path,"'","/",rep2,"_fc.bed",sep=""))
    system(paste("bigWigToBedGraph ","'",path,"'","/", rep1, ".pval.signal.bigwig ","'",path,"'","/",rep1,"_pval.bed",sep=""))
    system(paste("bigWigToBedGraph ","'",path,"'","/", rep2, ".pval.signal.bigwig ","'",path,"'","/",rep2,"_pval.bed",sep=""))
    
    
    system(paste("rm ","'",path,"'","/rep1.fc.signal.bigwig",sep=""))
    system(paste("rm ","'",path,"'","/rep2.fc.signal.bigwig",sep=""))
    system(paste("rm ","'",path,"'","/rep1.pval.signal.bigwig",sep=""))
    system(paste("rm ","'",path,"'","/rep2.pval.signal.bigwig",sep=""))
    system(paste("rm ","'",path,"'","/rep1",sep=""))
    system(paste("rm ","'",path,"'","/rep2",sep=""))
    
    
    if(args$signal=="fc"){
      
      system(paste("bedtools intersect -a ", "'",path,"'", "/",rep1,"_fc.bed -b ", "'",path,"'", "/",rep2,"_fc.bed -sorted | awk '{print $1}'>" ,"'",path,"'", "/intersect_file_chr_fc.bed",sep=""))
      system(paste("bedtools intersect -a ", "'",path,"'", "/",rep1,"_fc.bed -b ", "'",path,"'", "/",rep2,"_fc.bed -sorted | awk '{print $2,$3,$4}'>" ,"'",path,"'", "/intersect_file_rep1_fc.bed",sep=""))
      system(paste("bedtools intersect -a ", "'",path,"'", "/",rep2,"_fc.bed -b ", "'",path,"'", "/",rep1,"_fc.bed -sorted | awk '{print $4}'>" ,"'",path,"'", "/intersect_file_rep2_fc.bed",sep=""))
      

    }else if(args$signal=="pval"){
      
      system(paste("bedtools intersect -a ", "'",path,"'", "/",rep1,"_pval.bed -b ", "'",path,"'", "/",rep2,"_pval.bed -sorted | awk '{print $1}'>" ,"'",path,"'", "/intersect_file_chr_pval.bed",sep=""))
      system(paste("bedtools intersect -a ", "'",path,"'", "/",rep1,"_pval.bed -b ", "'",path,"'", "/",rep2,"_pval.bed -sorted | awk '{print $2,$3,$4}'>" ,"'",path,"'", "/intersect_file_rep1_pval.bed",sep=""))
      system(paste("bedtools intersect -a ", "'",path,"'", "/",rep2,"_pval.bed -b ", "'",path,"'", "/",rep1,"_pval.bed -sorted | awk '{print $4}'>" ,"'",path,"'", "/intersect_file_rep2_pval.bed",sep=""))
      
    }else if(args$signal=="both"){
      
      system(paste("bedtools intersect -a ", "'",path,"'", "/",rep1,"_fc.bed -b ", "'",path,"'", "/",rep2,"_fc.bed -sorted | awk '{print $1}'>" ,"'",path,"'", "/intersect_file_chr_fc.bed",sep=""))
      system(paste("bedtools intersect -a ", "'",path,"'", "/",rep1,"_fc.bed -b ", "'",path,"'", "/",rep2,"_fc.bed -sorted | awk '{print $2,$3,$4}'>" ,"'",path,"'", "/intersect_file_rep1_fc.bed",sep=""))
      system(paste("bedtools intersect -a ", "'",path,"'", "/",rep2,"_fc.bed -b ", "'",path,"'", "/",rep1,"_fc.bed -sorted | awk '{print $4}'>" ,"'",path,"'", "/intersect_file_rep2_fc.bed",sep=""))
      
      
      system(paste("bedtools intersect -a ", "'",path,"'", "/",rep1,"_pval.bed -b ", "'",path,"'", "/",rep2,"_pval.bed -sorted | awk '{print $1}'>" ,"'",path,"'", "/intersect_file_chr_pval.bed",sep=""))
      system(paste("bedtools intersect -a ", "'",path,"'", "/",rep1,"_pval.bed -b ", "'",path,"'", "/",rep2,"_pval.bed -sorted | awk '{print $2,$3,$4}'>" ,"'",path,"'", "/intersect_file_rep1_pval.bed",sep=""))
      system(paste("bedtools intersect -a ", "'",path,"'", "/",rep2,"_pval.bed -b ", "'",path,"'", "/",rep1,"_pval.bed -sorted | awk '{print $4}'>" ,"'",path,"'", "/intersect_file_rep2_pval.bed",sep=""))
      
      
    }
    
  }else if(getExtension(args$tas1)=="bed"){
    
    system(paste("cat ","'",args$tas1,"'",">","'",path,"'","/",rep1,".bed",sep=""))
    system(paste("cat ","'",args$tas2,"'",">","'",path,"'","/",rep2,".bed",sep=""))
    
    
    system(paste("bedtools intersect -a ", "'",path,"'", "/",rep1,".bed -b ", "'",path,"'", "/",rep2,".bed -sorted | awk '{print $1}'>" ,"'",path,"'", "/intersect_file_rep1_chr.bed",sep=""))
    system(paste("bedtools intersect -a ", "'",path,"'", "/",rep1,".bed -b ", "'",path,"'", "/",rep2,".bed -sorted | awk '{print $2,$3,$4}'>" ,"'",path,"'", "/intersect_file_rep1.bed",sep=""))
    system(paste("bedtools intersect -a ", "'",path,"'", "/",rep2,".bed -b ", "'",path,"'", "/",rep1,".bed -sorted | awk '{print $4}'>" ,"'",path,"'", "/intersect_file_rep2.bed",sep=""))
    
    
  }else if(getExtension(args$tas1)=="bedGraph"){
    
    system(paste("cat ","'",args$tas1,"'",">","'",path,"'","/",rep1,".bed",sep=""))
    system(paste("cat ","'",args$tas2,"'",">","'",path,"'","/",rep2,".bed",sep=""))
    
    system(paste("bedtools intersect -a ", "'",path,"'", "/",rep1,".bed -b ", "'",path,"'", "/",rep2,".bed -sorted | awk '{print $1}'>" ,"'",path,"'", "/intersect_file_rep1_chr.bed",sep=""))
    system(paste("bedtools intersect -a ", "'",path,"'", "/",rep1,".bed -b ", "'",path,"'", "/",rep2,".bed -sorted | awk '{print $2,$3,$4}'>" ,"'",path,"'", "/intersect_file_rep1.bed",sep=""))
    system(paste("bedtools intersect -a ", "'",path,"'", "/",rep2,".bed -b ", "'",path,"'", "/",rep1,".bed -sorted | awk '{print $4}'>" ,"'",path,"'", "/intersect_file_rep2.bed",sep=""))
    
    
  }else if(getExtension(args$tas1)=="bigWig"){
    
    system(paste("bigWigToBedGraph " ,"'",args$tas1,"'" ," ","'",path,"'","/",rep1,".bed",sep=""))
    system(paste("bigWigToBedGraph " ,"'",args$tas2,"'" ," ","'",path,"'","/",rep2,".bed",sep=""))
    
    system(paste("bedtools intersect -a ", "'",path,"'", "/",rep1,".bed -b ", "'",path,"'", "/",rep2,".bed -sorted | awk '{print $1}'>" ,"'",path,"'", "/intersect_file_rep1_chr.bed",sep=""))
    system(paste("bedtools intersect -a ", "'",path,"'", "/",rep1,".bed -b ", "'",path,"'", "/",rep2,".bed -sorted | awk '{print $2,$3,$4}'>" ,"'",path,"'", "/intersect_file_rep1.bed",sep=""))
    system(paste("bedtools intersect -a ", "'",path,"'", "/",rep2,".bed -b ", "'",path,"'", "/",rep1,".bed -sorted | awk '{print $4}'>" ,"'",path,"'", "/intersect_file_rep2.bed",sep=""))
    
  
  }else if(getExtension(args$tas1)=="bam" & args$signal=="raw"){
    
    system(paste("bedtools genomecov -ibam " ,"'",args$tas1,"'", " -bga >","'",path,"'", "/",rep1,".bed",sep=""))
    system(paste("bedtools genomecov -ibam " ,"'",args$tas2,"'", " -bga >","'",path,"'", "/",rep2,".bed",sep=""))
    
    system(paste("bedtools intersect -a ", "'",path,"'", "/",rep1,".bed -b ", "'",path,"'", "/",rep2,".bed -sorted | awk '{print $1}'>" ,"'",path,"'", "/intersect_file_rep1_chr.bed",sep=""))
    system(paste("bedtools intersect -a ", "'",path,"'", "/",rep1,".bed -b ", "'",path,"'", "/",rep2,".bed -sorted | awk '{print $2,$3,$4}'>" ,"'",path,"'", "/intersect_file_rep1.bed",sep=""))
    system(paste("bedtools intersect -a ", "'",path,"'", "/",rep2,".bed -b ", "'",path,"'", "/",rep1,".bed -sorted | awk '{print $4}'>" ,"'",path,"'", "/intersect_file_rep2.bed",sep=""))
    
    
  }
}


###################################
## VSS_train
## 
## Identifies the mean-variance relationship from the user provided replicates.
## See "parse_arguments" for the arguments passed to this function.
## 
## 
###################################


VSS_train <- function(args){
  
  start_time <- Sys.time()
  dir.create(args$traindir)
  train_path <- args$traindir
  path <- args$inputdir

  
  if(getExtension(args$tas1)=="bam" & args$signal=="fc"){
    
    intersect_intervals_chr <- read.table(paste(path,"/intersect_file_chr_fc.bed",sep=""))
    intersect_intervals <- read.table(paste(path,"/intersect_file_rep1_fc.bed",sep=""))
    intersect_intervals_rep2 <- read.table(paste(path,"/intersect_file_rep2_fc.bed",sep=""))
    intersect_intervals <- cbind(intersect_intervals_chr,intersect_intervals,intersect_intervals_rep2)
    
    ## uncomment the folllowing if working on a subset of genome
    
    # intersect_intervals <- intersect_intervals[intersect_intervals[,1]=="chr21",]
    
    intersect_intervals <- intersect_intervals[,-1]
    
    
    
  }else if(getExtension(args$tas1)=="bam" & args$signal=="pval"){
  
    intersect_intervals_chr <- read.table(paste(path,"/intersect_file_chr_pval.bed",sep=""))
    intersect_intervals <- read.table(paste(path,"/intersect_file_rep1_pval.bed",sep=""))
    intersect_intervals_rep2 <- read.table(paste(path,"/intersect_file_rep2_pval.bed",sep=""))
    intersect_intervals <- cbind(intersect_intervals_chr,intersect_intervals,intersect_intervals_rep2)
    
    # intersect_intervals <- intersect_intervals[intersect_intervals[,1]=="chr21",]
    
    intersect_intervals <- intersect_intervals[,-1]
    
  }else if(getExtension(args$tas1)=="bam" & args$signal=="raw"){
    
    intersect_intervals_chr <- read.table(paste(path,"/intersect_file_chr.bed",sep=""))
    intersect_intervals <- read.table(paste(path,"/intersect_file_rep1.bed",sep=""))
    intersect_intervals_rep2 <- read.table(paste(path,"/intersect_file_rep2.bed",sep=""))
    intersect_intervals <- cbind(intersect_intervals_chr,intersect_intervals,intersect_intervals_rep2)
    
    # intersect_intervals <- intersect_intervals[intersect_intervals[,1]=="chr21",]
    
    intersect_intervals <- intersect_intervals[,-1]
    
  }else{
    
    intersect_intervals_chr <- read.table(paste(path,"/intersect_file_chr.bed",sep=""))
    intersect_intervals <- read.table(paste(path,"/intersect_file_rep1.bed",sep=""))
    intersect_intervals_rep2 <- read.table(paste(path,"/intersect_file_rep2.bed",sep=""))
    intersect_intervals <- cbind(intersect_intervals_chr,intersect_intervals,intersect_intervals_rep2)
    
    # intersect_intervals <- intersect_intervals[intersect_intervals[,1]=="chr21",]
    
    intersect_intervals <- intersect_intervals[,-1]
    
  }
  
  ## Concatenates two input replicates
  scores <- big.matrix(nrow=(2 * length(intersect_intervals[,1])),ncol=3)
  scores[1:length(intersect_intervals[,1]),1]  <- intersect_intervals[,3]
  scores[1:length(intersect_intervals[,1]),2]  <- intersect_intervals[,4]
  scores[1:length(intersect_intervals[,1]),3]  <- intersect_intervals[,2]-intersect_intervals[,1]
  scores[(length(intersect_intervals[,1])+1):(2 * length(intersect_intervals[,1])),1]  <- intersect_intervals[,4]
  scores[(length(intersect_intervals[,1])+1):(2 * length(intersect_intervals[,1])),2]  <- intersect_intervals[,3]
  scores[(length(intersect_intervals[,1])+1):(2 * length(intersect_intervals[,1])),3]  <- intersect_intervals[,2]-intersect_intervals[,1]
  keep_rows <- morder(scores, 1, decreasing = FALSE)
  ordered_scores <- scores[keep_rows,]
  rm(scores)
  gc()
  ordered_scores <- data.frame(ordered_scores)
  names(ordered_scores)[1:3] <- c("Base_score","Aux_score","Aux_interval")
  L <- sum(as.numeric(ordered_scores$Aux_interval))
  zero_signals <- subset(ordered_scores,ordered_scores$Base_score==0)
  threshold <- 0.25 * L
  
  ## Identifying the mean-varaince trend has two modes:
  ## 1) Input replicates are either Fold enrichment or raw signals (Zero-inflated):
  ##    We define seperate bin for signals less than average signals.
  ## 2) Input replicates are represented as signal p-value (Non-zero inflated):
  ##    All signals are considered in mean-varaince trend identification process.
  
  if(sum(as.numeric(zero_signals$Aux_interval)) >= threshold){
    
    print("Zero inflated mode activated")
    
    ## mean-var curve smoothing parameters 
    beta_values <- 1000
    bin_sizes <- 100000
    alpha_values <- 2 ^ (1/as.integer(beta_values))
    width_values <- ceiling(log(1/0.01)/log(alpha_values))
    
    ## We define seperate bin for signals less than average signals
    
    mean_ordered <- sum(as.numeric(ordered_scores$Base_score*ordered_scores$Aux_interval))/sum(as.numeric(ordered_scores$Aux_interval))
    less_than_mean <- subset(ordered_scores,ordered_scores$Base_score<=mean_ordered)
    higher_than_mean <- subset(ordered_scores,ordered_scores$Base_score>mean_ordered)
    gc()
    Mean_1_over_sigma <- Weighted_Mean_1_over_sigma(width_values,as.integer(bin_sizes),alpha_values,higher_than_mean)
    rm(higher_than_mean)
    gc()
    Mean_1_over_sigma <- data.frame(Mean_1_over_sigma)
    mean_less_than_mean <- sum(as.numeric(less_than_mean$Aux_score*less_than_mean$Aux_interval))/sum(as.numeric(less_than_mean$Aux_interval))
    Mean_1_over_sigma[(nrow(Mean_1_over_sigma)+1),1] <- mean_less_than_mean
    sd_less_than_mean <- (sqrt(sum(as.numeric(((less_than_mean$Aux_score-mean_less_than_mean)^2)*less_than_mean$Aux_interval))/sum(as.numeric(less_than_mean$Aux_interval))))
    Mean_1_over_sigma[(nrow(Mean_1_over_sigma)),2] <- 1/sd_less_than_mean
    Mean_1_over_sigma[(nrow(Mean_1_over_sigma)+1),1] <- 0
    Mean_1_over_sigma[(nrow(Mean_1_over_sigma)),2] <- 1/sd_less_than_mean
    Mean_1_over_sigma <- data.frame(Mean_1_over_sigma)
    
    ## In special cases where IQR(data)=0, we use differrent smoothing parameters
    
    if(IQR(Mean_1_over_sigma[,1]) == 0){
      
      beta_values <- 1e+07
      bin_sizes <- 100000
      alpha_values <- 2 ^ (1/as.integer(beta_values))
      width_values <- ceiling(log(1/0.01)/log(alpha_values))
      
      mean_ordered <- sum(as.numeric(ordered_scores$Base_score*ordered_scores$Aux_interval))/sum(as.numeric(ordered_scores$Aux_interval))
      less_than_mean <- subset(ordered_scores,ordered_scores$Base_score<=mean_ordered)
      higher_than_mean <- subset(ordered_scores,ordered_scores$Base_score>mean_ordered)
      rm(ordered_scores)
      gc()
      Mean_1_over_sigma <- Weighted_Mean_1_over_sigma(width_values,as.integer(bin_sizes),alpha_values,higher_than_mean)
      rm(higher_than_mean)
      gc()
      Mean_1_over_sigma <- data.frame(Mean_1_over_sigma)
      mean_less_than_mean <- sum(as.numeric(less_than_mean$Aux_score*less_than_mean$Aux_interval))/sum(as.numeric(less_than_mean$Aux_interval))
      Mean_1_over_sigma[(nrow(Mean_1_over_sigma)+1),1] <- mean_less_than_mean
      sd_less_than_mean <- (sqrt(sum(as.numeric(((less_than_mean$Aux_score-mean_less_than_mean)^2)*less_than_mean$Aux_interval))/sum(as.numeric(less_than_mean$Aux_interval))))
      Mean_1_over_sigma[(nrow(Mean_1_over_sigma)),2] <- 1/sd_less_than_mean
      Mean_1_over_sigma[(nrow(Mean_1_over_sigma)+1),1] <- 0
      Mean_1_over_sigma[(nrow(Mean_1_over_sigma)),2] <- 1/sd_less_than_mean
      Mean_1_over_sigma <- data.frame(Mean_1_over_sigma)
    }
    
  }else if(sum(as.numeric(zero_signals$Aux_interval)) < threshold){
    
    print("Non-zero inflated mode activated")
    
    beta_values <- 1e+07
    bin_sizes <- 100000
    alpha_values <- 2^(1/as.integer(beta_values))
    width_values <- ceiling(log(1/0.01)/log(alpha_values))
    
    Mean_1_over_sigma <- Weighted_Mean_1_over_sigma(width_values,as.integer(bin_sizes),alpha_values,ordered_scores)
    Mean_1_over_sigma <- data.frame(Mean_1_over_sigma)
    
  }
  save(Mean_1_over_sigma,file=paste(train_path,"/Trained_mean_variance_model.Rdata",sep=""))
  end_time <- Sys.time()
  print(end_time - start_time)
}


###################################
## Weighted_Mean_1_over_sigma
## INTERNAL FUNCTION
## 
## 
## This function is called in "VSS_train" function. 
## It calculates the weighted (smoothed) mean-variance relationship
## 
###################################


Weighted_Mean_1_over_sigma <- function(distance,bin,alpha,ordered_scores){
  
  ordered_scores[,1] <- NULL
  ordered_scores[,3] <- ordered_scores$Aux_score*ordered_scores$Aux_interval
  ordered_scores[,4] <- (ordered_scores$Aux_score^2)*ordered_scores$Aux_interval
  names(ordered_scores)[3:4] <- c("sum","squared")
  L <- ceiling(sum(ordered_scores$Aux_interval)/bin)
  l <- bin
  ordered_scores[,5] <- floor(ordered_scores$Aux_interval/l)
  names(ordered_scores)[5] <- c("portion")
  Y <- matrix(0,nrow=L,ncol = 1)
  a <- matrix(0,nrow=L,ncol = 1)
  counter_bin <- 1
  z <- matrix(nrow=((2*distance)+1),ncol = 1)
  weighted_z <- matrix(nrow=((2*distance)+1),ncol=1)
  n <- (2*distance)+1
  mean_1_over_sigma <- matrix(nrow=L,ncol=2)
  y_value <- rep(0,nrow(Y))
  y_bins_counter <- rep(l,nrow(Y))
  a_value <- rep(0,nrow(Y))
  temp_interval <- c()
  temp_interval <- ordered_scores$Aux_interval
  temp_sum <- c()
  temp_sum <- ordered_scores$sum
  temp_squared <- c()
  temp_squared <- ordered_scores$squared
  temp_portion <- c()
  temp_portion <- ordered_scores$portion
  lll <- length(ordered_scores$Aux_score)
  rm(ordered_scores)
  gc()
  
  ## Reads the ordered_scores line by line to calculate "Y <- y_value" and "a <- a_value"
  
  for(i in 1:lll)
  {
    main_interval <- temp_interval[i]
    if(temp_portion[i]==0 && temp_interval[i]<=y_bins_counter[counter_bin]){
      y_value[counter_bin] <- y_value[counter_bin] + temp_sum[i]
      a_value[counter_bin] <- a_value[counter_bin] + temp_squared[i]
      y_bins_counter[counter_bin] <- y_bins_counter[counter_bin] - temp_interval[i]
      temp_interval[i] <- 0
      if(y_bins_counter[counter_bin]==0){
        counter_bin <- counter_bin+1
      }
    }else if(temp_portion[i]==0 && temp_interval[i]>y_bins_counter[counter_bin]){
      
      ## Fills the previous bin
      
      y_value[counter_bin] <- y_value[counter_bin] + (y_bins_counter[counter_bin]/main_interval)*temp_sum[i]
      a_value[counter_bin] <- a_value[counter_bin] + (y_bins_counter[counter_bin]/main_interval)*temp_squared[i]
      temp_interval[i] <- temp_interval[i] - y_bins_counter[counter_bin]
      y_bins_counter[counter_bin] <- 0
      counter_bin <- counter_bin+1
      
      ## Fills the current bin
      
      y_value[counter_bin] <- y_value[counter_bin] + (temp_interval[i]/main_interval)*temp_sum[i]
      a_value[counter_bin] <- a_value[counter_bin] + (temp_interval[i]/main_interval)*temp_squared[i]
      y_bins_counter[counter_bin] <- y_bins_counter[counter_bin] - temp_interval[i]
      temp_interval[i] <- 0
    }else if(temp_portion[i]!=0){
      left <- y_bins_counter[counter_bin]
      while(temp_interval[i]>left){
        y_value[counter_bin] <- y_value[counter_bin] + (left/main_interval)*temp_sum[i]
        a_value[counter_bin] <- a_value[counter_bin] + (left/main_interval)*temp_squared[i]
        y_bins_counter[counter_bin] <- 0
        temp_interval[i] <- temp_interval[i] - left
        counter_bin <- counter_bin + 1
        left <- y_bins_counter[counter_bin]
      }
      if(temp_interval[i]==left){
        y_value[counter_bin] <- y_value[counter_bin] + (left/main_interval)*temp_sum[i]
        a_value[counter_bin] <- a_value[counter_bin] + (left/main_interval)*temp_squared[i]
        y_bins_counter[counter_bin] <- 0
        temp_interval[i] <- 0
        counter_bin <- counter_bin + 1
      }else if(temp_interval[i]<left){
        y_value[counter_bin] <- y_value[counter_bin] + (temp_interval[i]/main_interval)*temp_sum[i]
        a_value[counter_bin] <- a_value[counter_bin] + (temp_interval[i]/main_interval)*temp_squared[i]
        y_bins_counter[counter_bin] <- y_bins_counter[counter_bin] - temp_interval[i]
        temp_interval[i] <- 0
      }
    }
  }
  Y[,1] <- y_value
  a[,1] <- a_value
  rm(temp_interval)
  rm(temp_portion)
  rm(temp_squared)
  rm(temp_sum)
  gc()
  
  ## Calculates the weight for bins contributing to the mean and variance of the current bin
  
  for (j in 1:n){
    z[j,1] <- (1/(alpha^(abs(distance+1-j)*l)))
    weighted_z[j,1] <- l*(1/(alpha^(abs(distance+1-j)*l)))
  }
  
  for(i in 1:L){
    output <- Smoothed_Weighted_Mean_1_over_sigma(Y,z,weighted_z,a,i,distance)
    mean_1_over_sigma[i,1] <- output[1]
    mean_1_over_sigma[i,2] <- output[2]
  }
  mean_1_over_sigma
}


###################################
## Smoothed_Weighted_Mean_1_over_sigma
## INTERNAL FUNCTION
## 
## Applies the bandwidth parameter (which measures the number of 
## adjacent bins to be considered in calculating current bin's mean) to 
## the raw mean and variance of the bins.
## 
## 
###################################

Smoothed_Weighted_Mean_1_over_sigma <- function(weighted_sum,weights,sum_of_weights,sum_square,center_bin,width){
  
  epsilon <- 0.0001
  y <- weighted_sum
  W <- weighted_sum
  L <- length(y)
  z <- weights
  weighted_z <- sum_of_weights
  k <- center_bin
  w <- width
  a <- sum_square
  A <- sum_square
  
  lower <- max(1,(k-w))
  upper <- min(L,(k+w))
  xx <- upper-lower+1
  center_z <- ceiling(length(z)/2)
  lower_z <- (center_z-(k-lower))
  counter <- lower_z
  upper_z <- (center_z+(upper-k))
  for(i in lower:upper){
    
    y[i,1] <- y[i,1]*z[counter,1]
    counter <- counter+1
  }
  weighted_mean <- sum(y[lower:upper])/sum(weighted_z[lower_z:upper_z])
  counter <- lower_z
  for(i in lower:upper){
    
    a[i,1] <- a[i,1]*z[counter,1]
    counter <- counter+1
  }
  weighted_var <- sum(a[lower:upper])/sum(weighted_z[lower_z:upper_z])
  weighted_var <- weighted_var -(weighted_mean^2)
  weighted_var <- weighted_var+epsilon
  mean_1_over_sigma_values <- c(weighted_mean,1/sqrt(weighted_var))
  mean_1_over_sigma_values
}


###################################
## VSS_transform
## 
## Takes the unstabilized input replicate (args$tas1) and transformes it using the
## identifed mean-variance relationship saved in args$traindir (Trained_mean_variance_model.Rdata)
## 
##
## Outputs variance-stabilized signals ("Variance_stabilized_signals.bed") at args$transformdir.
## 
###################################

VSS_transform <- function(args){
  
  dir.create(args$transformdir)
  transform_path <- args$transformdir
  train_path <- args$traindir
  input_path <- args$inputdir
  rep1 <- getBasename(args$tas1)
  if(getExtension(args$tas1)=="bam" & args$signal=="fc"){
    
    instabilized_replicate_signals <- read.table(paste(input_path,"/",rep1,"_fc.bed",sep=""))
    instabilized_replicate_signals <- data.frame(instabilized_replicate_signals)
    
  }else if(getExtension(args$tas1)=="bam" & args$signal=="pval"){
    
    instabilized_replicate_signals <- read.table(paste(input_path,"/",rep1,"_pval.bed",sep=""))
    instabilized_replicate_signals <- data.frame(instabilized_replicate_signals)
    
  }else if(getExtension(args$tas1)=="bam" & args$signal=="raw"){
    
    instabilized_replicate_signals <- read.table(paste(input_path,"/",rep1,".bed",sep=""))
    instabilized_replicate_signals <- data.frame(instabilized_replicate_signals)
  }else{
    
    instabilized_replicate_signals <- read.table(paste(input_path,"/",rep1,".bed",sep=""))
    instabilized_replicate_signals <- data.frame(instabilized_replicate_signals)
  }
  
  load(paste(train_path,"/Trained_mean_variance_model.Rdata",sep=""))
  Mean_1_over_sigma <- data.frame(Mean_1_over_sigma)
  start_t <- Sys.time()
  Replicate1_scores <- Calculating_vss_signals(Mean_1_over_sigma,instabilized_replicate_signals)
  end_t <- Sys.time()
  print(end_t - start_t)
  write.table(Replicate1_scores, file=paste(transform_path,"/Variance_stabilized_signals.bedGraph",sep=""), quote=F, sep="\t", row.names=F, col.names=F)
  
}

###################################
## Calculating_vss_signals
## INTERNAL FUNCTION
## 
## Calculates variance-stabilized signals
## 
##
###################################

Calculating_vss_signals <- function(Mean_1_over_sigma,rep1){
  
  replicate1_scores <- rep1
  ord <- order(replicate1_scores[,4])
  replicate1_scores <- replicate1_scores[ord,]
  replicate1_scores <- data.frame(replicate1_scores)
  unstab_signals <- replicate1_scores[,4]
  unique_numbers <- as.numeric(names(table(unstab_signals)))
  unique_unstab_signals <- unique_numbers
  frequency_unique_numbers <- as.numeric(table(unstab_signals))
  Mean_vs_sigma <- Mean_1_over_sigma
  Mean_vs_sigma <- data.frame(Mean_vs_sigma)
  Mean_vs_sigma[,1] <- Mean_1_over_sigma[,1]
  Mean_vs_sigma[,2] <- 1/Mean_1_over_sigma[,2]
  fit_1 <- smooth.spline(Mean_vs_sigma[,1], (Mean_vs_sigma[,2]), cv=TRUE)
  predicted_sigma <- stats:::predict.smooth.spline(fit_1,unique_unstab_signals)$y
  predicted_score <- matrix(nrow=length(predicted_sigma),ncol=2)
  predicted_score[,1] <- unique_unstab_signals
  predicted_score[,2] <- predicted_sigma
  if(min(predicted_score[,2])<0){
    less <- subset(predicted_score,predicted_score[,1]<min(Mean_vs_sigma[,1]))
    if(min(less[,2]) < 0){
      less[,2] <- min(Mean_vs_sigma[,2])
    }
    middle <- subset(predicted_score,predicted_score[,1] >= min(Mean_vs_sigma[,1]) & predicted_score[,1] <= max(Mean_vs_sigma[,1]))
    high <- subset(predicted_score,predicted_score[,1] > max(Mean_vs_sigma[,1]))
    if(min(high[,2]) < 0){
      high[,2] <- max(Mean_vs_sigma[,2])
    }
    predicted_score <- rbind(less,middle,high)
  }
  predicted_score[predicted_score[,2] < 0, 2] <- mean(Mean_vs_sigma[,2])
  predicted_score[,2] <- 1/predicted_score[,2]
  x <- predicted_score[,1]
  y <- predicted_score[,2]
  ordered_data <- predicted_score[order(x,y),]
  ordered_data <- data.frame(ordered_data)
  
  integral <- function(x)
  {
    sub <- subset(ordered_data, ordered_data[,1]<=x)
    score <-trapz(sub[,1],sub[,2])
    score
  }
  unique_stabilized_signals <- sapply(unique_unstab_signals,integral)
  stabilized_signals <- c()
  counter <- 1
  for(i in 1:length(unique_unstab_signals))
  {
    stabilized_signals[counter:(frequency_unique_numbers[i]+counter-1)] <- unique_stabilized_signals[i]
    counter <- frequency_unique_numbers[i]+counter
  }
  replicate1_scores[,5] <- stabilized_signals
  replicate1_scores[,4] <- NULL
  replicate1_scores <- replicate1_scores[order(ord),]
  replicate1_scores
}







main <- function(args){
  
  if(args[1]=="load_inputs"){
    args <- parse_arguments()
    load_inputs(args)
  }else if(args[1]=="train"){
    args <- parse_arguments()
    VSS_train(args)
  }else if(args[1]=="transform"){
    args <- parse_arguments_transform()
    VSS_transform(args)
  }
}

args <- commandArgs(trailingOnly=TRUE)
main(args)





