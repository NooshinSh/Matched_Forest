##################
library(randomForest)
library(ggplot2)
##################
#Inputs
dir_main = paste('C:/.../',sep='')
dir_input = paste(dir_main,'data/',sep='')
dir_output = paste(dir_main,'output/',sep='')
data_name='Childhood_Acute_Lymphoblastic_Leukemia_Data'
############################################################
df<-read.table(paste(dir_input,data_name,'.csv',sep=''),header = TRUE)
############################################################
Exp_indx=c(3:ncol(df))
Match_indx<-c()
Outcome_indx=1
Match<-length(Match_indx)
Exp<-length(Exp_indx)
case<-df[df[,Outcome_indx]==1,Exp_indx]
control<-df[df[,Outcome_indx]==0,Exp_indx]
matching<-df[df[,Outcome_indx]==0,Match_indx]
#############################################################
#1. Transform to Supervised Learning
case_star<-rbind(case,control)
control_star<-rbind(control,case)
d_star<-case_star-control_star
matching_plus<-rbind(matching,matching)
label<-rep(c(0,1),each=nrow(case))

if (is.null(Match_indx)){
  df_transformed<-cbind(label,case_star,control_star,d_star)
}else{
  df_transformed<-cbind(label,matching_plus,case_star,control_star,d_star)
}

###########################################################
#2. RF with matched variable importance
Run=10
ntree=5
M=dim(df_transformed)[2]-1
mtry=floor(sqrt(M))
sampsize=dim(df_transformed)[1]
maxnodes=NULL

#2.1. RF
imp<-matrix(nrow=M,ncol=Run)
for (i in c(1:Run)){
  print(paste('Run=',i))
  MF <- randomForest(df_transformed[,2:(ncol(df_transformed))], y =  as.factor(df_transformed[,1]), 
                     mtry=mtry, sampsize=sampsize,
                     replace = TRUE, importance=TRUE, ntree = ntree, 
                     maxnodes = maxnodes, keep.forest = TRUE)
  imp[,i]<-importance(MF,type=1)
}

#2.2. Matched variable importance

#2.2.1. Matching Variables
if (is.null(Match_indx)==FALSE){
  impM<-imp[1:Match,]
}

#2.2.2. Exposure variables
imp_case<-imp[c((Match+1):(Match+Exp)),]
imp_control<-imp[c((Match+1+Exp):(Match+2*Exp)),]
imp_d<-imp[c((Match+2*Exp+1):(Match+3*Exp)),]

impE<-imp_case+imp_control+imp_d
