exp.imp <- function(imp,Match,Exp){
    if (Match==0){
      df1=data.frame()
    }else{
      df1=data.frame(imp[1:Match,])
    }
    
    for (j in 1:Exp){
      a=imp[(Match+2*j-1),]+imp[(Match+2*j),]+imp[(Match+2*Exp+j),]
      df1=rbind(df1,a)
    }
  return(df1)
}

#pvalues
pVal <- function(obs, null.dist) {
  (sum(null.dist >= obs) + 1) / (length(null.dist) + 1)
}

##################
library(randomForest)
##################
#Inputs
dir_main = paste('C:/.../',sep='')
dir_input = paste(dir_main,'data/',sep='')
dir_output = paste(dir_main,'output/',sep='')
dir_input='C:/Users/Nooshin Shomalzadeh/Documents/Research/Matched Forest/Paper/Manuscript/code/'
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
p.iter=100
ntree=50000
M=dim(df_transformed)[2]-1
mtry=floor(sqrt(M))
N=dim(df_transformed)[1]
maxnodes=NULL

#2.1. RF
imp_null<-matrix(0,nrow=M,ncol=p.iter)
imp_obs<-matrix(0,nrow=M,ncol=Run)
impE_null<-matrix(0,nrow=Match+3*Exp,ncol=p.iter)
impE_obs<-matrix(0,nrow=Match+3*Exp,ncol=Run)
pval.mat <- matrix(0,nrow=Match+3*Exp,ncol=1)
#observed scores
X=df_transformed[,2:(ncol(df_transformed))]
y=as.factor(df_transformed[,1])
for (i in c(1:Run)){
  print(paste('Run=',i))
  MF <- randomForest(x=X, y = y, 
                     mtry=mtry, sampsize=N,
                     replace = TRUE, importance=TRUE, ntree = ntree, 
                     maxnodes = maxnodes, keep.forest = TRUE)
  imp_obs[,i]<-importance(MF,type=1)
}
#null distribution for scores

for (j in c(1:p.iter)){
  print(paste('permute = ',j,sep=''))
  XP=as.data.frame(matrix(0,nrow=N,ncol=M))
  if(Match>0){
    XP[,(1:Match)]=X[,c(1:Match)]
  }
  for(k in c(1:2)){
    cols=c(Match+2*k-1,Match+2*k)
    dat=X[1:(N/2),cols]
    datp<-apply(dat,1, function(z){
      sample(z,2,replace=FALSE)
    })
    #d*
    XP[1:(N/2),(Match+2*Exp+k)]=datp[1,]-datp[2,]
    XP[(N/2+1):N,(Match+2*Exp+k)]=datp[2,]-datp[1,]
    #case*, control*
    XP[1:(N/2),cols]<-t(datp)
    XP[(N/2+1):N,cols]<-cbind(datp[2,],datp[1,])
  }
  MFP <- randomForest(x=XP, y=y, 
                     mtry=mtry, sampsize=N,
                     replace = TRUE, importance=TRUE, ntree = ntree, 
                     maxnodes = maxnodes, keep.forest = TRUE)
  imp_null[,j]<-importance(MFP,type=1)
}





#2.2. Matched variable importance
impE_null<-exp.imp(imp_null,Match,Exp)
impE_obs<-exp.imp(imp_obs,Match,Exp)

#2.3. Variable selection

impExp.obs.mean.mat=apply(impE_obs,1,mean)

pval_exp = sapply(1:(Match+Exp),function(x){
  obs=impExp.obs.mean.mat[x]
  null.dist=impE_null[x,]
  pVal(obs, null.dist)
})

pval.mat=as.matrix(pval_exp)

#save output
row.names(impE_obs)=row.names(pval.mat)<-colnames(df)[c(Match_indx,Exp_indx)]
colnames(impE_obs)<-seq(1,Run,1)
Write.csv(impE_obs,paste(dir_output,'Variable_MFI_score.csv',sep=''))
Write.csv(pval.mat,paste(dir_output,'Variable_pvalue.csv',sep=''))

