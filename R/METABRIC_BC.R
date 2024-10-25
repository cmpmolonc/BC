####METABRIC expression data was obtained following approval from EGA

library(sva)

##create df of batch values where samples in Discovery set are batch 1 and samples in validation set are batch 2
batch<-data.frame("id"=colnames(METAdiscSet),batch="Discovery","PAM50"=METAdiscSet$Pam50Subtype)
batch2<-data.frame("id"=colnames(METAvalidSet),batch="Validation","PAM50"=METAvalidSet$Pam50Subtype)
batchfull<-rbind(batch,batch2)
rownames(batchfull)<-batchfull$id

###merge Discovery and validation set expression values into single matrix
mfull<-merge(as.data.frame(exprs(discfilt)),as.data.frame(exprs(validfilt)),by=0)
rownames(mfull)<-mfull$Row.names
mfull$Row.names<-NULL

#Parametric correction
METABRIC_full_par<-ComBat(dat=mfull, batch=batchfull$batch, par.prior=TRUE)

#Non-parametric correction
METABRIC_full_nonpar<-ComBat(dat=mfull, batch=batchfull$batch, par.prior=FALSE)

#Batch correction incorporating PAM50 as a covariate. Use PAM50cov vector of samples names and PAM50 subtype.
modPAM50 = model.matrix(~PAM50, data=PAM50cov)
METABRIC_full_PAM50cov<-ComBat(dat=mfull, batch=batchfull$batch, mod=modPAM50, par.prior=TRUE)

###Stratified correction of PAM50 molecular subtypes 
stratified_batchcorrect<-function(METAdiscovery,METAvalidation,subtype)
{
  merged_expression<-merge(as.data.frame(exprs(METAdiscovery[,METAdiscovery$Pam50Subtype %in% subtype])),as.data.frame(exprs(METAvalidation[,METAvalidation$Pam50Subtype %in% subtype])),by=0)
  rownames(merged_expression)<-merged_expression$Row.names
  merged_expression$Row.names<-NULL
  ##create subtype-specific df of batch values
  batch<-rbind(data.frame("id"=colnames(METAdiscovery[,METAdiscovery$Pam50Subtype %in% subtype]),batch=1),data.frame("id"=colnames(METAvalidation[,METAvalidation$Pam50Subtype %in% subtype]),batch=2))
  bc<-ComBat(dat=merged_expression, batch=batch$batch, par.prior=TRUE)
  return(bc)
}

METABRIC_full_basal<-stratified_batchcorrect(METAdiscSet,METAvalidSet,"Basal")
METABRIC_full_HER2<-stratified_batchcorrect(METAdiscSet,METAvalidSet,"HER2")
METABRIC_full_LUMA<-stratified_batchcorrect(METAdiscSet,METAvalidSet,"LumA")
METABRIC_full_LUMB<-stratified_batchcorrect(METAdiscSet,METAvalidSet,"LumB")
METABRIC_full_NORM<-stratified_batchcorrect(METAdiscSet,METAvalidSet,"Normal")

##Bind all Pam50 stratified corrections into a single matrix
METABRIC_full_stratified<-cbind(METABRIC_full_basal,METABRIC_full_HER2,METABRIC_full_LUMA,METABRIC_full_LUMB,METABRIC_full_NORM)

