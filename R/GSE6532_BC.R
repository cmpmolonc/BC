library(GEOquery)
library(reshape)
library(sva)
library(M3C)

gse6532 <- getGEO('GSE6532',GSEMatrix=TRUE)
batchfull<-readRDS("GSE6532_batch.rds")
batchfull<-batchfull[colnames(gds),]
PAM50.df<-readRDS("GSE6532_PAM50.df")
PAM50.df<-PAM50.df[colnames(gse),]

gse <- gse6532[[2]]
gse.exp<-na.omit(exprs(gse))
gse$PAM50<-PAM50.df
gse$batch<-batchfull

####plot UMAP showing batches
M3C::umap(gse.exp,labels=as.factor(batchfull),controlscale=TRUE,scale=3,dots=1,seed = 123)

##Parametric batch correction
gse_full_par<-ComBat(dat=gse.exp, batch=batchfull, par.prior=TRUE)

##Non-Parametric batch correction
gse_full_nonpar<-ComBat(dat=gse.exp, batch=batchfull, par.prior=FALSE)

##PAM50 cofactor batch correction
modPAM50 = model.matrix(~as.factor(gds$PAM50), data=as.data.frame(exprs(gse)))
combat_edata_cov = ComBat(dat=gse.exp, batch=batchfull, mod=modPAM50, par.prior=TRUE)

##PAM50 stratified batch correction
###Batch correct all BASAL samples together
gsefiltbasal<-gse[,gse$PAM50 %in% "Basal"] 
gse_basal_par<-ComBat(dat=exprs(gsefiltbasal), batch=gsefiltbasal$batch, par.prior=TRUE)

##Batch correct all Luminal A samples together
gsefiltLumA<-gse[,gse$PAM50 %in% "LumA"] 
gse_LumA_par<-ComBat(dat=exprs(gsefiltLumA), batch=gsefiltLumA$batch, par.prior=TRUE)

##Batch correct all Luminal B samples together
gsefiltLumB<-gse[,gse$PAM50 %in% "LumB"] 
gse_LumB_par<-ComBat(dat=exprs(gsefiltLumB), batch=gsefiltLumB$batch, par.prior=TRUE)

##Batch correct all Her2 samples together
gsefiltHer2<-gse[,gse$PAM50 %in% "Her2"] 
gse_Her2_par<-ComBat(dat=exprs(gsefiltHer2), batch=gsefiltHer2$batch, par.prior=TRUE)

##Batch correct all Normal samples together
gsefiltNormal<-gse[,gse$PAM50 %in% "Normal"] 
gse_Normal_par<-ComBat(dat=exprs(gsefiltNormal), batch=gsefiltNormal$batch, par.prior=TRUE)

##recombine subtype stratified corrections together
full_gse_stratified<-cbind(gse_basal_par,gse_Normal_par,gse_Her2_par,gse_LumA_par,gse_LumB_par)
