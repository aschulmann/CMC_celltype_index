### variables defined prior to this
#cmc_counts_f: counts of 603 individuals, for whom both RNA and clinical metadata exists (kept only 1 of 2 technical replicates for 10 individuals).
#metadata_clinical_f, metadata_rna_f: clinical and RNA metadata for the same 603 ndividuals 
#idx_mt: index of mitochondrial genes


#run limma logCPM transform
library(edgeR)
v=DGEList(counts = cmc_counts_f[-idx_mt,], genes = gene_data[-idx_mt,])
v=calcNormFactors(v)
v=voom(v,plot=T)
v_csi=data.frame(v$E)
v_csi$gene_name=gene_data$symbol[match(rownames(v_csi),gene_data$Geneid)]

# parse relevant subject variables
metadata_f=cbind(metadata_clinical_f,metadata_rna_f)
metadata_f$Age_of_Death=sub('90\\+',90,metadata_f$Age_of_Death) #90+ is binned, convert to 90
SubjectFactorVariables<-metadata_f[,c('Dx','Institution','Gender','Ethnicity','Cause_of_Death','DLPFC_Brodmann_Area','DLPFC_Hemisphere')]
SubjectContinuousVariables<-data.matrix(metadata_f[,c('Age_of_Death','PMI_hrs','pH','DLPFC_RNA_isolation_RIN','DLPFC_RNA_Sequencing_Percent_Aligned')])

colSums(is.na(SubjectFactorVariables))
table(rowSums(is.na(SubjectContinuousVariables)))
table(rowSums(is.na(SubjectFactorVariables))==0 & rowSums(is.na(SubjectContinuousVariables))==0)

idx_NoNA_clinvar=rowSums(is.na(SubjectFactorVariables))==0 & rowSums(is.na(SubjectContinuousVariables))==0
idx_NoNA_clinvar2=rowSums(is.na(SubjectFactorVariables))==0 & rowSums(is.na(SubjectContinuousVariables))==0 & metadata_rna_f$DLPFC_Brodmann_Area==9
idx_NoNA_clinvar3=rowSums(is.na(SubjectFactorVariables))==0 & rowSums(is.na(SubjectContinuousVariables))==0 & metadata_clinical_f$Age_of_Death!="90+" & metadata_rna_f$DLPFC_Brodmann_Area==9

#organize metadata
SubjectFactorVariables=droplevels(SubjectFactorVariables[idx_NoNA_clinvar2,])
SubjectContinuousVariables=SubjectContinuousVariables[idx_NoNA_clinvar2,]

v_csi_noNA=data.frame(v$E[,idx_NoNA_clinvar])
v_csi_noNA$gene_name=gene_data$symbol[match(rownames(v_csi_noNA),gene_data$Geneid)]

#run BrainInABlender
library(BrainInABlender)
cmc_csi_noNA=Sir_UnMixALot(userInput = v_csi_noNA,dataColumns = 1:443,geneColumn = 444,species = 'human')

#####
#PCA

pca.all<-prcomp(t(v$E))

PC1<-pca.all$x[,1]
PC2<-pca.all$x[,2]
PC3<-pca.all$x[,3]
PC4<-pca.all$x[,4]

plot(PC1~cmc_csi$AveragePrimary_CellTypeIndex["Neuron_All",],xlab="Neuron (all)")
lm.pc1.neurons=lm(PC1~cmc_csi$AveragePrimary_CellTypeIndex["Neuron_All",])
summary.lm.pc1.neurons=summary(lm.pc1.neurons)
abline(lm.pc1.neurons, lwd=3, col=ifelse(summary.lm.pc1.neurons$coefficients[2,1]<0,'blue','red'))
title(main=paste("r-squared =",format(summary.lm.pc1.neurons$r.squared, digits = 2),
                 "\np-value =",format(summary.lm.pc1.neurons$coefficients[2,4], digits = 3)))

lmPlot=function(x, y, xlab,ylab){
  plot(y~x,xlab=xlab,ylab=ylab)
  lm=lm(y~x)
  s.lm=summary(lm)
  abline(lm, lwd=3, col=ifelse(s.lm$coefficients[2,1]<0,'blue','red'))
  title(main=paste("r-squared =",format(s.lm$r.squared, digits = 2),
                   "\np-value =",format(s.lm$coefficients[2,4], digits = 3)))
}

pdf('pc1-4_5celltypes.pdf',h=15,w=12)
par(mfrow=c(5,4))
for (i in rownames(cmc_csi$AveragePrimary_CellTypeIndex)[c(1:2,5,7:8)]){
  lmPlot(cmc_csi$AveragePrimary_CellTypeIndex[i,],PC1,xlab=i,ylab="PC1")
  lmPlot(cmc_csi$AveragePrimary_CellTypeIndex[i,],PC2,xlab=i,ylab="PC2")
  lmPlot(cmc_csi$AveragePrimary_CellTypeIndex[i,],PC3,xlab=i,ylab="PC3")
  lmPlot(cmc_csi$AveragePrimary_CellTypeIndex[i,],PC4,xlab=i,ylab="PC4")
}
dev.off()

###
SubjectFactorVariables$Dx=factor(SubjectFactorVariables$Dx,levels=c('Control','BP','SCZ'))
SubjectFactorVariables$Dx=factor(SubjectFactorVariables$Dx,levels=c('Control','BP','SCZ'))
SubjectFactorVariables$Gender=factor(SubjectFactorVariables$Gender,levels=c('Male','Female'))
SubjectFactorVariables$Ethnicity=relevel(SubjectFactorVariables$Ethnicity,ref = 'Caucasian') # 80% are Caucasian

## Cause_of_Death:
#1 (Cardiovascular), 2 (Non-brain cancer), 3 (Infection and parasitic disease), 4 (COPD), 5 (Other), -1 (Missing)
## age of death binned to 90 when 90+

metadata_noNA1=data.frame(SubjectFactorVariables,SubjectContinuousVariables)
metadata_noNA2=metadata_noNA1[,-c(2,4:7,11:12)]
head(colnames(metadata_noNA2))

lm_csi_noNA2=lm(t(cmc_csi$AveragePrimary_CellTypeIndex[,idx_NoNA_clinvar2])~., data = metadata_noNA2)
s.lm_csi_noNA2=summary(lm_csi_noNA2)
pval.lm_csi_noNA2=list()
for(i in rownames(cmc_csi$AveragePrimary_CellTypeIndex)){pval.lm_csi_noNA2[[i]]=s.lm_csi_noNA2[[paste("Response",i)]]$coefficients[,4]}
pval.lm_csi_noNA2=data.frame(pval.lm_csi_noNA2)

coef.lm_csi_noNA2=list()
for(i in rownames(cmc_csi$AveragePrimary_CellTypeIndex)){coef.lm_csi_noNA2[[i]]=s.lm_csi_noNA2[[paste("Response",i)]]$coefficients[-1,1]}
coef.lm_csi_noNA2=data.frame(coef.lm_csi_noNA2)

###


metadata_f=cbind(metadata_clinical_f,metadata_rna_f)
metadata_f$Age_of_Death=sub('90\\+',90,metadata_f$Age_of_Death)
SubjectFactorVariables<-metadata_f[,c('Dx','Institution','Gender','Cause_of_Death')]
SubjectContinuousVariables<-data.frame(data.matrix(metadata_f[,c('Age_of_Death','PMI_hrs','pH','DLPFC_RNA_isolation_RIN')]))

SubjectFactorVariables$Cause_of_Death[is.na(SubjectFactorVariables$Cause_of_Death)]="(missing)" # cause of death: 74 NAs
colSums(is.na(SubjectContinuousVariables)) # pH: 87 NAs
colSums(is.na(SubjectFactorVariables)) 

idx_NoNA_clinvar=rowSums(is.na(SubjectFactorVariables))==0 & rowSums(is.na(SubjectContinuousVariables))==0 # 516 samples
#idx_NoNA_clinvar3=rowSums(is.na(SubjectFactorVariables))==0 & rowSums(is.na(SubjectContinuousVariables))==0 & metadata_clinical_f$Age_of_Death!="90+" & metadata_rna_f$DLPFC_Brodmann_Area==9

#1 (Cardiovascular), 2 (Non-brain cancer), 3 (Infection and parasitic disease), 4 (COPD), 5 (Other), -1 (Missing)
## age of death binned to 90 when 90+

SubjectFactorVariables_noNA=droplevels(SubjectFactorVariables[idx_NoNA_clinvar,])
SubjectContinuousVariables_noNA=SubjectContinuousVariables[idx_NoNA_clinvar,]

#reset reference variable to most common var
SubjectFactorVariables$Dx=factor(SubjectFactorVariables$Dx,levels=c('Control','SCZ','BP','AFF'))
SubjectFactorVariables$Gender=factor(SubjectFactorVariables$Gender,levels=c('Male','Female'))
library(plyr)
SubjectFactorVariables$Cause_of_Death=factor(mapvalues(SubjectFactorVariables$Cause_of_Death, from = 1:5, to = c("(Cardiovascular)","(Non-brain cancer)", "(Infection and parasitic disease)", "(COPD)", "(Other)")))

metadata_all=data.frame(SubjectFactorVariables,SubjectContinuousVariables)
metadata_nopH=metadata_all[,-7]
metadata_nopH_nocause=metadata_all[,-c(5,7)]


#v_csi_noNA=data.frame(v$E[,idx_NoNA_clinvar])
#v_csi_noNA$gene_name=gene_data$symbol[match(rownames(v_csi_noNA),gene_data$Geneid)]

lm_csi_all=lm(t(cmc_csi$AveragePrimary_CellTypeIndex[,idx_NoNA_clinvar])~., data = metadata_all[idx_NoNA_clinvar,])
s.lm_csi_all=summary(lm_csi_all)
pval.lm_csi_all=list()
for(i in rownames(cmc_csi$AveragePrimary_CellTypeIndex)){pval.lm_csi_all[[i]]=s.lm_csi_all[[paste("Response",i)]]$coefficients[,4]}
pval.lm_csi_all=data.frame(pval.lm_csi_all)
coef.lm_csi_all=list()
for(i in rownames(cmc_csi$AveragePrimary_CellTypeIndex)){coef.lm_csi_all[[i]]=s.lm_csi_all[[paste("Response",i)]]$coefficients[-1,1]}
coef.lm_csi_all=data.frame(coef.lm_csi_all)


for (i in names(s.lm_csi_all)){
  write.csv(s.lm_csi_all[[i]]$coefficients, paste("../coefficients/",sub("Response ","coef_",i),".csv", sep = ""))
}

###
# linear models design using cell type index 

design_noPH=model.matrix(~Dx + Institution + Gender + Cause_of_Death + Age_of_Death + PMI_hrs + DLPFC_RNA_isolation_RIN, data=metadata_all)
design_all=model.matrix(~Dx + Institution + Gender + Cause_of_Death + Age_of_Death + PMI_hrs + pH + DLPFC_RNA_isolation_RIN, data=metadata_all[idx_NoNA_clinvar,])

metadata_all.csi=data.frame(metadata_all,t(cmc_csi$AveragePrimary_CellTypeIndex))
design.m1=model.matrix(~Dx, data=metadata_all.csi[idx_NoNA_clinvar,])
design.m2=model.matrix(~Dx + Institution + Gender + Cause_of_Death + Age_of_Death + PMI_hrs + pH + DLPFC_RNA_isolation_RIN, data=metadata_all.csi[idx_NoNA_clinvar,])
design.m3=model.matrix(~Dx + Astrocyte + Microglia + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=metadata_all.csi[idx_NoNA_clinvar,])
design.m4=model.matrix(~Dx + Institution + Gender + Cause_of_Death + Age_of_Death + PMI_hrs + pH + DLPFC_RNA_isolation_RIN + Astrocyte + Microglia + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=metadata_all.csi[idx_NoNA_clinvar,])
design.m5=model.matrix(~Dx + Institution + Gender + Cause_of_Death + Age_of_Death + PMI_hrs + pH + DLPFC_RNA_isolation_RIN + Astrocyte + Microglia + Neuron_Interneuron + Neuron_Projection + Oligodendrocyte + Endothelial + Mural + Neuron_All + Oligodendrocyte_Immature + RBC, data=metadata_all.csi[idx_NoNA_clinvar,])


#differential gene expression analysis
getDEgenes=function(voom=v.noNA, design, coef=2){
  lfit=lmFit(voom, design = design)
  ebfit=eBayes(lfit)
  topTable(ebfit, coef = coef, n=nrow(ebfit$genes))
}

library(edgeR)
d.noNA=DGEList(counts = cmc_counts_f[-idx_mt,idx_NoNA_clinvar], genes = gene_data[-idx_mt,])
d.noNA=calcNormFactors(d.noNA)
v.noNA=voom(d.noNA,plot=T)

cpm.noNA=cpm(d.noNA)
keep.x=rowSums(cpm.noNA>1)>=50

d.noNA.f=DGEList(counts = d.noNA$counts[keep.x,], genes = d.noNA$genes[keep.x,])
d.noNA.f=calcNormFactors(d.noNA.f)
v.noNA.f=voom(d.noNA.f,plot=T)

top.m1=getDEgenes(voom = v.noNA, design = design.m1)
top.m2=getDEgenes(voom = v.noNA, design = design.m2)
top.m3=getDEgenes(voom = v.noNA, design = design.m3)
top.m4=getDEgenes(voom = v.noNA, design = design.m4)
top.m5=getDEgenes(voom = v.noNA, design = design.m5)

top.m1.f=getDEgenes(voom = v.noNA.f, design = design.m1)
top.m2.f=getDEgenes(voom = v.noNA.f, design = design.m2)
top.m3.f=getDEgenes(voom = v.noNA.f, design = design.m3)
top.m4.f=getDEgenes(voom = v.noNA.f, design = design.m4)
top.m5.f=getDEgenes(voom = v.noNA.f, design = design.m5)

top.bp.m1=getDEgenes(voom = v.noNA, design = design.m1, coef = 3)
top.bp.m2=getDEgenes(voom = v.noNA, design = design.m2, coef = 3)
top.bp.m3=getDEgenes(voom = v.noNA, design = design.m3, coef = 3)
top.bp.m4=getDEgenes(voom = v.noNA, design = design.m4, coef = 3)
top.bp.m5=getDEgenes(voom = v.noNA, design = design.m5, coef = 3)

top.bp.m1.f=getDEgenes(voom = v.noNA.f, design = design.m1, coef = 3)
top.bp.m2.f=getDEgenes(voom = v.noNA.f, design = design.m2, coef = 3)
top.bp.m3.f=getDEgenes(voom = v.noNA.f, design = design.m3, coef = 3)
top.bp.m4.f=getDEgenes(voom = v.noNA.f, design = design.m4, coef = 3)
top.bp.m5.f=getDEgenes(voom = v.noNA.f, design = design.m5, coef = 3)

m.compare.aic.f=selectModel(v.noNA.f, designlist = list(m1=design.m1,m2=design.m2,m3=design.m3,m4=design.m4,m5=design.m5))
m.compare.bic.f=selectModel(v.noNA.f, designlist = list(m1=design.m1,m2=design.m2,m3=design.m3,m4=design.m4,m5=design.m5), criterion = "bic")

top.m1=getDEgenes(design = design.m1)
top.m2=getDEgenes(design = design.m2)
top.m3=getDEgenes(design = design.m3)
top.m4=getDEgenes(design = design.m4)
top.m5=getDEgenes(design = design.m5)

m.compare.aic=selectModel(v.noNA, designlist = list(m1=design.m1,m2=design.m2,m3=design.m3,m4=design.m4,m5=design.m5))
m.compare.bic=selectModel(v.noNA, designlist = list(m1=design.m1,m2=design.m2,m3=design.m3,m4=design.m4,m5=design.m5), criterion = "bic")

# calculate AIC and BIC
aic=data.frame(model=1:5, IC=colMeans(m.compare.aic$IC), sem=apply(m.compare.aic$IC,2,function(x){sd(x)/sqrt(length(x))}))
bic=data.frame(model=1:5, IC=colMeans(m.compare.bic$IC), sem=apply(m.compare.bic$IC,2,function(x){sd(x)/sqrt(length(x))}))
abic=data.frame(rbind(aic,bic),criterion=c(rep("AIC",5),rep("BIC",5)))

aic.f=data.frame(model=1:5, IC=colMeans(m.compare.aic.f$IC), sem=apply(m.compare.aic.f$IC,2,function(x){sd(x)/sqrt(length(x))}))
bic.f=data.frame(model=1:5, IC=colMeans(m.compare.bic.f$IC), sem=apply(m.compare.bic.f$IC,2,function(x){sd(x)/sqrt(length(x))}))
abic.f=data.frame(rbind(aic.f,bic.f),criterion=c(rep("AIC",5),rep("BIC",5)))

library(ggplot2)
classic_theme=theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = 'black',fill=NA), legend.key = element_rect(fill = "white"),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))

pdf("models1-5_aic_bic_filtgenes.pdf",w=6,h=4)
ggplot(abic.f, aes(x=model, y=IC, colour=criterion)) +
  geom_errorbar(aes(ymin=IC-sem, ymax=IC+sem), width=.1, size=.7, colour="black") + geom_line(size=1) + classic_theme
dev.off()

# write out top genes for Dx differential gene expression
scz_de=list(m1=top.m1,m2=top.m2,m3=top.m3,m4=top.m4,m5=top.m5)
scz_de.f=list(m1=top.m1.f,m2=top.m2.f,m3=top.m3.f,m4=top.m4.f,m5=top.m5.f)

bp_de=list(m1=top.bp.m1,m2=top.bp.m2,m3=top.bp.m3,m4=top.bp.m4,m5=top.bp.m5)
bp_de.f=list(m1=top.bp.m1.f,m2=top.bp.m2.f,m3=top.bp.m3.f,m4=top.bp.m4.f,m5=top.bp.m5.f)


for (i in 1:5){
  write.csv(scz_de.f[[i]][,-(1:6)], paste("diffexp_scz/top_scz_filtgenes_m",i,".csv",sep=""))
}

for (i in 1:5){
  write.csv(bp_de[[i]][,-(1:6)], paste("diffexp_bp/top_bp_allgenes_m",i,".csv",sep=""))
}
