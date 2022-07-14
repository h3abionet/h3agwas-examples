library(data.table)

### read clump result
dataclump<-read.table('pheno_1_sub/clump_sim1_pheno1.clumped',header=T)
databim_array<-read.table('KGPH3abionet/geno_all/KGPH3abionet.bim')
DataFreq<-read.table('KGPH3abionet/geno_all/sim1_pheno1.frq', header=T)
databim<-fread('KGPH3abionet_imputed/KGPH3abionet_imputed_sub.bim')
listrs<-c(as.character(dataclump$SNP))
### choice of best solutions
dataclumpbest<-dataclump[order(dataclump$TOTAL, decreasing=T),]
#listrs<-c(as.character(dataclump$SNP))
#
around=10**5
listrsimp<-c()
for(CmtWind in 1:nrow(dataclumpbest)){
Chr<-dataclumpbest[CmtWind,'CHR'];BP<-dataclumpbest[CmtWind,'BP']
listrsimp<-unique(c(listrs,databim$V2[databim$V1==Chr & databim$V4>=BP-around & databim$V4<=BP+around]))
}
postokeep<-databim[databim$V2 %in% listrsimp ,c("V1", "V4", "V4", "V2")]

##
DataFreqX<-DataFreq[DataFreq$CHR==23,]
listrsx<-c(sample(as.character(DataFreqX$SNP[DataFreqX$MAF>0.05]), 1500))

DataFreq$SNP<-as.character(DataFreq$SNP)
listrsarray<-c(listrsx,sample(DataFreq[(DataFreq$CHR %in% 1:22)  & !(DataFreq$SNP %in% listrsimp),'SNP'], 50000 - length(c(listrsx))))
databim_arraypos<-databim_array[databim_array$V2 %in% listrsarray,c("V1", "V4", "V4", "V2")]

names(postokeep)<-c("CHR", "BEGIN", "END", "RSID")
names(databim_arraypos)<-c("CHR", "BEGIN", "END", "RSID")
postokeep<-rbind(postokeep,databim_arraypos)

write.table(postokeep,file='utils/listsnparray_pheno1.tokeep', sep='\t', row.names=F, col.names=T, quote=F)

## data 





