

indbefore_qc<-read.table('data/array_plk/array.fam')
indafter_qc<-read.table('data/imputed/imput_data.fam')
indafter_qc2<-read.table('KGPH3abionetsub_qc/KGPH3abionet_qc.fam')

dataphenoqc<-read.table('KGPH3abionet/simul_pheno/qual_pheno/KGPH3abionet_ql.pheno', header=T)
dataphenogwas<-read.table('KGPH3abionetsub_pheno/simul_pheno/quant_pheno/KGPH3abionetsub_pheno_qt.pheno', header=T)


dataqc<-read.table('utils/subpop_ind.tsv', header=T)
dataqc<- dataqc[,c('FID','IID','pop','super_pop','gender','sex_plinkf','notgoodsex','Family.ID')]
names(dataqc)<-c('FID','IID','pop','super_pop','gender1000G','Sex','SexChange','Family.ID')

dataphenoqc2<-dataphenoqc[,c('FID', 'IID', 'pheno_1')];names(dataphenoqc2)[3]<-'phenoqc_ql'

dataqc2<-merge(dataqc,dataphenoqc2, by=c('FID', 'IID'))

## 
dataphenogwas<-dataphenogwas[,c('FID', 'IID','pheno_1','pheno_2')];names(dataphenogwas)<-c('FID', 'IID', 'pheno_qt1', 'pheno_qt2')
dataqc3<-merge(dataqc2,dataphenogwas,by=c('FID', 'IID'), all=T)
getrandomval<-function(x){
y<-x[!is.na(x)]
x[is.na(x)]<-rnorm(length(which(is.na(x))), mean(y), sd(y))
x
}

dataqc3$pheno_qt1<-getrandomval(dataqc3$pheno_qt1)
dataqc3$pheno_qt2<-getrandomval(dataqc3$pheno_qt2)


dataqc3$gc10<-rnorm(nrow(dataqc3),0.8,1)
idtochange<-as.character(indafter_qc2$V1[!(indafter_qc2$V1 %in% dataphenogwas$FID)])
dataqc3$gc10[dataqc3$FID %in% idtochange]<-rnorm(length(idtochange),0.1,0.2)
dataqc3$batch<-sample(1:3, nrow(dataqc3), replace=T)


write.table(dataqc3, row.names=F, quote=F, col.names=T, sep='\t', file="data/pheno/pheno_test.all")

for(Pop in unique(dataqc3$super_pop)){
write.table(dataqc3[dataqc3$super_pop==Pop,], row.names=F, quote=F, col.names=T, sep='\t', file=paste("data/pheno/pheno_test.", Pop,sep=''))
}


