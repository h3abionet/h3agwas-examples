Data<-read.table('h3agwas-examples/data/pheno/pheno_test.all',header=T)
Data$pheno_qt1_miss<-Data$pheno_qt1
nbrow<-nrow(Data)
nbmissing<-round(nbrow*0.05)
Data$pheno_qt1_miss[sample(1:nbrow, nbmissing)]<-NA
write.table(Data, file='h3agwas-examples/data/pheno/pheno_test_missing.all',sep='\t', quote=F, row.names=F, col.names=T)

