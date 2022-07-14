library(data.table)
set.seed(255) 
## inforamtion relative to the sex from pipeline
dataind<-read.csv('KGPH3abionet/geno_all/KGPH3abionet_sexsex_change_plk.ind')
## select individual where has been change
dataind$notgoodsex<- (dataind$gender=='female' & dataind$sex==1) | (dataind$gender=='male' & dataind$sex==2) 
## we read initial datasex after qc
datai_fam<-read.table('KGPH3abionet_qc/KGPH3abionet_qc.fam')
## we selected in total 550 individual
### 15 individuals with secx
indprobsex<-sample(dataind$sample[dataind$notgoodsex],15)
#select other individual has been qc 
all<-merge(dataind,datai_fam,by.y=1, by.x='sample',all=T)
all$deleted_duringqc<-is.na(all$V2)
## 35 individuals has been deleted in preivous dataset
#indotherqc<-sample(all$sample[!all$notgoodsex & is.na(all$V2)],35)
## select 500 individuals
## utils/superPop.info
info1000G<-read.table('utils/20130606_g1k.ped',sep='\t', header=T)
all<-merge(all, info1000G, by.x='sample', by.y='Individual.ID', all.x=T)
#   sample X pop super_pop gender sex_plinkf notgoodsex   V2 V3 V4 V5 V6
#1 HG00096 1 GBR       EUR   male          2       TRUE <NA> NA NA NA NA
#  deleted_duringqc   qc Family.ID Paternal.ID Maternal.ID Gender Phenotype
#1             TRUE TRUE   HG00096           0           0      1         0
#  Population Relationship Siblings Second.Order Third.Order Other.Comments
all2<-all[,c("sample", "pop", "super_pop", "gender", "sex_plinkf", "notgoodsex","deleted_duringqc", "Family.ID")]

tbfmid<-table(all2$Family.ID)
listfamilydup<-names(tbfmid[tbfmid>1])

nbindbypop<-round(table(all$super_pop)/length(all$super_pop)*500)
nbindbypopkeepqc<-round(table(all$super_pop)/length(all$super_pop)*50)


familyrelatd<-all2[all2$Family.ID %in% listfamilydup,c('super_pop', 'Family.ID', 'deleted_duringqc')]
nbdelbyfamilyid<-aggregate(deleted_duringqc~super_pop+Family.ID, data=familyrelatd, function(x)length(which(x))>0 & length(which(!x))>0)
familyrelatd<-nbdelbyfamilyid[nbdelbyfamilyid$deleted_duringqc,]
all2$keep<-F
## we keep some individual with sex probleme
all2$keep[all2$sample %in% indprobsex]<-T
all2$sample<-as.character(all2$sample)
#all2$relatd[all2$sample %in% familyrelatd]<-T
for(pop in unique(all2$super_pop)){
 ##sselect some individual with some duplicate id
 listfamid<-sample(familyrelatd[familyrelatd$super_pop==pop,'Family.ID'], 2)
 listsample<-as.vector(unlist(sapply(listfamid, function(x)c(all2$sample[(all2$Family.ID %in% x) & all2$deleted_duringqc & !all2$keep][1], all2$sample[(all2$Family.ID %in% x) & !all2$deleted_duringqc & !all2$keep][1]) )))
 all2$keep[all2$Family.ID %in% listsample]<-T
 nbkeep<-length(which(!all2$deleted_duringqc & all2$super_pop==pop & all2$keep))
 nbdeleted<-length(which(all2$deleted_duringqc & all2$super_pop==pop & all2$keep))
 nbsample<-nbindbypop[names(nbindbypop)==pop] - nbkeep
 nbsampleqc<-nbindbypopkeepqc[names(nbindbypopkeepqc)==pop] - nbdeleted
 listidgoodpop<-sample(as.character(all2$sample[!all2$deleted_duringqc & all2$super_pop==pop & !all2$keep]) , nbsample)
 if(nbsampleqc>0)listidqcpop<-sample(as.character(all2$sample[all2$deleted_duringqc & all2$super_pop==pop & !all2$keep]) , nbsampleqc) else listidqcpop<-c()
 ## select individual   
 all2$keep[all2$sample %in% c(listidgoodpop, listidqcpop)]<-T
}
head<-c("sample", "sample","pop", "super_pop", "gender", "sex_plinkf", "notgoodsex","deleted_duringqc", "Family.ID", "keep")
all3<-all2[,head]
names(all3)<-c("FID", "IID", "pop", "super_pop", "gender", "sex_plinkf", "notgoodsex","deleted_duringqc", "Family.ID", "keep")

write.table(all3[all3$keep,], file='utils/subpop_ind.tsv', sep='\t', quote=F, col.names=T, row.names=F)
write.table(all3[all3$keep & !all2$deleted_duringqc,], file='utils/subpopnoqc_ind.tsv', sep='\t', quote=F, col.names=T, row.names=F)






#allindividual<-c(indprobsex,indotherqc,indnofilter)
#all$deleted_duringqc<-is.na(all$V2)
#sampleind<-all[all$sample %in% allindividual,]
#sampleind$sample2<-sampleind$sample
#sampleind<-sampleind[,c('sample', )]
##write.table(all[,])
#
### write a infosex
#if(!file.exists('infosample_sub.tsv'))write.table(AllPhenoSub, file='infosample_sub.tsv', sep='\t',quote=F, row.names=F)
#AllPhenoSub<-read.table('infosample_sub.tsv', header=T)
#
### sub select individual from imputed
#if(!file.exists("sub_sim_s1_s.bed"))system("plink -bfile /dataF/synthetic/raw/20k/sim1/imputed_sim/imp_plk_newid --keep  infosample_sub.tsv --out sub_sim_s1_s --make-bed --keep-allele-order --threads 5 --maf 0.05")
### sub select individual from array : keep X
#if(!file.exists("sub_sim_array_s.bed"))system("plink -bfile sim1_s_sex --keep  infosample_sub.tsv --out sub_sim_array_s --make-bed --keep-allele-order --threads 5 --maf 0.05")
### simulate phenotype
#system("/home/jeantristan/Travail/git/h3agwas/utils/build_example_data/simul-assoc_gcta.nf --input_pat sub_sim_s1_s --input_dir ./ -resume --simu_k 0.5 --simu_hsq 0.99  --nb_snp 5")
### run gwas on all 
#system("/home/jeantristan/Travail/git/h3agwas/assoc/main.nf --input_pat sub_sim_s1_s_gc --input_dir output/simul_pheno/pheno_format/ -resume --data output/simul_pheno/quant_pheno/out_qt.pheno --pheno pheno_1 --gemma 1 --sample_snps_rel 0 --gemma_multi 0 ")
#
## extracted positions of interrest from gws
#DataRes<-fread('output/gemma/sub_sim_s1_s_gc-pheno-1.assoc.txt')
#range(DataRes$p_wald);table(DataRes$p_wald<5*10**-8)
#DataPlink<-DataRes[,c('chr', 'rs','ps', 'allele1', 'allele0', 'p_wald')]
#names(DataPlink)<-c('CHR', 'SNP','BP', 'A1','A2','P')
#write.table(DataPlink, file='sim1_pheno1.assoc', row.names=F, col.names=T, quote=F,sep='\t')
#
### clump result
#if(!file.exists('sub_sim_s1_s'))system('plink -bfile sub_sim_s1_s --clump sim1_pheno1.assoc -out sim1_pheno1 --clump-p1 0.00000005 --clump-r2 0.1 --clump-kb 500 --maf 0.05')
### read clump result
#DataClp<-read.table('sim1_pheno1.clumped',header=T)
### choice of best solutions
#DataClpbest<-DataClp[order(DataClp$TOTAL, decreasing=T),][1:3,]
##listrs<-c(as.character(DataClp$SNP))
#
#around=10**5
#databim<-fread('sub_sim_s1_s.bim')
#listrsimp<-c()
#for(CmtWind in 1:nrow(DataClpbest)){
#Chr<-DataClpbest[CmtWind,'CHR'];BP<-DataClpbest[CmtWind,'BP']
#listrsimp<-unique(c(listrs,databim$V2[databim$V1==Chr & databim$V4>=BP-around & databim$V4<=BP+around]))
#}
### write 
#writeLines(listrs, con='list_rsimputed')
### extract rs from   imputed
#system('plink -bfile sub_sim_s1_s --keep-allele-order --extract list_rsimputed -out sub_imputed_s1_s --make-bed')
#
#
### computed frequencie for  
#if(!file.exists('sim1_pheno1.frq'))system('plink -bfile sub_sim_array_s --freq -out sim1_pheno1 --keep infosample_sub.tsv')
#DataFreq<-read.table('sim1_pheno1.frq', header=T)
#DataFreqX<-DataFreq[DataFreq$CHR==23,]
#listrsx<-c(sample(as.character(DataFreqX$SNP[DataFreqX$MAF>0.05]), 1500))
#
###
##listrs<-c(listrs,sample(DataFreq[(DataFreq$CHR %in% 1:22)  & DataFreq$MAF>0.05 & !(DataFreq$SNP %in% listrs),'SNP'], 20000))
#listrsarray<-c(listrsx,sample(DataFreq[(DataFreq$CHR %in% 1:22)  & !(DataFreq$SNP %in% listrsimp),'SNP'], 50000 - length(c(listrsimp,listrsx))))
#
#writeLines(listrsarray, con='listrsarray.tokeep')
#
#system('plink -bfile sim1_s_sex --keep infosample_sub.tsv --extract listrsarray.tokeep --make-bed -out sub_array_s1_s --update-sex infosample_sub.tsv 1')
#system('plink -bfile sub_array_s1_s -bmerge sub_imputed_s1_s --make-bed -out ../../data2/exampledata2 --update-sex infosample_sub.tsv 1')
#system('plink -bfile ../../data2/exampledata2 --out exampledata2 --K 3 --cluster')
#system('plink -bfile ../../data2/exampledata2 --out exampledata2 --pca')
#
#AllPhenoSub<-read.table('infosample_sub.tsv', header=T)
#DataPca<-read.table('exampledata2.eigenvec');names(DataPca)<-c('FID','IID', paste('Pcs_',1:(ncol(DataPca)-2),sep=''))
#DataCluster<-read.table('exampledata2.cluster2');names(DataCluster)<-c('FID','IID', 'ethnicity_predicted')
#DataPcsClus<-merge(DataPca, DataCluster)
#plot(DataPcsClus$Pcs_1, DataPcsClus$Pcs_2,col=DataPcsClus$Cluster+1)
#
#Pheno<-read.table('output/simul_pheno/quant_pheno/out_qt.pheno' ,header=T)[,c(1,2,3,4,5)]
#PhenoQT<-read.table('output/simul_pheno/qual_pheno/out_ql.pheno' ,header=T)[,c(1,2,3,4)];names(PhenoQT)<-c('FID','IID','CaseControl1','CaseControl2')
#
#names(DataCluster)<-c('FID','IID', 'Pop')
#DataPheno2<-merge(merge(merge(AllPhenoSub,Pheno,by=c(1,2),all=T),DataCluster,by=c(1,2),all=T),PhenoQT,by=c(1,2),all=T)
#write.table(DataPheno2, file='../../data2/exampledata2.pheno', row.names=F, col.names=T, quote=F,sep='\t')
#write.table(DataPheno2[DataPheno2$Pop==0,], file='../../data2/exampledata2.pop0.pheno', row.names=F, col.names=T, quote=F,sep='\t')
#write.table(DataPheno2[DataPheno2$Pop==1,], file='../../data2/exampledata2.pop1.pheno', row.names=F, col.names=T, quote=F,sep='\t')
#write.table(DataPheno2[DataPheno2$Pop==2,], file='../../data2/exampledata2.pop2.pheno', row.names=F, col.names=T, quote=F,sep='\t')
#
#a
