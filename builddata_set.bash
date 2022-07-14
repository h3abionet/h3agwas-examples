
dirh3agwas=h3abionet
dirh3agwas=~/Travail/git
if [ ! -f h3abionet.37.pos ]
then
## download H3Africa_2017_20021485_A3.zip on https://drive.google.com/file/d/14QBlpXjw2U4rxhDk6va4mG6mXGRdTRNX/view
# unzip and format
mkdir -p utils_data/
mv H3Africa_2017_20021485_A3.zip utils_data/
unzip utils_data/H3Africa_2017_20021485_A3.zip
awk -F"," '{if($9==37){
 if($10!='XY' && $10!='0' && $10!='MT')print $10"\t"$11}
}' utils_data/H3Africa_2017_20021485_A3.csv|sort|uniq > utils_data/h3abionet.37.pos
fi

# run 

ls /spaces/jeantristan/Data/1000Geno/AllVCF_Genome/*.vcf.gz|grep chr > listvcf
~/nextflow run $dirh3agwas/h3agwas/utils/build_example_data/main.nf --pos_allgeno h3abionet.37.pos -resume -profile slurmSingularity --list_chro 1-22,X --list_chro_pheno 1-22 --output_dir KGPH3abionet --output KGPH3abionet  --nb_snp 3 --simu_hsq 0.7 --list_vcf listvcf --ftp1000genome 0

# nexflow qc of genetics data
~/nextflow run $dirh3agwas/h3agwas/qc/main.nf -profile slurmSingularity -c utils/params_1000G_h3abionet_qc.params

#format output  in vcf for imputation
cd ./utils_data/
wget -c http://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
wget -c ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz  
cd ../
nextflow $dirh3agwas/h3agwas/formatdata/plk_in_vcf_imp.nf --input_dir   KGPH3abionet_qc --input_pat KGPH3abionet_qc --reffasta utils_data/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz --output_pat KGPH3abionet --output_dir KGPH3abionet_vcf -profile slurmSingularity --file_ref_gzip utils_data/All_20180423.vcf.gz  -resume

## tranfert data for imputation 

## after imputation format in plink dl in folder : KGP.vcfs/
ls KGP.vcfs//*.vcf.gz > listvcf_imputed 
~/nextflow run $dirh3agwas/h3agwas/formatdata/vcf_in_plink.nf --file_listvcf listvcf_imputed --output_pat KGPH3abionet_imputed --output_dir KGPH3abionet_imputed --reffasta utils_data/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz -profile slurmSingularity  --max_plink_cores 30

## sub sample of individual
```
##R script to selected 500 individuals 
Rscript utils/extractsubind_forqc1000genome.r
```

##build a new pheno
cd KGPH3abionet_imputed/
cp KGPH3abionet_imputed.fam  KGPH3abionet_imputed.save.fam
awk -F"[ _]" '{print $2" "$3" "$4" "$5" "$6" "$7}' KGPH3abionet_imputed.save.fam > KGPH3abionet_imputed.fam
cd ../
~/nextflow $dirh3agwas/h3agwas/utils/build_example_data/simul-assoc_gcta.nf --input_dir KGPH3abionet_imputed/  --input_pat KGPH3abionet_imputed --cut_maf 0.05 --keep utils/subpopnoqc_ind.tsv --simu_k 0.5 --simu_hsq 0.99  --nb_snp 5 -profile slurmSingularity  --output_dir KGPH3abionetsub_pheno --output KGPH3abionetsub_pheno



## we run 
plink --keep utils/subpopnoqc_ind.tsv -bfile KGPH3abionet_imputed/KGPH3abionet_imputed_sub --keep-allele-order -out KGPH3abionet_imputed/KGPH3abionet_imputed_sub -maf 0.05 -make-bed

~/nextflow $dirh3agwas/h3agwas/assoc/main.nf --data KGPH3abionetsub_pheno/simul_pheno/quant_pheno/KGPH3abionetsub_pheno_qt.pheno --pheno pheno_1   -profile slurmSingularity --gemma 1 --sample_snps_rel 1 --input_dir KGPH3abionet_imputed/ --input_pat KGPH3abionet_imputed_sub --gemma_multi 1 --output_dir 

## format for plink
echo -e "CHR\tSNP\tBP\tA1\tA2\tP" > pheno_1_sub/gemma/KGPH3abionet_imputed_sub-pheno-1.plink
sed '1d' pheno_1_sub/gemma/KGPH3abionet_imputed_sub-pheno-1.gemma |awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$12}'  >> pheno_1_sub/gemma/KGPH3abionet_imputed_sub-pheno-1.plink
##computed clump of best 
plink -bfile KGPH3abionet_imputed/KGPH3abionet_imputed_sub --clump output/gemma/KGPH3abionet_imputed_sub-pheno-1.plink -out output/clump_sim1_pheno1 --clump-p1 0.00000005 --clump-r2 0.1 --clump-kb 500 --maf 0.05

plink -bfile KGPH3abionet/geno_all/KGPH3abionet --freq -out KGPH3abionet/geno_all/sim1_pheno1 --keep utils/subpopnoqc_ind.tsv 


##  we choose positions for gwas and around 50000
Rscript utils/extractsubsnp_forqc1000genome.r

## build dataset final 
plink --extract range utils/listsnparray_pheno1.tokeep --keep utils/subpop_ind.tsv -bfile KGPH3abionet/geno_all/KGPH3abionet  -make-bed -out data/array_plk/array

## for gwas
plink --extract range utils/listsnparray_pheno1.tokeep --keep utils/subpopnoqc_ind.tsv -bfile KGPH3abionet_imputed/KGPH3abionet_imputed  -make-bed -out data//imputed/imput_data

## Extraction of vcf 
DirVcf=KGP.vcfs

plkinfo=data//imputed/imput_data
module load vcftools
awk '{print $1"\t"$4}' $plkinfo".bim" > listpos
awk '{print $1"_"$2}' $plkinfo".fam" > listind
DirOut=data//imputed/vcf/
mkdir -p $DirOut
for File in `ls $DirVcf/*.vcf.gz`
do
fileF=$DirOut/`basename $File`
vcftools  --gzvcf $File --positions listpos  --keep listind --recode --recode-INFO-all  --stdout | bgzip -c > $fileF &
done
wait

## qc on sub
~/nextflow run $dirh3agwas/h3agwas/qc/main.nf -profile slurmSingularity -c utils/params_1000G_h3abionet_qc_subtest.params 

## finalise phenotype 
Rscript utils/buildpheno_test.r

#building summary statistics


#pheno_test.AFR  pheno_test.all  pheno_test.AMR  pheno_test.EAS  pheno_test.EUR  pheno_test.SAS
for head in AFR all AMR EAS EUR SAS
do
filepheno=data/pheno/pheno_test.$head
~/nextflow $dirh3agwas/h3agwas/assoc/main.nf --data $filepheno --pheno pheno_qt1 -profile slurmSingularity --gemma 1 --sample_snps_rel 1 --input_dir data/imputed/ --input_pat imput_data  --output_dir gwas/$head -resume 
cp gwas/$head/gemma/imput_data-pheno-qt1.assoc.txt  data/summarystat/$head"_pheno.gemma"
done

head=all
filepheno=data/pheno/pheno_test.$head
~/nextflow $dirh3agwas/h3agwas/assoc/main.nf --data $filepheno --pheno pheno_qt2 -profile slurmSingularity --gemma 1 --sample_snps_rel 1 --input_dir data/imputed/ --input_pat imput_data  --output_dir gwas/$head
cp gwas/$head/gemma/imput_data-pheno-qt2.assoc.txt  data/summarystat/$head"_phenoq2.gemma"


poshead_chro_inforef=0
poshead_bp_inforef=1
poshead_rs_inforef=2
poshead_a1_inforef=3
poshead_a2_inforef=4
cat  data/imputed/imput_data.bim data/array_plk/array.bim|awk '{print $1"\t"$4"\t"$5"\t"$6}' > temp_pos
rsinfo=/home/jeantristan/Data/human_9606/All_20180423.vcf.gz
outrs=all_rsinfo
zcat $rsinfo | /home/jeantristan/Travail/git/h3agwas/formatdata/bin/extractrsid_bypos.py --file_chrbp temp_pos --out_file $outrs  --ref_file stdin --chro_ps ${poshead_chro_inforef} --bp_ps ${poshead_bp_inforef} --rs_ps ${poshead_rs_inforef} --a1_ps ${poshead_a1_inforef}  --a2_ps ${poshead_a2_inforef}

ls data/imputed/vcf/*.vcf.gz > listvcf
~/nextflow ~/Travail/git/h3agwas/formatdata/vcf_in_bgen_merge.nf --file_listvcf listvcf --output_pat imput_data  --output_dir ./data/imputed/bgen -profile slurmSingularity -resume --bgen_type bgen_v1.2 



