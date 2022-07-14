# Descriptif to build dataset for h3gwas using 1000 Genome

## build dataset as array
### extraction and format data positions
we download positions of [h3abionet array](https://www.h3abionet.org/h3africa-chip) using [manifest](https://drive.google.com/file/d/14QBlpXjw2U4rxhDk6va4mG6mXGRdTRNX/view)
```
unzip H3Africa_2017_20021485_A3.zip
awk -F"," '{if($9==37){
 if($10!='XY' && $10!='0' && $10!='MT')print $10"\t"$11}
}' H3Africa_2017_20021485_A3.csv|sort|uniq > h3abionet.37.pos
```

## Extraction of 1000 Genome positions
* we used [build_example_data/main.nf](https://github.com/h3abionet/h3agwas/tree/master/utils/build_example_data) script to build dataset
 * positions are diretly donwload from [1000 Genome project low coverage](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/)
```
nextflow run h3abionet/h3agwas/utils/build_example_data/main.nf --pos_allgeno h3abionet.37.pos -resume -profile slurmSingularity --list_chro 1-22,X --list_chro_pheno 1-22 --output_dir KGPH3abionet --output KGPH3abionet --nb_snp 3 --simu_hsq 0.7
```

 * if you have already download kgp, you can used `list_vcf` parameter with `ftp1000genome 0` :

```
ls folderkgp/*.vcf.gz|grep chr > listvcf
~/nextflow run $dirh3agwas/h3agwas/utils/build_example_data/main.nf --pos_allgeno h3abionet.37.pos -resume -profile slurmSingularity --list_chro 1-22,X --list_chro_pheno 1-22 --output_dir KGPH3abionet --output KGPH3abionet  --nb_snp 3 --simu_hsq 0.7 --list_vcf listvcf --ftp1000genome 0

```

###Output of pipeline:
 * `geno_all` : contains final  genotype
  * `geno_all/KGPH3abionet[.bed/fam/bim]` : plink file contains positions from array
   * `geno_all/KGPH3abionet_sex` : sex of individual after random  
   * `KGPH3abionet_sexsex_change_plk.ind` : list of individuals, where sex has been randomly changed
 * `gwascat` : folder contains gwas catalog download :
  * `KGPH3abionet/gwascat/gwascat_format_all.csv` : all positions from gwas cat 
  * `KGPH3abionet/gwascat/gwascat_format.[pos,bed]`  : positions associated to phenotype diabete
  * `KGPH3abionet/gwascat/gwascat_format_resume.csv` : information relative to diabetes
 * `simul_pheno` :
  * `KGPH3abionet/simul_pheno/datai/KGPH3abionet[.bed,plink,fam]` : plink contains all position extracted from gwas catalo link to phenotype
  * `KGPH3abionet/simul_pheno/data_format` :
   * positions clumped to avoid LD between positions extracted of phenotype :
   *  `KGPH3abionet.effect.rs_clump.clumped` : positions clumped to avoid LD between positions extracted of phenotype 
   *  `KGPH3abionet.effect.rs.assoc`  KGPH3abionet.effect.rs_clump.log      KGPH3abionet.effect.rs.infors
  * `KGPH3abionet.effect.rs` : effect of each position to build phenotype
  * `KGPH3abionet.effect.rs.infors` : information relative of positions to build phenotype 

## qc for imputation 
 *  we apply a Quality control before imputation 
```
~/nextflow run h3abionet/h3agwas/utils/qc/main.nf -profile slurmSingularity -c utils/params_1000G_h3abionet_qc.params
```
* report can be obtained [here](out_example/qc/KGPH3abionet_qc.pdf)

* `KGPH3abionet_qc/KGPH3abionet_qc.[bed,fam,bim]` : final file 
* `KGPH3abionet_qc/KGPH3abionet_qc.pdf` : report of qc


## Prepared data for imputation
* algorithm :
  * pipeline
```
cd ./utils_data/
wget -c http://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
wget -c ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz
cd ../
nextflow $dirh3agwas/h3agwas/formatdata/plk_in_vcf_imp.nf --input_dir   KGPH3abionet_qc --input_pat KGPH3abionet_qc --reffasta utils_data/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz --output KGPH3abionet --output_dir KGPH3abionet_vcf -profile slurmSingularity --file_ref_gzip utils_data/All_20180423.vcf.gz
```

* `KGPH3abionet_vcf/vcf/KGPH3abionet.vcf.gz` : data format in vcf to used for imputation
* pipeline do a post checking of final imputation file  :
 * [fixref plugin of bcftools]((https://samtools.github.io/bcftools/howtos/plugin.fixref.html) :
   * KGPH3abionet.checkbcf.err : contains position with error, if some positions has issue must be deleted
 * `KGPH3abionet_vcf/check/CheckVCF/` :
  * [CheckVCF.py](https://genome.sph.umich.edu/wiki/CheckVCF.py)
   * [out].check.af  
   * .check.dup  
   * .check.geno  
   * .check.log  
   * .check.mono 
   * .check.nonSnp  
   * .check.ref
  
 
## Imputation 
* Imputation done using sanger institute [see more explanation](documemtation/imputation_sangerinstitute.pdf)

## format vcf imputed in plink

```
nextflow run $dirh3agwas/h3agwas/formatdata/vcf_in_plink.nf --file_listvcf listvcf_imputed --output_pat KGPH3abionet_imputed --output_dir KGPH3abionet_imputed --reffasta utils_data/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz -profile slurmSingularity
```

* pipeline output :
 * bcffile : `ls KGPH3abionet_imputed/bcf_filter/*.bcf` :  vcf file in bcf
 * report in pdf on frequency distribution and score imputation [see output](documentation/)
 * plink file merged : `KGPH3abionet_imputed/KGPH3abionet_imputed.[bim/bed/fam]`
 * `KGPH3abionet_imputed/KGPH3abionet_imputed_report.pdf` : pdf report

##Sub selection of dataset
### individual selection
* script vselecting individual in function of randomisation of sex, super pop of 1000 genome, and individual deleted during first qc
```
##R script to selected 500 individuals  + 50 can be QC
Rscript utils/extractsubind_forqc1000genome.r
```

### Simulation of new phenotype 

```
plink --keep utils/subpopnoqc_ind.tsv -bfile KGPH3abionet_imputed/KGPH3abionet_imputed_sub --keep-allele-order -out KGPH3abionet_imputed/KGPH3abionet_imputed_sub -maf 0.05 -make-bed
```

### GWAS and clump : run a gwas on all phenotype 
* we ran a GWAS on phenotype selection 
```
plink --keep utils/subpopnoqc_ind.tsv -bfile KGPH3abionet_imputed/KGPH3abionet_imputed_sub --keep-allele-order -out KGPH3abionet_imputed/KGPH3abionet_imputed_sub -maf 0.05 -make-bed

~/nextflow $dirh3agwas/h3agwas/assoc/main.nf --data KGPH3abionetsub_pheno/simul_pheno/quant_pheno/KGPH3abionetsub_pheno_qt.pheno --pheno pheno_1   -profile slurmSingularity --gemma 1 --sample_snps_rel 1 --input_dir KGPH3abionet_imputed/ --input_pat KGPH3abionet_imputed_sub

echo -e "CHR\tSNP\tBP\tA1\tA2\tP" > pheno_1_sub/gemma/KGPH3abionet_imputed_sub-pheno-1.plink
sed '1d' pheno_1_sub/gemma/KGPH3abionet_imputed_sub-pheno-1.gemma |awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$12}'  >> pheno_1_sub/gemma/KGPH3abionet_imputed_sub-pheno-1.plink
##computed clump of best 
plink -bfile KGPH3abionet_imputed/KGPH3abionet_imputed_sub --clump output/gemma/KGPH3abionet_imputed_sub-pheno-1.plink -out output/clump_sim1_pheno1 --clump-p1 0.00000005 --clump-r2 0.1 --clump-kb 500 --maf 0.05

plink -bfile KGPH3abionet/geno_all/KGPH3abionet --freq -out KGPH3abionet/geno_all/sim1_pheno1 --keep utils/subpopnoqc_ind.tsv

```

### selection of positions test 

```
Rscript utils/extractsubsnp_forqc1000genome.r
```

### build dataset final 

```
## build dataset final  : plink array for qc 
plink --extract range utils/listsnparray_pheno1.tokeep --keep utils/subpop_ind.tsv -bfile KGPH3abionet/geno_all/KGPH3abionet  -make-bed -out data/array_plk/array

## build dataset final : plink for gwas
plink --extract range utils/listsnparray_pheno1.tokeep --keep utils/subpopnoqc_ind.tsv -bfile KGPH3abionet_imputed/KGPH3abionet_imputed  -make-bed -out data//imputed/imput_data
## Extraction of vcf 
DirVcf=KGP.vcfs

## vcf extracted for imputation
plkinfo=data//imputed/imput_data
module load vcftools
awk '{print $1"\t"$4}' $plkinfo".bim" > listpos
awk '{print $1"\_"$2}' $plkinfo".fam" > listind
DirOut=data//imputed/vcf/
mkdir -p $DirOut
for File in `ls $DirVcf/*.vcf.gz`
do
fileF=$DirOut/`basename $File`
vcftools  --gzvcf $File --positions listpos  --keep listind --recode --recode-INFO-all  --stdout | bgzip -c > $fileF &
done
wait
```

### for phenotype 

we resume all information to build phenotype 

```
utils/buildpheno_test.r
```

### doing summary statistic


```
for head in AFR all AMR EAS EUR SAS
do
filepheno=data/pheno/pheno_test.$head
~/nextflow $dirh3agwas/h3agwas/assoc/main.nf --data $filepheno --pheno pheno_qt1 -profile slurmSingularity --gemma 1 --sample_snps_rel 1 --input_dir data/imputed/ --input_pat imput_data  --output_dir gwas/$head -resume
cp gwas/$head/gemma/imput_data-pheno-qt1.assoc.txt  data/summarystat/$head"_pheno.gemma"
done

head=all
filepheno=data/pheno/pheno_test.$head
~/nextflow $dirh3agwas/h3agwas/assoc/main.nf --data $filepheno --pheno pheno_qt2 -profile slurmSingularity --gemma 1 --sample_snps_rel 1 --input_dir data/imputed/ --input_pat imput_data  --output_dir gwas/$head
```

### bgen file 
```
~/nextflow ~/Travail/git/h3agwas/formatdata/vcf_in_bgen_merge.nf --file_listvcf listvcf --output_pat imput_data  --output_dir ./data/imputed/bgen -profile slurmSingularity -resume --bgen_type bgen_v1.2  --other_opt  \"-bgen-bits 8\"
```
### other file 
we build a subsample file for test data
```
poshead_chro_inforef=0
poshead_bp_inforef=1
poshead_rs_inforef=2
poshead_a1_inforef=3
poshead_a2_inforef=4
cat  data/imputed/imput_data.bim data/array_plk/array.bim|awk '{print $1"\t"$4"\t"$5"\t"$6}' > temp_pos
rsinfo=/home/jeantristan/Data/human_9606/All_20180423.vcf.gz
outrs=all_rsinfo
zcat $rsinfo | /home/jeantristan/Travail/git/h3agwas/formatdata/bin/extractrsid_bypos.py --file_chrbp temp_pos --out_file $outrs  --ref_file stdin --chro_ps ${poshead_chro_inforef} --bp_ps ${poshead_bp_inforef} --rs_ps ${poshead_rs_inforef} --a1_ps ${poshead_a1_inforef}  --a2_ps ${poshead_a2_inforef}
```


