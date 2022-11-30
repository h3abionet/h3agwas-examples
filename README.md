# H3AGWAS : dataset and command line example 

Example and command line to run [h3agwas pipeline](https://github.com/h3abionet/h3agwas)

These are very simple examples to show working of the pipeline and get you started.
These examples show the basic working of some of the key sub-workflows only. We can't show all options here -- see the main documentation for more.


# 1. Set-up


A dataset has been build using 1000 Genomes Project data, selecting  50,000 SNPs found in the H3Africa Custom Chip Array and 500 individuals. Description and explanation can be found [here](README_buildataset.md) or bash script [here](builddata_set.bash)

Information to run locally can be found [here](runlocal/ubuntu/README.md)

## 1.1  Data directory :
* The data set can be found in [data folder](data), subfolder :
 * `data/array_plk/`: contains file for qc and prepared for imputation in plink format
   * data/array_plk/array.{.bed,bim,fam} : file to test qc
   * data/array_plk/array_qc.{.bed,bim,fam} : file to test format plink in vcf to prepare imputation
 * `imputed` : contains file for gwas in plink format
 * `pheno/pheno_test.[pop]` : contains phenotype to perform gwas for each superpopulation of KGP :  AFR,AMR,EAS,EUR,SAS and all
   * FID : Family ID 
   * IID : Individual ID 
   * `super_pop` : super-population used to split population
   * Sex : sex after randomisation
   * sex1000G : sex of individual 1000 Genome
   * SexChange : False : mean sex phenotype is sex genotype T / TRUE : mean sex phenotype is different than sex Genotype. pipeline had change sex phenotype to created example for qc
   * phenoqc\_ql : binary phenotype used 
   * pheno\_qt[12] : quantitative phenotype used for GWAS
 * `summarystat/[pop]_pheno.gemma` : result of GWAS using GEMMA for phenotype 1
 * `summarystat/all_phenoq2.gemma` : result of GWAS using gwas for phenotype 2 and all
 * `utils/all_rsinfo.init.gz` : contains information relative to rsid / positions, subsample of [files here](shorturl.at/HJLN0)


## 1.2  To run these examples

1. Fetch or update the workflow itself `nextflow pull h3abionet/h3agwas`
2. Clone the example directory: `git clone h3abionet/h3agwas-examples`
3. Change your working directory so you are in the new directory `cd h3agwas-examples`
4. You need to have have Singularity or Docker installed (or install all the software tools yourself)


# 2. Quality control 

The first example takes raw PLINK data as input (that is data that comes after genotype calling) and performs quality control. The `qc` workflow is used. The key parameters are
* input / output :
 * `--input_dir` specifies directory where data can be found
 * `--input_pat` specifies the base of the PLINK bed, bim and fam file name
 * `--output_dir` corresponding at the directory output directory 


In this example, we run QC on the raw data found in `data/array_plk/array`

```
nextflow run h3abionet/h3agwas/qc/main.nf --input_dir data/array_plk  --input_pat array --output_dir qc  --output array_qc \
 --phenotype data/pheno/pheno_test.all --pheno_col phenoqc_ql \
 --case_control data/pheno/pheno_test.all --case_control_col Sex \
 --batch data/pheno/pheno_test.all --batch_col batch \
 -profile singularity -resume
```

The `-profile` option asks for the use of Singularity. The first time you run the workflow, Nextflow will automatically download all the singularity images -- and this may take a few minutes (this only happens the first time and there's nothing special you have to do other than be patient). If you want to use Docker instead of Singularity you would say `-profile docker`. If you want to use both Singularity and a scheduler like slurm you say in this order `-profile singularity,slurm`.


In the QC we look for batch (broadly understood) effects in 3 ways (you don't have to do such a complex analysis). We check for batch effects of the `phenoqc_ql` phenotype, sex of individuals, and the batch from which the sample came. Default QC cut-offs for missingness found in the `nextflow.config` file are used.

Output goes to the `qc/` directory with `kgexample` as base for file names.. 

* output :
  * `$dir_output/$output.[bed/fam/bim]` : PLINK files with resulting QC data
  * `${dir_output}/$output.pdf` : [report of qc, for example here](out_example/qc/kgpexample.pdf)
  * `$dir_output/pca/` : contains figure and plink file used 
  * `$dir_output/samples/` : contains figures, files relative to qc apply to sample
  * `$dir_output/snps/` : contains figures, files relative to qc apply to SNPs
  * `$dir_output/phase1/` : contains figures, files relative to qc apply at phase 1 of qc

The PDF file gives the report of the QC. Note that because this is such a small and artificial file some of the pictures will look very odd.


# 3. Association pipeline

The pipeline offers multiple ways of testing for association. Only some examples:

 * pipeline needs at least the PLINK  (see `--input_dir` and `--input_pat`), phenotype file (see `--data`) and 1 or more phenotype (`--phenotype`)

Here's a simple example, assuming you've done the QC example above and haven't changed working directory.  We take the QCed data as genotype input, and the data in `data/pheno/pheno_test.all` as the source of the phenotype data. There are two columns in the phenotype file (`pheno_qt1`, `pheno_qt2`) and we'll test for associations against both. We specify which directory and base file name to use for output. We do linear regression. We use singularity.

```
nextflow run  h3abionet/h3agwas/assoc/main.nf --input_dir qc/    --input_pat kgpexample  --data data/pheno/pheno_test.all   --pheno pheno_qt1,pheno_qt2  --output_dir assoc1 --output qt12 --linear 1 -profile  singularity
```

Now we do a more complex example, adding on to this
  * We use imputed data as input
  * We do linear regression and GEMMA
  * `--sample_snps_rel 1` : builds a relatedness matrix sub-sampling SNPs
  * We used dosage data in BGEN format

```
nextflow run h3abionet/h3agwas/assoc/main.nf \
 --input_dir data/imputed/ --input_pat imput_data \
 --data data/pheno/pheno_test.all --pheno pheno_qt1,pheno_qt2 \
 --output_dir assoc --output assoc \
 --gemma 1 --sample_snps_rel 1 --linear 1 \ 
  -profile singularity \
 --bgen data/imputed/bgen/out.bgen --bgen_sample data/imputed/bgen/out.sample
```

Output :
  * [in the report, we combined best results, manhantan plot and qq plot](out_example/assoc-report.pdf) 
  * Each software as in own folder with output of software


In the next example, we have separate BGEN files as input -- we can analyse them use the `list_bgen` option.  Here we also use BOLTLMM, FASTGWA and RGENIE. There are some advanced options given for these tools
```
ls data/imputed/bgen_chro/*.bgen > listbgen
nextflow run h3abionet/h3agwas/assoc --input_dir data/imputed/ --input_pat input_data \
 --data data/pheno/pheno_test.all --pheno pheno_qt1,pheno_qt2 \
 --output_dir assoc_listbgen --output assoc \
 --boltlmm 1 --sample_snps_rel 1 --regenie 1 --fastgwa 1 --grm_nbpart 2\
  -profile singularity \
 --list_bgen listbgen --bgen_sample data/imputed/bgen/out.sample --saige 1 -resume
```

## Input types for different tools


The table below shows the different data types, the information they store and how referenced.

|     | plink | vcf | bgen | impute 2 |
| | genotype |  dosage | dosage | dosage |
| --- | --- | --- | --- | --- |
| Option |  `--input_dir`/`--input_pat`| `--list_vcf` | `--list_bgen`/ `--bgen`/`--bgen_sample` | `bolt_impute2filelist`/`bolt_impute2fidiid` |



The table below shows different data types used as input for the supported tools and the command used to activate


| Software | plink | vcf | bgen | impute 2 | option to activate |
| --- | --- | --- | --- | --- | --- |
| gemma | yes |  no | no | no | --gemma  / --gemma\_gxe |
| plink association | yes |  no | no | no | --[linear,logistic,assoc",fisher] 1 / --plink\_gxe |
| gcta/fastGWA | yes |  no | yes | no | --fastGWA 1 |
| saige | yes |  yes | yes | no | --saige 1 |
| bolt-LMM | yes |  no | yes | yes | --boltlmm 1|
| fast-lmm | yes |  no | no | no | --fastlmm 1 |
| regenie | yes |  no | yes | no | --regenie 1 |
| --- | --- | --- | --- | --- |

# 4.  Meta Analysis

### Build  an input file 
* a csv file is need to described each input, with an appropriate header

There are some pre-computed GEMMA output files we'll use in this example

```
ls data/summarystat/*.gemma
```

We create a file `input_meta.csv` with the file information

```
echo "rsID,Chro,Pos,A1,A2,Beta,Se,Pval,N,freqA1,direction,Imputed,Sep,File,Ncount" > utils/input_meta.csv
for File in `ls data/summarystat/{AFR,AMR,EAS,EUR,SAS}*.gemma`; do
   echo "rs,chr,ps,allele0,allele1,beta,se,p_wald,NA,af,NA,NA,TAB,$File,2500" >>  utils/input_meta.csv
done
```

### Run the meta-analysis pipeline

* input :  user can choose software that they want to run : metal (`--metal 1`), gwama (`--gwama 1`), metasoft (` --metasoft 1`) MrMega (`--mrmega 1`) and plink (`--plink 1`)

```

nextflow run h3abionet/h3agwas/meta/meta-assoc.nf   --metal 1 --gwama 1 --metasoft 1 --mrmega 1 --plink  1  --file_config utils/input_meta.csv -resume -profile singularity --output_dir meta

```

* output :
  * Each software in oits wn folder 
  * [Same report than association is generated](out_example/meta_report.pdf)

### Options for the meta-analysis software

| Software | `ma_genomic_cont` | `ma_inv_var_weigth` | `ma_overlap_sample` | `ma_random_effect` |
| --- | --- | --- | --- | --- |
| explanation | genomic control | invert variance weight | sample overlap | Random Effect |
| default | 0 | 0 | 0 | 0 |
| metal | yes | yes  | yes | no |
| gwama | yes | default  no | yes |
| Mr Mega| yes | no | no | no |
| plink | no | yes | no | no |
| Metasoft | no | no | no | default |
| --- | --- | --- | --- | --- |

1 'weighted-z' requests weighted Z-score-based p-values (as computed by the Abecasis Lab's METAL software)


# 5. Fine-mapping

Fine-mapping can be run on full summary statistics, or specific windows using two different scripta. There is also a  script which just does GCTA.

There are a number of general options.

First, we expect the summary statistics file(s) to have column headers. There are options that associated the column  headers with the the input values that the workflow expects


* header of summary statistics :
  * `--head_pval`  _p_-value of SNP
  * `--head_bp` position where SNP is
  * `--head_chr` Chromosome
  * `--head_rs` rs id of phenotype file
  * `--head_beta` effect on file 
  * `--head_se` standard error 
  * `--head_A1` effect of A1 allele
  * `--head_A2` effect of other allele
* Software :
  * by default all software will be run
* phenotype of gwas catalog:
  * `--list_phenogc`
  * `gwascat`
* significant threshold -- `--threshold_p` : by default to $5\times 10^{-8}$


## 5.1  Full summary statistics

The pipeline will apply clump to defined significant positions and run for each windows various software of fine mapping 

```
nextflow run  h3abionet/h3agwas/finemapping/main.nf --head_pval p_wald --head_bp ps --head_chr chr --head_rs rs --head_beta beta --head_se se --head_A1 allele1 --head_A2 allele0 --list_phenogc "Type 2 diabetes" --input_dir  data/imputed/  --input_pat imput_data --file_gwas data/summarystat/all_pheno.gemma  --output_dir finemapping_pheno1 --output finemapping_pheno1 -resume  -profile singularity
```

The output folder contains for each independant SNPs a folder with result :
  * figure with all resume data as locus zoom
  * file contains all positiosn with different probability of each software


## 5.2 Specific windows
Algorithm is same than previous, but user must specify the  chromosome (`--chro`), and range position (`--begin_seq`, `--end_seq`)

```
nextflow run h3abionet/h3agwas/finemapping/finemap_region.nf  --head_pval p_wald --head_bp ps --head_chr chr --head_rs rs --head_beta beta --head_se se --head_A1 allele1 --head_A2 allele0 --list_pheno "Type 2 diabetes" --input_dir  data/imputed/  --input_pat imput_data --file_gwas data/summarystat/all_pheno.gemma  --output_dir finemapping_pheno1_wind --output finemapping_pheno1 -resume  -profile singularity --begin_seq 112178657 --end_seq 113178657 --chro 10
```

## 5.3  GCTA-COJO

The `cojo-assoc.nf` pipeline performs conditional and joint analysis using summary data pipeline of cojo using two type of input :

* plink file + phenotype with `--input_dir`, `--input_pat`, `--data` and `--pheno`
* summary statistics :  `--file_gwas` one or more, and `--head_[lhead]`
* command line :
    * run cojo on each chromosome 
    * for each chromosome take best snp (`--cojo_top_snps_chro 1`)

The command to execute 

```
nextflow run h3abionet/h3agwas/finemapping/cojo-assoc.nf --head_pval p_wald --head_bp ps --head_chr chr --head_rs rs --head_beta beta --head_se se --head_A1 allele1 --head_A2 allele0 --input_dir data/imputed/  --input_pat imput_data  --output_dir cojo --data data/pheno/pheno_test.all --pheno pheno_qt1 --file_gwas data/summarystat/all_pheno.gemma  -resume   -profile singularity --cojo_top_snps_chro 1
```



# 6. Annotating data 


**Warning**
 aws : This workflow does not work on AWS at the moment


The following inputs are required
 * plink file
 * summary statistisc
 * phenotype and data file
 * rs from gwas (one or more, separate with comma)
 * Annotation used annovar file, if you don't `list_file_annot` and `info_file_annot`, data will be downloaded

```
nextflow run  h3abionet/h3agwas/utils/annotation/main.nf --head_pval p_wald --head_bp ps --head_chr chr --head_rs rs --head_beta beta --head_se se --head_A1 allele1 --head_A2 allele0 --input_dir data/imputed/  --input_pat imput_data --file_gwas data/summarystat/all_pheno.gemma  --output_dir annotation --list_rs "2:45832137:A:G,1:117539108:G:T" --data data/pheno/pheno_test.all --pheno pheno_qt1  -resume  -profile singularity --loczm_bin  "<full path to your locusZoom>"
```

# 7. Simulation of phenotype and build phenotype

The algorithm will use  the GWAS Catalog (default) and lead position to define position to build phenotype and gcta. The input required is a list of position for the genotype positions.

## 7.1 Basic operation
Key options
* `--list_vcf`:  If you provie a value for the `--list_vcf` option, the genotype data will be extracted from the given VCF file. If this option is not provided, then the data will be generated from 1000 Genomes Data
* `nb_snp` : snp number used to build phenotype


For our example, we construct a list of positions. For real work your data should be in a smilar format

```
awk '{print $1"\t"$4}' data/array_plk/array.bim  > utils/list_posarray
```


We can then execute the workflow:

```
nextflow run h3abionet/h3agwas/utils/build_example_data/main.nf -profile singularity   --pos_allgeno utils/list_posarray -resume --nb_snp 3 --output_dir simul_gcta_main
```

### Output of pipeline
* `geno_all` : contains final  genotype
  * `geno_all/$output[.bed/fam/bim]` : plink file contains positions from array
  * `geno_all/$output_sex` : sex of individual after random
  * `$output_sexsex_change_plk.ind` : list of individuals, where sex has been randomly changed
* `gwascat` : folder contains gwas catalog download :
  * `$output/gwascat/gwascat_format_all.csv` : all positions from gwas cat
  * `$output/gwascat/gwascat_format.[pos,bed]`  : positions associated to phenotype diabete
  * `$output/gwascat/gwascat_format_resume.csv` : information relative to diabetes
* `simul_pheno` :
  * `$output/simul_pheno/datai/$output[.bed,plink,fam]` : plink contains all position extracted from gwas catalo link to phenotype
  * `$output/simul_pheno/data_format` :
  * positions clumped to avoid LD between positions extracted of phenotype :
    *  `$output.effect.rs_clump.clumped` : positions clumped to avoid LD between positions extracted of phenotype
    *  `$output.effect.rs.assoc`  $output.effect.rs_clump.log      $output.effect.rs.infors
  * `$output.effect.rs` : effect of each position to build phenotype
  * `$output.effect.rs.infors` : information relative of positions to build phenotype

## 7.2  other example 

* see [page build example](README_buildataset.md)


## 7.3  Simulation using gcta and data


The algorithm is similar to the above but need genetic data in plink format

```
nextflow run h3abionet/h3agwas/utils/build_example_data/simul-assoc_gcta.nf -profile singularity  --input_dir data/imputed/  --input_pat  imput_data --output_dir simul_gcta

```

## 7.4 Simulation using random position  and phenosim

_phenosim_ will select random position in genome, defined effect of positions and run gemma or boltmm
* need plink file using `--input_dir` and `--input_pat`
* `--gemma 1` will run gemma after simulation of phenotype on genotype and phenotype


```
nextflow run h3abionet/h3agwas/utils/build_example_data/simul-assoc_phenosim.nf -profile singularity  --ph_normalise 0 --input_dir data/imputed/ --input_pat  imput_data --gemma 1
```

output :
* `output/gemma/res_gemma.stat` :  contains FP and TP comparison of position 
* `output/simul/` : contains phenotype and positions used for simulation

# 8. Conditional analysis using gemma


To perform a conditional analysis using gemma, we use come SNPs (positions) as covariate of the phenotype and the workflow seems is the specified positions (`--pos_ref`) are linked or indepependant.

The arguments are the same as gwas using gemma but you need to add :
 *  pipeline will run a raw gwas using phenotype and covariable and after performed gwas using genotypes of each position from `pos_cond` as covariable
 * `chro_cond` : chro where `pos_cond` and `pos_ref`
 * `pos_ref` : pos will be tested for independance or not to `pos_cond`
 * `pos_cond` : list of position as conditional to verify independance with `pos_ref`, include as covariate
The pipeline will compute LD 


```
nextflow run  h3abionet/h3agwas/finemapping/cond-assoc.nf --input_dir data/imputed/  --input_pat imput_data --output_dir cond --data data/pheno/pheno_test.all --pheno pheno_qt1  -profile singularity  -resume  --chro_cond 17  --pos_ref 78562921 --pos_cond 78494969,78502076,78534842  --sample_snps_rel=1
```

output :
* `cond/res/fig/` : contains figure and merge resultat :
  * out.csv   :

```
"rscond","chr","rs","ps","n_miss","allele1","allele0","af","beta","se","logl_H1","l_remle","p_wald","Num","CHR_A","BP_A","R2"
"Init",17,"17:78562921:G:A",78562921,0,"A","G",0.17,8.987846,1.361653,-2108.394,15.37693,1.050086e-10,1,NA,NA,1
"17:78494969:C:T",17,"17:78562921:G:A",78562921,0,"A","G",0.17,8.512503,2.45803,-2104.635,15.16254,0.0005800203,2,17,78494969,0.700696
"17:78534842:C:T",17,"17:78562921:G:A",78562921,0,"A","G",0.17,9.080245,1.380124,-2104.504,15.93441,1.200983e-10,3,17,78534842,0.00910659
"17:78502076:T:G",17,"17:78562921:G:A",78562921,0,"A","G",0.17,8.972201,1.362844,-2104.497,15.95574,1.170588e-10,4,17,78502076,0.0037796
"Merge",17,"17:78562921:G:A",78562921,0,"A","G",0.17,8.696802,2.570401,-2096.868,16.26486,0.0007726084,5,NA,NA,NA
```
    * rscond : rs used  for condition  
     * Init : gemma without pos cond
     *  Merge : all position together used as cond 
  * R2 : ld with posistion ref

* figures :
 * `out_ld.pdf` : ld between positions 
 * `out.pdf` : barplot of -log10(p\_wald), condition or not 


# 9. Format summary statistics, prepared data for imputation, or format data after imputation

## 9.1 format summary statistics

The `format_gwasfile.nf' script formats summary statistics, replaces header information

*Input* :
 *  new and old header, will be replaced
 * `input_dir` and `input_pat` : plink file contains used to rename rsid using informations
 * `--file_ref_gzip` : used to check ref alternatif or rsid
   * we used`utils/all_rsinfo.init.gz` : contains information relative to rsid / positions, subsample of files [here](shorturl.at/HJLN0)
 

**change header in file**

```
nextflow run  h3abionet/h3agwas/formatdata/format_gwasfile.nf --head_pval p_wald --head_bp ps --head_chr chr --head_rs rs --head_beta beta --head_se se --head_A1 allele1 --head_A2 allele0 --file_gwas data/summarystat/all_pheno.gemma  --output_dir format_assoc   -resume --headnew_pval p --headnew_bp bp --headnew_chr CHR --headnew_rs SNP --headnew_beta beta --headnew_se se --headnew_A1 allele1 --headnew_A2 allele0 --file_ref_gzip data/utils/all_rsinfo.init.gz --input_dir data/imputed/ --input_pat imput_data -profile singularity
```


**Add rsid to summary statistics**
 * if there is not rsid header in input, `format_gwasfile.nf` will add rsid extracted from file `--file_ref_gzip` based on chromosome, positions and alt/ref

```
nextflow run  h3abionet/h3agwas/formatdata/format_gwasfile.nf --head_pval p_wald  --head_chr chr --head_bp ps --head_beta beta --head_se se --head_A1 allele1 --head_A2 allele0 --file_gwas data/summarystat/all_pheno.gemma  --output_dir  ./ -resume  --file_ref_gzip data/utils/all_rsinfo.init.gz -profile singularity --output all_pheno.gemma.addrs -r dev
```

**Add chromosome and position in function of rsid**
 * if there is not chromosome / bp header in input, `format_gwasfile.nf` will add chromosome/positions extracted from file `--file_ref_gzip` based on rsid

```
nextflow run  h3abionet/h3agwas/formatdata/format_gwasfile.nf --head_pval p_wald   --head_beta beta --head_se se --head_A1 allele1 --head_A2 allele0 --file_gwas data/summarystat/all_pheno.gemma.addrs --output_dir  ./ -resume  --file_ref_gzip data/utils/all_rsinfo.init.gz --input_dir data/imputed/ --input_pat imput_data -profile slurmSingularity --head_out out_chrbp.gemma --head_rs rs --output all_pheno.gemma.addchrps -r dev
```

## 9.2 prepare data for imputation 

`plk_in_vcf_imp.nf` script take in input a plink file and prepared data for imputation

*Input :
 * `utils/all_rsinfo.init.gz` : contains information relative to rsid / positions, subsample of files [here](shorturl.at/HJLN0)

 * `input_dir` and `input_pat` plink file to convert in vcf 
 * `output_dir` and `output` : output dir and output header for the vcf ifle  
 * `file_ref_gzip` : contains information to check alternatif and reference
 * reffasta : fasta file of species


```
mkdir -p utils_data/
cd utils_data/
wget -c http://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
cd ../
nextflow run h3abionet/h3agwas/formatdata/plk_in_vcf_imp.nf -profile slurm  --input_dir='data/array_plk/' --input_pat='array_qc' --output_dir='vcf_output' --output='data_qc' --file_ref_gzip="data/utils/all_rsinfo.init.gz" --reffasta="utils_data/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
```


## 9.3 Format vcf file after imputation to  _plink_ format

* format file vcf come from to imputaion and produce a report with distribution of score and frequency
* options of interest: 
 * `--file_listvcf` : file contains list of VCF files
 * `--output_pat` : output pattern of plink file 
 * `--output_dir` : output direction where you file 
 * `--reffasta` : reference in fasta format
 * `--score_imp` : score imputation, depend of your imputation software [default: INFO]
 * `--min_scoreinfo` : minimum score [default : 0.6]
 
### 9.2.1 Example from Sanger imputation panel

```
 ls data/imputed/vcf/*.vcf.gz  > utils/listvcf
 mkdir -p utils_data/
 cd utils_data/
 wget -c http://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
 cd ../
 nextflow run h3abionet/h3agwas/formatdata/vcf_in_plink.nf --file_listvcf utils/listvcf --output_pat  kgp_imputed --output_dir plink_imputed/   --reffasta utils_data/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz  -profile singularity
```


### 9.2.2 Example from Michigan Imputation with genetic map


* file are zip with a password
* score imputation _R2_
* if you want to add genetics map in plink option you can download matp [here](https://zenodo.org/record/6422542/files/genetic_map_hg19.txt?download=1)

```
nextflow run h3abionet/h3agwas/formatdata/vcf_in_plink.nf --file_listvcf $ListVCF --min_scoreinfo 0.4 -profile slurm --output_pat 1000G_imp_mich --output_dir 1000G_imp_mich --genetic_maps genetic_map_hg19.txt --plink_mem_req 10GB -resume --big_time 1000h  --score_imp R2 --other_mem_req 20GB --reffasta hs37d5.fa.gz --unzip_zip 1 --unzip_password "xxx" --statfreq_vcf "%AF"
```

Output :
  * plink file 
  * report contains distribution of frequency and score info [example of report](out_example/KGPH3abionet_imputed_report.pdf)

Look at this  [other example](README_buildataset.md) 


## 9.3 Format VCF file in bimbam format

bimbam files contain dosage data used, for instance in gemma. Our pipeline does not use bimbam in this case, but some users may find this useful. The pipeline  uses _qctools (v2)_ and bcftools for format


  * `--file_listvcf` : file contains list of file vcf
  * `--output_pat` : output pattern of plink file
  * `--output_dir` : output direction where you file
  * `--score_imp` : score imputation, depend of your imputation software [default: INFO]
  * `--min_scoreinfo` : minimum score [default : 0.6]




```
nextflow run h3abionet/h3agwas/formatdata/vcf_in_bimbam.nf --file_listvcf utils/listvcf  --output_pat  kgp_imputed --output_dir bimbam/   -profile singularity
```

## 9.4 Format  VCF  file in  BGEN format

The BGEN Format is a dosage format used by _fastgwa_, _bolt-lmm_ or _regenie_ 

*Input option*:
  * `--file_listvcf` : file contains list of file vcf
  * `--output_pat` : output pattern of plink file
  * `--output_dir` : output direction where you file
  * `--score_imp` : score imputation, depend of your imputation software [default: INFO]
  * `--min_scoreinfo` : minimum score [default : 0.6]
* bcftools are used to filter vcf and qctools (v2) for format in bgen


### 9.4.1 Each vcf files independently

Each VCF files is  filtered independently; they are then merged and formatted in BGEN formay

```
nextflow run h3abionet/h3agwas/formatdata/vcf_in_bgen.nf --file_listvcf utils/listvcf --output_pat  exampledata2_imp --output_dir ./bgen -resume -profile singularity
```


### 9.4.2 Merge VCF files and reformat

Each vcf is filtered independantly; they are merged and the  final file is refromatten

```
nextflow run h3abionet/h3agwas/formatdata/vcf_in_bgen_merge.nf --output_dir bgen_v3 --output all  -profile singularity --file_listvcf listvcf -resume
```

### 9.4.3 Format first, and reformat

Each vcf is filtered reformatted in BGEN; all bgen are then merged


```
nextflow run h3abionet/h3agwas/formatdata/vcf_in_bgen_merge.nf --output_dir bgen_v3 --output all  -profile singularity --file_listvcf listvcf -resume
```



### 9.5 Format vcf in IMPUTE2 format

Dosage format for `bolt-lmm`

** Input option** :
  * `--file_listvcf` : file contains list of file vcf
  * `--output_pat` : output pattern of plink file
  * `--output_dir` : output direction where you file
  * `--score_imp` : score imputation, depend of your imputation software [default: INFO]
  * `--min_scoreinfo` : minimum score [default : 0.6]



```
nextflow h3abionet/h3agwas/formatdata/vcf_in_impute2.nf --file_listvcf listvcf --output_pat  utils/exampledata2_imp --output_dir ./impute2 -profile singularity -resume
```


# 10 Heritability and co-heritability estimation

The heritability pipeline can use summary statistics or genetics data for heritability  and co-heritatbility estimation: 

* With genetics data and phenotye :
    * `--input_dir`, `--input_pat`, `--data`, `--pheno`  `covariates`
* summary statistics :  `--file_gwas` one or more separated by a comma, and `--head_[lhead]`
* software options (ldsc or gemma):
  * `--ldsc_h2 1` : performed ldsc 
  * `--ldsc_h2_multi 1` : to do co heritability (must be have more than one file_gwas)
  * `--gemma_h2_pval 1` : performed gemma using summary statistics
  * `--Nind` : number of individuals (if not in summary statitics, if more than 1 studies, separated by a comma)
* genetics and phenotype :
  * `--gemma_h2 1` : using gemma
  * `--gcta_h2 1` : using gcta
  * `--gcta_h2_multi` : gcta and computed co - heritability
  * `--bolt_h2` : _bolt-lmm_
  * `--bolt_h2_multi` : _bolt-lmm_ coheritability between phenotype

## Example witthout ldsc

```
nextflow run h3abionet/h3agwas/heritabilities/main.nf \
  --input_dir data/imputed/  --input_pat imput_data --data data/pheno/pheno_test.all --pheno pheno_qt1,pheno_qt2 \
  --file_gwas data/summarystat/all_pheno.gemma,data/summarystat/all_phenoq2.gemma   --head_pval  "p_wald"  --head_freq  "af" --head_bp  "bp" --head_chr  "chr" --head_rs  "rs" --head_beta "beta" --head_se "se" --head_A1 "allele1" --head_A2 "allele0" --Nind 500,500 \
  --ldsc_h2 0 --ldsc_h2_multi 0 --bolt_h2 1 --bolt_h2_multi 1 --gcta_h2 0 --gcta_h2_imp 0 --gcta_h2_multi 0 --gemma_h2 1 --gemma_h2_pval 1 -resume --output_dir heritability/ -profile singularity --grm_cutoff 0.5
```


## Example wit ldsch

```
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
tar xvjf eur_w_ld_chr.tar.bz2
~/nextflow $DirNex/h3agwas/heritabilities/main.nf \
  --input_dir data/imputed/  --input_pat imput_data --data data/pheno/pheno_test.all --pheno pheno_qt1,pheno_qt2 \
  --file_gwas data/summarystat/all_pheno.gemma,data/summarystat/all_phenoq2.gemma   --head_pval  "p_wald"  --head_freq  "af" --head_bp  "ps" --head_chr  "chr" --head_rs  "rs" --head_beta "beta" --head_se "se" --head_A1 "allele1" --head_A2 "allele0" --Nind 500,500 \
  --ldsc_h2 1 --ldsc_h2_multi 1 --bolt_h2 1 --bolt_h2_multi 1 --gcta_h2 1 --gcta_h2_imp 0 --gcta_h2_multi 0 --gemma_h2 1 --gemma_h2_pval 1 -resume --output_dir heritability/ $profile --grm_cutoff 0.5 --dir_ref_ld_chr eur_w_ld_chr 

```

###Output :
* each log output by folder
* file contains all information extracted and merge 
* figure to compared each heritability

# 11. Multi-Trait Analysis of GWAS

Multi-trait analysis of GWAS (MTAG), a method for joint analysis of summary statistics from genome-wide association studies (GWAS) of different traits, possibly from overlapping samples. 

**Input** : 
* multi-trait analysis of GWAS (MTAG), a method for joint analysis of summary statistics from genome-wide association studies (GWAS) of different traits, possibly from overlapping samples. 
* input : 
 * list of summary statistic `file_gwas` and header from gwas file: `-head_[name]`
 * also you can give nformation relative to : ` --input_dir data/imputed/ --input_pat imput_data --pheno pheno_qt1,pheno_qt2 --data data/pheno/pheno_test.all `, can add N value to each summary statistic

**Output**:
 * each log contains heritability as output order by folder
 * 1 file contains all heritability extracted and merge 
 * figure to compared each heritability


```
nextflow run h3abionet/h3agwas/meta/mtag-assoc.nf --head_freq af --head_pval p_wald --head_bp ps --head_chr chr --head_rs rs --head_beta beta --head_se se --head_A1 allele1 --head_A2 allele0 --input_dir data/imputed/ --input_pat imput_data --file_gwas  data/summarystat/all_pheno.gemma,data/summarystat/all_phenoq2.gemma --pheno pheno_qt1,pheno_qt2 --data data/pheno/pheno_test.all -resume   -profile singularity
```

# 12 Conversion of positions  between build


## 12.1 GWAS catalog
Objective is to used two way to convert positions build using rs value and crossmap using chromosome positions

The objective is to convert build positions  using rs value and _crossmap_ using chromosome positions

Example : download gwas catalog in hg38, and convert in 19 / 37 position 
```
nextflow run h3abionet/h3agwas/formatdata/convert_posversiongenome.nf -profile singularity -resume
```

## 12.2 convert Summary statistics between hg38 and 19

Example : we used output of `h3abionet/h3agwas/formatdata/format_gwasfile.nf` to convert hg19 in 38 

```
nextflow run h3abionet/h3agwas/formatdata/convert_posversiongenome.nf -profile singularity -resume --file_toconvert  data/summarystat/assoc_testwithrs.sumstat  link_rs_info=ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz --link_data_crossmap=http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz --head_chr CHR --head_bp bp --head_rs SNP --sep SPACE
```


# 13 requirement 


