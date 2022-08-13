testdone=$1
h3agwasdir=$2
profile=$3
nextflowbin=$4
otheroption=$5
if [ $testdone == "installnf" ]
then 
wget -qO- https://get.nextflow.io | bash
chmod +x nextflow
fi
if [ $testdone == "qc" ]
then
$nextflowbin run $h3agwasdir//h3agwas/qc/main.nf --input_dir data/array_plk  --input_pat array --output_dir qc  --output kgpexample \
 --phenotype data/pheno/pheno_test.all --pheno_col phenoqc_ql \
 --case_control data/pheno/pheno_test.all --case_control_col Sex \
 --batch data/pheno/pheno_test.all --batch_col batch \
 -profile $profile -resume $otheroption
fi

if [ "$testdone" == "assoc" ]
then
awk '{print $1"\t"$4"\t"$4"\t"$2}'  data/array_plk/array.bim > list_pos
ls bgen_v3/bgen_chro/*.bgen > listbgen
$nextflowbin run $h3agwasdir/h3agwas/assoc --input_dir data/imputed/ --input_pat imput_data \
 --data data/pheno/pheno_test.all --pheno pheno_qt1,pheno_qt2 \
 --output_dir assoc_listbgen --output assoc \
 --boltlmm 1 --sample_snps_rel 1 --regenie 1 --fastgwa 1 --grm_nbpart 2\
  -profile $profile \
  --gemma_num_cores  6\
 --bgen data/imputed/bgen/out.bgen --bgen_sample data/imputed/bgen/out.sample --saige 1 -resume 
fi

if [ "$testdone" == "assocnobgen" ]
then
awk '{print $1"\t"$4"\t"$4"\t"$2}'  data/array_plk/array.bim > list_pos
$nextflowbin run $h3agwasdir/h3agwas/assoc --input_dir data/imputed/ --input_pat imput_data \
 --data data/pheno/pheno_test.all --pheno pheno_qt1,pheno_qt2 \
 --output_dir assoc_listbgen --output assoc \
 --boltlmm 1 --sample_snps_rel 1 --regenie 1 --fastgwa 1 --grm_nbpart 2\
  -profile $profile \
  --gemma_num_cores  6\
 --saige 1 -resume 
fi
if [ "$testdone" == "finemapping" ]
then
$nextflowbin run  $h3agwasdir/h3agwas/finemapping/main.nf --head_pval p_wald --head_bp ps --head_chr chr --head_rs rs --head_beta beta --head_se se --head_A1 allele1 --head_A2 allele0 --list_phenogc "Type 2 diabetes" --input_dir  data/imputed/  --input_pat imput_data --file_gwas data/summarystat/all_pheno.gemma  --output_dir finemapping_pheno1 -resume  -profile $profile --size_wind_kb 500 --data  data/pheno/pheno_test.all --pheno pheno_qt1
fi

if [ "$testdone" == "finemap_region" ]
then
$nextflowbin run  $h3agwasdir/h3agwas/finemapping/finemap_region.nf  --head_pval p_wald --head_bp ps --head_chr chr --head_rs rs --head_beta beta --head_se se --head_A1 allele1 --head_A2 allele0 --list_phenogc "Type 2 diabetes"   --file_gwas data/summarystat/all_pheno.gemma  --output_dir finemapping_pheno1_wind --output finemapping_pheno1 -resume  -profile $profile --begin_seq 112178657 --end_seq 113178657 --chro 10 
fi

if [ "$testdone" == "cojo" ]
then
$nextflowbin run  $h3agwasdir/h3agwas/finemapping/cojo-assoc.nf --head_pval p_wald --head_bp ps --head_chr chr --head_rs rs --head_beta beta --head_se se --head_A1 allele1 --head_A2 allele0 --input_dir data/imputed/  --input_pat imput_data  --output_dir cojo --data data/pheno/pheno_test.all --pheno pheno_qt1 --file_gwas data/summarystat/all_pheno.gemma  -resume   -profile $profile --cojo_top_snps_chro 5 
fi



if [ "$testdone" == "cond" ]
then
$nextflowbin run  $h3agwasdir/h3agwas/finemapping/cond-assoc.nf --input_dir data/imputed/  --input_pat imput_data --output_dir cond --data data/pheno/pheno_test.all --pheno pheno_qt1  -profile $profile  -resume  --chro_cond 17  --pos_ref 78562921 --pos_cond 78494969,78502076,78534842  --sample_snps_rel=1
fi

if [ "$testdone" == "meta" ]
then
$nextflowbin run $h3agwasdir/h3agwas/meta/meta-assoc.nf   --metal 1 --gwama 1 --metasoft 1 --mrmega 1 --plink  1  --file_config utils/input_meta.csv -resume -profile singularity --output_dir meta
fi
if [ "$testdone" == "formatgwasfile" ]
then
$nextflowbin run  $h3agwasdir/h3agwas/formatdata/format_gwasfile.nf --head_pval p_wald --head_bp ps --head_chr chr --head_rs rs --head_beta beta --head_se se --head_A1 allele1 --head_A2 allele0 --file_gwas data/summarystat/all_pheno.gemma  --output_dir format_assoc   -resume --headnew_pval p --headnew_bp bp --headnew_chr CHR --headnew_rs SNP --headnew_beta beta --headnew_se se --headnew_A1 allele1 --headnew_A2 allele0 --file_ref_gzip data/utils/all_rsinfo.init.gz --input_dir data/imputed/ --input_pat imput_data -profile $profile
fi

#ls data/imputed/vcf/*.gz > listvcf
#$nextflowbin run /home/jeantristan/Travail/git/h3agwas/formatdata/vcf_in_bgen_merge_chro.nf --output_dir bgen_v3 --output all  -profile $profile --file_listvcf listvcf -resume

if [ "$testdone" == "convertpos2" ]
then
$nextflowbin run $h3agwasdir/h3agwas/formatdata/convert_posversiongenome.nf -profile $profile -resume --file_toconvert  data/summarystat/assoc_testwithrs.sumstat  link_rs_info=ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz --link_data_crossmap=http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz --head_chr CHR --head_bp bp --head_rs SNP --sep SPACE
fi

if [ "$testdone" == "heritability" ]
then
$nextflowbin $h3agwasdir/h3agwas/heritabilities/main.nf \
  --input_dir data/imputed/  --input_pat imput_data --data data/pheno/pheno_test.all --pheno pheno_qt1,pheno_qt2 \
  --file_gwas data/summarystat/all_pheno.gemma,data/summarystat/all_phenoq2.gemma   --head_pval  "p_wald"  --head_freq  "af" --head_bp  "ps" --head_chr  "chr" --head_rs  "rs" --head_beta "beta" --head_se "se" --head_A1 "allele1" --head_A2 "allele0" --Nind 500,500 \
  --ldsc_h2 0 --ldsc_h2_multi 0 --bolt_h2 1 --bolt_h2_multi 1 --gcta_h2 1 --gcta_h2_imp 0 --gcta_h2_multi 0 --gemma_h2 1 --gemma_h2_pval 1 -resume --output_dir heritability/ -profile $profile --grm_cutoff 0.5
fi

if [ "$testdone" == "heritabilityldsc" ]
then
#wget -c https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_baselineLD_ldscores.tgz
#tar -xzf 1000G_Phase3_baselineLD_ldscores.tgz
#wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
#tar xvjf eur_w_ld_chr.tar.bz2

$nextflowbin $h3agwasdir/h3agwas/heritabilities/main.nf \
  --input_dir data/imputed/  --input_pat imput_data --data data/pheno/pheno_test.all --pheno pheno_qt1,pheno_qt2 \
  --file_gwas data/summarystat/all_pheno.gemma,data/summarystat/all_phenoq2.gemma   --head_pval  "p_wald"  --head_freq  "af" --head_bp  "ps" --head_chr  "chr" --head_rs  "rs" --head_beta "beta" --head_se "se" --head_A1 "allele1" --head_A2 "allele0" --Nind 500,500 \
  --ldsc_h2 1 --ldsc_h2_multi 1 --bolt_h2 1 --bolt_h2_multi 1 --gcta_h2 1 --gcta_h2_imp 0 --gcta_h2_multi 0 --gemma_h2 1 --gemma_h2_pval 1 -resume --output_dir heritability/ -profile $profile --grm_cutoff 0.5 --dir_ref_ld_chr eur_w_ld_chr -with-trace
fi
