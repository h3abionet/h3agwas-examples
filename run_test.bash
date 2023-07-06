testdone=$1
h3agwasdir=$2
profile=$3
nextflowbin=$4
otheroption=$5
maxmem=10GB
maxcpu=4
if [ $testdone == "installnf" ]
then 
wget -qO- https://get.nextflow.io | bash
chmod +x nextflow
fi
if [ $testdone == "qc" ]
then
$nextflowbin run $h3agwasdir/h3agwas/qc/main.nf --input_dir data/array_plk  --input_pat array --output_dir qc  --output array_qc \
 --phenotype data/pheno/pheno_test.all --pheno_col phenoqc_ql \
 --case_control data/pheno/pheno_test.all --case_control_col Sex \
 --batch data/pheno/pheno_test.all --batch_col batch \
 -profile $profile -resume $otheroption
fi

if [ "$testdone" == "assoc" ]
then
awk '{print $1"\t"$4"\t"$4"\t"$2}'  data/array_plk/array.bim > list_pos
ls data/imputed/bgen_chro/*.bgen > listbgen
echo "$nextflowbin run $h3agwasdir/h3agwas/assoc/main.nf --input_dir data/imputed/ --input_pat imput_data \
 --data data/pheno/pheno_test.all --pheno pheno_qt1,pheno_qt2 \
 --output_dir assoc_bgen --output assoc \
 --boltlmm 1 --sample_snps_rel 1 --regenie 1 --fastgwa 1 --grm_nbpart 2\
  -profile $profile \
 --gemma_num_cores $maxcpu  --plink_num_cores $maxcpu  --other_num_cores $maxcpu --bolt_num_cores $maxcpu --saige_num_cores $maxcpu --regenie_num_cores $maxcpu --fastgwa_num_cores $maxcpu \
 --gemma_mem_req $maxmem  --plink_mem_req $maxmem  --other_mem_req $maxmem --bolt_mem_req $maxmem --saige_mem_req $maxmem --regenie_mem_req $maxmem --fastgwa_mem_req $maxmem \
 --bgen data/imputed/bgen/out.bgen --bgen_sample data/imputed/bgen/out.sample --saige 1 -resume --assoc 1" > run_"$testdone".bash
bash run_"$testdone".bash

fi

if [ "$testdone" == "assocnobgen" ]
then
awk '{print $1"\t"$4"\t"$4"\t"$2}'  data/array_plk/array.bim > list_pos
$nextflowbin run $h3agwasdir/h3agwas/assoc/main.nf --input_dir data/imputed/ --input_pat imput_data \
 --data data/pheno/pheno_test.all --pheno pheno_qt1,pheno_qt2 \
 --output_dir assoc_nobgen --output assoc_nobgen \
 --boltlmm 1 --sample_snps_rel 1 --regenie 1 --fastgwa 1 --grm_nbpart 2 \
  -profile $profile \
 --gemma_num_cores $maxcpu  --plink_num_cores $maxcpu  --other_num_cores $maxcpu --bolt_num_cores $maxcpu --saige_num_cores $maxcpu --regenie_num_cores $maxcpu --fastgwa_num_cores $maxcpu \
 --gemma_mem_req $maxmem  --plink_mem_req $maxmem  --other_mem_req $maxmem --bolt_mem_req $maxmem --saige_mem_req $maxmem --regenie_mem_req $maxmem --fastgwa_mem_req $maxmem \
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

if [ "$testdone" == "metapval" ]
then
$nextflowbin run $h3agwasdir/h3agwas/meta/meta-assoc.nf   --metal 1 --gwama 1 --metasoft 1 --mrmega 1 --plink  1  --file_config utils/input_meta.csv -resume -profile singularity --output_dir metapval --used_pval_z 1
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




if [ "$testdone" == "prepareforimp" ]
then
if [ ! -f utils_data/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz ]
then
mkdir -p utils_data/
cd utils_data/
wget -c http://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
cd ../
fi
$nextflowbin $h3agwasdir/h3agwas/formatdata/plk_in_vcf_imp.nf -profile $profile --input_dir='data/array_plk/' --input_pat='array_qc' --output_dir='vcf_forimp' --output='data_qc' --file_ref_gzip="data/utils/all_rsinfo.init.gz" --reffasta="utils_data/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
fi

if [ "$testdone" == "prepareforimp_michigan" ]
then
if [ ! -f utils_data/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz ]
then
mkdir -p utils_data/
cd utils_data/
wget -c http://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
cd ../
fi
$nextflowbin $h3agwasdir/h3agwas/formatdata/plk_in_vcf_imp.nf -profile $profile --input_dir='data/array_plk/' --input_pat='array_qc' --output_dir='vcf_forimp_michigan' --output='data_qc' --file_ref_gzip="data/utils/all_rsinfo.init.gz" --reffasta="utils_data/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz" --michigan_qc 1 -resume
fi


if [ "$testdone" == "vcfinplink" ]
then
if [ -f utils_data/Homo_sapiens.GRCh37.dna.primary_assembly.fa_clean.fa.gz ]
then
 fasta=utils_data/Homo_sapiens.GRCh37.dna.primary_assembly.fa_clean.fa.gz
else
 if [ ! -f utils_data/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz ]
 then
 mkdir -p utils_data/
 cd utils_data/
 wget -c http://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
 cd ../
 fasta=utils_data/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
 fi
fi
nextflow run $h3agwasdir/h3agwas/formatdata/vcf_in_plink.nf --file_listvcf utils/listvcf --output_pat  kgp_imputed --output_dir plink_imputed/   --reffasta $fasta  -profile $profile -resume  --file_ref_gzip="data/utils/all_rsinfo.init.gz" --data data/pheno/pheno_test.all 
fi

if [ "$testdone" == "vcfinbimbam" ]
then
nextflow run $h3agwasdir/h3agwas/formatdata/vcf_in_bimbam.nf --file_listvcf utils/listvcf --output_pat kgp_imputed --output_dir bimbam_imputed/     -profile $profile -resume  --cut_hwe 0.0001
fi


if [ "$testdone" == "vcfinbimbamnohwe" ]
then
nextflow run $h3agwasdir/h3agwas/formatdata/vcf_in_bimbam.nf --file_listvcf utils/listvcf --output_pat kgp_imputed --output_dir bimbam_imputed_nohwe/     -profile $profile -resume
fi



if [ "$testdone" == "vcfinbgen" ]
then
nextflow run $h3agwasdir/h3agwas/formatdata/vcf_in_bgen.nf --file_listvcf utils/listvcf --output_pat kgp_imputed --output_dir bgen_imputed/     -profile $profile -resume  --cut_hwe 0.0001
fi


if [ "$testdone" == "vcfinbgennohwe" ]
then
nextflow run $h3agwasdir/h3agwas/formatdata/vcf_in_bgen.nf --file_listvcf utils/listvcf --output_pat kgp_imputed --output_dir bgen_imputed_nohwe/     -profile $profile -resume
fi


if [ "$testdone" == "vcfinbgenmerge" ]
then
nextflow run $h3agwasdir/h3agwas/formatdata/vcf_in_bgen_merge.nf --file_listvcf utils/listvcf --output_pat kgp_imputed --output_dir bgenmerge_imputed/     -profile $profile -resume  --cut_hwe 0.0001
fi

if [ "$testdone" == "vcfinbgenmergenohwe" ]
then
nextflow run $h3agwasdir/h3agwas/formatdata/vcf_in_bgen_merge.nf --file_listvcf utils/listvcf --output_pat kgp_imputed --output_dir bgenmerge_imputed_nohwe/     -profile $profile -resume  
fi


if [ "$testdone" == "vcfinbgenchromerge" ]
then
nextflow run $h3agwasdir/h3agwas/formatdata/vcf_in_bgen_merge_chro.nf --file_listvcf utils/listvcf --output_pat kgp_imputed --output_dir bgenmergechro_imputed/     -profile $profile -resume  --cut_hwe 0.0001
fi

if [ "$testdone" == "vcfinbgenmergechronohwe" ]
then
nextflow run $h3agwasdir/h3agwas/formatdata/vcf_in_bgen_merge_chro.nf --file_listvcf utils/listvcf --output_pat kgp_imputed --output_dir bgenmergechro_imputed_nohwe/     -profile $profile -resume  
fi


if [ "$testdone" == "vcfinimpute2chromerge" ]
then
nextflow run $h3agwasdir/h3agwas/formatdata/vcf_in_impute2.nf --file_listvcf utils/listvcf --output_pat kgp_imputed --output_dir impute2mergechro_imputed/     -profile $profile -resume  --cut_hwe 0.0001
fi

if [ "$testdone" == "vcfinimpute2mergechronohwe" ]
then
nextflow run $h3agwasdir/h3agwas/formatdata/vcf_in_impute2.nf --file_listvcf utils/listvcf --output_pat kgp_imputed --output_dir impute2mergechro_imputed_nohwe/     -profile $profile -resume
fi




if [ "$testdone" == "replication_gc" ]
then
$nextflowbin $h3agwasdir/h3agwas/replication/gwascat/main.nf  --head_pval p_wald --head_bp ps --head_chr chr --head_rs rs --head_beta beta --head_se se --head_A1 allele1 --head_A2 allele0  --file_gwas data/summarystat/all_pheno.gemma  --output_dir replication_gc -profile $profile -resume  --input_dir data/imputed/ --input_pat imput_data --head_af af --pheno_gc "Type 2 diabetes"
fi

if [ "$testdone" == "replication_ss" ]
then
#allowed_params_head_sumstat1 = ["file_gwas_sumstat1","head_pval_sumstat1", "head_freq_sumstat1", "head_bp_sumstat1", "head_chr_sumstat1", "head_rs_sumstat1", "head_beta_sumstat1", "head_se_sumstat1", "head_A1_sumstat1", "head_A2_sumstat1", "head_n_sumstat1", "n_count1",'head_z_sumstat1']
#allowed_params_head_sumstat2 = ["file_gwas_sumstat2","head_pval_sumstat2", "head_freq_sumstat2", "head_bp_sumstat2", "head_chr_sumstat2", "head_rs_sumstat2", "head_beta_sumstat2", "head_se_sumstat2", "head_A1_sumstat2", "head_A2_sumstat2", "head_n_sumstat2", 'n_count2','head_z_sumstat2']
$nextflowbin $h3agwasdir/h3agwas/replication/fullsumstat/main.nf  \
 --head_pval_sumstat1 p_wald --head_bp_sumstat1 ps --head_chr_sumstat1 chr --head_rs_sumstat1 rs --head_beta_sumstat1 beta --head_se_sumstat1 se --head_A1_sumstat1 allele1 --head_A2_sumstat1 allele0 --file_gwas_sumstat1 data/summarystat/all_pheno.gemma --head_frq_sumstat1 af --n_count1 500 \
 --head_pval_sumstat2 p_wald --head_bp_sumstat2 ps --head_chr_sumstat2 chr --head_rs_sumstat2 rs --head_beta_sumstat2 beta --head_se_sumstat2 se --head_A1_sumstat2 allele1 --head_A2_sumstat2 allele0 --file_gwas_sumstat2 data/summarystat/all_phenoq2.gemma --head_frq_sumstat1 af --n_count2 500\
  --output_dir replication_ss -profile $profile -resume  --input_dir data/imputed/ --input_pat imput_data 
fi

