# General instalation 


Installation recommended has been tested using docker image 

Binary used to test pipeline can be find :
 * binary/

## latex 

```
sudo apt-get install texlive-latex-base texlive-fonts-recommended texlive-fonts-extra texlive-latex-extra
```

## plink

plink version in pipeline tested was version of plink 1.9  (1.90~b6.16-200217-1) from repository ubuntu, plink1.9 binary must be renamed plink (cp or ln)

```
sudo apt-get update && sudo apt-get install plink1.9
sudo cp /bin/plink1.9 /bin/plink && chmod +x /bin/plink
```

## python

we tested with python 3.8 (should be work with python 3.6), it they is not in your repository

```
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:deadsnakes/ppa
```

```
sudo apt-get install python3.8 python3-pip
### upgrade pip
pip3 install --upgrade pip
```

## nextlow 

### install java, nextflow required java 11 or 18

we used tutorial [her](https://computingforgeeks.com/install-oracle-java-18-on-ubuntu-debian/)

```
## ensure that you havge get and curl
sudo apt -y install wget curl
## download debian lib
wget https://download.oracle.com/java/18/latest/jdk-18_linux-x64_bin.deb
## install
sudo apt install ./jdk-18_linux-x64_bin.deb -y
## Configure Java environment.
cat <<EOF | sudo tee /etc/profile.d/jdk18.sh
export JAVA_HOME=/usr/lib/jvm/jdk-18
export PATH=\$PATH:\$JAVA_HOME/bin
```

Install in your current directory and pull h3agwas/ha3bionet

```
wget -qO- https://get.nextflow.io | bash
nextflow pull h3abionet/h3agwas
```



## SAIGE
 * [installation](https://saigegit.github.io/SAIGE-doc/)
 * Version doesn't work on ubuntu 22.04
 * package : saige (release)
  * installation R, and devtools 3.6.3
  * cmake
  * libopenblas-base
  * python3

```
sudo apt-get update && \
    apt-get install -y \
    r-base \
    r-cran-devtools \
    build-essential \
    cmake \
    libopenblas-base \
    python3-pip \
    r-cran-devtools \
    git
pip3 install cget
```

```
git clone https://github.com/saigegit/SAIGE.git
cd SAIGE
Rscript extdata/install_packages.R
R CMD INSTALL .

mv step1_fitNULLGLMM.R step2_SPAtests.R step3_LDmat.R createSparseGRM.R /usr/local/bin/

chmod a+x /usr/local/bin/step1_fitNULLGLMM.R
chmod a+x /usr/local/bin/step2_SPAtests.R
chmod a+x /usr/local/bin/step3_LDmat.R
chmod a+x /usr/local/bin/createSparseGRM.R

```

## softwares

### gemma 

various version of GEMMA exist :
 * github you can do a compulation
 * conda
 * compilation, pipeline has been tested with :
   * [gemma : 98.4](https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.4/gemma-0.98.4-linux-static-AMD64.gz) with some bug
   * [gemma : 98.6, prerelease, fixed some error ](https://biogems.info/download/gemma-0.98.6-pre1.gz)

```
wget https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.4/gemma-0.98.4-linux-static-AMD64.gz
gunzip gemma-0.98.4-linux-static-AMD64.gz
## move gemma and rename in binary
mv gemma-0.98.4-linux-static-AMD64 /usr/local/bin/gemma 
chmod +x /usr/local/bin/gemma

```

### gcta

[Pipeline has been tested with version 1.93.beta](https://raw.githubusercontent.com/h3abionet/h3agwas/master/utils/bin/gcta_1.93.2beta.zip)

```
wget https://raw.githubusercontent.com/h3abionet/h3agwas/master/utils/bin/gcta_1.93.2beta.zip && unzip gcta_1.93.2beta.zip && cp gcta_1.93.2beta/gcta64 . && rm -rf gcta_1.93.2beta && rm -rf gcta_1.93.2beta.zip  && mv gcta64 /usr/local/bin/
```

### Bolt-LMM
* pipeline has been tested used using binary version 2.4
 * on ubuntu bolt need libiomp5.so, install libomp-dev
```
apt-get install -y libomp-dev
wget https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/downloads/BOLT-LMM_v2.4.tar.gz && tar -xzf  BOLT-LMM_v2.4.tar.gz && cp BOLT-LMM_v2.4/bolt bin/ && mv BOLT-LMM_v2.4.tar.gz src/ && rm -rf BOLT-LMM_v2.4/
```

### Regenie

 * version tested 3.1.3


```
wget https://github.com/rgcgithub/regenie/releases/download/v3.1.3/regenie_v3.1.3.gz_x86_64_Linux.zip && unzip regenie_v3.1.3.gz_x86_64_Linux.zip && mv regenie_v3.1.3.gz_x86_64_Linux /usr/local/bin/regenie  && chmod +x /usr/local/bin/regenie
```
## software utils
### bgenix

```
wget http://code.enkre.net/bgen/tarball/release/bgen.tgz && tar -xzf bgen.tgz && cd bgen.tgz/ && ./waf configure && ./waf && cp build/apps/bgenix /usr/bin/
```

### tabix 

```
sudo apt-get install tabix -y
```

## run data-set 

### download

```
### install git 
sudo apt-get install git 
## clone repository 
git clone https://github.com/h3abionet/h3agwas-examples
##build a folder, link data
mkdir testgit_qc && cd testgit_qc/ && ln -s ../h3agwas-examples/run_test.bash . && ln -s ../h3agwas-examples/data .
```

### run a test

a small test has been developped to run on test data set and differen part of pipeline. 
Script taked 5 arguments :
 * testdone : script to test qc, assoc
 * h3agwasdir : where find h3agwas repositotu 
 * profile slurm ? SIngularity?
 * nextflowbin : where find nextflow
 * other option  (optional)

```
./run_test.bash qc h3abionet batch standard nextflow
```


## Quality Control

test on ubuntu 20.04

need :
 * latex 
 * python
 * plink

## GWAS 

test on ubuntu 20.04

need :
 * latex 
 * python 
 * plink 1.9
optional :
 * gemma 
 * Regenie 
 * saige :
  * R
 * gcta
 * bolt-lmm
 * saige :
  * tabix for vcf input
  * bgenix for vcf input



Doesn't need other requierment that previous see [docker image](Docker/qc/)


