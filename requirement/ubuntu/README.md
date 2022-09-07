# General instalation using local machine

Installation recommended has been tested using docker image 

## nextlow 

### install dependencie of nextflow

[source](https://computingforgeeks.com/install-oracle-java-18-on-ubuntu-debian/)

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

### install nextlow 

Install nextlow in your current directory 

```
wget -qO- https://get.nextflow.io | bash
```

### pull  h3agwas 

```
nextflow pull h3abionet/h3agwas
```

## dependencie and softwares used by pipeline

### latex 

latex need to build report

```
sudo apt-get install texlive-latex-base texlive-fonts-recommended texlive-fonts-extra texlive-latex-extra
```

### plink

plink version in pipeline tested was version of plink 1.9  (1.90~b6.16-200217-1) from repository ubuntu, plink1.9 binary must be renamed plink (cp or ln)

```
sudo apt-get update && sudo apt-get install plink1.9
sudo cp /bin/plink1.9 /bin/plink && chmod +x /bin/plink
```

### python

we tested with python 3.8 and python 3.6, it they is not in your repository

```
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:deadsnakes/ppa
```

```
sudo apt-get install python3 python3-pip
### upgrade pip
pip3 install --upgrade pip
```

### python library

```

```

### SAIGE
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
pip3 install --upgrade pip
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

### gcta / fastgwa

[Pipeline has been tested with version 1.93.beta](https://raw.githubusercontent.com/h3abionet/h3agwas/master/utils/bin/gcta_1.93.2beta.zip)

```
wget https://raw.githubusercontent.com/h3abionet/h3agwas/master/utils/bin/gcta_1.93.2beta.zip && unzip gcta_1.93.2beta.zip && cp gcta_1.93.2beta/gcta64 . && rm -rf gcta_1.93.2beta && rm -rf gcta_1.93.2beta.zip  && mv gcta64 /usr/local/bin/
```

### Bolt-LMM
* pipeline has been tested used using binary version 2.4
 * ubuntu 20.04, libiomp5.so, install libomp-dev

```
apt-get install -y libomp-dev
wget https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/downloads/BOLT-LMM_v2.4.tar.gz 
tar -xzf  BOLT-LMM_v2.4.tar.gz 
cp BOLT-LMM_v2.4/bolt /usr/local/bin 
chmod +x /usr/local/bin/bolt
```

### Regenie
 * version tested 3.1.3
 * regenie need bgenix when run using bgen

```
wget https://github.com/rgcgithub/regenie/releases/download/v3.1.3/regenie_v3.1.3.gz_x86_64_Linux.zip && unzip regenie_v3.1.3.gz_x86_64_Linux.zip && mv regenie_v3.1.3.gz_x86_64_Linux /usr/local/bin/regenie  && chmod +x /usr/local/bin/regenie
```

### bgenix

```
wget http://code.enkre.net/bgen/tarball/release/bgen.tgz && tar -xzf bgen.tgz && cd bgen.tgz/ && ./waf configure && ./waf && cp build/apps/bgenix /usr/bin/
```

### tabix 

```
sudo apt-get install tabix -y
```

## How to run test dataset

### download and prepared folder

```
### install git 
sudo apt-get install git 
## clone repository 
git clone https://github.com/h3abionet/h3agwas-examples
##build a folder, link data
mkdir testgit_qc && cd testgit_qc/ && ln -s ../h3agwas-examples/run_test.bash . && ln -s ../h3agwas-examples/data .
```

### run a test

you can used `run_test.bash`, script has parameter to run the test data set  
Script taked 5 arguments :
 * testdone : script to test :
   * qc : perform qc  
 * h3agwasdir : where find h3agwas repositorie
 * profile slurm ? SIngularity? standart
 * nextflowbin : where find nextflow
 * other option  (optional)

```
./run_test.bash qc h3abionet batch standard nextflow
```


## Quality Control

test on ubuntu 22.04

need :
 * latex 
 * python
 * plink

* [docker image using for test](Docker/qc/)

## GWAS 

test on ubuntu 20.04

need :
 * latex 
 * python 3.6 / 3.8
 * plink 1.9

optional :
 * gemma 
 * Regenie 
 * saige :
  * bgenix for bgen input
  * R
 * gcta
 * bolt-lmm
 * saige :
  * tabix for vcf input
  * bgenix for bgen input





