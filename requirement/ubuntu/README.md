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

## prepared data set test

```
### install git 
sudo apt-get install git 
## clone repository 
git clone https://github.com/h3abionet/h3agwas-examples
##build a folder, link data
mkdir testgit_qc && cd testgit_qc/ && ln -s ../h3agwas-examples/run_test.bash . && ln -s ../h3agwas-examples/data .
```

## R 

```
```

## SAIGE
 * [installation here](https://saigegit.github.io/SAIGE-doc/)
 * installation R ( default 4.1.2)
 * gcc10 
 * cmake
 * cget using pip 3

```
apt-get install r-base cmake gcc-10 -y
pip3 install cget==0.2.0
```
 
 * download : saige (release 1.0.0)
 * untar
 * install package need
 * for library we need to install libssl-dev libcurl4-openssl-dev
```
apt-get install libssl-dev libcurl4-openssl-dev -y
```

```
wget https://github.com/saigegit/SAIGE/releases/download/1.0.0/SAIGE_1.0.0.tar.gz && tar -xzf SAIGE_1.0.0.tar.gz &&  Rscript ./SAIGE/extdata/install_packages.R

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

* version  2.4 of binary 
```
wget https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/downloads/BOLT-LMM_v2.4.tar.gz && tar -xzf  BOLT-LMM_v2.4.tar.gz && cp BOLT-LMM_v2.4/bolt bin/ && mv BOLT-LMM_v2.4.tar.gz src/ && rm -rf BOLT-LMM_v2.4/
```


```
wget https://github.com/rgcgithub/regenie/releases/download/v3.1.3/regenie_v3.1.3.gz_x86_64_Linux.zip && unzip regenie_v3.1.3.gz_x86_64_Linux.zip && mv regenie_v3.1.3.gz_x86_64_Linux /usr/local/bin/regenie  && chmod +x /usr/local/bin/regenie
```
## software utils

###bgenix

* work on bgen
```
wget http://code.enkre.net/bgen/tarball/release/bgen.tgz && tar -xzf bgen.tgz && cd bgen.tgz/ && ./waf configure && ./waf && cp build/apps/bgenix /usr/bin/
```

##run a test

script run test taking 4 arguments :
 * testdone : script to test qc..
 * h3agwasdir : where find h3agwas repositotu 
 * profile slurm ? SIngularity?
 * nextflowbin : where find nextflow
 * otheroption 

```
./run_test.bash qc h3abionet batch standard nextflow
```


## Quality Control

need :
 * latex 
 * python
 * plink

## GWAS 

need :
 * latex 
 * python
 * plink
 * gemma


Doesn't need other requierment that previous see [docker image](Docker/qc/)


