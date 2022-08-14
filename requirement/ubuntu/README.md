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

wget https://biogems.info/download/gemma-0.98.6-pre1.gz
gunzip gemma-0.98.6-pre1.gz
```

### gcta

[Pipeline has been tested with version 1.93.beta](https://raw.githubusercontent.com/h3abionet/h3agwas/master/utils/bin/gcta_1.93.2beta.zip)


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


