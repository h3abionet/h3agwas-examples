# General instalation 


Installation recommended has tested using docker image

##latex :

```
sudo apt-get install texlive-latex-base texlive-fonts-recommended texlive-fonts-extra texlive-latex-extra
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

# Quality Control

Doesn't need other requierment that previous see [docker image](Docker/qc/)


