FROM ubuntu:20.04
LABEL maintainer="jean-tristan.brandenburg@wits.ac.za"
# image to test qc
RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install software-properties-common -y
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install texlive-latex-base texlive-fonts-recommended texlive-fonts-extra texlive-latex-extra  python3.8 plink1.9 git wget python3-pip -y
RUN pip3 install --upgrade pip
RUN cp /usr/bin/plink1.9 /usr/bin/plink && chmod +x /usr/bin/plink
RUN pip3 install  matplotlib==3.5.3
RUN pip3 install  pandas==1.4.3
RUN pip3 install  scipy==1.9.0

RUN wget https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.4/gemma-0.98.4-linux-static-AMD64.gz && gunzip gemma-0.98.4-linux-static-AMD64.gz && mv gemma-0.98.4-linux-static-AMD64 gemma && chmod +x gemma && mv gemma /usr/local/bin/
RUN wget https://raw.githubusercontent.com/h3abionet/h3agwas/master/utils/bin/gcta_1.93.2beta.zip && unzip gcta_1.93.2beta.zip && cp gcta_1.93.2beta/gcta64 . && rm -rf gcta_1.93.2beta && rm -rf gcta_1.93.2beta.zip  && mv gcta64 /usr/local/bin/
RUN wget https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/downloads/BOLT-LMM_v2.4.tar.gz && tar -xzf  BOLT-LMM_v2.4.tar.gz && cp BOLT-LMM_v2.4/bolt /usr/local/bin/ && rm -rf BOLT-LMM_v2.4.tar.gz  && rm -rf BOLT-LMM_v2.4/


## for saige 
RUN apt-get install r-base \
    libopenblas-base \
    python3-pip \
    r-cran-devtools -y

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    cmake \
    libopenblas-base \
    python3-pip \
    r-cran-devtools

RUN pip3 install cget
RUN apt-get install git -y
RUN git clone https://github.com/saigegit/SAIGE

WORKDIR /SAIGE
COPY . .

RUN Rscript extdata/install_packages.R

# Force step_2 to use 1 single thread. More threads are ineffective
ENV OMP_NUM_THREADS=1

RUN R CMD INSTALL .

RUN mv extdata/step1_fitNULLGLMM.R extdata/step2_SPAtests.R extdata/step3_LDmat.R extdata/createSparseGRM.R /usr/local/bin/

RUN chmod a+x /usr/local/bin/step1_fitNULLGLMM.R
RUN chmod a+x /usr/local/bin/step2_SPAtests.R
RUN chmod a+x /usr/local/bin/step3_LDmat.R
RUN chmod a+x /usr/local/bin/createSparseGRM.R

RUN createSparseGRM.R  --help
RUN step1_fitNULLGLMM.R --help
RUN step2_SPAtests.R --help
RUN step3_LDmat.R --help

RUN apt-get update
RUN apt-get install time

RUN  wget http://code.enkre.net/bgen/tarball/release/bgen.tgz && tar -xzf bgen.tgz && cd bgen.tgz/ && ./waf configure && ./waf && cp build/apps/bgenix /usr/bin/

RUN wget https://github.com/rgcgithub/regenie/releases/download/v3.1.3/regenie_v3.1.3.gz_x86_64_Linux.zip && unzip regenie_v3.1.3.gz_x86_64_Linux.zip && mv regenie_v3.1.3.gz_x86_64_Linux /usr/local/bin/regenie  && chmod +x /usr/local/bin/regenie


#
#
#RUN wget https://download.oracle.com/java/18/latest/jdk-18_linux-x64_bin.deb && apt-get install ./jdk-18_linux-x64_bin.deb -y && rm -rf ./jdk-18_linux-x64_bin.deb
#RUN cat <<EOF | tee /etc/profile.d/jdk18.sh 
#ENV JAVA_HOME=/usr/lib/jvm/jdk-18
#ENV PATH=$PATH:$JAVA_HOME/bin
#RUN wget -qO- https://get.nextflow.io | bash && cp ./nextflow /usr/local/bin/
#RUN nextflow pull h3abionet/h3agwas
### test final
#RUN cd && git clone https://github.com/h3abionet/h3agwas-examples.git 
#RUN cd && mkdir testgit_qc && cd testgit_qc/ && ln -s ../h3agwas-examples/run_test.bash . && ln -s ../h3agwas-examples/data . && chmod +x ./run_test.bash
#RUN cd /root/testgit_qc && ./run_test.bash qc h3abionet standard nextflow
