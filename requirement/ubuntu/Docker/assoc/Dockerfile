FROM ubuntu:22.04
# image to test qc
RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install software-properties-common -y
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install texlive-latex-base texlive-fonts-recommended texlive-fonts-extra texlive-latex-extra  python3.8 plink1.9 git wget python3-pip -y
RUN pip3 install --upgrade pip
RUN cp /usr/bin/plink1.9 /usr/bin/plink && chmod +x /usr/bin/plink
RUN pip3 install  matplotlib==3.5.3
RUN pip3 install  pandas==1.4.3
RUN pip3 install  scipy==1.9.0

RUN apt-get install curl -y 
RUN wget https://download.oracle.com/java/18/latest/jdk-18_linux-x64_bin.deb && apt-get install ./jdk-18_linux-x64_bin.deb -y && rm -rf ./jdk-18_linux-x64_bin.deb
RUN cat <<EOF | tee /etc/profile.d/jdk18.sh 
ENV JAVA_HOME=/usr/lib/jvm/jdk-18
ENV PATH=$PATH:$JAVA_HOME/bin
RUN wget -qO- https://get.nextflow.io | bash && cp ./nextflow ./bin/
RUN nextflow pull h3abionet/h3agwas
## test final
RUN cd && git clone https://github.com/h3abionet/h3agwas-examples.git 
RUN cd && mkdir testgit_qc && cd testgit_qc/ && ln -s ../h3agwas-examples/run_test.bash . && ln -s ../h3agwas-examples/data . && chmod +x ./run_test.bash
RUN cd /root/testgit_qc && ./run_test.bash qc h3abionet standard nextflow
