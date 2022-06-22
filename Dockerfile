FROM phoenixbioinformatics/bioperl-centos
# Updating  packages
#RUN yum update -y
# Adding wget and bzip2
RUN yum install -y wget bzip2 git make gcc-c++ zlib-devel deltarpm
# Anaconda installing
RUN wget https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh
RUN bash Anaconda3-5.0.1-Linux-x86_64.sh -b -p /opt/anaconda3 
RUN rm Anaconda3-5.0.1-Linux-x86_64.sh
# Set path to conda
ENV PATH /opt/anaconda3/bin:$PATH
RUN conda update --all
#RUN conda install -c conda-forge vim 
#
WORKDIR /usr/local
RUN git clone https://github.com/lh3/bwa.git
WORKDIR /usr/local/bwa
RUN make
#
#COPY blast-2.2.26/ /usr/local/blast-2.2.26/
#WORKDIR /usr/local/blast-2.2.26
WORKDIR /usr/local
COPY externals/blast-2.2.26-x64-linux.tar.gz /usr/local/blast-2.2.26-x64-linux.tar.gz
RUN tar -xzf blast-2.2.26-x64-linux.tar.gz
ENV BLASTDIR=/usr/local/blast-2.2.26/bin
#
#COPY censor-4.2.29/ /usr/local/censor-4.2.29/
WORKDIR /usr/local
COPY externals/censor-4.2.29.tar.gz /usr/local/censor-4.2.29.tar.gz
RUN tar -xzf censor-4.2.29.tar.gz
WORKDIR /usr/local/censor-4.2.29
RUN sh ./configure
RUN make && make install && make clean
#
RUN pip install pysam
#
WORKDIR /usr/local
ADD https://api.github.com/repos/edwardslab-wustl/TE-detective/git/refs/heads/main version.json
RUN git clone https://github.com/edwardslab-wustl/TE-detective.git
WORKDIR /usr/local/TE-detective 
RUN pip install .
#RUN /opt/anaconda3/bin/python3 -m pip install git+https://github.com/edwardslab-wustl/TE-detective.git
#ENV PATH=/opt/conda/bin:/usr/local/utils/:$BLASTDIR:/usr/local/bwa/$PATH
