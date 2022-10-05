FROM phoenixbioinformatics/bioperl-centos

# Updating and install necessary linux packages
#RUN yum update -y
RUN yum install -y wget bzip2 git make gcc-c++ zlib-devel deltarpm

# Install Anaconda 
RUN wget https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh
RUN bash Anaconda3-5.0.1-Linux-x86_64.sh -b -p /opt/anaconda3 
RUN rm Anaconda3-5.0.1-Linux-x86_64.sh

# Set path to conda
ENV PATH /opt/anaconda3/bin:$PATH
RUN conda update --all
#RUN conda install -c conda-forge vim 

# Install bwa 
WORKDIR /usr/local
RUN git clone https://github.com/lh3/bwa.git
WORKDIR /usr/local/bwa
RUN make

# Install blast from externals
WORKDIR /usr/local
COPY externals/blast-2.2.26-x64-linux.tar.gz /usr/local/blast-2.2.26-x64-linux.tar.gz
RUN tar -xzf blast-2.2.26-x64-linux.tar.gz
ENV BLASTDIR=/usr/local/blast-2.2.26/bin

# Install Censor from externals
WORKDIR /usr/local
COPY externals/censor-4.2.29.tar.gz /usr/local/censor-4.2.29.tar.gz
RUN tar -xzf censor-4.2.29.tar.gz
WORKDIR /usr/local/censor-4.2.29
RUN sh ./configure
RUN make && make install && make clean

# Install pysam
RUN pip install pysam

# Install TE-detective
WORKDIR /usr/local
ADD https://api.github.com/repos/edwardslab-wustl/TE-detective/git/refs/heads/main version.json
RUN git clone https://github.com/edwardslab-wustl/TE-detective.git
#COPY /Users/EdwardsJ/source_code/TE-detective/ /usr/local/TE-detective/
WORKDIR /usr/local/TE-detective 
RUN pip install --use-feature=in-tree-build .

# set up PATH
ENV PATH=/opt/conda/bin:/usr/local/utils/:$BLASTDIR:/usr/local/bwa/:$PATH
