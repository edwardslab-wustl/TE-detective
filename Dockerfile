FROM continuumio/anaconda3:latest

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    git python-biopython make \
    gcc libc6-dev zlib1g-dev g++ \
    bioperl \
    && apt-get clean
#RUN conda install -c conda-forge biopython
#RUN conda config --add channels r
#RUN conda install bwa

WORKDIR /usr/local
RUN git clone https://github.com/lh3/bwa.git
WORKDIR /usr/local/bwa
RUN make 

COPY blast-2.2.26/ /usr/local/blast-2.2.26/
WORKDIR /usr/local/blast-2.2.26
ENV BLASTDIR=/usr/local/blast-2.2.26/bin

COPY censor-4.2.29/ /usr/local/censor-4.2.29/
WORKDIR /usr/local/censor-4.2.29
RUN sh ./configure
RUN make && make install && make clean

#RUN conda config --add channels bioconda
#RUN conda install -c anaconda libgcc-ng
#RUN conda update anaconda
#RUN conda update conda
#RUN conda install libgcc-ng
#RUN conda install -y -c bioconda pysam
RUN pip install pysam

WORKDIR /usr/local
RUN git clone https://github.com/edwardslab-wustl/TE-detective.git
RUN /opt/conda/bin/python3 -m pip install git+https://github.com/edwardslab-wustl/TE-detective.git

#ENV PATH=/opt/conda/bin:/usr/local/utils/$PATH

