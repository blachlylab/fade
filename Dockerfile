FROM ubuntu:hirsute

RUN apt-get update  # Ensure the package list is up to date
RUN apt-get install -y wget autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libphobos2-ldc-shared94 bzip2

ENV ZIP=htslib-1.12.tar.bz2
ENV URL=https://github.com/samtools/htslib/releases/download/1.12
ENV FOLDER=htslib-1.12
ENV DST=/tmp
ENV LD_LIBRARY_PATH=/usr/local/lib

RUN wget $URL/$ZIP -O $DST/$ZIP && \
    tar xvf $DST/$ZIP -C $DST && \
    rm $DST/$ZIP && \
    cd $DST/$FOLDER && \
    make && \
    make install && \
    cd / && \
    rm -rf $DST/$FOLDER

RUN wget https://github.com/jeffdaily/parasail/releases/download/v2.4.3/parasail-2.4.3-manylinux1_x86_64.tar.gz && \
    tar -xzf parasail-2.4.3-manylinux1_x86_64.tar.gz && \
    cd parasail-2.4.3-manylinux1_x86_64 && \
    cd lib/ && \
    cp *.* /usr/local/lib/ && \
    cd / && \
    rm -rf parasail-2.4.3-manylinux1_x86_64

RUN wget https://github.com/blachlylab/fade/releases/download/v0.3.1/fade && \
    chmod +x fade && \
    mv fade /usr/local/bin

# DEFINE DEFAULT COMMAND
ENTRYPOINT ["/usr/local/bin/fade"]
