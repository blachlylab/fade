FROM alpine

RUN apk add --no-cache \
    wget autoconf automake make gcc g++ perl git cmake \
    llvm clang musl musl-dev \
    ldc ldc-static dub \
    zlib zlib-dev zlib-static \
    xz xz-dev \
    binutils-gold \
    bzip2 bzip2-dev bzip2-static \
    llvm-libunwind-static

WORKDIR /home/

# apparently musl-c allocation is bad with multiple threads
# suggestion was to use mimalloc
# RUN git clone https://github.com/microsoft/mimalloc.git
# WORKDIR /home/mimalloc
# RUN git checkout v2.0.3
# RUN mkdir build 
# WORKDIR /home/mimalloc/build
# RUN cmake ..
# RUN make -j 8
# RUN make install
# WORKDIR /home/
# RUN rm -r mimalloc

# build libdeflate
# note we are using the same libdeflate that rust-htslib is using
ENV DEFLATE_VER=1.8
RUN wget https://github.com/ebiggers/libdeflate/archive/refs/tags/v$DEFLATE_VER.tar.gz
RUN tar -xf v$DEFLATE_VER.tar.gz
WORKDIR libdeflate-$DEFLATE_VER
RUN make install -j 8
WORKDIR /home/
RUN rm -r libdeflate-$DEFLATE_VER

# build openssl
ENV SSL_VER=3.0.0
RUN wget https://www.openssl.org/source/openssl-$SSL_VER.tar.gz
RUN tar -xf openssl-$SSL_VER.tar.gz
WORKDIR openssl-$SSL_VER
RUN ./Configure --prefix=/usr/local linux-x86_64 no-async no-engine -DOPENSSL_NO_SECURE_MEMORY
RUN make -j 8
RUN make install_sw -j 8
WORKDIR /home/
RUN rm -r openssl-$SSL_VER

# build libssh2
ENV SSH2_VER=1.10.0
RUN wget https://github.com/libssh2/libssh2/releases/download/libssh2-$SSH2_VER/libssh2-$SSH2_VER.tar.gz
RUN tar -xf libssh2-$SSH2_VER.tar.gz
WORKDIR libssh2-$SSH2_VER
RUN CFLAGS="-I/usr/local/include" LDFLAGS="-L/usr/local/lib64" ./configure --prefix=/usr/local
RUN make -j 8
RUN make install
WORKDIR /home/
RUN rm -r libssh2-$SSH2_VER

# build libcurl
ENV CURL_VER=7.80.0
RUN wget https://github.com/curl/curl/releases/download/curl-7_80_0/curl-$CURL_VER.tar.gz
RUN tar -xf curl-$CURL_VER.tar.gz
WORKDIR curl-$CURL_VER
RUN ./configure --with-libssh2 --with-ssl --prefix=/usr/local --enable-get-easy-options --enable-ftp
RUN make -j 8
RUN make install -j 8
WORKDIR /home/
RUN rm -r curl-$CURL_VER

# build htslib
ENV HTS_VER=1.14
RUN wget https://github.com/samtools/htslib/releases/download/$HTS_VER/htslib-$HTS_VER.tar.bz2
RUN tar xf htslib-$HTS_VER.tar.bz2
WORKDIR htslib-$HTS_VER
RUN ./configure
RUN make -j 8
RUN make install
WORKDIR /home/
RUN rm -r htslib-$HTS_VER

# download parasail
ENV PARASAIL_VER=2.4.3
RUN wget https://github.com/jeffdaily/parasail/releases/download/v$PARASAIL_VER/parasail-$PARASAIL_VER.tar.gz
RUN tar -xf parasail-$PARASAIL_VER.tar.gz
WORKDIR parasail-$PARASAIL_VER
RUN ./configure
RUN make -j 8
RUN make install

WORKDIR /home/
RUN rm -r parasail-$PARASAIL_VER

ADD . /home/fade

WORKDIR /home/fade

RUN dub build --compiler ldc2 -c static-alpine -b release

RUN cp fade /usr/local/bin

ENTRYPOINT ["/usr/local/bin/fade"]
