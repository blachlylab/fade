FROM alpine

RUN apk add --no-cache \
    wget autoconf automake make gcc perl git \
    llvm clang musl musl-dev \
    ldc ldc-static dub \
    zlib zlib-dev zlib-static \
    xz xz-dev \
    binutils-gold \
    bzip2 bzip2-dev bzip2-static \
    llvm-libunwind-static

WORKDIR /home/

# build libdeflate
# note we are using the same libdeflate that rust-htslib is using
ENV DEFLATE_VER=1.8
RUN wget https://github.com/ebiggers/libdeflate/archive/refs/tags/v$DEFLATE_VER.tar.gz
RUN tar -xf v$DEFLATE_VER.tar.gz
WORKDIR libdeflate-$DEFLATE_VER
RUN make install -j 8 CFLAGS="-fPIC -ffreestanding -DFREESTANDING"
WORKDIR /home/

# build openssl
ENV SSL_VER=3.0.0
RUN wget https://www.openssl.org/source/openssl-$SSL_VER.tar.gz
RUN tar -xf openssl-$SSL_VER.tar.gz
WORKDIR openssl-$SSL_VER
RUN ./Configure --prefix=/usr/local linux-x86_64 no-async no-engine -DOPENSSL_NO_SECURE_MEMORY
RUN make -j 8
RUN make install -j 8
WORKDIR /home/

# build libssh2
ENV SSH2_VER=1.10.0
RUN wget https://github.com/libssh2/libssh2/releases/download/libssh2-$SSH2_VER/libssh2-$SSH2_VER.tar.gz
RUN tar -xf libssh2-$SSH2_VER.tar.gz
WORKDIR libssh2-$SSH2_VER
RUN CFLAGS="-I/usr/local/include" LDFLAGS="-L/usr/local/lib64" ./configure --prefix=/usr/local
RUN make -j 8
RUN make install
WORKDIR /home/

# build libcurl
ENV CURL_VER=7.80.0
RUN wget https://github.com/curl/curl/releases/download/curl-7_80_0/curl-$CURL_VER.tar.gz
RUN tar -xf curl-$CURL_VER.tar.gz
WORKDIR curl-$CURL_VER
RUN ./configure --with-libssh2 --with-ssl --prefix=/usr/local --enable-get-easy-options --enable-ftp
RUN make -j 8
RUN make install -j 8
WORKDIR /home/

# build htslib
ENV HTS_VER=1.14
RUN wget https://github.com/samtools/htslib/releases/download/$HTS_VER/htslib-$HTS_VER.tar.bz2
RUN tar xf htslib-$HTS_VER.tar.bz2
WORKDIR htslib-$HTS_VER
RUN ./configure --with-libdeflate --enable-libcurl
RUN make -j 8
RUN make install

WORKDIR /home/

# download parasail
ENV PARASAIL_VER=2.4.3
RUN wget https://github.com/jeffdaily/parasail/releases/download/v$PARASAIL_VER/parasail-$PARASAIL_VER-manylinux1_x86_64.tar.gz
RUN tar -xf parasail-$PARASAIL_VER-manylinux1_x86_64.tar.gz
RUN cp -r parasail-$PARASAIL_VER-manylinux1_x86_64/* /usr/local

# set dflags for ldc static building
ENV DFLAGS="-link-defaultlib-shared=false -static --linker=gold"
ENV DFLAGS="$DFLAGS -L-L/lib -L-L/usr/local/lib -L-L/usr/local/lib64 -L-L/usr/lib"
ENV DFLAGS="$DFLAGS -L-lz -L-lbz2 -L-ldeflate -L-llzma"
ENV DFLAGS="$DFLAGS -L-lcurl -L-lssl -L-lssh2 -L-lcrypto -L-lunwind"

ADD . /home/fade

WORKDIR /home/fade
RUN dub build --compiler ldc2

RUN cp fade /usr/local/bin

ENTRYPOINT ["/usr/local/bin/fade"]
