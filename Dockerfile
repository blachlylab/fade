FROM blachlylab/dlang-htslib-static

# download parasail
ENV PARASAIL_VER=2.4.3
RUN wget https://github.com/jeffdaily/parasail/releases/download/v$PARASAIL_VER/parasail-$PARASAIL_VER.tar.gz \
    && tar -xf parasail-$PARASAIL_VER.tar.gz \
    && cd parasail-$PARASAIL_VER \
    && ./configure \
    && make -j 8 \
    && make install \
    && cd /home/ \
    && rm -r parasail-$PARASAIL_VER

ADD . /home/fade

WORKDIR /home/fade

ARG FADE_VER=$(git describe --tags)

RUN dub build --compiler ldc2 -c static-alpine -b release

RUN cp fade /usr/local/bin

ENTRYPOINT ["/usr/local/bin/fade"]
