FROM ubuntu:26.04
LABEL org.opencontainers.image.source="https://github.com/gi-bielefeld/plast"
#Install build essentials, python and unzip
RUN apt-get update \
	&& apt-get install -y build-essential cmake zlib1g-dev git \
	&& apt-get install -y --no-install-recommends python3 unzip \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*
#Install Bifrost and adjust library paths
RUN git clone https://github.com/pmelsted/bifrost.git \
	&& cd bifrost \
	&& mkdir build \
	&& cd build \
	&& cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCOMPILATION_ARCH=OFF -DMAX_KMER_SIZE=64 .. \
	&& make \
	&& make install
ENV C_INCLUDE_PATH=/usr/local/include/
ENV CPLUS_INCLUDE_PATH=/usr/local/include/
ENV LD_LIBRARY_PATH=/usr/local/lib
ENV LIBRARY_PATH=/usr/local/lib
#Install PLAST
RUN git clone https://github.com/gi-bielefeld/plast.git \
	&& cd plast/src \
	&& sed -i 's/march=native/DMAX\_KMER\_SIZE=64/g' makefile \
	&& make \
	&& cd /usr/local/bin \
	&& ln -s /plast/src/PLAST
