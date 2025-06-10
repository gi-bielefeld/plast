FROM ubuntu
#Install build essentials
RUN apt-get update && apt-get install -y build-essential cmake zlib1g-dev git && apt-get clean
#Install Bifrost and adjust library paths
#Notes:
#	-All version of Bifrost >1.2.1 currently seem to suffer from a bug. Thus, we use version 1.2.1 here
#	-Here we build Bifrost for a maximal k-mer size of 64 which seems sufficient at this point...
RUN git clone https://github.com/pmelsted/bifrost.git && cd bifrost && git checkout v1.2.1 && mkdir build && cd build && cmake -DMAX_KMER_SIZE=64 .. && make -j && make install
ENV C_INCLUDE_PATH=/usr/local/include/
ENV CPLUS_INCLUDE_PATH=/usr/local/include/
ENV LD_LIBRARY_PATH=/usr/local/lib
ENV LIBRARY_PATH=/usr/local/lib
#Install PLAST
RUN git clone https://github.com/gi-bielefeld/plast.git && cd plast/src && sed -i 's/O3/O3 -DMAX\_KMER\_SIZE=64/g' makefile && make
ENTRYPOINT ["/bin/bash", "-c", "/plast/src/PLAST"]