CC = g++
CFLAGS = -DMAX_KMER_SIZE=64 -march=native -lbifrost -lroaring -pthread -lz -lpthread -std=c++11 -Wall -O3 -mcmodel=medium

all:
	$(CC) PLAST.cpp $(CFLAGS) -o PLAST
