#ifndef INDEX_HPP
#define INDEX_HPP

#include "Smer.h"
#include "Seedlist.h"

// //The number of DNA base pairs which can be encoded using 1 byte
// #define BASES_PER_BYTE 4
// //The number of bits needed to store 1 DNA base pair
// #define BITS_PER_BASE 2

//This function writes PLAST's q-gram profile to disk
void saveProfile(ofstream &fStr, const uint32_t &profSize, uint32_t* qProfile);

//This function loads the q-gram profile from a file and returns it
uint32_t* loadQProf(ifstream &iFile, const uint32_t size);

// //This function takes a DNA sequence and compresses it s.t. each base is stored using only 2 bit
// inline const char* cmprSeq(const string seq, const uint8_t &size){
// 	char c = 0;
// 	uint8_t pos = 0;
// 	//Allocate memory
// 	char *cSeq = (char*) malloc(size);

// 	//Iterate over memory allocated for the compressed sequence
// 	for(int16_t i = 0; i  < size; ++i){
// 		//Compress next four bases into a char
// 		for(uint16_t j = 0; j < BASES_PER_BYTE; ++j, ++pos){
// 			//Shift previously stored base
// 			c = c << BITS_PER_BASE;

// 			//Check whether we have reached seq's end
// 			if(!seq[pos]){
// 				//Decrement pos to avoid reading outside of the string
// 				--pos;	
// 			}

// 			//Add the next base
// 			c |= r(seq[pos]);
// 		}

// 		//Add compressed bases to sequence
// 		cSeq[i] = c;
// 	}

// 	return cSeq;
// }

// //This function takes a compressed DNA sequence and decompesses it
// inline const char* decmprSeq(const char* cSeq, const uint8_t seqLen){
// 	unsigned char c;
// 	uint16_t i = 0, j;
// 	//Allocate memory
// 	char *seq = (char*) malloc(seqLen);

// 	//Check whether we have decompressed the whole sequence
// 	while(i < seqLen){
// 		j = 0;

// 		//Consider each compressed character separately
// 		while(j < BASES_PER_BYTE){
// 			//Get rid of information in front of current character to decompress
// 			c = *cSeq << j * BITS_PER_BASE;
// 			//Shift character to end of memory section
// 			c = c >> 6;
// 			//Recalculate actual base pair
// 			seq[i] = rr(c);
			
// 			//Check whether the whole sequence has been decompressed
// 			if(++i == seqLen) break;

// 			++j;
// 		}

// 		++cSeq;
// 	}

// 	return seq;
// }

//This function saves PLAST's index data structures as binary
void saveIndexBin(const char *filename, const uint32_t &profSize, uint32_t* qProfile, const uint32_t numUni, UnitigColorMap<seedlist>* uniArr, struct s_mer_pos*& pos, size_t &seedNum, const int32_t &k);

//This function loads the indices from a binary file
void loadIndexesBin(const char* fName, const uint32_t &qProfSize, uint32_t *qProf, ColoredCDBG<seedlist>& graph, UnitigColorMap<seedlist>*& uniArr, size_t &posSize, struct s_mer_pos*& pos);

//This function loads the indexes from a file
void loadIndexes(const string fName, const uint32_t qSize, uint32_t*& prof, CompactedDBG<seedlist> &cdbg, uint32_t &seedNum, struct s_mer_pos*& lArr);

//This function builds all necessary indexes
void buildIndex(ColoredCDBG<seedlist> &cdbg, const int32_t &minSeedLength, uint32_t* qProf, size_t &numSmers, const uint32_t &profSize, struct s_mer_pos*& linkArr, UnitigColorMap<seedlist>*& uniArray);

#endif