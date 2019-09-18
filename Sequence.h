#ifndef SEQUENCE_HPP
#define SEQUENCE_HPP

//DNA nucleotide alphabet size
#define SIGMAR 4
//Ranks for characters in nucleotide alphabet
#define CHAR_RANK_A 0
#define CHAR_RANK_C 1
#define CHAR_RANK_G 2
#define CHAR_RANK_T 3
//Character constants used for nucleotide bases
#define NUCL_BASE_A 'A'
#define NUCL_BASE_C 'C'
#define NUCL_BASE_G 'G'
#define NUCL_BASE_T 'T'
//Complements of nucleotide bases
#define CMPL_BASE_A 'T'
#define CMPL_BASE_C 'G'
#define CMPL_BASE_G 'C'
#define CMPL_BASE_T 'A'
//The number of DNA base pairs which can be encoded using 1 byte
#define BASES_PER_BYTE 4
//The number of bits needed to store 1 DNA base pair
#define BITS_PER_BASE 2

class ColSet {

	public:

		//Default constructor
		ColSet (uint32_t pos, vector<string> cols) { endPos = pos; colNames = cols; }

		//Alignment position where color set ends
		uint32_t endPos;
		vector<string> colNames;

};

//Ranking function on nucleotides
inline uint16_t r(char c){
	switch (c){
		case NUCL_BASE_A: return CHAR_RANK_A;
		case NUCL_BASE_C: return CHAR_RANK_C;
		case NUCL_BASE_G: return CHAR_RANK_G;
		case NUCL_BASE_T: return CHAR_RANK_T;
		default: cerr << "Invalid nucleotide character detected\nRank of " << c << " cannot be calculated correctly" << endl;
			 return CHAR_RANK_A;
	}
}

//Inverse ranking function on nucleotides
inline char rr(const char &c){
	switch(c){
		case CHAR_RANK_A: return NUCL_BASE_A;
		case CHAR_RANK_C: return NUCL_BASE_C;
		case CHAR_RANK_G: return NUCL_BASE_G;
		case CHAR_RANK_T: return NUCL_BASE_T;
		default: cerr << "Invalid rank detected\nNucleotide of rank " << c << " cannot be calculated correctly" << endl;
			return NUCL_BASE_A;
	}
}

//This function takes a DNA sequence and compresses it s.t. each base is stored using only 2 bit
inline const char* cmprSeq(const string seq, const uint8_t &size){
	char c = 0;
	uint8_t pos = 0;
	//Allocate memory
	char *cSeq = (char*) malloc(size);

	//Iterate over memory allocated for the compressed sequence
	for(int16_t i = 0; i  < size; ++i){
		//Compress next four bases into a char
		for(uint16_t j = 0; j < BASES_PER_BYTE; ++j, ++pos){
			//Shift previously stored base
			c = c << BITS_PER_BASE;

			//Check whether we have reached seq's end
			if(!seq[pos]){
				//Decrement pos to avoid reading outside of the string
				--pos;	
			}

			//Add the next base
			c |= r(seq[pos]);
		}

		//Add compressed bases to sequence
		cSeq[i] = c;
	}

	return cSeq;
}

//This function takes a compressed DNA sequence and decompesses it
inline const char* decmprSeq(const char* cSeq, const uint8_t seqLen){
	unsigned char c;
	uint16_t i = 0, j;
	//Allocate memory
	char *seq = (char*) malloc(seqLen);

	//Check whether we have decompressed the whole sequence
	while(i < seqLen){
		j = 0;

		//Consider each compressed character separately
		while(j < BASES_PER_BYTE){
			//Get rid of information in front of current character to decompress
			c = *cSeq << j * BITS_PER_BASE;
			//Shift character to end of memory section
			c = c >> 6;
			//Recalculate actual base pair
			seq[i] = rr(c);
			
			//Check whether the whole sequence has been decompressed
			if(++i == seqLen) break;

			++j;
		}

		++cSeq;
	}

	return seq;
}

//Function to compute a k-mer's (falling) rank
const int32_t compRank(const string k, const int32_t &prevRank);

//This function calculates the reverse complement of a DNA sequence
string revComp(const string &seq);

#endif