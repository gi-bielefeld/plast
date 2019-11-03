#include "Sequence.h"

//Function to compute a k-mer's (falling) rank
const int32_t compRank(const string k, const int32_t &prevRank){
	//Rank to be computed
	int32_t rank;

	if(prevRank >= 0){
		return ((prevRank) % (int32_t) pow(SIGMAR, k.length() - 1)) * SIGMAR + r(k.back());
	} else{
		rank = 0;

		for(uint32_t i = 1; i <= k.length(); ++i){
			rank += r(k[i-1]) * (int32_t) pow(SIGMAR, k.length() - i);
		}
	}

	return rank;
}

//This function calculates the reverse complement of a DNA sequence
string revComp(const string &seq){
	//The result string
	string revSeq;

	//Go through the query from the end to the beginning
	for(int32_t i = seq.length() - 1; i >= 0; --i){
		//Check which base we are dealing with and append its complement
		switch(seq[i]){
			case NUCL_BASE_A:
				revSeq += CMPL_BASE_A;
				break;
			case NUCL_BASE_C:
				revSeq += CMPL_BASE_C;
				break;
			case NUCL_BASE_G:
				revSeq += CMPL_BASE_G;
				break;
			case NUCL_BASE_T:
				revSeq += CMPL_BASE_T;
				break;
			default:
				cerr << "ERROR: Unknown nucleotide base detected in query. Only A,C,G and T are supported" << endl;
				exit(EXIT_FAILURE);
		}
	}

	return revSeq;
}