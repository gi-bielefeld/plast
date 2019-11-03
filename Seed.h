#ifndef SEED_HPP
#define SEED_HPP

//Data structure to store a seed
struct seed{
	//Offset in the unitig
	uint32_t offsetU;
	//Offset in the query
	uint32_t offsetQ;
	//The seed's length
	uint32_t len;
	//The seed's score (A score of 0 means a seed has not been extended yet)
	uint32_t score;
	//Pointer to the next seed of unitig's seed list
	struct seed *nextSeed;
};

#endif
