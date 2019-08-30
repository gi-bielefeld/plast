#ifndef SEED_HPP
#define SEED_HPP

//Data structure to store a seed
struct seed{//TODO Think about a better name for this
	//Offset in the unitig
	uint32_t offsetU;//TODO The default maximum length of a unitig is 100000, so this should be completely sufficient. However, keep in mind that the default might change!
	//Offset in the query
	uint32_t offsetQ;//TODO Keep in mind that the query might be longer in general!
	//The seed's length
	uint32_t len;//TODO If the data type for offsetU has to be changed, this one has to be changed too!
	//The seed's score (A score of 0 means a seed has not been extended yet)//TODO Do we actually still need the score?
	uint32_t score;
	//Pointer to the next seed of unitig's seed list
	struct seed *nextSeed;

	//Constructor to initalize *nextSeed TODO Why is this not working? - Probably the constructor is never called...
	//seed(): nextSeed(NULL){}
};

#endif
