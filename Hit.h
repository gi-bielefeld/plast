#ifndef HIT_HPP
#define HIT_HPP

#include "Seedlist.h"

struct hit{
	//The hit's length
	uint32_t length;
	//The hit's score
	uint32_t score;
	//The hit's left border offsets
	uint32_t lSeedUoff, lSeedQoff;
	//The hit's right border offsets
	uint32_t rSeedUoff, rSeedQoff;
	//Unitigs where the borders are located
	UnitigColorMap<seedlist> lUnitig, rUnitig;
	//Pointer to the next hit
	struct hit *nextHit;
};

//This function inserts a hit into the hit list where all hits are order increasingly by their left border's q index
void insertHit(struct hit*& hitToIns, struct hit* hitList);

//This function scans through a hit list decreasingly ordered by score and replaces the worst entry by the given hit if it is not the worst one itself
void replWorseRes(list<struct hit*> &hList, struct hit* hit);

//This function scans through a hit list and inserts the given hit right in front of the first entry which score is smaller or equal to the given hit's score
void insRes(list<struct hit*> &hList, struct hit* hit);

//This function retreats an ungapped extension of a seed on the reverse complementary strand if its right border ends within the k-1 overlap in a unitig sequence's end
void makeValidBorders(struct hit *hit, const string &q);

#endif