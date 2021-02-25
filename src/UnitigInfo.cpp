#include "UnitigInfo.h"

//This function searches for the closest seed to reach during an extension and returns it if it exists (otherwise NULL)
struct Seed* searchRightNeighbor(struct Seed* sLSeed, const uint32_t &iniExtSeedOffsQ, const uint32_t &extLen, const uint32_t &curUPos, struct Seed*& prevSeed){
	struct Seed *lastSeed = NULL, *nearestSeed = NULL;
	
	//Reset prevSeed
	prevSeed = NULL;

	//Find the nearest seed that we might be able to reach during our extension
	while(sLSeed != NULL){
		//Since our seed lists are decreasingly ordered by their offset in q we only need to search for a reachable seed as long as the current seed's q offset is not smaller than the q offset of the seed we are just trying to extend
		if(iniExtSeedOffsQ >= sLSeed->offsetQ){
			return nearestSeed;
		} else{
			//Check whether we might be able to reach this seed
			if(sLSeed->offsetQ == iniExtSeedOffsQ + extLen - (curUPos - sLSeed->offsetU)){
				//Store a pointer to the previous seed
				prevSeed = lastSeed;
				//Store the reachable seed
				nearestSeed = sLSeed;
			}
		}

		lastSeed = sLSeed;
		sLSeed = sLSeed->nextSeed;
	}

	return nearestSeed;
}

//This function searches for the closed seed to reach during an extension and returns it if it exists (otherwise NULL)
struct Seed* searchLeftNeighbor(struct Seed *seedList, const uint32_t &qOff, const uint32_t &uOff, struct Seed*& prevSeed){
	struct Seed *nearestSeed = NULL;
	struct Seed *lastSeed = NULL;

	//Iterate over the seed list
	while(seedList != NULL){
		//Check whether we could reach the current seed
		if(seedList->offsetQ <= qOff && seedList->offsetU <= uOff && qOff - seedList->offsetQ == uOff - seedList->offsetU){
			//Save the seed as the nearest
			nearestSeed = seedList;
			//Save last iteration's seed as prevSeed
			prevSeed = lastSeed;
			break;
		}
		
		//Update lastSeed and seedList
		lastSeed = seedList;
		seedList = seedList->nextSeed;
	}

	return nearestSeed;
}