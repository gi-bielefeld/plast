#include "Seedlist.h"

// //This function deletes a seed that has been used during an extension over multiple unitigs
// void delExtPtrSeed(const bool& inSeedList, struct seed*& seedToDel, seedlist* seedList){
// 	//A pointer to save the predecessor of the seed we want to delete
// 	struct seed *prev = NULL;
	
// 	//Check whether our seed is already part of a seed list
// 	if(inSeedList){
// 		//Check whether our seed is the first element in its seed list
// 		if(seedToDel == seedList->getLastProcSeed()){
// 			//Exclude the seed from the seed list
// 			seedList->setLastProcSeed(seedToDel->nextSeed);
// 			//Delete the seed
// 			extractProcSeed(seedList, seedToDel, prev);
// 		} else{
// 			//Go through the seed list to find the seed we are looking for
// 			prev = seedList->getLastProcSeed();

// 			while(prev->nextSeed != seedToDel){ prev = prev->nextSeed; }

// 			//Delete the seed
// 			extractProcSeed(seedList, seedToDel, prev);
// 		}
// 	} else{
// 		//Delete the seed right away
// 		delSeed(seedToDel, prev);
// 	}
// }

// //This function excludes a seed from its seed list
// void excludeSeed(seedlist* sl, const struct seed*& seedToExcl, struct seed*& prevSeed){
// 	//Check if the seed has a predecessor in its seed list
// 	if(prevSeed != NULL){
// 		//Link seed's predecessor and successor
// 		prevSeed->nextSeed = seedToExcl->nextSeed;

// 		// //Delete the extracted seed
// 		// delSeed(seedToExtr, prevSeed);
// 		// //Note that we should look for another seed
// 		// return true;
// 	} else{
// 		//If prevSeed is NULL this either means that the seed to delete was the last remaining seed in the seed list or just the first one. The procedure, however, stays the same
// 		//Set seed's successor as head of the seed list
// 		sl->setSeed(seedToExcl->nextSeed);

// 		// //Delete extracted seed
// 		// delSeed(seedToExtr, prevSeed);
// 		// //There is no further seed to reach
// 		// seedToExtr = NULL;
// 		// //We do not need to look for another seed to reach
// 		// return false;
// 	}
// }

// //This function extracts a seed from a seed list of processed seeds. Returns true if there might be another seed to reach
// bool extractProcSeed(seedlist* sl, struct seed*& seedToExtr, struct seed*& prevSeed){
// 	//First things first - check if we should search for another seed once we are done here
// 	if(prevSeed != NULL){
// 		//Delete the extracted seed
// 		delSeed(seedToExtr, prevSeed);
// 		//Note that we should look for another seed
// 		return true;
// 	} else{
// 		//If prevSeed is NULL this either means that the seed to delete was the last remaining seed in the seed list or just the first one. The procedure, however, stays the same
// 		//Link second seed list's seed to the unitig
// 		sl->setLastProcSeed(seedToExtr->nextSeed);
// 		//Delete extracted seed
// 		delSeed(seedToExtr, prevSeed);
// 		//There is no further seed to reach
// 		seedToExtr = NULL;
// 		//We do not need to look for another seed to reach
// 		return false;
// 	}
// }

//This function searches for the closest seed to reach during an extension and returns it if it exists (otherwise NULL)
//TODO Each time a seed is searched the seed list is gone through from the beginning again. If the last closest seed was not the first one that was found it might be possible not to start from the beginning again!
//TODO The idea behind ordering all seeds in a seed list decreasingly by their q offset was to avoid to iterate over the complete seed list during the search for the closest seed to reach. Unfortunately, since new seeds are now inserted into the seed lists during extension this does not work anymore now. However, if we are able to exchange the score of a seed by a bool to indicate whether it has already been extended or not we might be able to incorporate this idea again. For now we have to skip it <- We can forget about this anyways, because with incorporating the reverse, complementary sequence the ordering is always wrong for one of both lists...
struct seed* searchRightNeighbor(struct seed* sLSeed, const uint32_t &iniExtSeedOffsQ, const uint32_t &extLen, const uint32_t &curUPos, struct seed*& prevSeed){
	struct seed *lastSeed = NULL, *nearestSeed = NULL;
	
	//Reset prevSeed
	prevSeed = NULL;

	// //If there is no seed list there is no nearerst seed
	// if(sLSeed == NULL){
	// 	//Testing
	// 	//cout << "1 Option 1" << endl;

	// 	return nearestSeed;
	// }

	//Testing
	//cout << "1 Option 2" << endl;

	//Find the nearest seed that we might be able to reach during our extension
	while(sLSeed != NULL){
		//Since our seed lists are decreasingly ordered by their offset in q we only need to search for a reachable seed as long as the current seed's q offset is not smaller than the q offset of the seed we are just trying to extend
		if(iniExtSeedOffsQ >= sLSeed->offsetQ){
			//Testing
			//cout << "2 Option 1" << endl;

			return nearestSeed;
		} else{
			//Testing
			//cout << "2 Option 2" << endl;
			//cout << "3 Option ";

			//Check whether we might be able to reach this seed
			if(sLSeed->offsetQ == iniExtSeedOffsQ + extLen - (curUPos - sLSeed->offsetU)){
				//Testing
				//cout << "not ";

				//Store a pointer to the previous seed
				prevSeed = lastSeed;
				//Store the reachable seed
				nearestSeed = sLSeed;
			}

			//Testing
			//cout << "2" << endl;
		}

		lastSeed = sLSeed;
		sLSeed = sLSeed->nextSeed;
	}

	//Testing
	//cout << "4" << endl;

	return nearestSeed;
}

//This function searches for the closed seed to reach during an extension and returns it if it exists (otherwise NULL)
struct seed* searchLeftNeighbor(struct seed *seedList, const uint32_t &qOff, const uint32_t &uOff, struct seed*& prevSeed){
	struct seed *nearestSeed = NULL;
	struct seed *lastSeed = NULL;

	//Testing
	// cout << "Start of searchLeftNeighbor" << endl;

	//Check if we are dealing with an empty seed list
	// if(seedList == NULL){
	// 	//Testing
	// 	//cout << "1 Option 1" << endl;

	// 	return nearestSeed;
	// }

	//Testing
	//cout << "1 Option 2" << endl;
	// if(report){
	// 	cerr << "searchLeftNeighbor" << endl;
	// }
	// cout << "Seedlist is not empty" << endl;

	//Iterate over the seed list
	while(seedList != NULL){
		//Testing
		// cout << "Iterate over seedlist" << endl;
		// cout << "Seed is offsetU: " << seedList->offsetU << " offsetQ: " << seedList->offsetQ << " len: " << seedList->len << " nextSeed is " << (seedList->nextSeed == NULL ? "" : "not ") << "NULL" << endl;

		//Check whether we could reach the current seed
		if(seedList->offsetQ <= qOff && seedList->offsetU <= uOff && qOff - seedList->offsetQ == uOff - seedList->offsetU){
			//Testing
			//cout << "2 Option 1" << endl;

			//Save the seed as the nearest
			nearestSeed = seedList;
			//Save last iteration's seed as prevSeed
			prevSeed = lastSeed;
			break;
		}
		
		//Testing
		//cout << "2 Option 2" << endl;

		//Update lastSeed and seedList
		lastSeed = seedList;
		seedList = seedList->nextSeed;
	}

	//Testing
	//cout << "3" << endl;
	// cout << "Done" << endl;

	return nearestSeed;
}