#include "Search.h"
// #include "Extension.h"
#include "Extension.cpp"
#include "GappedAlignment.cpp"
#include "Sequence.cpp"
#include "Hit.cpp"
#include "Statistics.h"

// //This function performs the seed detection if no search color set is given
// void detectSeeds(const int32_t &k, const int32_t &minSeedLength, const size_t &numSmers, const uint32_t &quorum, const uint32_t &profileSize, const uint32_t *qProfile, const string &q, const UnitigColorMap<seedlist> *uArr, const struct s_mer_pos *posArray, hit *hitArr, const bool isRefSeq){
// 	uint32_t occnum;
// 	int32_t pRank = -1;
// 	size_t uniLen;
// 	struct seed *newSeed, *lastSeed;

// 	//Go through the query and check for each s-mer whether it appears in the graph
// 	for(uint32_t i = 0; i < q.length() - minSeedLength + 1; ++i){
// 		//Once we are iterating over q anyways initialize hitArr NOTE: Initialized hits have a length of 0//TODO: This shouldn't be in here!
// 		hitArr[i].length = 0;
// 		//Calculate current s-mer's rank
// 		pRank = compRank(q.substr(i, minSeedLength), pRank);

// 		//Testing
// 		//cout << "Do we get here as well?" << endl;

// 		//To avoid border effects check whether we are at the end of the q-gram profile
// 		if(static_cast<unsigned>(pRank) < profileSize - 1 && qProfile[pRank] < numSmers){
// 			//Testing
// 			//cout << "Testcase: Not at the end of s-mer index" << endl;

// 			occnum = qProfile[pRank + 1] - qProfile[pRank];//if occnum is zero the current s-mer doesn't exist
// 		} else{
// 			//Testing
// 			// cout << "Testcase: At the end of s-mer index" << endl;
// 			// cout << "compRank(q.substr(i, minSeedLength), pRank): " << compRank(q.substr(i, minSeedLength), pRank) << endl;

// 			occnum = numSmers - qProfile[pRank];//if occnum is zero the current s-mer doesn't exist
// 		}

// 		//Process all occurrences of a certain s-mer
// 		for(int32_t j = occnum - 1; j >= 0; --j){//for(uint32_t j = 0; j < occnum; ++j){
// 			//Get current unitig's length
// 			uniLen = uArr[posArray[qProfile[pRank] + j].unitig_id].size;

// 			//Check if the seed list is empty
// 			if(uArr[posArray[qProfile[pRank] + j].unitig_id].getData()->getData(uArr[posArray[qProfile[pRank] + j].unitig_id])->getSeed(isRefSeq) == NULL){
// 				//Check whether our quorum is fulfilled move on otherwise
// 				if(checkSearchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], quorum, posArray[qProfile[pRank] + j].offset) == 0){
// 					//Testing
// 					//cout << "Testcase: S-mer has to be inserted as the first element in a seed list and quorum is not fulfilled" << endl;

// 					continue;
// 				}

// 				//Testing
// 				// if(i == 489 && !strcmp(uArr[posArray[qProfile[pRank] + j].unitig_id].referenceUnitigToString().c_str(), s)){
// 				// 	cout << "j:" << j << endl;
// 				// 	exit(0);
// 				// }
// 				//cout << "Testcase: S-mer has to be inserted as the first element in a seed list and quorum is fulfilled" << endl;

// 				//Create a new seed and...
// 				newSeed = (struct seed*) malloc(sizeof(struct seed));
// 				newSeed->offsetU = compOffset(posArray[qProfile[pRank] + j].offset, minSeedLength, uniLen, isRefSeq);
// 				newSeed->offsetQ = compOffset(i, minSeedLength, q.length(), isRefSeq);
// 				newSeed->len = minSeedLength;
// 				newSeed->score = 0;
// 				newSeed->nextSeed = NULL;

// 				//Testing
// 				// if(newSeed->offsetU == 20 && newSeed->offsetQ == 8024 && !isRefSeq){
// 				// 	cerr << "Seed uoff:" << newSeed->offsetU << " qoff:" << newSeed->offsetQ << " has been created" << endl << "Unitig is " << uArr[posArray[qProfile[pRank] + j].unitig_id].mappedSequenceToString() << endl;
// 				// }

// 				//...make it the first element of the seed list
// 				uArr[posArray[qProfile[pRank] + j].unitig_id].getData()->getData(uArr[posArray[qProfile[pRank] + j].unitig_id])->pushSeed(newSeed, isRefSeq);
// 			} else{
// 				//Testing
// 				//cout << "Testcase: S-mer has not to be inserted as the first element in a seed list and quorum is " << (checkSearchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], quorum, posArray[qProfile[pRank] + j].offset) != 0 ? "" : "not ") << "fulfilled" << endl;
// 				//Check on which strand we are
// 				if(isRefSeq){
// 					//Get the last inserted seed
// 					lastSeed = uArr[posArray[qProfile[pRank] + j].unitig_id].getData()->getData(uArr[posArray[qProfile[pRank] + j].unitig_id])->getSeed(isRefSeq);

// 					//Check if we can extend the last seed inserted seed in our seed list
// 					if(lastSeed->offsetU + lastSeed->len - (minSeedLength - 1) == posArray[qProfile[pRank] + j].offset && lastSeed->offsetQ + lastSeed->len - (minSeedLength - 1) == i){
// 						//Check whether we have to check our quorum at the current position
// 						if(lastSeed->len % k == 0){
// 							//Check whether our quorum is fulfilled move on otherwise
// 							if(checkSearchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], quorum, posArray[qProfile[pRank] + j].offset) == 0){ continue; }
// 						}

// 						//Extending simply means to increment the seed's length
// 						++lastSeed->len;
// 					} else {
// 						//Check whether our quorum is fulfilled move on otherwise
// 						if(checkSearchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], quorum, posArray[qProfile[pRank] + j].offset) == 0){ continue; }

// 						//If an extension is not possible create a new seed
// 						newSeed = (struct seed*) malloc(sizeof(struct seed));
// 						newSeed->offsetU = posArray[qProfile[pRank] + j].offset;
// 						newSeed->offsetQ = i;
// 						newSeed->len = minSeedLength;
// 						newSeed->score = 0;
// 						newSeed->nextSeed = NULL;
// 						//Insert the new seed into the seed list
// 						uArr[posArray[qProfile[pRank] + j].unitig_id].getData()->getData(uArr[posArray[qProfile[pRank] + j].unitig_id])->pushSeed(newSeed, isRefSeq);
// 					}
// 				} else{
// 					//Get the last inserted seed
// 					lastSeed = uArr[posArray[qProfile[pRank] + j].unitig_id].getData()->getData(uArr[posArray[qProfile[pRank] + j].unitig_id])->getLastInSd();

// 					//Check if we can extend the last seed inserted seed in our seed list
// 					if(lastSeed->offsetU > 0 && lastSeed->offsetU - 1 == compOffset(posArray[qProfile[pRank] + j].offset, minSeedLength, uniLen, isRefSeq) && lastSeed->offsetQ > 0 && lastSeed->offsetQ - 1 == compOffset(i, minSeedLength, q.length(), isRefSeq)){
// 						//Check whether we have to check our quorum at the current position
// 						if(lastSeed->len % k == 0){
// 							//Check whether our quorum is fulfilled move on otherwise
// 							if(checkSearchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], quorum, posArray[qProfile[pRank] + j].offset) == 0){ continue; }
// 						}

// 						//Extending means to increment the seed's length and to adjust the offsets
// 						++lastSeed->len;
// 						--lastSeed->offsetU;
// 						--lastSeed->offsetQ;
// 					} else{
// 						//Check whether our quorum is fulfilled move on otherwise
// 						if(checkSearchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], quorum, posArray[qProfile[pRank] + j].offset) == 0){ continue; }

// 						//If an extension is not possible create a new seed
// 						newSeed = (struct seed*) malloc(sizeof(struct seed));
// 						newSeed->offsetU = compOffset(posArray[qProfile[pRank] + j].offset, minSeedLength, uniLen, isRefSeq);
// 						newSeed->offsetQ = compOffset(i, minSeedLength, q.length(), isRefSeq);
// 						newSeed->len = minSeedLength;
// 						newSeed->score = 0;
// 						newSeed->nextSeed = NULL;

// 						//Testing
// 						// if(newSeed->offsetU == 20 && newSeed->offsetQ == 8024 && !isRefSeq){
// 						// 	cerr << "Seed uoff:" << newSeed->offsetU << " qoff:" << newSeed->offsetQ << " has been created" << endl << "Unitig is " << uArr[posArray[qProfile[pRank] + j].unitig_id].mappedSequenceToString() << endl;
// 						// 	cerr << "Real offsets were u:" << posArray[qProfile[pRank] + j].offset << " q:" << i << endl << "unitig length:" << uArr[posArray[qProfile[pRank] + j].unitig_id].size << endl;
// 						// 	//exit(0);
// 						// }

// 						//Insert the new seed into the seed list
// 						uArr[posArray[qProfile[pRank] + j].unitig_id].getData()->getData(uArr[posArray[qProfile[pRank] + j].unitig_id])->pushSeed(newSeed, isRefSeq);
// 					}
// 				}
// 			}

// 			//Testing
// 			// if(newSeed->offsetQ > 9999){
// 			// 	cout << "Something is wrong here" << endl << "isRefSeq:" << isRefSeq << endl << "i:" << i << "minSeedLength:" << minSeedLength << "uniLen:" << uniLen << endl;
// 			// 	exit(0);
// 			// }
// 		}
// 	}
// }

//This function checks if the search criteria are fulfilled for a region on a unitig starting at some offset and having a certain length
bool vfySrchCrit(UnitigColorMap<seedlist> unitig, const uint32_t &off, const int32_t &len, const uint32_t &quorum, const list<pair<string, size_t>> &srchColSet){
	//Check if search criteria have already been checked for this unitig and check if necessary
	if(!unitig.getData()->getData(unitig)->srchCritChckd()){
		//Testing
		// cout << "1 Option 2" << endl;

		calcSrchCritBrds(unitig, quorum, srchColSet);
	}// else{
	// 	//Testing
		// cout << "1 Option 1" << endl;
	// }

	//Check if search criteria are not fulfilled at all on this unitig
	if(unitig.getData()->getData(unitig)->getlBrd() < 0 && unitig.getData()->getData(unitig)->getrBrd() == (int32_t) unitig.size){
		//Testing
		// cout << "2" << endl;

		return false;
	}

	//Check if a seed falls completely behind the left or in front of the right border (rBrd is never negative)
	if((int32_t) off + len <= unitig.getData()->getData(unitig)->getlBrd() || off > (uint32_t) unitig.getData()->getData(unitig)->getrBrd()){
		//Testing
		// if(off + len <= (uint32_t) unitig.getData()->getData(unitig)->getlBrd()){
		// 	cout << "3 Option 1" << endl;
		// } else{
		// 	cout << "3 Option 2" << endl;
		// }

		return true;
	}

	//Testing
	// cout << "3 Option 3" << endl;
	// cout << "off: " << off << " len: " << len << " lBrd: " << unitig.getData()->getData(unitig)->getlBrd() << " rBrd: " << unitig.getData()->getData(unitig)->getrBrd() << endl;

	return false;
}

//This function performs the seed detection if a search color set is given
void detectSeeds(const int32_t &k, const int32_t &minSeedLength, const size_t &numSmers, const uint32_t &quorum, const uint32_t &profileSize, const uint32_t *qProfile, const string &q, const UnitigColorMap<seedlist> *uArr, const struct s_mer_pos *posArray, hit *hitArr, const list<pair<string, size_t>> &searchSet, const bool isRefSeq){
	uint32_t occnum;
	int32_t pRank = -1;
	size_t uniLen;
	struct seed *newSeed, *lastSeed;

	//Go through the query and check for each s-mer whether it appears in the graph
	for(uint32_t i = 0; i < q.length() - minSeedLength + 1; ++i){
		//Once we are iterating over q anyways initialize hitArr NOTE: Initialized hits have a length of 0//TODO: This shouldn't be in here!
		hitArr[i].length = 0;

		//Testing
		// cout << "Query: " << q << endl;
		// cout << "i: " << i << endl;

		//Calculate current s-mer's rank
		pRank = compRank(q.substr(i, minSeedLength), pRank);

		//To avoid border effects check whether we are at the end of the q-gram profile
		if(static_cast<unsigned>(pRank) < profileSize - 1 && qProfile[pRank] < numSmers){
			occnum = qProfile[pRank + 1] - qProfile[pRank];//if occnum is zero the current s-mer doesn't exist
		} else{
			occnum = numSmers - qProfile[pRank];//if occnum is zero the current s-mer doesn't exist
		}

		//Testing
		//cout << "Testcase: S-mer is " << (occnum != 0 ? "" : "not ") << "found" << endl;

		//Process all occurrences of a certain s-mer
		for(int32_t j = occnum - 1; j >= 0; --j){//for(uint32_t j = 0; j < occnum; ++j){
			//Get current unitig's length
			uniLen = uArr[posArray[qProfile[pRank] + j].unitig_id].size;

			//Testing
			// cout << "Seed at " << uArr[posArray[qProfile[pRank] + j].unitig_id].mappedSequenceToString() << " offset " << posArray[qProfile[pRank] + j].offset << " found" << endl;

			//Find out if unitig fulfills search criteria
			if(!vfySrchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], posArray[qProfile[pRank] + j].offset, minSeedLength, quorum, searchSet)){
				continue;
			}

			//Testing
			// cout << "Search criteria are fulfilled" << endl;

			//Check if the seed list is empty
			if(uArr[posArray[qProfile[pRank] + j].unitig_id].getData()->getData(uArr[posArray[qProfile[pRank] + j].unitig_id])->getSeed(isRefSeq) == NULL){
				// //Check whether our quorum is fulfilled move on otherwise
				// if(checkSearchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], quorum, posArray[qProfile[pRank] + j].offset, searchSet) == 0){

				// //Check if a seed falls completely behind the left or in front of the right border
				// if(posArray[qProfile[pRank] + j].offset + minSeedLength > lOff || posArray[qProfile[pRank] + j].offset <= rOff){
				// 	//Testing
				// 	//cout << "Testcase: S-mer has to be inserted as the first element in a seed list and quorum is not fulfilled" << endl;

				// 	continue;
				// }

				//Testing
				// cout << "Testcase: S-mer has to be inserted as the first element in a seed list and quorum is fulfilled" << endl;

				//Create a new seed and...
				newSeed = (struct seed*) malloc(sizeof(struct seed));
				newSeed->offsetU = compOffset(posArray[qProfile[pRank] + j].offset, minSeedLength, uniLen, isRefSeq);
				newSeed->offsetQ = compOffset(i, minSeedLength, q.length(), isRefSeq);
				newSeed->len = minSeedLength;
				newSeed->score = 0;
				newSeed->nextSeed = NULL;

				//...make it the first element of the seed list
				uArr[posArray[qProfile[pRank] + j].unitig_id].getData()->getData(uArr[posArray[qProfile[pRank] + j].unitig_id])->pushSeed(newSeed, isRefSeq);
			} else{
				//Check on which strand we are
				if(isRefSeq){
					//Testing
					// cout << "6 Option 1" << endl;

					//Get the last inserted seed
					lastSeed = uArr[posArray[qProfile[pRank] + j].unitig_id].getData()->getData(uArr[posArray[qProfile[pRank] + j].unitig_id])->getSeed(isRefSeq);

					//Check if we can extend the last seed inserted in our seed list
					if(lastSeed->offsetU + lastSeed->len - (minSeedLength - 1) == posArray[qProfile[pRank] + j].offset && lastSeed->offsetQ + lastSeed->len - (minSeedLength - 1) == i){
						//Testing
						// cout << "7 Option ";

						// //Check whether we have to check our quorum at the current position
						// if(lastSeed->len % k == 0){
						// 	//Testing
						// 	// cout << "not ";

						// 	//Check whether our quorum is fulfilled move on otherwise
						// 	if(checkSearchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], quorum, posArray[qProfile[pRank] + j].offset, searchSet) == 0){
						// 		//Testing
						// 		// cout << "2" << endl;
						// 		// cout << "8 Option 2" << endl;

						// 		continue;
						// 	}

						// 	//Testing							
						// 	// cout << "8 Option 1" << endl;
						// }

						//Testing
						// cout << "2" << endl;

						//Extending simply means to increment the seed's length
						++lastSeed->len;
					} else {
						//Testing
						// cout << "7 Option 3" << endl;

						// //Check whether our quorum is fulfilled move on otherwise
						// if(checkSearchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], quorum, posArray[qProfile[pRank] + j].offset, searchSet) == 0){

						// //Check if a seed falls completely behind the left or in front of the right border
						// if(posArray[qProfile[pRank] + j].offset + minSeedLength > lOff || posArray[qProfile[pRank] + j].offset <= rOff){
						// 	//Testing
						// 	// cout << "9 Option 2" << endl;

						// 	continue;
						// }

						//Testing
						// cout << "9 Option 1" << endl;

						//If an extension is not possible create a new seed
						newSeed = (struct seed*) malloc(sizeof(struct seed));
						newSeed->offsetU = posArray[qProfile[pRank] + j].offset;
						newSeed->offsetQ = i;
						newSeed->len = minSeedLength;
						newSeed->score = 0;
						newSeed->nextSeed = NULL;
						//Insert the new seed into the seed list
						uArr[posArray[qProfile[pRank] + j].unitig_id].getData()->getData(uArr[posArray[qProfile[pRank] + j].unitig_id])->pushSeed(newSeed, isRefSeq);
					}
				} else{
					//Testing
					// cout << "6 Option 2" << endl;

					//Get the last inserted seed
					lastSeed = uArr[posArray[qProfile[pRank] + j].unitig_id].getData()->getData(uArr[posArray[qProfile[pRank] + j].unitig_id])->getLastInSd();

					//Check if we can extend the last seed inserted seed in our seed list
					if(lastSeed->offsetU > 0 && lastSeed->offsetU - 1 == compOffset(posArray[qProfile[pRank] + j].offset, minSeedLength, uniLen, isRefSeq) && lastSeed->offsetQ > 0 && lastSeed->offsetQ - 1 == compOffset(i, minSeedLength, q.length(), isRefSeq)){
						//Testing
						// cout << "10 Option ";

						// //Check whether we have to check our quorum at the current position
						// if(lastSeed->len % k == 0){
						// 	//Testing
						// 	// cout << "not ";

						// 	//Check whether our quorum is fulfilled move on otherwise
						// 	if(checkSearchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], quorum, posArray[qProfile[pRank] + j].offset, searchSet) == 0){
						// 		//Testing
						// 		// cout << "2" << endl;
						// 		// cout << "11 Option 2" << endl;

						// 		continue;
						// 	}

						// 	//Testing
						// 	// cout << "11 Option 1" << endl;
						// }

						//Testing
						// cout << "2" << endl;

						//Extending means to increment the seed's length and to adjust the offsets
						++lastSeed->len;
						--lastSeed->offsetU;
						--lastSeed->offsetQ;
					} else{
						//Testing
						// cout << "10 Option 3" << endl;

						// //Check whether our quorum is fulfilled move on otherwise
						// if(checkSearchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], quorum, posArray[qProfile[pRank] + j].offset, searchSet) == 0){
						// 	//Testing
						// 	// cout << "12 Option 2" << endl;

						// 	continue;
						// }

						//Testing
						// cout << "12 Option 1" << endl;

						//If an extension is not possible create a new seed
						newSeed = (struct seed*) malloc(sizeof(struct seed));
						newSeed->offsetU = compOffset(posArray[qProfile[pRank] + j].offset, minSeedLength, uniLen, isRefSeq);
						newSeed->offsetQ = compOffset(i, minSeedLength, q.length(), isRefSeq);
						newSeed->len = minSeedLength;
						newSeed->score = 0;
						newSeed->nextSeed = NULL;
						//Insert the new seed into the seed list
						uArr[posArray[qProfile[pRank] + j].unitig_id].getData()->getData(uArr[posArray[qProfile[pRank] + j].unitig_id])->pushSeed(newSeed, isRefSeq);
					}
				}
			}
		}
	}
}

// //This function extends all seeds found on the queries reference strand considering a quorum
// void extendRefSeeds(ColoredCDBG<seedlist> &cdbg, const string &q, const int16_t &X, hit *hitArr, const uint32_t &quorum){
// 	struct seed *currSeed;

// 	//struct seed *delSeed = NULL;

// 	hit newHit;
// 	UnitigColorMap<seedlist> currUni;

// 	//Iterate over all seeds of all unitigs
// 	for(ColoredCDBG<seedlist>::iterator i = cdbg.begin(); i != cdbg.end(); ++i){
// 		//Get current unitig
// 		currUni = *i;
// 		//Get the first seed
// 		currSeed = currUni.getData()->getData(currUni)->getSeed(currUni.strand);

// 		//Iterate over all seeds of a unitig
// 		while(currSeed != NULL){
// 			//TODO: This procedure might be worth to be put into an external (inline) function!
// 			//Extract seed from the list of unextended seeds
// 			currUni.getData()->getData(currUni)->setSeed(currSeed->nextSeed, currUni.strand);
// 			//Setup initial hit infos
// 			newHit.score = 0;
// 			newHit.length = currSeed->len;
// 			newHit.offU = currSeed->offsetU;
// 			newHit.offQ = currSeed->offsetQ;
// 			newHit.origUni = currUni;
// 			newHit.nextHit = NULL;

// 			//Testing
// 			//cout << "Start right extension" << endl;
// 			//cout << "Current seed is uoff:" << currSeed->offsetU << " qoff:" << currSeed->offsetQ << " len:" << currSeed->len << endl;
// 			// if(currSeed->offsetU == 0 && currSeed->offsetQ == 870 && currSeed->len == 32){
// 			// 	report = true;
// 			// }

// 			//Extend hit to the right
// 			newHit.rExt = startRightX_Drop(&newHit, q, X, quorum);
// 			/*TODO Ideas to make the seed extension faster:
// 				1. It might be faster to delete seeds not before a maximum scoring path over many untigs has been calculated. Many seeds could be incorporated into different extensions.
// 				2. If s < k s-mers in the beginning of two unitigs that are both successors of a third one are the same. If an extension is already pretty bad when these unitigs are explored s.t. it drops while comparing base pairs of the identical s-mers it is sufficient to do this once and not for each unitig separately. (But the X-drop never drops while exploring a reached seed does it?)*/

// 			//Testing
// 			// if(newHit.rSeedUoff == 5 && newHit.rSeedQoff == 9999 && newHit.lSeedUoff == 0 && newHit.lSeedQoff == 4545){
// 			// cerr << "After right extension: length: " << newHit.length << " score: " << newHit.score << endl;
// 			// 	//cerr << "Left seed: uoff:" << newHit.lSeedUoff << " qoff: " << newHit.lSeedQoff << endl;
// 			// }

// 			//Filter out some seeds; the second condition ensures that we do not miss seeds in the end of the query
// 			if(newHit.length - currSeed->len > 0 || currSeed->offsetQ + currSeed->len == q.length()){
// 				//Extend hit to the left
// 				newHit.lExt = startLeftX_Drop(&newHit, q, X, quorum);

// 				//Check whether we have created a hit for this query position already
// 				if(hitArr[newHit.offQ].length == 0){
// 					//Copy hit into the array
// 					hitArr[newHit.offQ] = newHit;

// 					//Use hit space in the array
// 					// hitArr[newHit.lSeedQoff].score = newHit.score;
// 					// hitArr[newHit.lSeedQoff].length = newHit.length;
// 					// hitArr[newHit.lSeedQoff].lSeedUoff = newHit.lSeedUoff;
// 					// hitArr[newHit.lSeedQoff].lSeedQoff = newHit.lSeedQoff;
// 					// hitArr[newHit.lSeedQoff].lUnitig = newHit.lUnitig;
// 					// hitArr[newHit.lSeedQoff].rSeedUoff = newHit.rSeedUoff;
// 					// hitArr[newHit.lSeedQoff].rSeedQoff = newHit.rSeedQoff;
// 					// hitArr[newHit.lSeedQoff].rUnitig = newHit.rUnitig;
// 					// hitArr[newHit.lSeedQoff].sdOffU = newHit.sdOffU;
// 					// hitArr[newHit.lSeedQoff].sdOffQ = newHit.sdOffQ;
// 					// hitArr[newHit.lSeedQoff].origUni = newHit.origUni;
// 					// hitArr[newHit.lSeedQoff].rExt = newHit.rExt;
// 					// hitArr[newHit.lSeedQoff].lExt = newHit.lExt;
// 					// hitArr[newHit.lSeedQoff].nextHit = newHit.nextHit;
// 				} else{
// 					//Save pointer to the hit list at the current q position
// 					newHit.nextHit = hitArr[newHit.offQ].nextHit;
// 					//Create a new hit
// 					hitArr[newHit.offQ].nextHit = new hit(newHit);
					
// 					//Copy information
// 					// hitArr[newHit.lSeedQoff].nextHit->score = newHit.score;
// 					// hitArr[newHit.lSeedQoff].nextHit->length = newHit.length;
// 					// hitArr[newHit.lSeedQoff].nextHit->lSeedUoff = newHit.lSeedUoff;
// 					// hitArr[newHit.lSeedQoff].nextHit->lSeedQoff = newHit.lSeedQoff;
// 					// hitArr[newHit.lSeedQoff].nextHit->lUnitig = newHit.lUnitig;
// 					// hitArr[newHit.lSeedQoff].nextHit->rSeedUoff = newHit.rSeedUoff;
// 					// hitArr[newHit.lSeedQoff].nextHit->rSeedQoff = newHit.rSeedQoff;
// 					// hitArr[newHit.lSeedQoff].nextHit->rUnitig = newHit.rUnitig;
// 					// hitArr[newHit.lSeedQoff].nextHit->sdOffU = newHit.sdOffU;
// 					// hitArr[newHit.lSeedQoff].nextHit->sdOffQ = newHit.sdOffQ;
// 					// hitArr[newHit.lSeedQoff].nextHit->origUni = newHit.origUni;
// 					// hitArr[newHit.lSeedQoff].nextHit->rExt = newHit.rExt;
// 					// hitArr[newHit.lSeedQoff].nextHit->lExt = newHit.lExt;

// 					//Testing
// 					// cout << "We get until here" << endl;

// 					hitArr[newHit.offQ].nextHit->gAlgn.aSeqG = "";//TODO Check if this is necessary!
// 					hitArr[newHit.offQ].nextHit->gAlgn.aSeqQ = "";

// 					//Testing
// 					// cout << "But not here" << endl;

// 					// hitArr[newHit.lSeedQoff].nextHit->nextHit = newHit.nextHit;
// 				}
// 			}

// 			//Delete the current seed
// 			free(currSeed);
// 			//Move on to the next seed
// 			currSeed = currUni.getData()->getData(currUni)->getSeed(currUni.strand);
// 		}
// 	}
// }

// //This function extends all seeds found on the queries reverse complement considering a quorum
// void extendRevCompSeeds(ColoredCDBG<seedlist> &cdbg, const string &q, const int16_t &X, hit *hitArr, const uint32_t &quorum){
// 	bool isValid = true;
// 	struct seed *currSeed;
// 	hit newHit;
// 	hit *hitIt;
// 	UnitigColorMap<seedlist> currUni;

// 	//Iterate over all seeds of all unitigs
// 	for(ColoredCDBG<seedlist>::iterator i = cdbg.begin(); i != cdbg.end(); ++i){
// 		//Get current unitig
// 		currUni = *i;
// 		//We are on the reverse complementary sequence
// 		currUni.strand = false;
// 		//Get the first seed
// 		currSeed = currUni.getData()->getData(currUni)->getSeed(currUni.strand);

// 		//Testing
// 		//if(currSeed != NULL) cout << "Extending seed: upos:" << currSeed->offsetU << " qpos:" << currSeed->offsetQ << endl;

// 		//Testing
// 		//cout << "Current unitig is " << i->mappedSequenceToString() << endl;

// 		//Iterate over all seeds of a unitig
// 		while(currSeed != NULL){
// 			//Testing
// 			//cout << "Current seed is u offset:" << currSeed->offsetU << " q offset:" << currSeed->offsetQ << " len:" << currSeed->len << " score:" << currSeed->score << endl;

// 			//TODO: This procedure might be worth to be put into an external (inline) function!
// 			//Extract seed from the list of unextended seeds
// 			currUni.getData()->getData(currUni)->setSeed(currSeed->nextSeed, currUni.strand);
// 			//Setup initial hit infos
// 			newHit.score = 0;
// 			newHit.length = currSeed->len;
// 			newHit.offU = currSeed->offsetU;
// 			newHit.offQ = currSeed->offsetQ;
// 			newHit.origUni = *i;
// 			newHit.nextHit = NULL;

// 			//Testing
// 			//cout << "Next we extend seed: upos:" << currSeed->offsetU << " qpos:" << currSeed->offsetQ << endl;	
// 			// if(currSeed->offsetQ == 9){
// 			// 	cout << "Extending our seed of interest" << endl;
// 			// }// else{
// 			// 	cout << "There is none" << endl;
// 			// }

// 			//Extend hit to the right
// 			startRightX_Drop_OnRevComp(&newHit, q, X, quorum);
// 			/*TODO Ideas to make the seed extension faster:
// 				1. It might be faster to delete seeds not before a maximum scoring path over many untigs has been calculated. Many seeds could be incorporated into different extensions.
// 				2. If s < k s-mers in the beginning of two unitigs that are both successors of a third one are the same. If an extension is already pretty bad when these unitigs are explored s.t. it drops while comparing base pairs of the identical s-mers it is sufficient to do this once and not for each unitig separately. (But the X-drop never drops while exploring a reached seed does it?)*/

// 			//Testing
// 			// if(currSeed->offsetQ == 9){
// 			// 	cout << "After right extension the hit's right border in q is " << newHit.rSeedQoff << endl;
// 			// }

// 			//Filter out some seeds; the second condition ensures that we do not miss seeds in the end of the query
// 			if(newHit.length - currSeed->len > 0 || currSeed->offsetQ + currSeed->len == q.length()){
// 				//Extend hit to the left
// 				startLeftX_Drop_OnRevComp(&newHit, q, X, quorum);

// 				//Check whether we have created a hit for this query position already
// 				if(hitArr[newHit.offQ].length != 0){
// 					//Get the first hit at this position
// 					hitIt = &hitArr[newHit.offQ];
// 					//Go through the list of all hits at this position
// 					while(hitIt != NULL){
// 						//Check if the current hit and the last extended hit might be cooptimal solutions (i.e. they have the same q offsets and score)
// 						if(hitIt->offQ == newHit.offQ && hitIt->score == newHit.score){
// 							//Mark hit as not interesting
// 							isValid = false;
// 							break;
// 						}

// 						//Move on to the next hit
// 						hitIt = hitIt->nextHit;
// 					}

// 					//Check if extended hit is worth to be kept
// 					if(!isValid){
// 						//Reset flag
// 						isValid = true;
// 					} else{
// 						//Add hit to the hit array//

// 						//Save pointer to the hit list at the current q position
// 						newHit.nextHit = hitArr[newHit.offQ].nextHit;
// 						//Create a new hit
// 						hitArr[newHit.offQ].nextHit = new hit(newHit);
// 					}
// 				} else{
// 					//Copy hit into the array
// 					hitArr[newHit.offQ] = newHit;
// 				}
// 			}
// 			//Testing
// 			// if(newHit.score == 4294967295){
// 			// 	cerr << "Seed was uoff:" << currSeed->offsetU << " qoff:" << currSeed->offsetQ << " len:" << currSeed->len << endl;
// 			// 	exit(0);
// 			// }
			
// 			//Delete the current seed
// 			free(currSeed);
// 			//Move on to the next seed
// 			currSeed = currUni.getData()->getData(currUni)->getSeed(currUni.strand);
// 		}
// 	}
// }

//This function extends all seeds found on the queries reference strand considering a quorum and a search color set
void extendRefSeeds(ColoredCDBG<seedlist> &cdbg, const string &q, const int16_t &X, hit *hitArr, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	struct seed *currSeed;
	hit newHit;
	UnitigColorMap<seedlist> currUni;

	//Testing
	// cout << "Extension on reference strand" << endl;

	//Iterate over all seeds of all unitigs
	for(ColoredCDBG<seedlist>::iterator i = cdbg.begin(); i != cdbg.end(); ++i){
		//Get current unitig
		currUni = *i;
		//Get the first seed
		currSeed = currUni.getData()->getData(currUni)->getSeed(currUni.strand);

		//Iterate over all seeds of a unitig
		while(currSeed != NULL){
			//TODO: This procedure might be worth to be put into an external (inline) function!
			//Extract seed from the list
			currUni.getData()->getData(currUni)->setSeed(currSeed->nextSeed, currUni.strand);
			//Setup initial hit infos
			newHit.score = 0;
			newHit.length = currSeed->len;
			newHit.offU = currSeed->offsetU;
			newHit.offQ = currSeed->offsetQ;
			newHit.origUni = currUni;
			newHit.nextHit = NULL;

			//Extend hit to the right
			startRightX_Drop(&newHit, q, X, quorum, searchSet);
			/*TODO Ideas to make the seed extension faster:
				1. It might be faster to delete seeds not before a maximum scoring path over many untigs has been calculated. Many seeds could be incorporated into different extensions.
				2. If s < k s-mers in the beginning of two unitigs that are both successors of a third one are the same. If an extension is already pretty bad when theses unitigs are explored s.t. it drops while comparing base pairs of the identical s-mers it is sufficient to do this once and not for each unitig separately. (But the X-drop never drops while exploring a reached seed does it?)*/			
			
			//Testing
			// if(currUni.mappedSequenceToString() == "TCTTTACGGCGAAGTTCAGCGCCCTCATAGCC" && currSeed->offsetU == 0 && currSeed->offsetQ == 19){
			// 	// cout << "Seed of the unitig:" << endl << "offsetU: " << currSeed->offsetU << " offsetQ: " << currSeed->offsetQ << " score: " << currSeed->score << " length: " << currSeed->len << endl;
			// 	cout << "extendRefSeeds: After right extension: offU: " << newHit.offU << " offQ: " << newHit.offQ << " score: " << newHit.score << " length: " << newHit.length << endl;
			// 	exit(0);
			// }
			// cout << "We do get here" << endl;
			// cout << "extendRefSeeds: Processing seed offsetU: " << currSeed->offsetU << " offsetQ: " << currSeed->offsetQ << " on unitig " << currUni.mappedSequenceToString() << endl; 

			//Filter out some seeds; the second condition ensures that we do not miss seeds in the end of the query
			//Note: What we do not consider here is that some seeds might not be extended to the right because search criteria are not fullfilled anymore. This is intended though. We should not miss too much, because a good hit should have more than one seed
			if(newHit.length - currSeed->len > 0 || currSeed->offsetQ + currSeed->len == q.length()){
				//Testing
				// cout << "We get here" << endl;

				//Extend hit to the left
				startLeftX_Drop(&newHit, q, X, quorum, searchSet);

				//Testing
				// if(newHit.origUni.mappedSequenceToString() == "TTTCACCGCATCGGCAATGGCGGGCAGGGCG"){
				// 	cerr << "extendRefSeeds: We have found the interesting unitig" << endl;
				// 	exit(0);
				// }

				//Check whether we have created a hit for this query position already
				if(hitArr[newHit.offQ].length == 0){
					//Copy hit into the array
					hitArr[newHit.offQ] = newHit;
				} else{
					//Testing
					//cout << "4 Option 2" << endl;

					//Save pointer to the hit list at the current q position
					newHit.nextHit = hitArr[newHit.offQ].nextHit;
					//Create a new hit
					hitArr[newHit.offQ].nextHit = new hit(newHit);

					hitArr[newHit.offQ].nextHit->gAlgn.aSeqG = "";//TODO Check if this is necessary!
					hitArr[newHit.offQ].nextHit->gAlgn.aSeqQ = "";
				}
			}

			//Delete the extended seed
			free(currSeed);
			//Move on to the next seed
			currSeed = currUni.getData()->getData(currUni)->getSeed(currUni.strand);
		}

		//Testing
		// if(currUni.mappedSequenceToString() == "TCTTTACGGCGAAGTTCAGCGCCCTCATAGCC"){
		// 	exit(0);
		// }
	}
}

//This function extends all seeds found on the queries reverse complement considering a quorum and a search color set
//Note: What we do not consider here is that some seeds might not be extended to the right because search criteria are not fullfilled anymore. This is intended though. We should not miss too much, because a good hit should have more than one seed
void extendRevCompSeeds(ColoredCDBG<seedlist> &cdbg, const string &q, const int16_t &X, hit *hitArr, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	bool isValid = true;
	struct seed *currSeed;
	hit newHit;
	hit *hitIt;
	UnitigColorMap<seedlist> currUni;

	//Testing
	// cout << "We are at the beginning of extendRevCompSeeds" << endl;

	//Iterate over all seeds of all unitigs
	for(ColoredCDBG<seedlist>::iterator i = cdbg.begin(); i != cdbg.end(); ++i){
		//Get current unitig
		currUni = *i;
		//We are on the reverse complementary sequence
		currUni.strand = false;
		//Get the first seed
		currSeed = currUni.getData()->getData(currUni)->getSeed(currUni.strand);

		//Testing
		// uint32_t curLen = 0, curOffU = 0, curOffQ = 0;
		// cout << "Process the first unitig" << endl;
		// cout << "Hit array: " << endl;
		// for(uint16_t j = 0; j < q.length(); ++j){
		// 	cout << hitArr[j].length << endl;
		// }

		//Iterate over all seeds of a unitig
		while(currSeed != NULL){
			//TODO: This procedure might be worth to be put into an external (inline) function!
			//Extract seed from the list
			currUni.getData()->getData(currUni)->setSeed(currSeed->nextSeed, currUni.strand);
			//Setup initial hit infos
			newHit.score = 0;
			newHit.length = currSeed->len;
			newHit.offU = currSeed->offsetU;
			newHit.offQ = currSeed->offsetQ;
			newHit.origUni = currUni;
			newHit.nextHit = NULL;

			//Testing
			if(newHit.origUni.mappedSequenceToString() == "CAGTACGGTATCGGCCCCCAACGCAATCATGCGCACAACGTCCAGACCGTTACGGATCCCGCTGTCTGCCAGAATGGTGATGTCGCCTTTCACCGCATCGGCAATGGCGGG"){
				cerr << "Initial seed is offsetU: " << currSeed->offsetU << " offsetQ: " << currSeed->offsetQ << " len: " << currSeed->len << endl;
				// exit(0);
			}

			//Extend hit to the right
			startRightX_Drop_OnRevComp(&newHit, q, X, quorum, searchSet);
			/*TODO Ideas to make the seed extension faster:
				1. It might be faster to delete seeds not before a maximum scoring path over many untigs has been calculated. Many seeds could be incorporated into different extensions.
				2. If s < k s-mers in the beginning of two unitigs that are both successors of a third one are the same. If an extension is already pretty bad when theses unitigs are explored s.t. it drops while comparing base pairs of the identical s-mers it is sufficient to do this once and not for each unitig separately. (But the X-drop never drops while exploring a reached seed does it?)*/			

			//Testing	
			// if(currUni.mappedSequenceToString() == "TTCCTCACAGAATTCAATTCGATATGGTTACACACGGCGTTTGCTCTGAGTTATC"){
			// 	cout << "Seed of the unitig:" << endl << "offsetU: " << currSeed->offsetU << " offsetQ: " << currSeed->offsetQ << " score: " << currSeed->score << " length: " << currSeed->len << endl;
			// 	cout << "extendRevCompSeeds: After right extension: offU: " << newHit.offU << " offQ: " << newHit.offQ << " score: " << newHit.score << " length: " << newHit.length << endl;
			// 	// exit(0);
			// }
			// curLen = currSeed->len;
			// curOffU = currSeed->offsetU;
			// curOffQ = currSeed->offsetQ;
			// cout << "Seed is offsetU: " << currSeed->offsetU << " offsetQ: " << currSeed->offsetQ << " on unitig " << currUni.mappedSequenceToString() << endl;
			// cout << "Seed is offsetU: " << currSeed->offsetU << " offsetQ: " << currSeed->offsetQ << endl;
			// if(currSeed->offsetU == 0 && currSeed->offsetQ == 19){
			// 	cout << "Seed of interest extended" << endl << " score: " << newHit.score << " length: " << newHit.length << " on unitig " << newHit.origUni.mappedSequenceToString() << endl;
			// 	exit(0);
			// }
			// cout << "Right extension done" << endl;
			// cout << "Hit array: " << endl;
			// for(uint16_t j = 0; j < q.length(); ++j){
			// 	cout << hitArr[j].length << endl;
			// }

			//Filter out some seeds; the second condition ensures that we do not miss seeds in the end of the query and the third that we do not miss seeds in the end of a unitig if the unitig does not have successors (i.e. predecessors on the reference strand)
			if(newHit.length - currSeed->len > 0 || currSeed->offsetQ + currSeed->len == q.length() || (!i->getPredecessors().hasPredecessors() && currSeed->offsetU + currSeed->len == i->size)){
				//Testing
				// if(newHit.length - currSeed->len > 0){
				// 	cout << "1 Option 1" << endl;	
				// }
				// if(currSeed->offsetQ + currSeed->len == q.length()){
				// 	// cout << "2 Option 1" << endl;
				// }
				// if(!i->getPredecessors().hasPredecessors() && currSeed->offsetU + currSeed->len == i->size){
				// 	cout << "3 Option 1" << endl;
				// }
				// cout << "We passed the filter" << endl;

				//Extend hit to the left
				startLeftX_Drop_OnRevComp(&newHit, q, X, quorum, searchSet);

				//Testing
				// cout << "Left extension done" << endl;TCTTTACGGCGAAGTTCAGCGCCCTCATAGCC
				// if(newHit.origUni.mappedSequenceToString() == "TTTGTTTTCAATTGCTGATGAATGGGGTATGAGTAAACTGAG"){
				// 	// cout << "Seed of the unitig:" << endl << "offsetU: " << currSeed->offsetU << " offsetQ: " << currSeed->offsetQ << " score: " << currSeed->score << " length: " << currSeed->len << endl;
				// 	// cout << "After left extension: offU: " << newHit.offU << " offQ: " << newHit.offQ << " score: " << newHit.score << " length: " << newHit.length << endl;
				// 	// cout << "origUni: " << currUni.mappedSequenceToString() << endl;
				// 	cout << "extendRevCompSeeds: Unitig found" << endl;
				// 	exit(0);
				// }
				// cout << "Left extension done" << endl;

				// //For now we skip seeds ending within the overlap at the end of a unitig sequence
				// if(currUni.size - newHit.offU < (uint32_t) cdbg.getK()){
				// 	//Delete the extended seed
				// 	free(currSeed);
				// 	//Move on to the next seed
				// 	currSeed = currUni.getData()->getData(currUni)->getSeed(currUni.strand);
				// 	continue;
				// }

				//Check whether we have created a hit for this query position already
				if(hitArr[newHit.offQ].length != 0){
					//Testing
					// cout << "4 Option 1" << endl;

					//Get the first hit at this position
					hitIt = &hitArr[newHit.offQ];
					//Go through the list of all hits at this position
					while(hitIt != NULL){
						//Check if the current hit and the last extended hit might be cooptimal solutions (i.e. they have the same q offsets and score)
						if(hitIt->length == newHit.length && hitIt->score == newHit.score){
							//Testing
							// cout << "This is triggered" << endl;

							//Mark hit as not interesting
							isValid = false;
							break;
						}

						//Move on to the next hit
						hitIt = hitIt->nextHit;
					}

					//Check if extended hit is worth to be kept
					if(!isValid){
						//Reset flag
						isValid = true;
					} else{
						//Testing
						// cout << "Make valid borders" << endl;

						//Adjust hit borders if necessary//TODO This should be done not before we select the best extended hits for the gapped alignment, but for the moment this is good enough. Apart from that we can just assume that a gapped extension will not be possible to extend our hits much further so that we can just take it out here and consider it later!
						//makeValidBorders(&newHit, q);

						//Testing
						// cout << "Valid borders made" << endl;

						//Add hit to the hit array//

						//Save pointer to the hit list at the current q position
						newHit.nextHit = hitArr[newHit.offQ].nextHit;
						//Create a new hit
						hitArr[newHit.offQ].nextHit = new hit(newHit);
						hitArr[newHit.offQ].nextHit->gAlgn.aSeqG = "";//TODO Check if this is necessary!
						hitArr[newHit.offQ].nextHit->gAlgn.aSeqQ = "";
					}
				} else{
					//Testing
					// cout << "4 Option 2" << endl;
					// cout << "Hit array: " << endl;
					// for(uint16_t j = 0; j < q.length(); ++j){
					// 	cout << hitArr[j].length << endl;
					// }
					// cout << "Make valid borders" << endl;

					//Adjust hit borders if necessary//TODO This should be done not before we select the best extended hits for the gapped alignment, but for the moment this is good enough. Apart from that we can just assume that a gapped extension will not be possible to extend our hits much further so that we can just take it out here and consider it later!
					//makeValidBorders(&newHit, q);

					//Copy hit into the array
					hitArr[newHit.offQ] = newHit;

					//Testing
					// cout << "Valid borders made" << endl;
				}
			}// else{
			// 	//Testing
			// 	// cout << "No left extension for this unitig" << endl;
			// }

			//Delete the extended seed
			free(currSeed);
			//Move on to the next seed
			currSeed = currUni.getData()->getData(currUni)->getSeed(currUni.strand);

			//Testing
			// cout << "Do we get here?" << endl;
		}

		//Testing (This is just a sanity test to make sure that there are no processed seeds which are saved anymore)
		// if(currUni.mappedSequenceToString() == "TCTTTACGGCGAAGTTCAGCGCCCTCATAGCC"){
		// 	exit(0);
		// }
		// if(newHit.length - curLen <= 0){
		// 	cout << "1 Option 2" << endl;
		// }
		// if(curOffQ + curLen != q.length()){
		// 	cout << "2 Option 2" << endl;
		// }
		// if(i->getSuccessors().hasSuccessors() || curOffU + curLen != i->size){
		// 	cout << "3 Option 2" << endl;
		// }
		// if(i->getData()->getData(*i)->getLastProcSeed() != NULL) cerr << "There is a list of processed seeds which is not empty!" << endl;
		// cout << "This is the end of the current iteration" << endl;
		// cout << "Hit array: " << endl;
		// for(uint16_t j = 0; j < q.length(); ++j){
		// 	cout << hitArr[j].length << endl;
		// }
	}
}

// //This function calculates a banded, semi-global, gapped alignment on a list of results considering a quorum and outputs the result if demanded
// void calcGappedAlignment(const list<hit*> &resList, const string &q, const int16_t &X, const bool &noOutput, const uint32_t &quorum){
// 	bool isDupl = false;
// 	uint32_t bandRadius;
// 	list<hit*>::const_iterator iter;

// 	//Testing
// 	// uint32_t reslistSize = 0;
// 	// uint32_t nb_duplicates = 0;

// 	//Go through result list
// 	for(list<hit*>::const_iterator it = resList.begin(); it != resList.end(); ++it){
// 		//Calculate the band width to be used during the gapped extension
// 		bandRadius = (*it)->length / GAP_RATIO;

// 		//Calculate gapped extension to the right
// 		startRightGappedAlignment(*it, q, X, bandRadius, quorum);

// 		//Testing
// 		// if(counter < 100){
// 		// 	cout << "Result of right extension: rSeed uoff:" << (*it)->rSeedUoff << " rSeed qoff:" << (*it)->rSeedQoff << " rUnitig:" << (*it)->rUnitig.mappedSequenceToString() << endl;
// 		// 	++counter;
// 		// }

// 		//Check whether this is not the first hit in the result list
// 		if(it != resList.begin()){
// 			iter = resList.begin();

// 			//Check for duplicates //TODO: If we filter for duplicates we could keep a statistic on how many duplicates of a specific result have been filtered out which could then be made a part of the result to inform the user about it!
// 			while(it != iter){
// 				//Testing
// 				//if(counter < 100) cout << "it != iter" << endl;

// 				//We assume to have a duplicate if right border offsets and unitigs are identical
// 				if((*it)->offU == (*iter)->offU && (*it)->offQ == (*iter)->offQ && (*it)->origUni == (*it)->origUni && (*iter)->score == (*it)->score && (*iter)->length == (*it)->length){
// 					//cerr << "It seems that the hit going spanning from q=" << (*iter)->lSeedQoff << " to q=" << (*iter)->rSeedQoff << " is a duplicate" << endl;
// 					isDupl = true;
// 					//++nb_duplicates;
// 					break;
// 				}

// 				++iter;
// 			}
// 		}

// 		//Testing
// 		// cout << "After right gapped alignment: Alignment: query: " << (*it)->gAlgn.aSeqQ << " graph: " << (*it)->gAlgn.aSeqG << endl;
// 		// cout << "Starting left gapped extension" << endl;

// 		//We want to avoid doing left extensions for duplicates
// 		if(!isDupl){
// 			startLeftGappedAlignment(*it, q, X, bandRadius, quorum);
// 		}

// 		//Testing
// 		// cout << "After left gapped alignment: Alignment: query: " << (*it)->gAlgn.aSeqQ << " graph: " << (*it)->gAlgn.aSeqG << endl;

// 		//Output result
// 		if(!noOutput && !isDupl){

// 			// cout << "Result info: score: " << (*it)->score << " left seed: offset q: " << (*it)->lSeedQoff << " offset u: " << (*it)->lSeedUoff << " unitig: " << (*it)->lUnitig.mappedSequenceToString() << " right seed: offset q: " << (*it)->rSeedQoff << " offset u: " << (*it)->rSeedUoff << " unitig: " << (*it)->rUnitig.mappedSequenceToString() << endl;

// 			//Output alignment
// 			repAlgn(*it);
// 		}

// 		//Reset duplicate flag
// 		isDupl = false;
// 	}

// 	//Testing
// 	//cerr << "Result list length:" << reslistSize << endl << "Number of duplicates:" << nb_duplicates << endl;
// }

//This function calculates a banded, semi-global, gapped alignment on a list of results considering a quorum and a search color set and outputs the result if demanded
void calcGappedAlignment(ColoredCDBG<seedlist> &cdbg, list<hit*> &resList, const string &q, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const double &lambda, const double &C){
	bool isDupl = false;
	uint32_t bandRadius;
	list<hit*>::const_iterator iter;

	//Testing
	// uint16_t hitCount = 0;
	// bool duplFound = false;
	// cout << "Do we get here?" << endl;

	//Go through result list
	for(list<hit*>::const_iterator it = resList.begin(); it != resList.end(); ++it){
		//Testing
		// cout << "Let's check out the next result" << endl;

		//Calculate the band width to be used during the gapped extension
		bandRadius = (*it)->length / GAP_RATIO;

		//Testing
		// if((*it)->origUni.mappedSequenceToString() == "GCAAGAGCCGCTGTTTCTTGAACAATATCTCG"){
		// 	cout << "Before right GappedAlignment: (*it)->offU: " << (*it)->offU << " (*it)->offQ: " << (*it)->offQ << endl << "aSeqG: " << (*it)->gAlgn.aSeqG << endl << "aSeqQ: " << (*it)->gAlgn.aSeqQ << endl;
		// 	report = true;
		// }

		//Calculate gapped extension to the right
		startRightGappedAlignment(*it, q, X, bandRadius, quorum, searchSet);

		//Testing
		// if((*it)->origUni.mappedSequenceToString() == "GATTTTGGCTTCACGTTTAAAAAGCAGCGAAA"){
		// 	// cout << "Seed of the unitig:" << endl << "offsetU: " << (*it)->offU << " offsetQ: " << (*it)->offQ << " score: " << (*it)->score << " length: " << (*it)->length << endl;
		// 	cout << "After right gapped extension: aSeqG: " << (*it)->gAlgn.aSeqG << endl << "aSeqQ: " << (*it)->gAlgn.aSeqQ << endl;
		// 	report = false;
		// }
		// if((*it)->origUni.mappedSequenceToString() == "TCTTTACGGCGAAGTTCAGCGCCCTCATAGCC"){
		// 	// cout << "Seed of the unitig:" << endl << "offsetU: " << currSeed->offsetU << " offsetQ: " << currSeed->offsetQ << " score: " << currSeed->score << " length: " << currSeed->len << endl;
		// 	cout << "After right GappedAlignment: offU: " << (*it)->offU << " offQ: " << (*it)->offQ << endl << "aSeqG: " << (*it)->gAlgn.aSeqG << endl << "aSeqQ: " << (*it)->gAlgn.aSeqQ << endl;
		// 	report = true;
		// 	// cout << "After right GappedAlignment: (*it)->offU: " << (*it)->offU << " (*it)->offQ: " << (*it)->offQ << endl;
		// 	// exit(0);
		// }
		// cout << "Completed right gapped extension" << endl;

		//Check whether this is not the first hit in the result list
		if(it != resList.begin()){
			iter = resList.begin();

			//Check for duplicates //TODO: If we filter for duplicates we could keep a statistic on how many duplicates of a specific result have been filtered out which could then be made a part of the result to inform the user about it!
			while(it != iter){
				//We assume to have a duplicate if right border offsets and unitigs are identical
				if((*it)->offU == (*iter)->offU && (*it)->offQ == (*iter)->offQ && (*it)->origUni == (*it)->origUni && (*iter)->score == (*it)->score && (*iter)->length == (*it)->length){
					// cerr << "It seems that the hit going spanning from q=" << (*iter)->lSeedQoff << " to q=" << (*iter)->rSeedQoff << " is a duplicate" << endl;

					//Testing
					// cout << "Duplicate detected" << endl;
					
					isDupl = true;
					break;
				}

				++iter;
			}
		}

		//Testing
		// cout << "origUni: " << (*it)->origUni.mappedSequenceToString() << endl;

		//We want to avoid doing left extensions for duplicates
		if(!isDupl){
			//Testing
			// cout << "Calculate left gapped alignment" << endl;

			startLeftGappedAlignment(*it, q, X, bandRadius, quorum, searchSet);
		} else{
			//Testing
			// cout << "2 Option 1" << endl;
			// duplFound = true;

			//Delete result
			resList.erase(it);
		}

		//Recalculate e-value
		(*it)->eval = calcEVal((*it)->score, lambda, C, q.length());

		//Testing
		// if((*it)->offQ == 2878){
		// 	cout << "Initial unitig was " << (*it)->origUni.mappedSequenceToString() << endl;
		// 	// cout << "Seed of the unitig:" << endl << "offsetU: " << currSeed->offsetU << " offsetQ: " << currSeed->offsetQ << " score: " << currSeed->score << " length: " << currSeed->len << endl;
		// 	// cout << "After left GappedAlignment: score: " << (*it)->score << endl << "aSeqG: " << (*it)->gAlgn.aSeqG << endl << "aSeqQ: " << (*it)->gAlgn.aSeqQ << endl;
		// 	// cout << "After right GappedAlignment: (*it)->offU: " << (*it)->offU << " (*it)->offQ: " << (*it)->offQ << endl;
		// 	// exit(0);
		// }
		// if((*it)->origUni.mappedSequenceToString() == "TTTGTTTTCAATTGCTGATGAATGGGGTATGAGTAAACTGAG"){
		// 	cout << "We have found the interesting unitig" << endl;
		// 	exit(0);
		// }

		//Testing
		// cout << "Inside calcGappedAlignment: Calculations done" << endl;

		//Output result
		// if(!noOutput && !isDupl){
		// 	// cout << "Result info: score: " << (*it)->score << " left seed: offset q: " << (*it)->offQ << " offset u: " << (*it)->offU << " unitig: " << (*it)->origUni.mappedSequenceToString() << " right seed: offset q: " << (*it)->offQ << " offset u: " << (*it)->offU << " unitig: " << (*it)->origUni.mappedSequenceToString() << endl;

		// 	//Testing
		// 	// cout << "Before repAlgn: (*it)->offU: " << (*it)->offU << " (*it)->offQ: " << (*it)->offQ << endl;

		// 	//Output alignment
		// 	repAlgn(*it);

		// 	//Testing
		// 	// cout << "We get here" << endl;
		// 	// cout << "Before outpColSets: (*it)->offU: " << (*it)->offU << " (*it)->offQ: " << (*it)->offQ << endl;

		// 	//Output color sets
		// 	// outpColSets(cdbg, *it);

		// 	// cout << "We get even out again" << endl;
		// }

		//Reset duplicate flag
		isDupl = false;
	}

	//Testing
	// cout << "We are done with calcGappedAlignment" << endl;
}

//This function performs the actual graph search for a query
void searchQuery(ColoredCDBG<seedlist> &cdbg, const int32_t &kMerLength, const int32_t &minSeedLength, const size_t &numSmers, const uint32_t &quorum, const uint32_t &profileSize, const uint32_t *qProfile, const string &q, const SrchStrd &strand, const UnitigColorMap<seedlist> *uArr, const struct s_mer_pos *posArray, const list<pair<string, size_t>> &searchColors, const int16_t &X, const bool &calcRT, uint16_t nRes, const double &lambda, const double &C, const double &eLim){
	//Staff we need to measure run times
	auto startTime = std::chrono::system_clock::now();
	auto endTime = std::chrono::system_clock::now();
	std::chrono::duration<double> tDiff;

	//No hit can start behind this point
	uint32_t hitArrSize = q.length() - minSeedLength + 1;
	//The query's reverse complement
	string revQ;
	
	//Query counter
	//uint32_t qCounter = 0;

	//Output which query we are working on
	//cout << "Query " << ++qCounter << endl;

	//Pointer needed to free hit array
	// hit *h, *delHit;
	//Variables needed for the seed extension
	hit hitArr[hitArrSize];
	//Variables needed for the gapped extension
	list<hit*> resList;
	list<hit*>::iterator iter;

	/*Searching for seeds*/
	cout << "Searching for seeds" << endl;

	// if(searchColors.empty()){//TODO It is not nice to copy all is! Is it maybe possible to circumvent this by giving a function as parameter or does it actually not even save time?
	// 	detectSeeds(kMerLength, minSeedLength, numSmers, quorum, profileSize, qProfile, q, uArr, posArray, hitArr);
	// } else{

	//Testing
	// if(strand == Plus){
	// 	cout << "1 Option 1" << endl;
	// }
	// if(strand == Minus){
	// 	cout << "1 Option 2" << endl;
	// }
	// if(strand == Both){
	// 	cout << "1 Option 3" << endl;
	// }

	//Detect seeds on the original query
	if(strand != Minus) detectSeeds(kMerLength, minSeedLength, numSmers, quorum, profileSize, qProfile, q, uArr, posArray, hitArr, searchColors, true);

	//Detect seeds on the reverse complementary query
	if(strand != Plus){
		//Calculate the query's reverse complement
		revQ = revComp(q);

		detectSeeds(kMerLength, minSeedLength, numSmers, quorum, profileSize, qProfile, revQ, uArr, posArray, hitArr, searchColors, false);
	}

	// }

	//Measure and output current runtime if demanded
	if(calcRT){
		//Get current time
		endTime = std::chrono::system_clock::now();
		//Calculate the time difference
		tDiff = endTime - startTime;
		//Output the measured time
		cout << "Seed detection took " << tDiff.count() << " s" << endl;
		//Update start time for next part
		startTime = std::chrono::system_clock::now();
	} else{
		//Testing
		// cout << "Let's see which seeds we have found" << endl;
		// for(ColoredCDBG<seedlist>::iterator i = cdbg.begin(); i != cdbg.end(); ++i){
		// 	struct seed *next = i->getData()->getData(*i)->getSeed(false);

		// 	while(next != NULL){
		// 		cout << i->referenceUnitigToString() << ": u:" << next->offsetU << " q:" << next->offsetQ << " len:" << next->len << " nextSeed is " << (next->nextSeed == NULL ? "" : "not ") << "NULL" << endl;
		// 			next = next->nextSeed;
		// 	}
		// }
	}

	/*Seed Extension*/
	cout << "Extending seeds" << endl;

	// //Initialize rest of hit array
	// for(uint32_t i = hitArrSize - minSeedLength + 1; i < hitArrSize; ++i){
	// 	hitArr[i].length = 0;
	// }

	//Check which function we have to use
	// if(searchColors.empty()){//TODO It is not nice to copy all is! Is it maybe possible to circumvent this by giving a function as parameter or does it actually not even save time?
	// 	extendSeeds(cdbg, q, X, hitArr, quorum);
	// } else{

	//Extend seeds lying on the reference strand if demanded
	if(strand != Minus) extendRefSeeds(cdbg, q, X, hitArr, quorum, searchColors);

	//Testing
	// cout << "Reference strand done" << endl;

	//Extend seeds lying on the reverse complementary strand if demanded
	if(strand != Plus) extendRevCompSeeds(cdbg, q, X, hitArr, quorum, searchColors);

	//Testing
	// cout << "Reverse complementary strand done" << endl;
	// exit(0);

	// }

	//Measure and output current runtime if demanded
	if(calcRT){
		//Get current time
		endTime = std::chrono::system_clock::now();
		//Calculate the time difference
		std::chrono::duration<double> tDiff = endTime - startTime;
		//Output the measured time
		cout << "Seed extension took " << tDiff.count() << " s" << endl;
		//Update start time for next part
		startTime = std::chrono::system_clock::now();
	} else{
		//Testing
		// cout << "Hit array:" << endl;
		// for(uint32_t i = 0; i < q.length() - minSeedLength + 1; ++i){
		// 	struct hit *hitList;

		// 	if(hitArr[i].length == 0){
		// 		hitList = NULL;
		// 	} else{
		// 		hitList	= &hitArr[i];
		// 	}

		// 	while(hitList != NULL){
		// 		cout << "Hit info: length: " << hitList->length << " score: " << hitList->score << " left seed: Uoff: " << hitList->lSeedUoff << " Qoff: " << hitList->lSeedQoff << " Unitig: " << hitList->lUnitig.mappedSequenceToString() << " right seed: Uoff: " << hitList->rSeedUoff << " Qoff: " << hitList->rSeedQoff << " Unitig: " << hitList->rUnitig.mappedSequenceToString() << endl;

		// 		hitList = hitList->nextHit;
		// 	}
		// }
	}

	/*Gapped extension*/
	cout << "Performing gapped extension" << endl;

	//Find best results
	for(uint32_t i = 0; i < hitArrSize; ++i){//TODO: Here, we should scan for duplicates already!//TODO This may be put into a function!
		hit *hitList;

		if(hitArr[i].length == 0){
			//Testing
			// cout << "2 Option 2" << endl;

			hitList = NULL;
		} else{
			//Testing
			// cout << "2 Option 1" << endl;

			hitList	= &hitArr[i];
		}

		//Go through all hits of a certain array position
		while(hitList != NULL){
			//Testing 
			// cout << "Let's make valid borders" << endl;
			// cout << "hitList->score: " << hitList->score << endl;

			// //Adjust right hit border if necessary
			// makeValidBorders(hitList, q);

			//Calculate e-value
			hitList->eval = calcEVal(hitList->score, lambda, C, q.length());

			//Testing
			// cout << "hitList->score: " << hitList->score << endl;

			//Make sure the hit should still be considered
			if(hitList->score != 0 && hitList->eval <= eLim){
				//Testing
				// cout << "4 Option 2" << endl << "5 Option 2" << endl;

				//Check if we can still just add new results to the result list or need to replace worse existing ones
				if(nRes == 0){
					//Testing
					// cout << "3 Option 1" << endl;

					//Try to replace a worse result
					replWorseRes(resList, hitList);//TODO: Here, we should start to free the memory of all hits not needed anymore!<-Why not doing this in the very end? This would avoid unnecessary pointer juggling
				} else{
					//Testing
					// cout << "3 Option 2" << endl;

					//Insert a new result
					insRes(resList, hitList);
					--nRes;
				}
			}

			//Testing
			// if(hitList->score == 0) cout << "4 Option 1" << endl;
			// if(hitList->eval > eLim) cout << "5 Option 1" << endl;

			hitList = hitList->nextHit;
		}
	}

	//Check which function we have to use
	// if(searchColors.empty()){//TODO It is not nice to copy all is! Is it maybe possible to circumvent this by giving a function as parameter or does it actually not even save time?
	// 	calcGappedAlignment(resList, q, X, calcRT, quorum);
	// } else{

	//Calculate gapped alignments
	calcGappedAlignment(cdbg, resList, q, X, quorum, searchColors, lambda, C);

	//Sort results
	resList.sort(compEvals);

	//Measure and output current runtime if demanded
	if(calcRT){
		//Get current time
		endTime = std::chrono::system_clock::now();
		//Calculate the time difference
		tDiff = endTime - startTime;
		//Output the measured time
		cout << "Gapped extension took " << tDiff.count() << " s" << endl;
	} else{
		//Iterate over alignments
		for(list<hit*>::const_iterator it = resList.begin(); it != resList.end(); ++it){
			//Output alignment
			repAlgn(*it);
			//Output color sets
			// outpColSets(cdbg, *it);
		}
	}

	//Free memory in hit array
	freeHitArray(hitArr, hitArrSize);

	//Testing
	// struct seed *s;
	// for(ColoredCDBG<seedlist>::iterator i = cdbg.begin(); i != cdbg.end(); ++i){
	// 	s = i->getData()->getData(*i)->getSeed(true);
	// 	if(s != NULL) cerr << "Detected remaining seed in seedlist after calculations" << endl;
	// 	s = i->getData()->getData(*i)->getSeed(false);
	// 	if(s != NULL) cerr << "Detected remaining seed in seedlist after calculations" << endl;
	// }
}

//This function frees all memory additionally allocated for an hit array
void freeHitArray(hit *arr, uint32_t &arrLen){
	bool frst = false;
	hit *h, *delHit;

	for(uint32_t i = 0; i < arrLen; ++i){
		//Check if hit array at current position is empty
		if(arr[i].length == 0){
			//Testing
			// cout << "1 Option 1" << endl;

			h = NULL;
		} else{
			//Testing
			// cout << "1 Option 2" << endl;

			//Get pointer of hit at current position
			h = &arr[i];
			//Mark hit as being the first in the hit array
			frst = true;
		}

		//Iterate over all hits of a certain array position
		while(h != NULL){
			//Check if we are dealing with the first hit at this position
			if(frst){
				//Testing
				// cout << "2 Option 1" << endl;

				//Next hit won't be the first anymore
				frst = false;
				//Move to the next hit
				h = h->nextHit;
			} else{
				//Testing
				// cout << "2 Option 2" << endl;

				//Save hit to delete
				delHit = h;
				//Get pointer to next hit
				h = h->nextHit;
				//Free current hit
				free(delHit);
			}
		}
	}
}