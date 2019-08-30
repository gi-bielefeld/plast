#include "Search.h"
#include "Extension.h"
#include "Extension.cpp"
#include "GappedAlignment.cpp"
#include "Sequence.cpp"
#include "Hit.cpp"

//This function performs the seed detection if no search color set is given
void detectSeeds(const int32_t &k, const int32_t &minSeedLength, const size_t &numSmers, const uint32_t &quorum, const uint32_t &profileSize, const uint32_t *qProfile, const string &q, const UnitigColorMap<seedlist> *uArr, const struct s_mer_pos *posArray, struct hit *hitArr, const bool isRefSeq){
	uint32_t occnum;
	int32_t pRank = -1;
	size_t uniLen;
	struct seed *newSeed, *lastSeed;

	//Testing
	// cout << "Do we get here?" << endl;

	//Go through the query and check for each s-mer whether it appears in the graph
	for(uint32_t i = 0; i < q.length() - minSeedLength + 1; ++i){
		//Once we are iterating over q anyways initialize hitArr NOTE: Initialized hits have a length of 0//TODO: This shouldn't be in here!
		hitArr[i].length = 0;
		//Calculate current s-mer's rank
		pRank = compRank(q.substr(i, minSeedLength), pRank);

		//Testing
		//cout << "Do we get here as well?" << endl;

		//To avoid border effects check whether we are at the end of the q-gram profile
		if(static_cast<unsigned>(pRank) < profileSize - 1 && qProfile[pRank] < numSmers){
			//Testing
			//cout << "Testcase: Not at the end of s-mer index" << endl;

			occnum = qProfile[pRank + 1] - qProfile[pRank];//if occnum is zero the current s-mer doesn't exist
		} else{
			//Testing
			// cout << "Testcase: At the end of s-mer index" << endl;
			// cout << "compRank(q.substr(i, minSeedLength), pRank): " << compRank(q.substr(i, minSeedLength), pRank) << endl;

			occnum = numSmers - qProfile[pRank];//if occnum is zero the current s-mer doesn't exist
		}

		//Process all occurrences of a certain s-mer
		for(int32_t j = occnum - 1; j >= 0; --j){//for(uint32_t j = 0; j < occnum; ++j){
			//Get current unitig's length
			uniLen = uArr[posArray[qProfile[pRank] + j].unitig_id].size;

			//Check if the seed list is empty
			if(uArr[posArray[qProfile[pRank] + j].unitig_id].getData()->getData(uArr[posArray[qProfile[pRank] + j].unitig_id])->getSeed(isRefSeq) == NULL){
				//Check whether our quorum is fulfilled move on otherwise
				if(checkSearchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], quorum, posArray[qProfile[pRank] + j].offset) == 0){
					//Testing
					//cout << "Testcase: S-mer has to be inserted as the first element in a seed list and quorum is not fulfilled" << endl;

					continue;
				}

				//Testing
				// if(i == 489 && !strcmp(uArr[posArray[qProfile[pRank] + j].unitig_id].referenceUnitigToString().c_str(), s)){
				// 	cout << "j:" << j << endl;
				// 	exit(0);
				// }
				//cout << "Testcase: S-mer has to be inserted as the first element in a seed list and quorum is fulfilled" << endl;

				//Create a new seed and...
				newSeed = (struct seed*) malloc(sizeof(struct seed));
				newSeed->offsetU = compOffset(posArray[qProfile[pRank] + j].offset, minSeedLength, uniLen, isRefSeq);
				newSeed->offsetQ = compOffset(i, minSeedLength, q.length(), isRefSeq);
				newSeed->len = minSeedLength;
				newSeed->score = 0;
				newSeed->nextSeed = NULL;

				//Testing
				// if(newSeed->offsetU == 20 && newSeed->offsetQ == 8024 && !isRefSeq){
				// 	cerr << "Seed uoff:" << newSeed->offsetU << " qoff:" << newSeed->offsetQ << " has been created" << endl << "Unitig is " << uArr[posArray[qProfile[pRank] + j].unitig_id].mappedSequenceToString() << endl;
				// }

				//...make it the first element of the seed list
				uArr[posArray[qProfile[pRank] + j].unitig_id].getData()->getData(uArr[posArray[qProfile[pRank] + j].unitig_id])->pushSeed(newSeed, isRefSeq);
			} else{
				//Testing
				//cout << "Testcase: S-mer has not to be inserted as the first element in a seed list and quorum is " << (checkSearchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], quorum, posArray[qProfile[pRank] + j].offset) != 0 ? "" : "not ") << "fulfilled" << endl;
				//Check on which strand we are
				if(isRefSeq){
					//Get the last inserted seed
					lastSeed = uArr[posArray[qProfile[pRank] + j].unitig_id].getData()->getData(uArr[posArray[qProfile[pRank] + j].unitig_id])->getSeed(isRefSeq);

					//Check if we can extend the last seed inserted seed in our seed list
					if(lastSeed->offsetU + lastSeed->len - (minSeedLength - 1) == posArray[qProfile[pRank] + j].offset && lastSeed->offsetQ + lastSeed->len - (minSeedLength - 1) == i){
						//Check whether we have to check our quorum at the current position
						if(lastSeed->len % k == 0){
							//Check whether our quorum is fulfilled move on otherwise
							if(checkSearchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], quorum, posArray[qProfile[pRank] + j].offset) == 0){ continue; }
						}

						//Extending simply means to increment the seed's length
						++lastSeed->len;
					} else {
						//Check whether our quorum is fulfilled move on otherwise
						if(checkSearchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], quorum, posArray[qProfile[pRank] + j].offset) == 0){ continue; }

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
					//Get the last inserted seed
					lastSeed = uArr[posArray[qProfile[pRank] + j].unitig_id].getData()->getData(uArr[posArray[qProfile[pRank] + j].unitig_id])->getLastInSd();

					//Check if we can extend the last seed inserted seed in our seed list
					if(lastSeed->offsetU > 0 && lastSeed->offsetU - 1 == compOffset(posArray[qProfile[pRank] + j].offset, minSeedLength, uniLen, isRefSeq) && lastSeed->offsetQ > 0 && lastSeed->offsetQ - 1 == compOffset(i, minSeedLength, q.length(), isRefSeq)){
						//Check whether we have to check our quorum at the current position
						if(lastSeed->len % k == 0){
							//Check whether our quorum is fulfilled move on otherwise
							if(checkSearchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], quorum, posArray[qProfile[pRank] + j].offset) == 0){ continue; }
						}

						//Extending means to increment the seed's length and to adjust the offsets
						++lastSeed->len;
						--lastSeed->offsetU;
						--lastSeed->offsetQ;
					} else{
						//Check whether our quorum is fulfilled move on otherwise
						if(checkSearchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], quorum, posArray[qProfile[pRank] + j].offset) == 0){ continue; }

						//If an extension is not possible create a new seed
						newSeed = (struct seed*) malloc(sizeof(struct seed));
						newSeed->offsetU = compOffset(posArray[qProfile[pRank] + j].offset, minSeedLength, uniLen, isRefSeq);
						newSeed->offsetQ = compOffset(i, minSeedLength, q.length(), isRefSeq);
						newSeed->len = minSeedLength;
						newSeed->score = 0;
						newSeed->nextSeed = NULL;

						//Testing
						// if(newSeed->offsetU == 20 && newSeed->offsetQ == 8024 && !isRefSeq){
						// 	cerr << "Seed uoff:" << newSeed->offsetU << " qoff:" << newSeed->offsetQ << " has been created" << endl << "Unitig is " << uArr[posArray[qProfile[pRank] + j].unitig_id].mappedSequenceToString() << endl;
						// 	cerr << "Real offsets were u:" << posArray[qProfile[pRank] + j].offset << " q:" << i << endl << "unitig length:" << uArr[posArray[qProfile[pRank] + j].unitig_id].size << endl;
						// 	//exit(0);
						// }

						//Insert the new seed into the seed list
						uArr[posArray[qProfile[pRank] + j].unitig_id].getData()->getData(uArr[posArray[qProfile[pRank] + j].unitig_id])->pushSeed(newSeed, isRefSeq);
					}
				}
			}

			//Testing
			// if(newSeed->offsetQ > 9999){
			// 	cout << "Something is wrong here" << endl << "isRefSeq:" << isRefSeq << endl << "i:" << i << "minSeedLength:" << minSeedLength << "uniLen:" << uniLen << endl;
			// 	exit(0);
			// }
		}
	}
}

//This function performs the seed detection if a search color set is given
void detectSeeds(const int32_t &k, const int32_t &minSeedLength, const size_t &numSmers, const uint32_t &quorum, const uint32_t &profileSize, const uint32_t *qProfile, const string &q, const UnitigColorMap<seedlist> *uArr, const struct s_mer_pos *posArray, struct hit *hitArr, const list<pair<string, size_t>> &searchSet, const bool isRefSeq){
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

			//Check if the seed list is empty
			if(uArr[posArray[qProfile[pRank] + j].unitig_id].getData()->getData(uArr[posArray[qProfile[pRank] + j].unitig_id])->getSeed(isRefSeq) == NULL){
				//Check whether our quorum is fulfilled move on otherwise
				if(checkSearchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], quorum, posArray[qProfile[pRank] + j].offset, searchSet) == 0){
					//Testing
					//cout << "Testcase: S-mer has to be inserted as the first element in a seed list and quorum is not fulfilled" << endl;

					continue;
				}

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

						//Check whether we have to check our quorum at the current position
						if(lastSeed->len % k == 0){
							//Testing
							// cout << "not ";

							//Check whether our quorum is fulfilled move on otherwise
							if(checkSearchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], quorum, posArray[qProfile[pRank] + j].offset, searchSet) == 0){
								//Testing
								// cout << "2" << endl;
								// cout << "8 Option 2" << endl;

								continue;
							}

							//Testing							
							// cout << "8 Option 1" << endl;
						}

						//Testing
						// cout << "2" << endl;

						//Extending simply means to increment the seed's length
						++lastSeed->len;
					} else {
						//Testing
						// cout << "7 Option 3" << endl;

						//Check whether our quorum is fulfilled move on otherwise
						if(checkSearchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], quorum, posArray[qProfile[pRank] + j].offset, searchSet) == 0){
							//Testing
							// cout << "9 Option 2" << endl;

							continue;
						}

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

						//Check whether we have to check our quorum at the current position
						if(lastSeed->len % k == 0){
							//Testing
							// cout << "not ";

							//Check whether our quorum is fulfilled move on otherwise
							if(checkSearchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], quorum, posArray[qProfile[pRank] + j].offset, searchSet) == 0){
								//Testing
								// cout << "2" << endl;
								// cout << "11 Option 2" << endl;

								continue;
							}

							//Testing
							// cout << "11 Option 1" << endl;
						}

						//Testing
						// cout << "2" << endl;

						//Extending means to increment the seed's length and to adjust the offsets
						++lastSeed->len;
						--lastSeed->offsetU;
						--lastSeed->offsetQ;
					} else{
						//Testing
						// cout << "10 Option 3" << endl;

						//Check whether our quorum is fulfilled move on otherwise
						if(checkSearchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], quorum, posArray[qProfile[pRank] + j].offset, searchSet) == 0){
							//Testing
							// cout << "12 Option 2" << endl;

							continue;
						}

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

//This function extends all seeds found on the queries reference strand considering a quorum
void extendRefSeeds(ColoredCDBG<seedlist> &cdbg, const string &q, const int16_t &X, struct hit *hitArr, const uint32_t &quorum){
	struct seed *currSeed;
	//struct seed *delSeed = NULL;
	struct hit newHit;
	UnitigColorMap<seedlist> currUni;

	//Testing
	//uint32_t nb_seeds = 0;

	//Iterate over all seeds of all unitigs
	for(ColoredCDBG<seedlist>::iterator i = cdbg.begin(); i != cdbg.end(); ++i){
		//Get current unitig
		currUni = *i;
		//Get the first seed
		currSeed = currUni.getData()->getData(currUni)->getSeed(currUni.strand);

		//Iterate over all seeds of a unitig
		while(currSeed != NULL){
			//TODO: This procedure might be worth to be put into an external (inline) function!
			//Extract seed from the list of unextended seeds
			currUni.getData()->getData(currUni)->setSeed(currSeed->nextSeed, currUni.strand);
			//Setup initial hit infos
			newHit.score = 0;
			newHit.length = currSeed->len;
			newHit.lSeedUoff = currSeed->offsetU;
			newHit.lSeedQoff = currSeed->offsetQ;
			newHit.lUnitig = currUni;
			newHit.rSeedUoff = currSeed->offsetU;
			newHit.rSeedQoff = currSeed->offsetQ;
			newHit.rUnitig = currUni;
			newHit.nextHit = NULL;

			//Testing
			//cout << "Start right extension" << endl;
			//cout << "Current seed is uoff:" << currSeed->offsetU << " qoff:" << currSeed->offsetQ << " len:" << currSeed->len << endl;
			// if(currSeed->offsetU == 0 && currSeed->offsetQ == 870 && currSeed->len == 32){
			// 	report = true;
			// }

			//Extend hit to the right
			startRightX_Drop(&newHit, q, X, quorum);
			/*TODO Ideas to make the seed extension faster:
				1. It might be faster to delete seeds not before a maximum scoring path over many untigs has been calculated. Many seeds could be incorporated into different extensions.
				2. If s < k s-mers in the beginning of two unitigs that are both successors of a third one are the same. If an extension is already pretty bad when these unitigs are explored s.t. it drops while comparing base pairs of the identical s-mers it is sufficient to do this once and not for each unitig separately. (But the X-drop never drops while exploring a reached seed does it?)*/

			//Testing
			// if(newHit.rSeedUoff == 5 && newHit.rSeedQoff == 9999 && newHit.lSeedUoff == 0 && newHit.lSeedQoff == 4545){
			// 	cerr << "After right extension: length: " << newHit.length << " score: " << newHit.score << endl;
			// 	//cerr << "Left seed: uoff:" << newHit.lSeedUoff << " qoff: " << newHit.lSeedQoff << endl;
			// }
			//cout << "Right extension done" << endl;
			// if((currSeed->offsetU == 0 && currSeed->offsetQ == 9492 && currSeed->len == 14) || report){
			// 	cerr << "Interesting seed found! lSeedUoff: " << newHit.lSeedUoff << " lSeedQoff: " << newHit.lSeedQoff << " rSeedUoff: " << newHit.rSeedUoff << " rSeedQoff: " << newHit.rSeedQoff << " score: " << newHit.score << endl;
			// }
			//++nb_seeds;
			//if(currSeed->offsetU == 0 && currSeed->offsetQ == 870 && currSeed->len == 32){
			//	cout << "After right extension: rSeed: uoff:" << newHit.rSeedUoff << " qoff:" << newHit.rSeedQoff << " length: " << newHit.length << " score: " << newHit.score << endl;
				//exit(0);
				//cerr << "Left seed: uoff:" << newHit.lSeedUoff << " qoff: " << newHit.lSeedQoff << endl;
			//}

			//Filter out some seeds; the second condition ensures that we do not miss seeds in the end of the query
			if(newHit.length - currSeed->len > 0 || currSeed->offsetQ + currSeed->len == q.length()){
				//Testing
				//cout << "We passed the filter" << endl;
				// cout << "After right extension: rSeed: uoff:" << newHit.rSeedUoff << " qoff:" << newHit.rSeedQoff << " length: " << newHit.length << " score: " << newHit.score << endl
				
				//Extend hit to the left
				startLeftX_Drop(&newHit, q, X, quorum);

				//Testing
				// if(newHit.length == 11 && newHit.score == 11){
				// 	cerr << "This should not happen" << endl << "Seed is uoff:" << currSeed->offsetU << " qoff:" << currSeed->offsetQ << " length:" << currSeed->len << endl;;
				// 	exit(0);
				// }

				//Check whether we have created a hit for this query position already
				if(hitArr[newHit.lSeedQoff].length == 0){
					//Use hit space in the array
					hitArr[newHit.lSeedQoff].score = newHit.score;
					hitArr[newHit.lSeedQoff].length = newHit.length;
					hitArr[newHit.lSeedQoff].lSeedUoff = newHit.lSeedUoff;
					hitArr[newHit.lSeedQoff].lSeedQoff = newHit.lSeedQoff;
					hitArr[newHit.lSeedQoff].lUnitig = newHit.lUnitig;
					hitArr[newHit.lSeedQoff].rSeedUoff = newHit.rSeedUoff;
					hitArr[newHit.lSeedQoff].rSeedQoff = newHit.rSeedQoff;
					hitArr[newHit.lSeedQoff].rUnitig = newHit.rUnitig;
					hitArr[newHit.lSeedQoff].nextHit = newHit.nextHit;
				} else{
					//Save pointer to the hit list at the current q position
					newHit.nextHit = hitArr[newHit.lSeedQoff].nextHit;
					//Create a new hit
					hitArr[newHit.lSeedQoff].nextHit = (struct hit*) malloc(sizeof(struct hit));
					//Copy information
					hitArr[newHit.lSeedQoff].nextHit->score = newHit.score;
					hitArr[newHit.lSeedQoff].nextHit->length = newHit.length;
					hitArr[newHit.lSeedQoff].nextHit->lSeedUoff = newHit.lSeedUoff;
					hitArr[newHit.lSeedQoff].nextHit->lSeedQoff = newHit.lSeedQoff;
					hitArr[newHit.lSeedQoff].nextHit->lUnitig = newHit.lUnitig;
					hitArr[newHit.lSeedQoff].nextHit->rSeedUoff = newHit.rSeedUoff;
					hitArr[newHit.lSeedQoff].nextHit->rSeedQoff = newHit.rSeedQoff;
					hitArr[newHit.lSeedQoff].nextHit->rUnitig = newHit.rUnitig;
					hitArr[newHit.lSeedQoff].nextHit->nextHit = newHit.nextHit;
				}
			}// else{
			//Testing
			// if(newHit.lSeedUoff == 12 && newHit.lSeedQoff == 0 && newHit.rSeedUoff == 2 && newHit.rSeedQoff == 9999){
			// 	cout << "Seed found: upos:" << currSeed->offsetU << " qpos:" << currSeed->offsetQ << " len:" << currSeed->len << endl;
			// 	exit(0);
			// }

			// 	//Mark that this seed can be deleted
			// 	delSeed = currSeed;
			// }

			//Delete the current seed
			free(currSeed);
			//Move on to the next seed
			currSeed = currUni.getData()->getData(currUni)->getSeed(currUni.strand);
		}
	}
}

//This function extends all seeds found on the queries reverse complement considering a quorum
void extendRevCompSeeds(ColoredCDBG<seedlist> &cdbg, const string &q, const int16_t &X, struct hit *hitArr, const uint32_t &quorum){
	bool isValid = true;
	struct seed *currSeed;
	struct hit newHit;
	struct hit *hitIt;
	UnitigColorMap<seedlist> currUni;

	//Iterate over all seeds of all unitigs
	for(ColoredCDBG<seedlist>::iterator i = cdbg.begin(); i != cdbg.end(); ++i){
		//Get current unitig
		currUni = *i;
		//We are on the reverse complementary sequence
		currUni.strand = false;
		//Get the first seed
		currSeed = currUni.getData()->getData(currUni)->getSeed(currUni.strand);

		//Testing
		//if(currSeed != NULL) cout << "Extending seed: upos:" << currSeed->offsetU << " qpos:" << currSeed->offsetQ << endl;

		//Testing
		//cout << "Current unitig is " << i->mappedSequenceToString() << endl;

		//Iterate over all seeds of a unitig
		while(currSeed != NULL){
			//Testing
			//cout << "Current seed is u offset:" << currSeed->offsetU << " q offset:" << currSeed->offsetQ << " len:" << currSeed->len << " score:" << currSeed->score << endl;

			//TODO: This procedure might be worth to be put into an external (inline) function!
			//Extract seed from the list of unextended seeds
			currUni.getData()->getData(currUni)->setSeed(currSeed->nextSeed, currUni.strand);
			//Setup initial hit infos
			newHit.score = 0;
			newHit.length = currSeed->len;
			newHit.lSeedUoff = currSeed->offsetU;
			newHit.lSeedQoff = currSeed->offsetQ;
			newHit.lUnitig = currUni;
			newHit.rSeedUoff = currSeed->offsetU;
			newHit.rSeedQoff = currSeed->offsetQ;
			newHit.rUnitig = currUni;
			newHit.nextHit = NULL;

			//Testing
			//cout << "Next we extend seed: upos:" << currSeed->offsetU << " qpos:" << currSeed->offsetQ << endl;	
			// if(currSeed->offsetQ == 9){
			// 	cout << "Extending our seed of interest" << endl;
			// }// else{
			// 	cout << "There is none" << endl;
			// }

			//Extend hit to the right
			startRightX_Drop_OnRevComp(&newHit, q, X, quorum);
			/*TODO Ideas to make the seed extension faster:
				1. It might be faster to delete seeds not before a maximum scoring path over many untigs has been calculated. Many seeds could be incorporated into different extensions.
				2. If s < k s-mers in the beginning of two unitigs that are both successors of a third one are the same. If an extension is already pretty bad when these unitigs are explored s.t. it drops while comparing base pairs of the identical s-mers it is sufficient to do this once and not for each unitig separately. (But the X-drop never drops while exploring a reached seed does it?)*/

			//Testing
			// if(currSeed->offsetQ == 9){
			// 	cout << "After right extension the hit's right border in q is " << newHit.rSeedQoff << endl;
			// }

			//Filter out some seeds; the second condition ensures that we do not miss seeds in the end of the query
			if(newHit.length - currSeed->len > 0 || currSeed->offsetQ + currSeed->len == q.length()){
				//Extend hit to the left
				startLeftX_Drop_OnRevComp(&newHit, q, X, quorum);

				//Testing
				// if(currSeed->offsetQ == 9){
				// 	cout << "After left extension our hit's q offset is " << newHit.lSeedQoff << endl;
				// }

				//Check whether we have created a hit for this query position already
				if(hitArr[newHit.lSeedQoff].length != 0){
					//Get the first hit at this position
					hitIt = &hitArr[newHit.lSeedQoff];
					//Go through the list of all hits at this position
					while(hitIt != NULL){
						//Check if the current hit and the last extended hit might be cooptimal solutions (i.e. they have the same q offsets and score)
						if(hitIt->lSeedQoff == newHit.lSeedQoff && hitIt->rSeedQoff == newHit.rSeedQoff && hitIt->score == newHit.score){
							//Testing
							//cerr << "Does this happen often?" << endl;

							//Mark hit as not interesting
							isValid = false;
							break;
						}

						//Move on to the next hit
						hitIt = hitIt->nextHit;
					}

					//Testing
					// if(currSeed->offsetU == 57 && currSeed->offsetQ == 1201 && currSeed->len == 12){
					// 	cerr << "Hit is lSeed: " << newHit.lUnitig.mappedSequenceToString() << " len:" << newHit.lUnitig.size << " uoff:" << newHit.lSeedUoff << " qoff:" << newHit.lSeedQoff << " rSeed: " << newHit.rUnitig.mappedSequenceToString() << " len:" << newHit.rUnitig.size << " uoff:" << newHit.rSeedUoff << " qoff:" << newHit.rSeedQoff << " score:" << newHit.score << endl;
					// 	exit(0);
					// }

					//Check if extended hit is worth to be kept
					if(!isValid){
						//Reset flag
						isValid = true;
						//continue;
					} else{
						//Adjust hit borders if necessary//TODO This should be done not before we select the best extended hits for the gapped alignment, but for the moment this is good enough. Apart from that we can just assume that a gapped extension will not be possible to extend our hits much further so that we can just take it out here and consider it later!
						//makeValidBorders(&newHit, q);

						//Add hit to the hit array//

						//Save pointer to the hit list at the current q position
						newHit.nextHit = hitArr[newHit.lSeedQoff].nextHit;
						//Create a new hit
						hitArr[newHit.lSeedQoff].nextHit = (struct hit*) malloc(sizeof(struct hit));
						//Copy information
						hitArr[newHit.lSeedQoff].nextHit->score = newHit.score;
						hitArr[newHit.lSeedQoff].nextHit->length = newHit.length;
						hitArr[newHit.lSeedQoff].nextHit->lSeedUoff = newHit.lSeedUoff;
						hitArr[newHit.lSeedQoff].nextHit->lSeedQoff = newHit.lSeedQoff;
						hitArr[newHit.lSeedQoff].nextHit->lUnitig = newHit.lUnitig;
						hitArr[newHit.lSeedQoff].nextHit->rSeedUoff = newHit.rSeedUoff;
						hitArr[newHit.lSeedQoff].nextHit->rSeedQoff = newHit.rSeedQoff;
						hitArr[newHit.lSeedQoff].nextHit->rUnitig = newHit.rUnitig;
						hitArr[newHit.lSeedQoff].nextHit->nextHit = newHit.nextHit;
					}
				} else{
					//Adjust hit borders if necessary//TODO This should be done not before we select the best extended hits for the gapped alignment, but for the moment this is good enough. Apart from that we can just assume that a gapped extension will not be possible to extend our hits much further so that we can just take it out here and consider it later!
					//makeValidBorders(&newHit, q);
					//Use hit space in the array
					hitArr[newHit.lSeedQoff].score = newHit.score;
					hitArr[newHit.lSeedQoff].length = newHit.length;
					hitArr[newHit.lSeedQoff].lSeedUoff = newHit.lSeedUoff;
					hitArr[newHit.lSeedQoff].lSeedQoff = newHit.lSeedQoff;
					hitArr[newHit.lSeedQoff].lUnitig = newHit.lUnitig;
					hitArr[newHit.lSeedQoff].rSeedUoff = newHit.rSeedUoff;
					hitArr[newHit.lSeedQoff].rSeedQoff = newHit.rSeedQoff;
					hitArr[newHit.lSeedQoff].rUnitig = newHit.rUnitig;
					hitArr[newHit.lSeedQoff].nextHit = newHit.nextHit;
				}
			}
			//Testing
			// if(newHit.score == 4294967295){
			// 	cerr << "Seed was uoff:" << currSeed->offsetU << " qoff:" << currSeed->offsetQ << " len:" << currSeed->len << endl;
			// 	exit(0);
			// }
			
			//Delete the current seed
			free(currSeed);
			//Move on to the next seed
			currSeed = currUni.getData()->getData(currUni)->getSeed(currUni.strand);
		}
	}
}

//This function extends all seeds found on the queries reference strand considering a quorum and a search color set
void extendRefSeeds(ColoredCDBG<seedlist> &cdbg, const string &q, const int16_t &X, struct hit *hitArr, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	struct seed *currSeed;
	struct hit newHit;
	UnitigColorMap<seedlist> currUni;

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
			newHit.lSeedUoff = currSeed->offsetU;
			newHit.lSeedQoff = currSeed->offsetQ;
			newHit.lUnitig = currUni;
			newHit.rSeedUoff = currSeed->offsetU;
			newHit.rSeedQoff = currSeed->offsetQ;
			newHit.rUnitig = currUni;
			newHit.nextHit = NULL;
			//Extend hit to the right
			startRightX_Drop(&newHit, q, X, quorum, searchSet);
			/*TODO Ideas to make the seed extension faster:
				1. It might be faster to delete seeds not before a maximum scoring path over many untigs has been calculated. Many seeds could be incorporated into different extensions.
				2. If s < k s-mers in the beginning of two unitigs that are both successors of a third one are the same. If an extension is already pretty bad when theses unitigs are explored s.t. it drops while comparing base pairs of the identical s-mers it is sufficient to do this once and not for each unitig separately. (But the X-drop never drops while exploring a reached seed does it?)*/			
			
			//Testing
			// cout << "We do get here" << endl;

			//Filter out some seeds; the second condition ensures that we do not miss seeds in the end of the query
			if(newHit.length - currSeed->len > 0 || currSeed->offsetQ + currSeed->len == q.length()){
				//Testing
				// cout << "We get here" << endl;

				//Extend hit to the left
				startLeftX_Drop(&newHit, q, X, quorum, searchSet);

				//Check whether we have created a hit for this query position already
				if(hitArr[newHit.lSeedQoff].length == 0){
					//Testing
					//cout << "4 Option 1" << endl;

					//Use hit space in the array
					hitArr[newHit.lSeedQoff].score = newHit.score;
					hitArr[newHit.lSeedQoff].length = newHit.length;
					hitArr[newHit.lSeedQoff].lSeedUoff = newHit.lSeedUoff;
					hitArr[newHit.lSeedQoff].lSeedQoff = newHit.lSeedQoff;
					hitArr[newHit.lSeedQoff].lUnitig = newHit.lUnitig;
					hitArr[newHit.lSeedQoff].rSeedUoff = newHit.rSeedUoff;
					hitArr[newHit.lSeedQoff].rSeedQoff = newHit.rSeedQoff;
					hitArr[newHit.lSeedQoff].rUnitig = newHit.rUnitig;
					hitArr[newHit.lSeedQoff].nextHit = newHit.nextHit;
				} else{
					//Testing
					//cout << "4 Option 2" << endl;

					//Save pointer to the hit list at the current q position
					newHit.nextHit = hitArr[newHit.lSeedQoff].nextHit;
					//Create a new hit
					hitArr[newHit.lSeedQoff].nextHit = (struct hit*) malloc(sizeof(struct hit));
					//Copy information
					hitArr[newHit.lSeedQoff].nextHit->score = newHit.score;
					hitArr[newHit.lSeedQoff].nextHit->length = newHit.length;
					hitArr[newHit.lSeedQoff].nextHit->lSeedUoff = newHit.lSeedUoff;
					hitArr[newHit.lSeedQoff].nextHit->lSeedQoff = newHit.lSeedQoff;
					hitArr[newHit.lSeedQoff].nextHit->lUnitig = newHit.lUnitig;
					hitArr[newHit.lSeedQoff].nextHit->rSeedUoff = newHit.rSeedUoff;
					hitArr[newHit.lSeedQoff].nextHit->rSeedQoff = newHit.rSeedQoff;
					hitArr[newHit.lSeedQoff].nextHit->rUnitig = newHit.rUnitig;
					hitArr[newHit.lSeedQoff].nextHit->nextHit = newHit.nextHit;
				}
			}

			//Delete the extended seed
			free(currSeed);
			//Move on to the next seed
			currSeed = currUni.getData()->getData(currUni)->getSeed(currUni.strand);
		}
	}
}

//This function extends all seeds found on the queries reverse complement considering a quorum and a search color set
void extendRevCompSeeds(ColoredCDBG<seedlist> &cdbg, const string &q, const int16_t &X, struct hit *hitArr, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	bool isValid = true;
	struct seed *currSeed;
	struct hit newHit;
	struct hit *hitIt;
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

		//Iterate over all seeds of a unitig
		while(currSeed != NULL){
			//TODO: This procedure might be worth to be put into an external (inline) function!
			//Extract seed from the list
			currUni.getData()->getData(currUni)->setSeed(currSeed->nextSeed, currUni.strand);
			//Setup initial hit infos
			newHit.score = 0;
			newHit.length = currSeed->len;
			newHit.lSeedUoff = currSeed->offsetU;
			newHit.lSeedQoff = currSeed->offsetQ;
			newHit.lUnitig = currUni;
			newHit.rSeedUoff = currSeed->offsetU;
			newHit.rSeedQoff = currSeed->offsetQ;
			newHit.rUnitig = currUni;
			newHit.nextHit = NULL;
			//Extend hit to the right
			startRightX_Drop_OnRevComp(&newHit, q, X, quorum, searchSet);
			/*TODO Ideas to make the seed extension faster:
				1. It might be faster to delete seeds not before a maximum scoring path over many untigs has been calculated. Many seeds could be incorporated into different extensions.
				2. If s < k s-mers in the beginning of two unitigs that are both successors of a third one are the same. If an extension is already pretty bad when theses unitigs are explored s.t. it drops while comparing base pairs of the identical s-mers it is sufficient to do this once and not for each unitig separately. (But the X-drop never drops while exploring a reached seed does it?)*/			

			//Testing
			// curLen = currSeed->len;
			// curOffU = currSeed->offsetU;
			// curOffQ = currSeed->offsetQ;
			// cout << "Seed is offsetU: " << currSeed->offsetU << " offsetQ: " << currSeed->offsetQ << endl;

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
				// cout << "Left extension done" << endl;

				//Check whether we have created a hit for this query position already
				if(hitArr[newHit.lSeedQoff].length != 0){
					//Testing
					// cout << "4 Option 1" << endl;

					//Get the first hit at this position
					hitIt = &hitArr[newHit.lSeedQoff];
					//Go through the list of all hits at this position
					while(hitIt != NULL){
						//Check if the current hit and the last extended hit might be cooptimal solutions (i.e. they have the same q offsets and score)
						if(hitIt->lSeedQoff == newHit.lSeedQoff && hitIt->rSeedQoff == newHit.rSeedQoff && hitIt->score == newHit.score){
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
						//Adjust hit borders if necessary//TODO This should be done not before we select the best extended hits for the gapped alignment, but for the moment this is good enough. Apart from that we can just assume that a gapped extension will not be possible to extend our hits much further so that we can just take it out here and consider it later!
						//makeValidBorders(&newHit, q);

						//Add hit to the hit array//

						//Save pointer to the hit list at the current q position
						newHit.nextHit = hitArr[newHit.lSeedQoff].nextHit;
						//Create a new hit
						hitArr[newHit.lSeedQoff].nextHit = (struct hit*) malloc(sizeof(struct hit));
						//Copy information
						hitArr[newHit.lSeedQoff].nextHit->score = newHit.score;
						hitArr[newHit.lSeedQoff].nextHit->length = newHit.length;
						hitArr[newHit.lSeedQoff].nextHit->lSeedUoff = newHit.lSeedUoff;
						hitArr[newHit.lSeedQoff].nextHit->lSeedQoff = newHit.lSeedQoff;
						hitArr[newHit.lSeedQoff].nextHit->lUnitig = newHit.lUnitig;
						hitArr[newHit.lSeedQoff].nextHit->rSeedUoff = newHit.rSeedUoff;
						hitArr[newHit.lSeedQoff].nextHit->rSeedQoff = newHit.rSeedQoff;
						hitArr[newHit.lSeedQoff].nextHit->rUnitig = newHit.rUnitig;
						hitArr[newHit.lSeedQoff].nextHit->nextHit = newHit.nextHit;
					}
				} else{
					//Testing
					// cout << "4 Option 2" << endl;

					//Adjust hit borders if necessary//TODO This should be done not before we select the best extended hits for the gapped alignment, but for the moment this is good enough. Apart from that we can just assume that a gapped extension will not be possible to extend our hits much further so that we can just take it out here and consider it later!
					//makeValidBorders(&newHit, q);
					//Use hit space in the array
					hitArr[newHit.lSeedQoff].score = newHit.score;
					hitArr[newHit.lSeedQoff].length = newHit.length;
					hitArr[newHit.lSeedQoff].lSeedUoff = newHit.lSeedUoff;
					hitArr[newHit.lSeedQoff].lSeedQoff = newHit.lSeedQoff;
					hitArr[newHit.lSeedQoff].lUnitig = newHit.lUnitig;
					hitArr[newHit.lSeedQoff].rSeedUoff = newHit.rSeedUoff;
					hitArr[newHit.lSeedQoff].rSeedQoff = newHit.rSeedQoff;
					hitArr[newHit.lSeedQoff].rUnitig = newHit.rUnitig;
					hitArr[newHit.lSeedQoff].nextHit = newHit.nextHit;

					//Testing
					// cout << "Do we get here?" << endl;
				}
			}

			//Delete the extended seed
			free(currSeed);
			//Move on to the next seed
			currSeed = currUni.getData()->getData(currUni)->getSeed(currUni.strand);

			//Testing
			// cout << "Do we get here?" << endl;
		}

		//Testing (This is just a sanity test to make sure that there are no processed seeds which are saved anymore)
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
	}
}

//This function calculates a banded, semi-global, gapped alignment on a list of results considering a quorum and outputs the result if demanded
void calcGappedAlignment(const list<struct hit*> &resList, const string &q, const int16_t &X, const bool &noOutput, const uint32_t &quorum){
	bool isDupl = false;
	uint32_t bandRadius;
	list<struct hit*>::const_iterator iter;

	//Testing
	// uint32_t reslistSize = 0;
	// uint32_t nb_duplicates = 0;

	//Go through result list
	for(list<struct hit*>::const_iterator it = resList.begin(); it != resList.end(); ++it){
		//Calculate the band width to be used during the gapped extension
		bandRadius = (*it)->length / GAP_RATIO;

		//Testing
		// if((*it)->score > 10000){
		// 	cerr << "Hit to extend: rSeed uoff:" << (*it)->rSeedUoff << " rSeed qoff:" << (*it)->rSeedQoff << " score" << (*it)->score << endl;
		// }
		//++reslistSize;

		//Calculate gapped extension to the right
		startRightGappedAlignment(*it, q, X, bandRadius, quorum);

		//Testing
		// if(counter < 100){
		// 	cout << "Result of right extension: rSeed uoff:" << (*it)->rSeedUoff << " rSeed qoff:" << (*it)->rSeedQoff << " rUnitig:" << (*it)->rUnitig.mappedSequenceToString() << endl;
		// 	++counter;
		// }

		//Check whether this is not the first hit in the result list
		if(it != resList.begin()){
			iter = resList.begin();

			//Check for duplicates //TODO: If we filter for duplicates we could keep a statistic on how many duplicates of a specific result have been filtered out which could then be made a part of the result to inform the user about it!
			while(it != iter){
				//Testing
				//if(counter < 100) cout << "it != iter" << endl;

				//We assume to have a duplicate if right border offsets and unitigs are identical
				if((*it)->lSeedQoff == (*iter)->lSeedQoff && (*it)->rSeedQoff == (*iter)->rSeedQoff && (*it)->score == (*iter)->score){//TODO What is a duplicate? To answer this question a comparison with BLAST should be done and the definition should be adjusted. Also it is questionable whether this is the correct place to filter out duplicated results!
					//cerr << "It seems that the hit going spanning from q=" << (*iter)->lSeedQoff << " to q=" << (*iter)->rSeedQoff << " is a duplicate" << endl;
					isDupl = true;
					//++nb_duplicates;
					break;
				}

				++iter;
			}
		}

		//We want to avoid doing left extensions for duplicates
		if(!isDupl){
			startLeftGappedAlignment(*it, q, X, bandRadius, quorum);
		}

		//Output result
		if(!noOutput && !isDupl){
			cout << "Result info: score: " << (*it)->score << " left seed: offset q: " << (*it)->lSeedQoff << " offset u: " << (*it)->lSeedUoff << " unitig: " << (*it)->lUnitig.mappedSequenceToString() << " right seed: offset q: " << (*it)->rSeedQoff << " offset u: " << (*it)->rSeedUoff << " unitig: " << (*it)->rUnitig.mappedSequenceToString() << endl;
		}

		//Reset duplicate flag
		isDupl = false;
	}

	//Testing
	//cerr << "Result list length:" << reslistSize << endl << "Number of duplicates:" << nb_duplicates << endl;
}

//This function calculates a banded, semi-global, gapped alignment on a list of results considering a quorum and a search color set and outputs the result if demanded
void calcGappedAlignment(const list<struct hit*> &resList, const string &q, const int16_t &X, const bool &noOutput, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	bool isDupl = false;
	uint32_t bandRadius;
	list<struct hit*>::const_iterator iter;

	//Testing
	// uint16_t hitCount = 0;
	// bool duplFound = false;
	// cout << "Do we get here?" << endl;

	//Go through result list
	for(list<struct hit*>::const_iterator it = resList.begin(); it != resList.end(); ++it){
		//Testing
		// ++hitCount;
		// cout << "Do we get here?" << endl;

		//Calculate the band width to be used during the gapped extension
		bandRadius = (*it)->length / GAP_RATIO;
		//Calculate gapped extension to the right
		startRightGappedAlignment(*it, q, X, bandRadius, quorum, searchSet);

		//Check whether this is not the first hit in the result list
		if(it != resList.begin()){
			iter = resList.begin();

			//Check for duplicates //TODO: If we filter for duplicates we could keep a statistic on how many duplicates of a specific result have been filtered out which could then be made a part of the result to inform the user about it!
			while(it != iter){
				//We assume to have a duplicate if right border offsets and unitigs are identical
				if((*it)->rSeedUoff == (*iter)->rSeedUoff && (*it)->rSeedQoff == (*iter)->rSeedQoff && (*it)->rUnitig == (*iter)->rUnitig){
					//cerr << "It seems that the hit going spanning from q=" << (*iter)->lSeedQoff << " to q=" << (*iter)->rSeedQoff << " is a duplicate" << endl;
					isDupl = true;
					break;
				}

				++iter;
			}
		}

		//We want to avoid doing left extensions for duplicates
		if(!isDupl){
			startLeftGappedAlignment(*it, q, X, bandRadius, quorum, searchSet);
		} else{
			//Testing
			// cout << "2 Option 1" << endl;
			// duplFound = true;
		}

		//Output result
		if(!noOutput && !isDupl){
			//Testing
			// cout << "3 Option 1" << endl;

			cout << "Result info: score: " << (*it)->score << " left seed: offset q: " << (*it)->lSeedQoff << " offset u: " << (*it)->lSeedUoff << " unitig: " << (*it)->lUnitig.mappedSequenceToString() << " right seed: offset q: " << (*it)->rSeedQoff << " offset u: " << (*it)->rSeedUoff << " unitig: " << (*it)->rUnitig.mappedSequenceToString() << endl;
		}// else{
		// 	//Testing
		// 	cout << "3 Option 2" << endl;
		// }

		//Reset duplicate flag
		isDupl = false;
	}
}

//This function performs the actual graph search for a query
void searchQuery(ColoredCDBG<seedlist> &cdbg, const int32_t &kMerLength, const int32_t &minSeedLength, const size_t &numSmers, const uint32_t &quorum, const uint32_t &profileSize, const uint32_t *qProfile, const string &q, const SrchStrd &strand, const UnitigColorMap<seedlist> *uArr, const struct s_mer_pos *posArray, const list<pair<string, size_t>> &searchColors, const int16_t &X, const bool &calcRT, uint16_t nRes){
	//Staff we need to measure run times
	auto startTime = std::chrono::system_clock::now();
	auto endTime = std::chrono::system_clock::now();
	std::chrono::duration<double> tDiff;


	//The query's reverse complement
	string revQ;
	
	//Query counter
	//uint32_t qCounter = 0;

	//Output which query we are working on
	//cout << "Query " << ++qCounter << endl;

	//Variables needed for the seed extension
	struct hit hitArr[q.length()];
	//Variables needed for the gapped extension
	list<struct hit*> resList;
	list<struct hit*>::iterator iter;	

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
	if(strand == Both){
		cout << "1 Option 3" << endl;
	}

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

	//Initialize rest of hit array
	for(uint32_t i = q.length() - minSeedLength + 1; i < q.length(); ++i){
		hitArr[i].length = 0;
	}

	//Check which function we have to use
	// if(searchColors.empty()){//TODO It is not nice to copy all is! Is it maybe possible to circumvent this by giving a function as parameter or does it actually not even save time?
	// 	extendSeeds(cdbg, q, X, hitArr, quorum);
	// } else{

	//Extend seeds lying on the reference strand if demanded
	if(strand != Minus) extendRefSeeds(cdbg, q, X, hitArr, quorum, searchColors);

	//Extend seeds lying on the reverse complementary strand if demanded
	if(strand != Plus) extendRevCompSeeds(cdbg, q, X, hitArr, quorum, searchColors);

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
	for(uint32_t i = 0; i < q.length() - minSeedLength + 1; ++i){//TODO: Here, we should scan for duplicates already!//TODO This may be put into a function!
		struct hit *hitList;

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

			//Adjust right hit border if necessary
			makeValidBorders(hitList, q);

			//Testing
			// cout << "hitList->score: " << hitList->score << endl;

			//Make sure the hit should still be considered
			if(hitList->score != 0){
				//Check if we can still just add new results to the result list or need to replace worse existing ones
				if(nRes == 0){
					//Testing
					// cout << "3 Option 1" << endl;

					//Try to replace a worse result
					replWorseRes(resList, hitList);
				} else{
					//Testing
					// cout << "3 Option 2" << endl;

					//Insert a new result
					insRes(resList, hitList);
					--nRes;
				}
			}

			hitList = hitList->nextHit;
		}
	}

	//Check which function we have to use
	// if(searchColors.empty()){//TODO It is not nice to copy all is! Is it maybe possible to circumvent this by giving a function as parameter or does it actually not even save time?
	// 	calcGappedAlignment(resList, q, X, calcRT, quorum);
	// } else{

	calcGappedAlignment(resList, q, X, calcRT, quorum, searchColors);

	// }

	//Measure and output current runtime if demanded
	if(calcRT){
		//Get current time
		endTime = std::chrono::system_clock::now();
		//Calculate the time difference
		tDiff = endTime - startTime;
		//Output the measured time
		cout << "Gapped extension took " << tDiff.count() << " s" << endl;
	}
}