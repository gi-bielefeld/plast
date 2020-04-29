#include "Search.h"
#include "Extension.cpp"
#include "GappedAlignment.cpp"
#include "Sequence.cpp"
#include "Hit.cpp"
#include "Statistics.h"

//This function checks if the search criteria are fulfilled for a region on a unitig starting at some offset and having a certain length
bool vfySrchCrit(UnitigColorMap<seedlist> unitig, const uint32_t &off, const int32_t &len, const uint32_t &quorum, const list<pair<string, size_t>> &srchColSet){
	//Check if search criteria have already been checked for this unitig and check if necessary
	if(!unitig.getData()->getData(unitig)->srchCritChckd()) calcSrchCritBrds(unitig, quorum, srchColSet);

	//Check if search criteria are not fulfilled at all on this unitig
	if(unitig.getData()->getData(unitig)->getlBrd() < 0 && unitig.getData()->getData(unitig)->getrBrd() == (int32_t) unitig.size) return false;

	//Check if a seed falls completely behind the left or in front of the right border (rBrd is never negative)
	if((int32_t) off + len <= unitig.getData()->getData(unitig)->getlBrd() || off > (uint32_t) unitig.getData()->getData(unitig)->getrBrd()) return true;

	return false;
}

//This function performs the seed detection if a search color set is given
void detectSeeds(const int32_t &k, const int32_t &minSeedLength, const size_t &numSmers, const uint32_t &quorum, const uint32_t &profileSize, const uint32_t *qProfile, const string &q, const UnitigColorMap<seedlist> *uArr, const struct s_mer_pos *posArray, Hit *hitArr, const list<pair<string, size_t>> &searchSet, const bool isRefSeq){
	uint32_t occnum;
	int32_t pRank = -1;
	size_t uniLen;
	struct Seed *newSeed, *lastSeed;

	//Go through the query and check for each s-mer whether it appears in the graph
	for(uint32_t i = 0; i < q.length() - minSeedLength + 1; ++i){
		//Once we are iterating over q anyways initialize hitArr NOTE: Initialized hits have a length of 0
		hitArr[i].length = 0;
		//Calculate current s-mer's rank
		pRank = compRank(q.substr(i, minSeedLength), pRank);

		//To avoid border effects check whether we are at the end of the q-gram profile
		if(static_cast<unsigned>(pRank) < profileSize - 1 && qProfile[pRank] < numSmers){
			occnum = qProfile[pRank + 1] - qProfile[pRank];//if occnum is zero the current s-mer doesn't exist
		} else{
			occnum = numSmers - qProfile[pRank];//if occnum is zero the current s-mer doesn't exist
		}

		//Process all occurrences of a certain s-mer
		for(int32_t j = occnum - 1; j >= 0; --j){
			//Get current unitig's length
			uniLen = uArr[posArray[qProfile[pRank] + j].unitig_id].size;

			//Find out if unitig fulfills search criteria
			if(!vfySrchCrit(uArr[posArray[qProfile[pRank] + j].unitig_id], posArray[qProfile[pRank] + j].offset, minSeedLength, quorum, searchSet)) continue;

			//Check if the seed list is empty
			if(uArr[posArray[qProfile[pRank] + j].unitig_id].getData()->getData(uArr[posArray[qProfile[pRank] + j].unitig_id])->getSeed(isRefSeq) == NULL){
				//Create a new seed and...
				newSeed = (struct Seed*) malloc(sizeof(struct Seed));
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
					//Get the last inserted seed
					lastSeed = uArr[posArray[qProfile[pRank] + j].unitig_id].getData()->getData(uArr[posArray[qProfile[pRank] + j].unitig_id])->getSeed(isRefSeq);

					//Check if we can extend the last seed inserted in our seed list
					if(lastSeed->offsetU + lastSeed->len - (minSeedLength - 1) == posArray[qProfile[pRank] + j].offset && lastSeed->offsetQ + lastSeed->len - (minSeedLength - 1) == i){
						//Extending simply means to increment the seed's length
						++lastSeed->len;
					} else {
						//If an extension is not possible create a new seed
						newSeed = (struct Seed*) malloc(sizeof(struct Seed));
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
						//Extending means to increment the seed's length and to adjust the offsets
						++lastSeed->len;
						--lastSeed->offsetU;
						--lastSeed->offsetQ;
					} else{
						//If an extension is not possible create a new seed
						newSeed = (struct Seed*) malloc(sizeof(struct Seed));
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

//This function extends all seeds found on the queries reference strand considering a quorum and a search color set
void extendRefSeeds(ColoredCDBG<seedlist> &cdbg, const string &q, const int16_t &X, Hit *hitArr, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	struct Seed *currSeed;
	Hit newHit;
	UnitigColorMap<seedlist> currUni;

	//Iterate over all seeds of all unitigs
	for(ColoredCDBG<seedlist>::iterator i = cdbg.begin(); i != cdbg.end(); ++i){
		//Get current unitig
		currUni = *i;
		//Get the first seed
		currSeed = currUni.getData()->getData(currUni)->getSeed(currUni.strand);

		//Iterate over all seeds of a unitig
		while(currSeed != NULL){
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
			
			//Filter out some seeds; the second condition ensures that we do not miss seeds in the end of the query.
			//Note: What we do not consider here is that some seeds might not be extended to the right because search criteria are not fullfilled anymore. This is intended though. We should not miss too much, because a good hit should have more than one seed
			if(newHit.length - currSeed->len > 0 || currSeed->offsetQ + currSeed->len == q.length()){
				//Extend hit to the left
				startLeftX_Drop(&newHit, q, X, quorum, searchSet);

				//Check whether we have created a hit for this query position already
				if(hitArr[newHit.offQ].length == 0){
					//Copy hit into the array
					hitArr[newHit.offQ] = newHit;
				} else{
					//Save pointer to the hit list at the current q position
					newHit.nextHit = hitArr[newHit.offQ].nextHit;
					//Create a new hit
					hitArr[newHit.offQ].nextHit = new Hit(newHit);

					hitArr[newHit.offQ].nextHit->gAlgn.aSeqG = "";
					hitArr[newHit.offQ].nextHit->gAlgn.aSeqQ = "";
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
//Note: What we do not consider here is that some seeds might not be extended to the right because search criteria are not fullfilled anymore. This is intended though. We should not miss too much, because a good hit should have more than one seed
void extendRevCompSeeds(ColoredCDBG<seedlist> &cdbg, const string &q, const int16_t &X, Hit *hitArr, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	bool isValid = true;
	struct Seed *currSeed;
	Hit newHit;
	Hit *hitIt;
	UnitigColorMap<seedlist> currUni;

	//Iterate over all seeds of all unitigs
	for(ColoredCDBG<seedlist>::iterator i = cdbg.begin(); i != cdbg.end(); ++i){
		//Get current unitig
		currUni = *i;
		//We are on the reverse complementary sequence
		currUni.strand = false;
		//Get the first seed
		currSeed = currUni.getData()->getData(currUni)->getSeed(currUni.strand);

		//Iterate over all seeds of a unitig
		while(currSeed != NULL){
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
			startRightX_Drop_OnRevComp(&newHit, q, X, quorum, searchSet);

			//Filter out some seeds; the second condition ensures that we do not miss seeds in the end of the query and the third that we do not miss seeds in the end of a unitig if the unitig does not have successors (i.e. predecessors on the reference strand)
			if(newHit.length - currSeed->len > 0 || currSeed->offsetQ + currSeed->len == q.length() || (!i->getPredecessors().hasPredecessors() && currSeed->offsetU + currSeed->len == i->size)){
				//Extend hit to the left
				startLeftX_Drop_OnRevComp(&newHit, q, X, quorum, searchSet);

				//Check whether we have created a hit for this query position already
				if(hitArr[newHit.offQ].length != 0){
					//Get the first hit at this position
					hitIt = &hitArr[newHit.offQ];
					//Go through the list of all hits at this position
					while(hitIt != NULL){
						//Check if the current hit and the last extended hit might be cooptimal solutions (i.e. they have the same q offsets and score)
						if(hitIt->length == newHit.length && hitIt->score == newHit.score){
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
						//Add hit to the hit array//

						//Save pointer to the hit list at the current q position
						newHit.nextHit = hitArr[newHit.offQ].nextHit;
						//Create a new hit
						hitArr[newHit.offQ].nextHit = new Hit(newHit);
						hitArr[newHit.offQ].nextHit->gAlgn.aSeqG = "";
						hitArr[newHit.offQ].nextHit->gAlgn.aSeqQ = "";
					}
				} else{
					//Copy hit into the array
					hitArr[newHit.offQ] = newHit;
				}
			}

			//Delete the extended seed
			free(currSeed);
			//Move on to the next seed
			currSeed = currUni.getData()->getData(currUni)->getSeed(currUni.strand);
		}
	}
}

//This function calculates a banded, semi-global, gapped alignment on a list of results considering a quorum and a search color set and outputs the result if demanded
void calcGappedAlignment(ColoredCDBG<seedlist> &cdbg, list<Hit*> &resList, const string &q, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const double &lambda, const double &C){
	bool isDupl = false;
	uint32_t bandRadius;
	list<Hit*>::const_iterator iter;

	//Go through result list
	for(list<Hit*>::const_iterator it = resList.begin(); it != resList.end(); ++it){
		//Calculate the band width to be used during the gapped extension
		bandRadius = (*it)->length / GAP_RATIO;
		//Calculate gapped extension to the right
		startRightGappedAlignment(*it, q, X, bandRadius, quorum, searchSet);

		//Check whether this is not the first hit in the result list
		if(it != resList.begin()){
			iter = resList.begin();

			//Check for duplicates
			while(it != iter){
				//We assume to have a duplicate if right border offsets and unitigs are identical
				if((*it)->offU == (*iter)->offU && (*it)->offQ == (*iter)->offQ && (*it)->origUni == (*it)->origUni && (*iter)->score == (*it)->score && (*iter)->length == (*it)->length){
					// cerr << "It seems that the hit going spanning from q=" << (*iter)->lSeedQoff << " to q=" << (*iter)->rSeedQoff << " is a duplicate" << endl;
					isDupl = true;
					break;
				}

				++iter;
			}
		}

		//We want to avoid doing left extensions for duplicates
		if(!isDupl){
			startLeftGappedAlignment(*it, q, X, bandRadius, quorum, searchSet);
			//Recalculate e-value
			(*it)->eval = calcEVal((*it)->score, lambda, C, q.length());
		} else{
			//Delete result
			it = resList.erase(it);
			//Make sure we do not miss a result
			--it;
		}

		//Reset duplicate flag
		isDupl = false;
	}
}

//This function performs the actual graph search for a query
void searchQuery(ColoredCDBG<seedlist> &cdbg, const int32_t &kMerLength, const int32_t &minSeedLength, const size_t &numSmers, const uint32_t &quorum, const uint32_t &profileSize, const uint32_t *qProfile, const string &q, const SrchStrd &strand, const UnitigColorMap<seedlist> *uArr, const struct s_mer_pos *posArray, const list<pair<string, size_t>> &searchColors, const int16_t &X, const bool &calcRT, uint16_t nRes, const double &lambda, const double &lambdaGap, const double &C, const double &Cgap, const double &eLim, const bool &colOut, const bool &isSim){
	//Staff we need to measure run times
	auto startTime = std::chrono::system_clock::now();
	auto endTime = std::chrono::system_clock::now();
	std::chrono::duration<double> tDiff;

	//No hit can start behind this point
	uint32_t hitArrSize = q.length() - minSeedLength + 1;
	//The query's reverse complement
	string revQ;
	//Variables needed for the seed extension
	Hit hitArr[hitArrSize];
	//Variables needed for the gapped extension
	list<Hit*> resList;

	/*Searching for seeds*/
	cout << "Searching for seeds" << endl;

	//Detect seeds on the original query
	if(strand != Minus) detectSeeds(kMerLength, minSeedLength, numSmers, quorum, profileSize, qProfile, q, uArr, posArray, hitArr, searchColors, true);

	//Detect seeds on the reverse complementary query
	if(strand != Plus){
		//Calculate the query's reverse complement
		revQ = revComp(q);

		detectSeeds(kMerLength, minSeedLength, numSmers, quorum, profileSize, qProfile, revQ, uArr, posArray, hitArr, searchColors, false);
	}

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
	}

	/*Seed Extension*/
	cout << "Extending seeds" << endl;

	//Extend seeds lying on the reference strand if demanded
	if(strand != Minus) extendRefSeeds(cdbg, q, X, hitArr, quorum, searchColors);

	//Extend seeds lying on the reverse complementary strand if demanded
	if(strand != Plus) extendRevCompSeeds(cdbg, q, X, hitArr, quorum, searchColors);

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
	}

	/*Gapped extension*/
	cout << "Performing gapped extension" << endl;

	//Find best results
	for(uint32_t i = 0; i < hitArrSize; ++i){
		Hit *hitList;

		if(hitArr[i].length == 0){
			hitList = NULL;
		} else{
			hitList	= &hitArr[i];
		}

		//Go through all hits of a certain array position
		while(hitList != NULL){
			//Calculate e-value
			hitList->eval = calcEVal(hitList->score, lambda, C, q.length());

			//Make sure the hit should still be considered
			if(hitList->score != 0 && hitList->eval <= eLim){
				//Check if we can still just add new results to the result list or need to replace worse existing ones
				if(nRes == 0){
					//Try to replace a worse result
					replWorseRes(resList, hitList);
				} else{
					//Insert a new result
					insRes(resList, hitList);
					--nRes;
				}
			}

			hitList = hitList->nextHit;
		}
	}

	//Check if this is a simulation run
	if(isSim){
		//Go through results and report ungapped scores
		for(list<Hit*>::iterator iter = resList.begin(); iter != resList.end(); ++iter) cout << "Score (ungapped): " << (*iter)->score << endl;
	}

	//Calculate gapped alignments
	calcGappedAlignment(cdbg, resList, q, X, quorum, searchColors, lambdaGap, Cgap);

	//Check if this is a simulation run
	if(isSim){
		//Go through results and report gapped scores
		for(list<Hit*>::iterator iter = resList.begin(); iter != resList.end(); ++iter) cout << "Score (gapped): " << (*iter)->score << endl;
		
		//Free memory in hit array
		freeHitArray(hitArr, hitArrSize);
		return;
	}

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
		for(list<Hit*>::const_iterator it = resList.begin(); it != resList.end(); ++it){
			//Check e-value once again
			if((*it)->eval <= eLim){
				//Output alignment
				repAlgn(*it);

				//Output color sets if demanded
				if(colOut) outpColSets(cdbg, *it);
			}
		}
	}

	//Free memory in hit array
	freeHitArray(hitArr, hitArrSize);
}

//This function frees all memory additionally allocated for an hit array
void freeHitArray(Hit *arr, uint32_t &arrLen){
	bool frst = false;
	Hit *h, *delHit;

	for(uint32_t i = 0; i < arrLen; ++i){
		//Check if hit array at current position is empty
		if(arr[i].length == 0){
			h = NULL;
		} else{
			//Get pointer of hit at current position
			h = &arr[i];
			//Mark hit as being the first in the hit array
			frst = true;
		}

		//Iterate over all hits of a certain array position
		while(h != NULL){
			//Check if we are dealing with the first hit at this position
			if(frst){
				//Next hit won't be the first anymore
				frst = false;
				//Move to the next hit
				h = h->nextHit;
			} else{
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