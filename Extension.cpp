#include "Extension.h"
#include "Hit.h"
#include "Search.h"
#include "Seedlist.cpp"

bool report = false;

// //This function searches for the closest seed to reach during an extension and returns it if it exists (otherwise NULL)
// //TODO Each time a seed is searched the seed list is gone through from the beginning again. If the last closest seed was not the first one that was found it might be possible not to start from the beginning again!
// //TODO The idea behind ordering all seeds in a seed list decreasingly by their q offset was to avoid to iterate over the complete seed list during the search for the closest seed to reach. Unfortunately, since new seeds are now inserted into the seed lists during extension this does not work anymore now. However, if we are able to exchange the score of a seed by a bool to indicate whether it has already been extended or not we might be able to incorporate this idea again. For now we have to skip it
// struct seed* searchRightNeighbor(struct seed* sLSeed, const uint32_t &iniExtSeedOffsQ, const uint32_t &extLen, const uint32_t &curUPos, struct seed*& prevSeed){
// 	struct seed *lastSeed = NULL, *nearestSeed = NULL;

// 	//If there is no seed list there is no nearerst seed
// 	if(sLSeed == NULL){ return nearestSeed;	}

// 	//Find the nearest seed that we might be able to reach during our extension
// 	do{
// 		//Since our seed lists are decreasingly ordered by their offset in q we only need to search for a reachable seed as long as the current seed's q offset is not smaller than the q offset of the seed we are just trying to extend
// 		if(iniExtSeedOffsQ >= sLSeed->offsetQ){
// 			return nearestSeed;
// 		} else{
// 			//Check whether we might be able to reach this seed
// 			if(sLSeed->offsetQ == iniExtSeedOffsQ + extLen - (curUPos - sLSeed->offsetU)){
// 				//Store a pointer to the previous seed
// 				prevSeed = lastSeed;
// 				//Store the reachable seed
// 				nearestSeed = sLSeed;

// 				//Testing
// 				//cout << "We should get here twice" << endl;
// 			}
// 		}

// 		lastSeed = sLSeed;
// 		sLSeed = sLSeed->nextSeed;
// 	} while(sLSeed != NULL);

// 	return nearestSeed;
// }

//This function goes through an extension list and incorporates all seeds that are not already part of a seed list into their corresponding seed lists
void processExtPtr(struct Ext_ptr*& extPtr){
	//Check whether there is an extension seed and whether it has to be incorporated into a seed list
	if(extPtr != NULL){
		//We do not want to have seeds in the processed seed list
		delExtPtr(extPtr);
		return; //We are done here

		// //Make currently first element of the respective seed list the successor of the new seed to incorporate
		// extPtr->extSeed->nextSeed = extPtr->origUnitig.getData()->getData(extPtr->origUnitig)->getLastProcSeed();

		// //Make the new seed to incorporate the corresponding seed list's first element
		// extPtr->origUnitig.getData()->getData(extPtr->origUnitig)->setLastProcSeed(extPtr->extSeed);
	}
}

// //This function initiates the extension on all successors of a unitig and returns the best one considering a quorum
// int32_t extendAtNextUnitig(const ForwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &iniUniPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum){
// 	uint16_t sucID;
// 	uint32_t tmpHitLen = 0;
// 	int32_t maxScore = 0, currScore;

// 	// struct Ext_ptr *lExtSeed = NULL;
// 	//ForwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> sucIter = lastUni.getSuccessors();

// 	//Check whether we have reached the maximum recursion depth of an extension
// 	if(++explCount > MAXRECURSIONDEPTH){
// 		//Report this incident//TODO: If this happens the condition that an extended seed cannot be reached anymore is violated. Implement a procedure after the extension that merges such seeds with other seeds that would have reached them!
// 		//cerr << "Maximum recursion depth reached during extension. Position in q: " << iniQoff << endl;
// 		//Terminate this extension
// 		return 0;
// 	}

// 	//Testing
// 	// if(report){
// 	// 	cerr << "extendAtNextUnitig" << endl;
// 	// }

// 	maxScore = 0;
// 	sucID = 0;

// 	//Iterate over successors
// 	for(neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> nI = sucIter.begin(); nI != sucIter.end(); ++nI){
// 		// //Temporary extension pointer to store the current extensions border seed and unitig
// 		// struct Ext_ptr *tmpExtSeed = NULL;

// 		//Temporary extention path
// 		list<uint16_t> tmpPth;
// 		//Note which successor we are on
// 		++sucID;

// 		//Calculate the score of an extension of a successor
// 		currScore = contRightX_Drop(nI, iniQoff, tmpHitLen, extLen, q, X, lastExtSeedTmpScore, iniUniPos, tmpPth, explCount, quorum);

// 		//Check whether the score of the current successors extension is the best one found so far
// 		if(currScore > maxScore){
// 			//Update maxScore
// 			maxScore = currScore;
// 			//Update hit length
// 			hitLen = tmpHitLen;

// 			// //Check whether we saved an extension list before already
// 			// if(lExtSeed != NULL){
// 			// 	//Delete old extension list
// 			// 	delExtPtr(lExtSeed);
// 			// }

// 			extPth = tmpPth;
// 			//Save which successor we have chosen
// 			extPth.push_front(sucID);
			
// 			// //Save the new extension list
// 			// lExtSeed = tmpExtSeed;
// 		}// else if(tmpExtSeed != NULL){
// 		// 	//Delete the new extension list
// 		// 	delExtPtr(tmpExtSeed);
// 		// }
// 	}

// 	// //Check whether we could find a good extension
// 	// if(maxScore > 0){ 
// 	// 	//if(record) cerr << "Successive unitigs visited" << endl;
// 	// 		return lExtSeed;
// 	// }

// 	// }

// 	//Testing
// 	//if(record) cerr << "Successive unitigs visited" << endl;

// 	//Nothing found
// 	return maxScore;
// }

// //This function initiates the extension on all successors of a unitig and returns the best one considering a quorum. This function is explicitly designed for seeds lying on the query's reverse complement (considering the overlap between unitigs in sequences' beginning)
// struct Ext_ptr* extendAtNextUnitig_OnRevComp(const ForwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &iniUniPos, uint32_t &explCount, const uint32_t &quorum){
// 	uint32_t tmpHitLen = 0;
// 	int32_t maxScore = 0, currScore;
// 	struct Ext_ptr *lExtSeed = NULL;

// 	//Check whether we have reached the maximum recursion depth of an extension
// 	if(++explCount > MAXRECURSIONDEPTH){
// 		//Terminate this extension
// 		return NULL;
// 	}

// 	//Iterate over successors
// 	for(neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> nI = sucIter.begin(); nI != sucIter.end(); ++nI){
// 		//Temporary extension pointer to store the current extensions border seed and unitig
// 		struct Ext_ptr *tmpExtSeed = NULL;

// 		//Calculate the score of an extension of a successor
// 		currScore = contRightX_Drop_OnRevComp(nI, iniQoff, tmpHitLen, extLen, q, X, lastExtSeedTmpScore, iniUniPos, explCount, quorum);

// 		//Check whether the score of the current successors extension is the best one found so far
// 		if(currScore > maxScore){
// 			//Update maxScore
// 			maxScore = currScore;
// 			//Update hit length
// 			hitLen = tmpHitLen;
// 			//Check whether we saved an extension list before already
// 			if(lExtSeed != NULL){
// 				//Delete old extension list
// 				delExtPtr(lExtSeed);
// 			}
// 			//Save the new extension list
// 			lExtSeed = tmpExtSeed;
// 		} else if(tmpExtSeed != NULL){
// 			//Delete the new extension list
// 			delExtPtr(tmpExtSeed);
// 		}
// 	}

// 	//Check whether we could find a good extension
// 	if(maxScore > 0){
// 		return lExtSeed;
// 	}

// 	//Nothing found
// 	return NULL;
// }

//This function initiates the extension on all successors of a unitig and returns the best one considering a quorum and a search color set
int32_t extendAtNextUnitig(const ForwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &uniPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	uint16_t sucID;
	uint32_t tmpHitLen = 0;
	int32_t maxScore = 0, currScore;
	
	// struct Ext_ptr *lExtSeed = NULL;

	//ForwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> sucIter = lastUni.getSuccessors();

	//Check whether we have reached the maximum recursion depth of an extension
	if(++explCount > MAXRECURSIONDEPTH){
		//Report this incident//TODO: If this happens the condition that an extended seed cannot be reached anymore is violated. Implement a procedure after the extension that merges such seeds with other seeds that would have reached them!
		//cerr << "Maximum recursion depth reached during extension. Position in q: " << iniQoff << endl;
		//Terminate this extension
		return 0;
	}

	//Testing
	// if(record){
	// 	cerr << "Traverse successive unitigs" << endl;
	// }

	maxScore = 0;
	sucID = 0;

	//Iterate over successors
	for(neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> nI = sucIter.begin(); nI != sucIter.end(); ++nI){
		// //Temporary extension pointer to store the current extensions border seed and unitig
		// struct Ext_ptr *tmpExtSeed = NULL;

		//Temporary extention path
		list<uint16_t> tmpPth;
		//Note which successor we are on
		++sucID;

		//Testing
		// cout << "Unitig is " << nI->mappedSequenceToString() << endl;

		//Calculate the score of an extension of a successor
		currScore = contRightX_Drop(nI, iniQoff, tmpHitLen, extLen, q, X, lastExtSeedTmpScore, uniPos, tmpPth, explCount, quorum, searchSet);

		//Check whether the score of the current successors extension is the best one found so far
		if(currScore > maxScore){
			//Testing
			//if(record) cerr << "Is new maximum" << endl;
			// cout << "New maximum found" << endl;

			//Update maxScore
			maxScore = currScore;
			//Update hit length (tmpHitLen does not have to be reseted, because it is only set to and never read from except here)
			hitLen = tmpHitLen;

			// //Check whether we saved an extension list before already
			// if(lExtSeed != NULL){
			// 	//Delete old extension list
			// 	delExtPtr(lExtSeed);
			// }

			extPth = tmpPth;
			//Save which successor we have chosen
			extPth.push_front(sucID);

			// //Save the new extension list
			// lExtSeed = tmpExtSeed;
		}// else if(tmpExtSeed != NULL){
		// 	//Delete the new extension list
		// 	delExtPtr(tmpExtSeed);
		// }
	}

	//Testing
	// if(record){
	// 	//cerr << "Successive unitigs visited" << endl;
	// 	// if(maxScore > 0){
	// 	// 	cerr << "Best unitig was " << lExtSeed->origUnitig.mappedSequenceToString() << endl;
	// 	// 	cerr << "Extension pointer is " << explCount << endl;
	// 	// 	cerr << "Score is " << maxScore << endl;
	// 	// 	cerr << "Right extension seed: offsetU:" << lExtSeed->extSeed->offsetU << " offsetQ:" << lExtSeed->extSeed->offsetQ << " len:" << lExtSeed->extSeed->len << " score:" << lExtSeed->extSeed->score << endl;
	// 	// } else{
	// 	// 	cerr << "No good extension found" << "Extension pointer is " << explCount << endl;
	// 	// }
	// }
		
	// //Check whether we could find a good extension
	// if(maxScore > 0){ return lExtSeed; }
	
	//}

	//Testing
	//if(record) cerr << "Successive unitigs visited" << endl;

	//Nothing found
	return maxScore;
}

//This function  initiates the extension on all successors of a unitig and returns the best one considering a quorum and a search color set. This function is explicitly designed for seeds lying on the query's reverse complement (considering the overlap between unitigs in sequences' beginning)
int32_t extendAtNextUnitig_OnRevComp(const ForwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const int16_t &X, const int32_t &lastExtSeedTmpScore, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	uint16_t sucID;
	uint32_t tmpHitLen = 0;
	int32_t maxScore = 0, currScore;
	
	// struct Ext_ptr *lExtSeed = NULL;

	//Check whether we have reached the maximum recursion depth of an extension
	if(++explCount > MAXRECURSIONDEPTH){
		//Testing
		// cout << "1 Option 1" << endl;

		//Terminate this extension
		return 0;
	}

	//Testing
	// cout << "1 Option 2" << endl;

	maxScore = 0;
	sucID = 0;

	//Iterate over successors
	for(neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> nI = sucIter.begin(); nI != sucIter.end(); ++nI){
		// //Temporary extension pointer to store the current extensions border seed and unitig
		// struct Ext_ptr *tmpExtSeed = NULL;

		//Temporary extention path
		list<uint16_t> tmpPth;
		//Note which successor we are on
		++sucID;
		//Calculate the score of an extension of a successor
		currScore = contRightX_Drop_OnRevComp(nI, iniQoff, tmpHitLen, extLen, q, X, lastExtSeedTmpScore, tmpPth, explCount, quorum, searchSet);

		//Check whether the score of the current successors extension is the best one found so far
		if(currScore > maxScore){
			//Testing
			// cout << "2 Option 1" << endl;
			// cout << "3 Option ";

			//Update maxScore
			maxScore = currScore;
			//Update hit length
			hitLen = tmpHitLen;

			// //Check whether we saved an extension list before already
			// if(lExtSeed != NULL){
			// 	//Testing
			// 	// cout << "not ";

			// 	//Delete old extension list
			// 	delExtPtr(lExtSeed);
			// }

			//Testing
			// cout << "1" << endl;

			extPth = tmpPth;
			//Save which successor we have chosen
			extPth.push_front(sucID);

			// //Save the new extension list
			// lExtSeed = tmpExtSeed;
		}// else if(tmpExtSeed != NULL){
		// 	//Testing
		// 	// cout << "4" << endl;

		// 	//Delete the new extension list
		// 	delExtPtr(tmpExtSeed);
		// }

		//Testing
		// if(currScore <= maxScore){
		// 	cout << "2 Option 2" << endl;
		// }
		// cout << "We are on unitig " << nI->mappedSequenceToString() << endl;
	}
		
	// //Check whether we could find a good extension
	// if(maxScore > 0){
	// 	//Testing
	// 	// cout << "5 Option 1" << endl;

	// 	return lExtSeed;
	// }

	//Testing
	// cout << "5 Option 2" << endl;

	//Nothing found
	return maxScore;
}

// //The good old X-drop algorithm (extension to the right) for seeds matching the query's reference strand considering quorum. Returns an extension pointer storing the extension path through the graph/*NOTE: For now the seed extension is a little bit redundant. Would be cool if this can be changed in the future! - Is the seed extension still redundant?*/
// struct ExtPth startRightX_Drop(hit* hit, const string &q, const int16_t &X, const uint32_t &quorum){
// 	//Initialization of auxiliary variables
// 	int32_t tmpScore = 0;
// 	uint32_t tmpSeedLen = hit->length;
// 	uint32_t endBuf = 0;
// 	uint32_t iniUniPos;
// 	uint32_t checkedPos = 0;
// 	size_t offset;
// 	//Counter to count tries to explore a further unitig
// 	uint32_t explCount

// 	string uSeq = hit->origUni.mappedSequenceToString();
// 	list<uint16_t> extPth;

// 	// //struct seed *nearestSeed, *prevSeed = NULL;
// 	// struct Ext_ptr *extPtr = NULL;
// 	// //Make a copy of our unitig
// 	// UnitigColorMap<seedlist> cpyUnitig = hit->rUnitig;

// 	//Calculate hit's initial score
// 	hit->score = hit->length;//TODO If we want to use an arbitrary cost function this has to be changed first!

// 	//Check how far we should explore the unitig's sequence
// 	if(hit->origUni.getSuccessors().hasSuccessors()){
// 		endBuf = hit->origUni.getGraph()->getK() - 1;
// 	}

// 	//We are done if we have reached the end of the query
// 	while(hit->offQ + tmpSeedLen < q.length()){//TODO Perhaps it would be better to outsource the X-drop Alg.!
// 		// //Check whether we have reached the next seed
// 		// if(nearestSeed != NULL && hit->rSeedQoff + tmpSeedLen == nearestSeed->offsetQ){
// 		// 	//Update temporary seed length
// 		// 	tmpSeedLen += nearestSeed->length;
// 		// 	//Update number of positions checked for quorum already
// 		// 	if(checkedPos >= nearestSeed->length - 1){
// 		// 		checkedPos -= (nearestSeed->length - 1);
// 		// 	} else{
// 		// 		checkedPos = 0;
// 		// 	}

// 		// 	//Check if the reached seed leads to a positive temporary score
// 		// 	if((tmpScore += nearestSeed->length) > 0){//TODO This needs to be changed as soon as we want to use a non-unit score!
// 		// 		//Update hit's length
// 		// 		hit->length = tmpSeedLen;
// 		// 		//Update hit's score
// 		// 		hit->score += tmpScore;
// 		// 		//Reset temporary score
// 		// 		tmpScore= 0;
// 		// 		//Delete the reached seed
// 		// 		if(extractSeed(cpyUnitig->getData()->getData(cpyUnitig)->getLastSeed(cpyUnitig.strand), nearestSeed, prevSeed)){<-nearestSeed should be NULL after this
// 		// 			//Search for the next neighbor
// 		// 			nearestSeed = searchRightNeighbor(cpyUnitig, hit->rSeedQoff, tmpSeedLen, hit->rSeedUoff, prevSeed);
// 		// 		}
// 		// 	}

// 		// //Check whether we have reached the end of the unitig's sequence
// 		// if(hit->rSeedUoff + tmpSeedLen < hit->rUnitig.size - endBuf){
// 	// while(hit->offQ + tmpSeedLen < q.length()){//TODO Perhaps it would be better to outsource the X-drop Alg.!
// 		// //Check how far we should explore the unitig's sequence
// 		// if(hit->origUni.getSuccessors().hasSuccessors()){
// 		// 	endBuf = hit->origUni.getGraph()->getK() - 1;
// 		// }

// 		//Check whether we have reached the end of the unitig's sequence
// 		if(hit->offU + tmpSeedLen < hit->origUni.size - endBuf){
// 			//Check whether quorum has to be checked
// 			if(checkedPos == 0){
// 				//Calculate offset
// 				offset = hit->offU + tmpSeedLen;

// 				//Check whether quorum is still fulfilled
// 				if((checkedPos = checkSearchCrit(hit->origUni, quorum, offset)) == 0) break;
// 			}

// 			//Testing
// 			//if(report) cout << "We do extend first on the unitig of our of interest" << endl;
// 			// cout << "Compare next bases" << endl << "tmpSeedLen:" << tmpSeedLen << endl << "Unitig length:" << cpyUnitig.size << endl << "tmpScore:" << tmpScore << endl << "hit->rSeedUoff:" << hit->rSeedUoff << endl;
// 			// cout << "Bases compared: " << uSeq[hit->rSeedUoff + tmpSeedLen] << " " << q[hit->rSeedQoff + tmpSeedLen] << endl;

// 			//Check whether the score of our extension is positive
// 			if((tmpScore += compUScore(uSeq[hit->offU + tmpSeedLen], q[hit->offQ + tmpSeedLen])) > 0){//TODO Change to calculate arbitrary score functions
// 				//Update the seed info
// 				hit->score += tmpScore;

// 				hit->length = tmpSeedLen + 1;//+1 because we haven't increased tmpSeedLen yet
// 				//Reset the temporary score
// 				tmpScore = 0;
// 			} else{
// 				//Check if the current extension is already too bad
// 				if(tmpScore < -X){ break; }
// 			}
// 			//Proceed with the next two positions
// 			++tmpSeedLen;
// 			--checkedPos;
// 		} else{
// 			//Check if the current unitig has successors
// 			if(endBuf != 0){
// 				//Calculate the unitig sequence position we have to start with in the successive unitig
// 				iniUniPos = hit->offU + tmpSeedLen - hit->origUni.size + endBuf;
// 				//Initialize explCount
// 				explCount = 0;
// 				//Explore unitig's successors
// 				hit->score += extendAtNextUnitig(hit->origUni.getSuccessors(), hit->offQ, hit->length, tmpSeedLen, q, X, tmpScore, iniUniPos, extPth, explCount, quorum);
// 			}

// 			//Testing
// 			//if(report) cout << "We have returned to startRightX_Drop and our new q offset is " << extPtr->extSeed->offsetQ << endl;

// 			break;
// 		}
// 	}

// 	// //Check whether our extension has been continued on successive unitigs
// 	// if(extPtr != NULL){
// 	// 	//Update hits score
// 	// 	hit->score += extPtr->extSeed->score;

// 	// 	//Recalculate border seed's offsets
// 	// 	extPtr->extSeed->offsetU += extPtr->extSeed->len - 1;

// 	// 	extPtr->extSeed->offsetQ += extPtr->extSeed->len - 1;

// 	// 	//Exchange hit's right border offsets
// 	// 	hit->offU = extPtr->extSeed->offsetU;
// 	// 	hit->offQ = extPtr->extSeed->offsetQ;

// 	// 	//Save link to the last reached unitig
// 	// 	hit->rUnitig = extPtr->origUnitig;
// 	// 	//Delete the extension pointer
// 	// 	delExtPtr(extPtr);
// 	// } else{
// 	// 	//Recalculate border seed's offsets
// 	// 	hit->rSeedUoff += hit->length - 1;
// 	// 	hit->rSeedQoff += hit->length - 1;
// 	// }

// 	// //Incorporate seed from extension ptr into seed list
// 	// processExtPtr(extPtr);

// 	return cmprExtPth(extPth);
// }

// //The good old X-drop algorithm (extension to the right) for seeds matching the query's reverse complement considering quorum/*NOTE: For now the seed extension is a little bit redundant. Would be cool if this can be changed in the future! - Is the seed extension still redundant?*/
// void startRightX_Drop_OnRevComp(struct hit* hit, const string &q, const int16_t &X, const uint32_t &quorum){
// 	//Initialization of auxiliary variables
// 	int32_t tmpScore = 0;
// 	uint32_t tmpSeedLen = hit->length;
// 	uint32_t iniUniPos;
// 	uint32_t checkedPos = 0;
// 	size_t offset;
// 	//Counter to count tries to explore a further unitig
// 	uint32_t explCount;
// 	string uSeq = hit->rUnitig.mappedSequenceToString();
// 	//struct seed *nearestSeed, *prevSeed = NULL;
// 	struct Ext_ptr *extPtr = NULL;
// 	//Make a copy of our unitig
// 	UnitigColorMap<seedlist> cpyUnitig = hit->rUnitig;

// 	//Calculate hit's initial score
// 	hit->score = hit->length;//TODO If we want to use an arbitrary cost function this has to be changed first

// 	//We are done if we have reached the end of the query
// 	while(hit->rSeedQoff + tmpSeedLen < q.length()){//TODO Perhaps it would be better to outsource the X-drop Alg.!
// 		//Check whether we have reached the end of the unitig's sequence
// 		if(hit->rSeedUoff + tmpSeedLen < cpyUnitig.size){
// 			//Check whether quorum has to be checked
// 			if(checkedPos == 0){
// 				//Calculate offset
// 				offset = hit->rSeedUoff + tmpSeedLen;

// 				//Check whether quorum is still fulfilled
// 				if((checkedPos = checkSearchCrit(cpyUnitig, quorum, offset)) == 0) break;
// 			}

// 			//Check whether the score of our extension is positive
// 			if((tmpScore += compUScore(uSeq[hit->rSeedUoff + tmpSeedLen], q[hit->rSeedQoff + tmpSeedLen])) > 0){//TODO Change to calculate arbitrary score functions
// 				//Update the seed info
// 				hit->score += tmpScore;

// 				hit->length = tmpSeedLen + 1;//+1 because we haven't increased tmpSeedLen yet
// 				//Reset the temporary score
// 				tmpScore = 0;
// 			} else{
// 				//Check if the current extension is already too bad
// 				if(tmpScore < -X){ break; }
// 			}
// 			//Proceed with the next two positions
// 			++tmpSeedLen;
// 			--checkedPos;
// 		} else{
// 			//Check if the current unitig has successors
// 			if(cpyUnitig.getSuccessors().hasSuccessors()){
// 				//Calculate the unitig sequence position we have to start with in the successive unitig
// 				iniUniPos = cpyUnitig.getGraph()->getK() - 1;
// 				//Initialize explCount
// 				explCount = 0;
// 				//Explore unitig's successors
// 				extPtr = extendAtNextUnitig_OnRevComp(cpyUnitig.getSuccessors(), hit->rSeedQoff, hit->length, tmpSeedLen, q, X, tmpScore, iniUniPos, explCount, quorum);
// 			}
			
// 			break;
// 		}
// 	}

// 	//Check whether our extension has been continued on successive unitigs
// 	if(extPtr != NULL){
// 		//Update hits score
// 		hit->score += extPtr->extSeed->score;
// 		//Recalculate border seed's offsets
// 		extPtr->extSeed->offsetU += extPtr->extSeed->len - 1;
// 		extPtr->extSeed->offsetQ += extPtr->extSeed->len - 1;
// 		//Exchange hit's right border offsets
// 		hit->rSeedUoff = extPtr->extSeed->offsetU;
// 		hit->rSeedQoff = extPtr->extSeed->offsetQ;
// 		//Save link to the last reached unitig
// 		hit->rUnitig = extPtr->origUnitig;
// 		//Delete the extension pointer
// 		delExtPtr(extPtr);
// 	} else{
// 		//Recalculate border seed's offsets
// 		hit->rSeedUoff += hit->length - 1;
// 		hit->rSeedQoff += hit->length - 1;
// 	}
// }

//The good old X-drop algorithm (extension to the right) for seeds matching the query's reference strand considering quorum and search color set. Returns an extension pointer storing the extension path through the graph/*NOTE: For now the seed extension is a little bit redundant. Would be cool if this can be changed in the future! - Is the seed extension still redundant?*/
void startRightX_Drop(hit* hit, const string &q, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	//Initialization of auxiliary variables
	int32_t tmpScore = 0;
	uint32_t tmpSeedLen = hit->length;
	uint32_t overlap = 0;
	uint32_t iniUniPos;
	int32_t checkedPos = getSrchCritCov(hit->origUni, quorum, searchSet, hit->length, true);
	//Counter to count tries to explore a further unitig
	uint32_t explCount;
	string uSeq = hit->origUni.mappedSequenceToString();
	list<uint16_t> extPth;

	// struct Ext_ptr *extPtr = NULL;

	//Testing
	// cout << "checkedPos: " << checkedPos << endl;
	
	// //Make a copy of our unitig
	// UnitigColorMap<seedlist> cpyUnitig = hit->rUnitig;

	//Calculate hit's initial score
	hit->score = hit->length;//TODO If we want to use an arbitrary cost function this has to be changed first!

	//Check how far we should explore the unitig's sequence
	if(hit->origUni.getSuccessors().hasSuccessors()) overlap = hit->origUni.getGraph()->getK() - 1;

	//We are done if we have reached the end of the query
	while(hit->offQ + tmpSeedLen < q.length()){//TODO Perhaps it would be better to outsource the X-drop Alg.!
		//Check whether we have reached the end of the unitig's sequence
		if(hit->offU + tmpSeedLen < hit->origUni.size - overlap){
			//Ensure that search criteria are still fulfilled
			if(checkedPos == 0){
				break;

				// //Calculate offset
				// offset = hit->offU + tmpSeedLen;

				//Testing
				// cout << "We need to check the quorum uPos " << hit->offU + tmpSeedLen << endl;

				// //Check whether quorum is still fulfilled
				// if((checkedPos = checkSearchCrit(hit->origUni, quorum, offset, searchSet)) == 0) break;

				//Testing
				// cout << "Done" << endl;
			}

			//Testing
			// cout << "Bases to compare:" << uSeq[hit->offU + tmpSeedLen] << " " << q[hit->offQ + tmpSeedLen] << endl;
			// cout << "Base comparisons leads to a";

			//Check whether the score of our extension is positive
			if((tmpScore += compUScore(uSeq[hit->offU + tmpSeedLen], q[hit->offQ + tmpSeedLen])) > 0){//TODO Change to calculate arbitrary score functions
				//Update the seed info
				hit->score += tmpScore;

				hit->length = tmpSeedLen + 1;//+1 because we haven't increased tmpSLen yet
				//Reset the temporary score
				tmpScore = 0;
			} else{
				//Check if the current extension is already too bad
				if(tmpScore < -X){
					//Testing
					//cout << " too negative score" << endl;

					break;
				}

				//Testing
				//cout << " negative score" << endl;
			}
			//Proceed with the next two positions
			++tmpSeedLen;
			--checkedPos;
		} else{
			//Check if the current unitig has successors
			if(overlap != 0){
				//Calculate the unitig sequence position we have to start with in the successive unitig
				iniUniPos = hit->offU + tmpSeedLen - hit->origUni.size + overlap;
				//Initialize explCount
				explCount = 0;
				//Explore unitig's successors
				hit->score += extendAtNextUnitig(hit->origUni.getSuccessors(), hit->offQ, hit->length, tmpSeedLen, q, X, tmpScore, iniUniPos, extPth, explCount, quorum, searchSet);
			}

			break;
		}
	}

	// //Check whether our extension has been continued on successive unitigs
	// if(extPtr != NULL){
	// 	//Update hits score
	// 	hit->score += extPtr->extSeed->score;

	// 	//Recalculate border seed's offsets
	// 	extPtr->extSeed->offsetU += extPtr->extSeed->len - 1;

	// 	extPtr->extSeed->offsetQ += extPtr->extSeed->len - 1;

	// 	//Exchange hit's right border offsets
	// 	hit->offU = extPtr->extSeed->offsetU;
	// 	hit->offQ = extPtr->extSeed->offsetQ;
	// } else{
	// 	//Recalculate border seed's offsets
	// 	hit->offU += hit->length - 1;
	// 	hit->offQ += hit->length - 1;
	// }

	//Testing
	// if(uSeq == "GCAAGAGCCGCTGTTTCTTGAACAATATCTCG"){
	// 	cout << "Inside startRightX_Drop: Final extension path:" << endl;
	// 	for(list<uint16_t>::const_iterator i = extPth.begin(); i != extPth.end(); ++i) cout << *i << " ";
	// 	cout << endl;
	// }

	hit->rExt = cmprExtPth(extPth);

	//Testing
	// if(uSeq == "GCAAGAGCCGCTGTTTCTTGAACAATATCTCG"){
	// 	extPth = list<uint16_t>();
	// 	extPth = decmprExtPth(hit->rExt);
	// 	cout << "Inside startRightX_Drop: Decompressed extension path:" << endl;
	// 	for(list<uint16_t>::const_iterator i = extPth.begin(); i != extPth.end(); ++i) cout << *i << " ";
	// 	cout << endl;
	// 	exit(0);
	// }
}

//The good old X-drop algorithm (extension to the right) for seeds matching the query's reverse complement considering quorum and search color set. Returns an extension pointer storing the extension path through the graph/*NOTE: For now the seed extension is a little bit redundant. Would be cool if this can be changed in the future! - Is the seed extension still redundant?*/
void startRightX_Drop_OnRevComp(hit* hit, const string &q, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	//Initialization of auxiliary variables
	int32_t tmpScore = 0;
	uint32_t tmpSeedLen = hit->length;

	// uint32_t iniUniPos;

	//Get the number of checked positions considering that we do not start in the sequence's very end
	int32_t checkedPos = getSrchCritCov(hit->origUni, quorum, searchSet, compOffset(hit->offU + hit->length, 1, hit->origUni.size, false), false);// - (hit->origUni.size - compOffset(hit->offU, hit->length, hit->origUni.size, false));

	// size_t offset;

	//Counter to count tries to explore a further unitig
	uint32_t explCount;
	string uSeq = hit->origUni.mappedSequenceToString();
	list<uint16_t> extPth;

	// struct Ext_ptr *extPtr = NULL;

	// //Make a copy of our unitig
	// UnitigColorMap<seedlist> cpyUnitig = hit->rUnitig;

	//Calculate hit's initial score
	hit->score = hit->length;//TODO If we want to use an arbitrary cost function this has to be changed first!

	//We are done if we have reached the end of the query
	while(hit->offQ + tmpSeedLen < q.length()){//TODO Perhaps it would be better to outsource the X-drop Alg.!
		//Testing
		// cout << "1 Option 2" << endl;

		//Check whether we have reached the end of the unitig's sequence
		if(hit->offU + tmpSeedLen < hit->origUni.size){
			//Testing
			// cout << "2 Option 2" << endl;
			// cout << "3 Option ";

			//Check whether the current position still fulfills the search criteria
			if(checkedPos == 0){
				//Testing
				// cout << "not ";

				// //Check if we are already in the overlap at the end of the unitig sequence
				// if(hit->offU + tmpSeedLen > hit->origUni.size - hit->origUni.getGraph()->getK()){
				// 	//We check the first k-mer on the reference strand
				// 	offset = 0;
				// } else{
				// 	//Move k-1 positions to the right and calculate the corresponding offset on the reference strand
				// 	offset = compOffset(hit->offU + tmpSeedLen, hit->origUni.getGraph()->getK(), hit->origUni.size, false);
				// }

				// // //Calculate offset
				// // offset = compOffset(hit->rSeedUoff + tmpSeedLen, cpyUnitig.getGraph()->getK(), cpyUnitig.size, false);

				// //Check whether quorum is still fulfilled
				// if((checkedPos = checkSearchCrit(hit->origUni, quorum, offset, searchSet)) == 0){
					//Testing
					// cout << "2" << endl;
					// cout << "4 Option 2" << endl;

				break;

				// }

				//Testing
				// cout << "4 Option 1" << endl;
			}

			//Testing
			// cout << "2" << endl;

			//Check whether the score of our extension is positive
			if((tmpScore += compUScore(uSeq[hit->offU + tmpSeedLen], q[hit->offQ + tmpSeedLen])) > 0){//TODO Change to calculate arbitrary score functions
				//Testing
				// cout << "5 Option 1" << endl;

				//Update the seed info
				hit->score += tmpScore;
				hit->length = tmpSeedLen + 1;//+1 because we haven't increased tmpSLen yet
				//Reset the temporary score
				tmpScore = 0;
			} else{
				//Check if the current extension is already too bad
				if(tmpScore < -X){
					//Testing
					// cout << "5 Option 3" << endl;

					break;
				}

				//Testing
				// cout << "5 Option 2" << endl;
			}
			//Proceed with the next two positions
			++tmpSeedLen;
			--checkedPos;
		} else{
			//Testing
			// cout << "2 Option 1" << endl;

			//Check if the current unitig has successors
			if(hit->origUni.getSuccessors().hasSuccessors()){
				//Testing
				// cout << "6 Option 1" << endl;

				//Initialize explCount
				explCount = 0;
				//Explore unitig's successors
				hit->score += extendAtNextUnitig_OnRevComp(hit->origUni.getSuccessors(), hit->offQ, hit->length, tmpSeedLen, q, X, tmpScore, extPth, explCount, quorum, searchSet);
			}// else{
				//Testing
				// cout << "6 Option 2" << endl;
			//}

			break;
		}
	}

	//Testing
	// cout << "1 Option 1" << endl;
	if(hit->origUni.mappedSequenceToString() == "CAGTACGGTATCGGCCCCCAACGCAATCATGCGCACAACGTCCAGACCGTTACGGATCCCGCTGTCTGCCAGAATGGTGATGTCGCCTTTCACCGCATCGGCAATGGCGGG"){
		cerr << "After right extension: length: " << hit->length << " score: " << hit->score << " right extension path: " << endl;
		for(list<uint16_t>::const_iterator i = extPth.begin(); i != extPth.end(); ++i){
			cerr << *i << " ";
		}
		cerr << endl;
	}

	// //Check whether our extension has been continued on successive unitigs
	// if(extPtr != NULL){
	// 	//Testing
	// 	// cout << "7 Option 1" << endl;

	// 	//Update hits score
	// 	hit->score += extPtr->extSeed->score;
	// 	//Recalculate border seed's offsets
	// 	extPtr->extSeed->offsetU += extPtr->extSeed->len - 1;
	// 	extPtr->extSeed->offsetQ += extPtr->extSeed->len - 1;
	// 	//Exchange hit's right border offsets
	// 	hit->rSeedUoff = extPtr->extSeed->offsetU;
	// 	hit->rSeedQoff = extPtr->extSeed->offsetQ;
	// 	//Save link to the last reached unitig
	// 	hit->rUnitig = extPtr->origUnitig;
	// 	//Delete the extension pointer
	// 	delExtPtr(extPtr);
	// } else{
	// 	//Testing
	// 	// cout << "7 Option 2" << endl;

	// 	//Recalculate border seed's offsets
	// 	hit->rSeedUoff += hit->length - 1;
	// 	hit->rSeedQoff += hit->length - 1;
	// }

	hit->rExt = cmprExtPth(extPth);
}

// //This function continues an extension to the right on a successive unitig of a seed lying on the query's reference strand considering a quorum. Returns the maximum score reached.
// int32_t contRightX_Drop(const neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> &sucUnitig, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t &extLen, const string &q, const int16_t &X, const int32_t &lastSeedTmpScore, uint32_t uniSeqPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum){
// 	int32_t tmpScore, progress, score = 0;
// 	uint32_t endBuf = 0;
// 	uint32_t iniSeqPos;
// 	uint32_t tmpSLen;
// 	uint32_t checkedPos = 0;
// 	size_t offset;

// 	//bool inSeedList = false;

// 	string sucUniSeq = sucUnitig->mappedSequenceToString();
// 	struct seed *nearestSeed, *prevSeed;//, *extSeed = NULL;

// 	//Save the initial offset in the current unitig which we need for all nearest neighbor calculations
// 	iniSeqPos = uniSeqPos;

// 	// //Check whether we are working on the reference strand of the curent unitig
// 	// if(sucUnitig->strand){

// 	//Testing
// 	// if(report){
// 	// 	cerr << "contRightX_Drop" << endl;
// 	// }

// 	//Find the nearest seed that we might be able to reach during our extension
// 	nearestSeed = searchRightNeighbor(sucUnitig->getData()->getData(*sucUnitig)->getSeed(sucUnitig->strand), iniQoff, extLen, iniSeqPos, prevSeed);

// 	// } else{
// 	// 	//So far we have not considered the reverse complementary strand for the seed detection //TODO: Change this!
// 	// 	nearestSeed = NULL;
// 	// }

// 	//Perform the X-drop algorithm on the successive unitig
// 	tmpScore = lastSeedTmpScore;
// 	tmpSLen = 0;

// 	//We are done if we have reached the end of the query
// 	while(iniQoff + extLen + tmpSLen < q.length()){
// 		//Testing
// 		// if(nearestSeed != NULL && nearestSeed->offsetU == 52 && nearestSeed->offsetQ == 5155 && nearestSeed->len == 17){
// 		// 	cout << "The seed has been reached before in contRightX_Drop" << endl;
// 		// 	exit(0);
// 		// }

// 		//Check whether we have reached the next seed
// 		if(nearestSeed != NULL && iniQoff + extLen + tmpSLen >= nearestSeed->offsetQ){
// 			//Calculate the gain we get by incorporating the reached seed
// 			progress = nearestSeed->offsetQ + nearestSeed->len - (iniQoff + extLen + tmpSLen);

// 			//Testing
// 			// if(sucUniSeq == "CGAGCGCCTGTACTGATCCACTCGACGCCTACGCGTCAGTGATCGCGCCTCAGCCCCGCCCCCGCC"){
// 			// 	cerr << "We have reached our interesting seed during right extension" << endl;
// 			// 	record = true;
// 			// }

// 			//Update temporary seed length
// 			tmpSLen += progress;

// 			//Update current position in the unitig sequence
// 			uniSeqPos += progress;

// 			//If we have reached a seed check if it suffices to get a positive tempScore
// 			if((tmpScore += progress) > 0){//TODO This needs to be changed as soon as we want to use a non-unit score!
// 				//Update hit's length
// 				hitLen = extLen + tmpSLen;

// 				// //Check weather we have already created an extension seed for this unitig
// 				// if(extSeed != NULL){
// 				// 	//Update the seed's score
// 				// 	extSeed->score += tmpScore;

// 				//Check if there is another seed to reach
// 				if(prevSeed != NULL){
// 					//Exclude the reached seed from its seed list
// 					prevSeed->nextSeed = nearestSeed->nextSeed;
// 					//Delete the reached seed
// 					free(nearestSeed);
// 					//Search for the next neighbor
// 					nearestSeed = searchRightNeighbor(sucUnitig->getData()->getData(*sucUnitig)->getSeed(sucUnitig->strand), iniQoff, extLen, iniSeqPos, prevSeed);
// 				} else{
// 					//Set seed's successor as head of the seed list
// 					sucUnitig->getData()->getData(*sucUnitig)->setSeed(nearestSeed->nextSeed, sucUnitig->strand);
// 					//Delete the reached see
// 					free(nearestSeed);
// 					//Reset nearestSeed
// 					nearestSeed = NULL;
// 				}

// 				score += tmpScore;
// 				//Reset the temporary score
// 				tmpScore = 0;
// 			} else{
// 				//Check if there is another seed to reach
// 				if(prevSeed != NULL){
// 					//Exclude the reached seed from its seed list
// 					prevSeed->nextSeed = nearestSeed->nextSeed;
// 					//Delete the reached seed
// 					free(nearestSeed);
// 					//Search for the next one
// 					nearestSeed = searchRightNeighbor(sucUnitig->getData()->getData(*sucUnitig)->getSeed(sucUnitig->strand), iniQoff, extLen, iniSeqPos, prevSeed);
// 				} else{
// 					//Set seed's successor as head of the seed list
// 					sucUnitig->getData()->getData(*sucUnitig)->setSeed(nearestSeed->nextSeed, sucUnitig->strand);
// 					//Delete the reached seed
// 					free(nearestSeed);
// 					//Reset nearestSeed
// 					nearestSeed = NULL;
// 				}
// 			}
// 		} else{
// 			//Check up to which point we have to compare the unitig sequence
// 			if(sucUnitig->getSuccessors().hasSuccessors()){
// 				endBuf = sucUnitig->getGraph()->getK() - 1;
// 			}

// 			//Check whether we have reached the end of the unitig's sequence
// 			if(uniSeqPos < sucUnitig->size - endBuf){
// 				//Check whether quorum has to be checked
// 				if(checkedPos == 0){
// 					//Calculate offset
// 					offset = uniSeqPos;

// 					//Check whether quorum is still fulfilled
// 					if((checkedPos = checkSearchCrit(*sucUnitig, quorum, offset)) == 0) break;
// 				}

// 				//Check whether the score of our extension is positive
// 				if((tmpScore += compUScore(sucUniSeq[uniSeqPos], q[iniQoff + extLen + tmpSLen])) > 0){
// 					// //Update the seed info (Here we must use uniSeqPos instead of tmpSLen which is okay since our extension seed always starts from the beginning of the unitig)
// 					// extSeed->len = uniSeqPos + 1;//+1 because the index starts with 0
// 					// extSeed->score += tmpScore;

// 					//Update score
// 					score += tmpScore;
// 					//Update hit's length
// 					hitLen = extLen + tmpSLen + 1;
// 					//Reset the temporary score
// 					tmpScore = 0;
// 				} else{
// 					//Check if the current extension is already too bad
// 					if(tmpScore < -X){ break; }
// 				}

// 				//Proceed with the next two positions
// 				++tmpSLen;
// 				++uniSeqPos;
// 				--checkedPos;
// 			} else{
// 				//Check if the current unitig has successors
// 				if(endBuf != 0){
// 					//Calculate the position in the next unitig's sequence we have to start with
// 					uniSeqPos = uniSeqPos - sucUnitig->size + endBuf;
// 					//Try to extend on successive unitigs
// 					score += extendAtNextUnitig(sucUnitig->getSuccessors(), iniQoff, hitLen, extLen + tmpSLen, q, X, tmpScore, uniSeqPos, extPth, explCount, quorum);
// 				}

// 				break;
// 			}
// 		}
// 	}

// 	// if(extPtr != NULL){ return extPtr->extSeed->score; }

// 	return score;
// }

// //This function continues an extension to the right on a successive unitig of a seed lying on the query's reverse complement considering a quorum. Returns the maximum score reached.
// int32_t contRightX_Drop_OnRevComp(const neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> &sucUnitig, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t &extLen, const string &q, const int16_t &X, struct Ext_ptr*& extPtr, const int32_t &lastSeedTmpScore, uint32_t uniSeqPos, uint32_t &explCount, const uint32_t &quorum){
// 	int32_t tmpScore, progress;
// 	uint32_t iniSeqPos;
// 	uint32_t tmpSLen;
// 	uint32_t checkedPos = 0;
// 	size_t offset;
// 	string sucUniSeq = sucUnitig->mappedSequenceToString();
// 	struct seed *nearestSeed, *prevSeed, *extSeed = NULL;

// 	//Testing
// 	//bool seedReached = true;

// 	//Save the initial offset in the current unitig which we need for all nearest neighbor calculations
// 	iniSeqPos = uniSeqPos;
// 	//Find the nearest seed that we might be able to reach during our extension
// 	nearestSeed = searchRightNeighbor(sucUnitig->getData()->getData(*sucUnitig)->getSeed(sucUnitig->strand), iniQoff, extLen, iniSeqPos, prevSeed);

// 	//Perform the X-drop algorithm on the successive unitig
// 	tmpScore = lastSeedTmpScore;
// 	tmpSLen = 0;

// 	//We are done if we have reached the end of the query
// 	while(iniQoff + extLen + tmpSLen < q.length()){
// 		//Testing
// 		//cout << "2 Option 2" << endl;
// 		//cout << "Do we get here once more?" << endl << "iniQoff + extLen + tmpSLen:" << iniQoff + extLen + tmpSLen << endl;
// 		// cout << "nearestSeed is " << (nearestSeed == NULL ? "" : "not ") << "NULL" << endl;
// 		// if(nearestSeed != NULL){
// 		// 	cout << "uoff:" << nearestSeed->offsetU << " qoff:" << nearestSeed->offsetQ << endl;
// 		// }

// 		//Check whether we have reached the next seed
// 		if(nearestSeed != NULL && iniQoff + extLen + tmpSLen >= nearestSeed->offsetQ){
// 			//Testing
// 			//seedReached = true;

// 			//Calculate the gain we get by incorporating the reached seed
// 			progress = nearestSeed->offsetQ + nearestSeed->len - (iniQoff + extLen + tmpSLen);
// 			//Update temporary seed length
// 			tmpSLen += progress;
// 			//Update current position in the unitig sequence
// 			uniSeqPos += progress;

// 			//Testing
// 			//cout << "We do get here" << endl;

// 			//If we have reached a seed check if it suffices to get a positive tempScore
// 			if((tmpScore += progress) > 0){//TODO This needs to be changed as soon as we want to use a non-unit score!//TODO Maybe this whole progress stuff can be skipped due to efficiency reasons! Instead, it would be sufficient to delete all seeds which could be reached after the extension on that unitig is finished and why do we actually need an extension seed per unitig again?
// 				//Testing
// 				//cout << "3 Option 1" << endl;

// 				//Update hit's length
// 				hitLen = extLen + tmpSLen;

// 				//Check weather we have already created an extension seed for this unitig
// 				if(extSeed != NULL){
// 					//Testing
// 					//cout /*<< "4 Option 1" << endl */<< "5";
// 					//cout << "prevSeed is " << (prevSeed == NULL ? "" : "not ") << "NULL" << endl;

// 					//Update the seed's score
// 					extSeed->score += tmpScore;

// 					//Check if there is another seed to reach
// 					if(prevSeed != NULL){
// 						//Exclude the reached seed from its seed list
// 						prevSeed->nextSeed = nearestSeed->nextSeed;
// 						//Delete the reached seed
// 						free(nearestSeed);
// 						//Search for the next neighbor
// 						nearestSeed = searchRightNeighbor(sucUnitig->getData()->getData(*sucUnitig)->getSeed(sucUnitig->strand), iniQoff, extLen, iniSeqPos, prevSeed);
// 					} else{
// 						//Set seed's successor as head of the seed list
// 						sucUnitig->getData()->getData(*sucUnitig)->setSeed(nearestSeed->nextSeed, sucUnitig->strand);
// 						//Delete the reached seed
// 						free(nearestSeed);
// 						//Reset nearestSeed
// 						nearestSeed = NULL;
// 					}
// 					//Testing
// 					//cout << " Option 2" << endl;
// 				} else{
// 					//Testing
// 					//cout << "4 Option 2" << endl;

// 					//Let the reached seed become the extension seed
// 					extSeed = nearestSeed;
// 					//Adjust seed's offsets, length and score
// 					extSeed->offsetU = 0;//TODO Sometimes we might start the sequence comparison not at position 0 in the unitig sequence (because there was a seed in the predecessive unitig reaching into the buffer area). In such a case the score for the extension seed is not bases on the comparison of all bases of the unitig but offsetU = 0 suggested that. Check whether this affects other parts of the code and change it!
// 					extSeed->offsetQ = iniQoff + extLen - iniSeqPos;
// 					//Set up the seed's score
// 					extSeed->score = tmpScore;
// 					//Reset nearestSeed
// 					nearestSeed = NULL;

// 					//Check whether our seed has a predecessor
// 					if(prevSeed != NULL){
// 						//Exclude the seed from its seed list, i.e. link its predecessor to its successor
// 						prevSeed->nextSeed = extSeed->nextSeed;
// 						nearestSeed = searchRightNeighbor(sucUnitig->getData()->getData(*sucUnitig)->getSeed(sucUnitig->strand), iniQoff, extLen, iniSeqPos, prevSeed);
// 					} else{
// 						//Testing
// 						//cout << "6 Option 2" << endl;

// 						//Exclude the seed from its current seed list
// 						sucUnitig->getData()->getData(*sucUnitig)->setSeed(extSeed->nextSeed, sucUnitig->strand);
// 					}
// 				}
// 				//Transfer seed length and score to the extension seed
// 				extSeed->len = uniSeqPos;
// 				//Reset the temporary score
// 				tmpScore = 0;
// 			} else{
// 				//Check if there is another seed to reach
// 				if(prevSeed != NULL){
// 					//Exclude the reached seed from its seed list
// 					prevSeed->nextSeed = nearestSeed->nextSeed;
// 					//Delete the reached seed
// 					free(nearestSeed);
// 					//Search for the next one
// 					nearestSeed = searchRightNeighbor(sucUnitig->getData()->getData(*sucUnitig)->getSeed(sucUnitig->strand), iniQoff, extLen, iniSeqPos, prevSeed);
// 				} else{
// 					//Set seed's successor as head of the seed list
// 					sucUnitig->getData()->getData(*sucUnitig)->setSeed(nearestSeed->nextSeed, sucUnitig->strand);
// 					//Delete the reached seed
// 					free(nearestSeed);
// 					//Reset nearestSeed
// 					nearestSeed = NULL;
// 				}
// 				//Testing
// 				//cout << " Option 2" << endl;
// 			}

// 			//Testing
// 			//if(nearestSeed != NULL) seedReached = false;
// 		} else{
// 			//Check whether we have reached the end of the unitig's sequence
// 			if(uniSeqPos < sucUnitig->size){
// 				//Check whether quorum has to be checked
// 				if(checkedPos == 0){
// 					//Testing
// 					//cout << "not";

// 					//Calculate offset
// 					offset = uniSeqPos;

// 					//Check whether quorum is still fulfilled
// 					if((checkedPos = checkSearchCrit(*sucUnitig, quorum, offset)) == 0) break;
// 				}
// 				//Testing
// 				//cout << " Option 2" << endl;
// 				//cout << "sucUniSeq[uniSeqPos]:" << sucUniSeq[uniSeqPos] << "q[iniQoff + extLen + tmpSLen]:" << q[iniQoff + extLen + tmpSLen]/*<< " iniQoff + extLen + tmpSLen:" << iniQoff + extLen + tmpSLen */<< endl;

// 				//Check whether the score of our extension is positive
// 				if((tmpScore += compUScore(sucUniSeq[uniSeqPos], q[iniQoff + extLen + tmpSLen])) > 0){//TODO It might be that we start on a new unitig and bring a positive tmpSco
// 					//Testing
// 					//cout << "Option 1" << endl;

// 					//Check weather we haven't created an extension seed for this unitig yet
// 					if(extSeed == NULL){
// 						//Create a new seed
// 						extSeed = (struct seed*) malloc(sizeof(struct seed));
// 						//Creating a new seed means we are right at the beginning of the unitig sequence
// 						extSeed->offsetU = 0;
// 						//In q we start where the predecessive unitig's seed ended
// 						extSeed->offsetQ = iniQoff + extLen - iniSeqPos;
// 						//Seed's length is the distance within the unitig's sequence we have processed so far
// 						//Seed's score is the score gained so far
// 						extSeed->score = 0;
// 					}
					
// 					//Update the seed info (Here we must use uniSeqPos instead of tmpSLen which is okay since our extension seed always starts from the beginning of the unitig)
// 					extSeed->len = uniSeqPos + 1;//+1 because the index starts with 0
// 					extSeed->score += tmpScore;
// 					//Update hit's length
// 					hitLen = extLen + tmpSLen + 1;
// 					//Reset the temporary score
// 					tmpScore = 0;
// 				} else{
// 					//Testing
// 					//cout /*<< "11 Option 2" << endl */<< "12 ";

// 					//Check if the current extension is already too bad
// 					if(tmpScore < -X){
// 						//Testing
// 						//cout << "12 Option 2" << endl;
					
// 						break;
// 					}
// 					//Testing
// 					//cout << "Option 1" << endl;
// 				}

// 				//Proceed with the next two positions
// 				++tmpSLen;
// 				++uniSeqPos;
// 				--checkedPos;
// 			} else{
// 				//Check if the current unitig has successors
// 				if(sucUnitig->getSuccessors().hasSuccessors()){
// 					//Calculate the position in the next unitig's sequence we have to start with
// 					uniSeqPos = sucUnitig->getGraph()->getK() - 1;
// 					//Try to extend on successive unitigs
// 					extPtr = extendAtNextUnitig_OnRevComp(sucUnitig->getSuccessors(), iniQoff, hitLen, extLen + tmpSLen, q, X, tmpScore, uniSeqPos, explCount, quorum);
// 				}

// 				break;
// 			}
// 		}
// 	}

// 	//Testing
// 	//cout << "2 Option 1" << endl;

// 	//Check whether we had a successful extension within this unitig
// 	if(extSeed != NULL){
// 		//Check whether there was a successful extension while exploring successive seeds
// 		if(extPtr == NULL){
// 			//Testing
// 			//cout << "14 Option 2" << endl;

// 			//Create a new extension pointer
// 			extPtr = (struct Ext_ptr*) malloc(sizeof(struct Ext_ptr));
// 			//Link this unitig to the extension pointer
// 			extPtr->origUnitig = *sucUnitig;
// 			//Link the extension seed
// 			extPtr->extSeed = extSeed;
// 		} else{
// 			//Testing
// 			//cout << "14 Option 1" << endl;

// 			//So far the score gain while exploring this unitig is not considered in the extension pointer seed
// 			extPtr->extSeed->score += extSeed->score;
// 			//Delete the extension seed
// 			free(extSeed);
// 		}
// 	}

// 	if(extPtr != NULL){
// 		//Testing
// 		//cout << "13 Option 2" << endl;

// 		return extPtr->extSeed->score;
// 	}
// 	//Testing
// 	//cout << "13 Option 1" << endl;

//This function continues an extension to the right on a successive unitig of a seed lying on the query's reference strand considering a quorum and a search color set. Returns the maximum score reached.
int32_t contRightX_Drop(const neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> &sucUnitig, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t &extLen, const string &q, const int16_t &X, const int32_t &lastSeedTmpScore, uint32_t uniSeqPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	int32_t tmpScore, progress, score = 0;
	int32_t overlap = sucUnitig->getGraph()->getK() - 1;
	uint32_t iniSeqPos;
	uint32_t tmpSLen;
	//Even if we are on the reverse complementary strand and no position is covered on this unitig, because we have checked the first k - 1 position on the last unitig already
	int32_t checkedPos = getSrchCritCov(*sucUnitig, quorum, searchSet, compOffset(uniSeqPos, 1, sucUnitig->size, sucUnitig->strand), sucUnitig->strand);// - uniSeqPos;

	//Testing
	// cout << "checkedPos: " << checkedPos << endl;

	//bool inSeedList = false;

	string sucUniSeq = sucUnitig->mappedSequenceToString();
	struct seed *nearestSeed, *prevSeed;//, *extSeed = NULL;

	//Save the initial offset in the current unitig which we need for all nearest neighbor calculations
	iniSeqPos = uniSeqPos;

	// //Check whether we are working on the reference strand of the curent unitig
	// if(sucUnitig->strand){

	// //If we are on the reverse complementary strand we have to consider the overlap for checkedPos
	// if(!sucUnitig->strand && checkedPos == 8) checkedPos -= overlap;

	//Find the nearest seed that we might be able to reach during our extension
	nearestSeed = searchRightNeighbor(sucUnitig->getData()->getData(*sucUnitig)->getSeed(sucUnitig->strand), iniQoff, extLen, iniSeqPos, prevSeed);

	// } else{
	// 	//So far we have not considered the reverse complementary strand for the seed detection //TODO: Change this!
	// 	nearestSeed = NULL;
	// }

	//Perform the X-drop algorithm on the successive unitig
	tmpScore = lastSeedTmpScore;
	tmpSLen = 0;

	//We are done if we have reached the end of the query
	while(iniQoff + extLen + tmpSLen < q.length()){
		//Check whether we have reached the next seed
		if(nearestSeed != NULL && iniQoff + extLen + tmpSLen >= nearestSeed->offsetQ){
			//Calculate the gain we get by incorporating the reached seed
			progress = nearestSeed->offsetQ + nearestSeed->len - (iniQoff + extLen + tmpSLen);
			//Update temporary seed length
			tmpSLen += progress;
			//Update current position in the unitig sequence
			uniSeqPos += progress;

			//Adjust number of remaining covered positions
			checkedPos -= progress;

			//If we have reached a seed check if it suffices to get a positive tempScore
			if((tmpScore += progress) > 0){//TODO This needs to be changed as soon as we want to use a non-unit score!
				//Update hit's length
				hitLen = extLen + tmpSLen;

				// //Check weather we have already created an extension seed for this unitig
				// if(extSeed != NULL){
				// 	//Update the seed's score
				// 	extSeed->score += tmpScore;

				//Check if there is another seed to reach
				if(prevSeed != NULL){
					//Exclude the reached seed from its seed list
					prevSeed->nextSeed = nearestSeed->nextSeed;
					//Delete the reached seed
					free(nearestSeed);
					//Search for the next neighbor
					nearestSeed = searchRightNeighbor(sucUnitig->getData()->getData(*sucUnitig)->getSeed(sucUnitig->strand), iniQoff, extLen, iniSeqPos, prevSeed);
				} else{
					//Set seed's successor as head of the seed list
					sucUnitig->getData()->getData(*sucUnitig)->setSeed(nearestSeed->nextSeed, sucUnitig->strand);
					//Delete the reached seed
					free(nearestSeed);
					//Reset nearestSeed
					nearestSeed = NULL;
				}

				// } else{
				// 	//Let the reached seed become the extension seed
				// 	extSeed = nearestSeed;
				// 	//Adjust seed's offsets, length and score
				// 	extSeed->offsetU = 0;
				// 	extSeed->offsetQ = iniQoff + extLen - iniSeqPos;
				// 	//Set up the seed's score (has to be != 0 to mark it as visited)
				// 	extSeed->score = tmpScore;

				// 	// //Mark that we recycle a seed instead of creating a new one for our extension
				// 	// inSeedList = true;

				// 	//Reset nearestSeed
				// 	nearestSeed = NULL;

				// 	//Check whether our seed has a predecessor
				// 	if(prevSeed != NULL){
				// 		//Exclude the seed from its current seed list
				// 		prevSeed->nextSeed = extSeed->nextSeed;
				// 		nearestSeed = searchRightNeighbor(sucUnitig->getData()->getData(*sucUnitig)->getSeed(sucUnitig->strand), iniQoff, extLen, iniSeqPos, prevSeed);
				// 	} else{
				// 		//Exclude the seed from its current seed list
				// 		sucUnitig->getData()->getData(*sucUnitig)->setSeed(extSeed->nextSeed, sucUnitig->strand);
				// 	}

				// 	// //Insert the seed into the seed list of processed seeds
				// 	// extSeed->nextSeed = sucUnitig->getData()->getData(*sucUnitig)->getLastProcSeed();
				// 	// sucUnitig->getData()->getData(*sucUnitig)->setLastProcSeed(extSeed);
				// }
				// //Transfer seed length and score to the extension seed
				// extSeed->len = uniSeqPos;

				//Update score
				score += tmpScore;
				//Reset the temporary score
				tmpScore = 0;
			} else{
				//Check if there is another seed to reach
				if(prevSeed != NULL){
					//Exclude the reached seed from its seed list
					prevSeed->nextSeed = nearestSeed->nextSeed;
					//Delete the reached seed
					free(nearestSeed);
					//Search for the next one
					nearestSeed = searchRightNeighbor(sucUnitig->getData()->getData(*sucUnitig)->getSeed(sucUnitig->strand), iniQoff, extLen, iniSeqPos, prevSeed);
				} else{
					//Set seed's successor as head of the seed list
					sucUnitig->getData()->getData(*sucUnitig)->setSeed(nearestSeed->nextSeed, sucUnitig->strand);
					//Delete the reached seed
					free(nearestSeed);
					//Reset nearestSeed
					nearestSeed = NULL;
				}
			}
		} else{
			//Check up to which point we have to compare the unitig sequence
			if(!sucUnitig->getSuccessors().hasSuccessors()){
				overlap = 0;
			}

			//Check whether we have reached the end of the unitig's sequence
			if(uniSeqPos < sucUnitig->size - overlap){
				//Are search criteria still fulfilled?
				if(checkedPos <= 0){
					break;

					// //Calculate offset
					// offset = compOffset(uniSeqPos, sucUnitig->getGraph()->getK(), sucUnitig->size, sucUnitig->strand);

					// //Check whether quorum is still fulfilled
					// if((checkedPos = checkSearchCrit(*sucUnitig, quorum, offset, searchSet)) == 0) break;
				}

				//Testing
				// cout << "Compare bases unitig: " << sucUniSeq[uniSeqPos] << " query: " << q[iniQoff + extLen + tmpSLen] << endl;

				//Check whether the score of our extension is positive
				if((tmpScore += compUScore(sucUniSeq[uniSeqPos], q[iniQoff + extLen + tmpSLen])) > 0){
					// //Check weather we haven't created an extension seed for this unitig yet
					// if(extSeed == NULL){
					// 	//Create a new seed
					// 	extSeed = (struct seed*) malloc(sizeof(struct seed));
					// 	//Creating a new seed means we are right at the beginning of the unitig sequence
					// 	extSeed->offsetU = 0;
					// 	//In q we start where the predecessive unitig's seed ended
					// 	extSeed->offsetQ = iniQoff + extLen - iniSeqPos;
					// 	//Seed's length is the distance within the unitig's sequence we have processed so far
					// 	//Seed's score is the score gained so far
					// 	extSeed->score = 0;
						
					// 	// //Mark that this seed has still to be added to the seed list
					// 	// inSeedList = false;
					// }
					
					// //Update the seed info (Here we must use uniSeqPos instead of tmpSLen which is okay since our extension seed always starts from the beginning of the unitig)
					// extSeed->len = uniSeqPos + 1;//+1 because the index starts with 0
					// extSeed->score += tmpScore;

					//Update score
					score += tmpScore;
					//Update hit's length
					hitLen = extLen + tmpSLen + 1;
					//Reset the temporary score
					tmpScore = 0;
				} else{
					//Check if the current extension is already too bad
					if(tmpScore < -X){ break; }
				}

				//Proceed with the next two positions
				++tmpSLen;
				++uniSeqPos;
				--checkedPos;
			} else{
				//Check if the current unitig has successors
				if(overlap != 0){
					//Calculate the position in the next unitig's sequence we have to start with
					uniSeqPos = uniSeqPos - sucUnitig->size + overlap;
					//Check out next unitig
					score += extendAtNextUnitig(sucUnitig->getSuccessors(), iniQoff, hitLen, extLen + tmpSLen, q, X, tmpScore, uniSeqPos, extPth, explCount, quorum, searchSet);
				}

				break;
			}
		}
	}

	// //Check whether we had a successful extension within this unitig
	// if(extSeed != NULL){
	// 	//Check whether there was a successful extension while exploring successive seeds
	// 	if(extPtr == NULL){
	// 		//Create a new extension pointer
	// 		extPtr = (struct Ext_ptr*) malloc(sizeof(struct Ext_ptr));
	// 		//Link this unitig to the extension pointer
	// 		extPtr->origUnitig = *sucUnitig;
	// 		//Link the extension seed
	// 		extPtr->extSeed = extSeed;

	// 		// //Save whether the linked seed still needs to be incorporated into a seed list
	// 		// extPtr->inSeedList = inSeedList;
	// 	} else{
	// 		//So far the score gain while exploring this unitig is not considered in the extension pointer seed
	// 		extPtr->extSeed->score += extSeed->score;
	// 		//Delete this extension seed
	// 		free(extSeed);

	// 		//delExtPtrSeed(inSeedList, extSeed, sucUnitig->getData()->getData(*sucUnitig));
	// 	}
	// }

	// if(extPtr != NULL){ return extPtr->extSeed->score; }

	return score;
}

//This function continues an extension to the right on a successive unitig of a seed lying on the query's reverse complement considering a quorum and a search color set. Returns the maximum score reached.
int32_t contRightX_Drop_OnRevComp(const neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> &sucUnitig, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t &extLen, const string &q, const int16_t &X, const int32_t &lastSeedTmpScore, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	int32_t tmpScore, progress, score = 0;
	//Since we consider the overlap in the sequence's beginning for seed's on the reverse complementary strand the initial unitig position is always the same
	uint32_t iniSeqPos = sucUnitig->getGraph()->getK() - 1;
	uint32_t uniSeqPos = iniSeqPos;
	uint32_t tmpSLen;
	//Get the number of covered positions
	int32_t checkedPos = getSrchCritCov(*sucUnitig, quorum, searchSet, compOffset(uniSeqPos, 1, sucUnitig->size, sucUnitig->strand), sucUnitig->strand);

	//Testing
	// cout << "checkedPos: " << checkedPos << endl;
	// cout << "(int32_t) sucUnitig->size - iniSeqPos: " << (int32_t) sucUnitig->size - iniSeqPos << endl;
	// cout << "(int32_t) sucUnitig->size - iniSeqPos > checkedPos: " << ((int32_t) (sucUnitig->size - iniSeqPos) > checkedPos) << endl;
	// cout << "sucUnitig: " << sucUnitig->mappedSequenceToString() << endl;

	// //Check if the number of checked positions from the right are enough to reach our current position
	// if((int32_t) (sucUnitig->size - iniSeqPos) > checkedPos){
	// 	//Testing
	// 	// cout << "Do we get here?" << endl;

	// 	checkedPos = getSrchCritCov(*sucUnitig, quorum, searchSet, !sucUnitig->strand);

	// 	//Check if we have enough positions from the left
	// 	if(checkedPos <= compOffset(iniSeqPos, 0, sucUnitig->size, sucUnitig->strand)){
	// 		checkedPos = 0;
	// 	} else{
	// 		checkedPos -= compOffset(iniSeqPos, 0, sucUnitig->size, sucUnitig->strand);
	// 	}
	// }

	// size_t offset;

	//Testing
	// cout << "checkedPos: " << checkedPos << endl;
	// cout << "iniSeqPos: " << iniSeqPos << endl;
	// cout << "iniQoff: " << iniQoff << " extLen: " << extLen << " tmpSLen: " << 0 << endl;
	// cout << "getSrchCritCov(*sucUnitig, quorum, searchSet, sucUnitig->strand): " << getSrchCritCov(*sucUnitig, quorum, searchSet, sucUnitig->strand) << endl;

	string sucUniSeq = sucUnitig->mappedSequenceToString();
	struct seed *nearestSeed, *prevSeed;//, *extSeed = NULL;

	//Find the nearest seed that we might be able to reach during our extension
	nearestSeed = searchRightNeighbor(sucUnitig->getData()->getData(*sucUnitig)->getSeed(sucUnitig->strand), iniQoff, extLen, iniSeqPos, prevSeed);
	//Perform the X-drop algorithm on the successive unitig
	tmpScore = lastSeedTmpScore;
	tmpSLen = 0;

	//Testing
	// cout << "nearestSeed is " << (nearestSeed == NULL ? "" : "not ") << "NULL" << endl;

	//We are done if we have reached the end of the query
	while(iniQoff + extLen + tmpSLen < q.length()){
		//Testing
		// cout << "1 Option 2" << endl;
		// cout << "uniSeqPos: " << uniSeqPos << endl;

		//Check whether we have reached the next seed
		if(nearestSeed != NULL && iniQoff + extLen + tmpSLen >= nearestSeed->offsetQ){
			//Testing
			// cout << "2 Option 1" << endl;

			//Calculate the gain we get by incorporating the reached seed
			progress = nearestSeed->offsetQ + nearestSeed->len - (iniQoff + extLen + tmpSLen);
			//Update temporary seed length
			tmpSLen += progress;
			//Update current position in the unitig sequence
			uniSeqPos += progress;
			//Adjust number of remaining covered positions
			checkedPos -= progress;

			//Testing
			// cout << "progress: " << progress << endl;

			//If we have reached a seed check if it suffices to get a positive tempScore
			if((tmpScore += progress) > 0){//TODO This needs to be changed as soon as we want to use a non-unit score!
				//Testing
				// cout << "3 Option 1" << endl;

				//Update hit's length
				hitLen = extLen + tmpSLen;

				// //Check weather we have already created an extension seed for this unitig
				// if(extSeed != NULL){
				// 	//Testing
				// 	// cout << "4 Option 1" << endl;

				// 	//Update the seed's score
				// 	extSeed->score += tmpScore;

				//Check if there is another seed to reach
				if(prevSeed != NULL){
					//Testing
					// cout << "5 Option 1" << endl;

					//Exclude the reached seed from its seed list
					prevSeed->nextSeed = nearestSeed->nextSeed;
					//Delete the reached seed
					free(nearestSeed);
					//Search for the next neighbor
					nearestSeed = searchRightNeighbor(sucUnitig->getData()->getData(*sucUnitig)->getSeed(sucUnitig->strand), iniQoff, extLen, iniSeqPos, prevSeed);
				} else{
					//Testing
					// cout << "5 Option 2" << endl;

					//Set seed's successor as head of the seed list
					sucUnitig->getData()->getData(*sucUnitig)->setSeed(nearestSeed->nextSeed, sucUnitig->strand);
					//Delete the reached seed
					free(nearestSeed);
					//Reset nearestSeed
					nearestSeed = NULL;
				}

				// } else{
				// 	//Testing
				// 	// cout << "4 Option 2" << endl;

				// 	//Let the reached seed become the extension seed
				// 	extSeed = nearestSeed;
				// 	//Adjust seed's offsets, length and score
				// 	extSeed->offsetU = iniSeqPos;
				// 	extSeed->offsetQ = iniQoff + extLen + tmpSLen - (uniSeqPos - iniSeqPos);

				// 	//Testing
				// 	// cout << "iniQoff: " << iniQoff << "extLen: " << extLen << "uniSeqPos: " << uniSeqPos << "iniSeqPos: " << iniSeqPos << "extSeed->offsetQ: " << extSeed->offsetQ << endl;

				// 	//Set up the seed's score (has to be != 0 to mark it as visited)
				// 	extSeed->score = tmpScore;
				// 	//Reset nearestSeed
				// 	nearestSeed = NULL;

				// 	//Check whether our seed has a predecessor
				// 	if(prevSeed != NULL){
				// 		//Testing
				// 		// cout << "6 Option 1" << endl;

				// 		//Exclude the seed from its current seed list
				// 		prevSeed->nextSeed = extSeed->nextSeed;
				// 		nearestSeed = searchRightNeighbor(sucUnitig->getData()->getData(*sucUnitig)->getSeed(sucUnitig->strand), iniQoff, extLen, iniSeqPos, prevSeed);
				// 	} else{
				// 		//Testing
				// 		// cout << "6 Option 2" << endl;

				// 		//Exclude the seed from its current seed list
				// 		sucUnitig->getData()->getData(*sucUnitig)->setSeed(extSeed->nextSeed, sucUnitig->strand);
				// 	}
				// }
				// //Transfer seed length and score to the extension seed
				// extSeed->len = uniSeqPos - iniSeqPos;

				//Update score
				score += tmpScore;
				//Reset the temporary score
				tmpScore = 0;
			} else{
				//Testing
				// cout << "3 Option 2" << endl;

				//Check if there is another seed to reach
				if(prevSeed != NULL){
					//Testing
					// cout << "7 Option 1" << endl;

					//Exclude the reached seed from its seed list
					prevSeed->nextSeed = nearestSeed->nextSeed;
					//Delete the reached seed
					free(nearestSeed);
					//Search for the next one
					nearestSeed = searchRightNeighbor(sucUnitig->getData()->getData(*sucUnitig)->getSeed(sucUnitig->strand), iniQoff, extLen, iniSeqPos, prevSeed);
				} else{
					//Testing
					// cout << "7 Option 2" << endl;

					//Set seed's successor as head of the seed list
					sucUnitig->getData()->getData(*sucUnitig)->setSeed(nearestSeed->nextSeed, sucUnitig->strand);
					//Delete the reached seed
					free(nearestSeed);
					//Reset nearestSeed
					nearestSeed = NULL;
				}
			}
		} else{
			//Testing
			// cout << "2 Option 2" << endl;
			// cout << "uniSeqPos: " << uniSeqPos << "sucUnitig->size: " << sucUnitig->size << endl;

			//Check whether we have reached the end of the unitig's sequence
			if(uniSeqPos < sucUnitig->size){
				//Testing
				// cout << "8 Option 2" << endl;
				// cout << "9 Option ";

				//Check if the current position still fulfills the search criteria
				if(checkedPos <= 0){
					//Testing
					// cout << "not ";

					// //Check on which strand we are
					// if(sucUnitig->strand){
					// 	//Calculate offset
					// 	offset = compOffset(uniSeqPos, sucUnitig->getGraph()->getK(), sucUnitig->size, true);
					// } else{
					// 	//On the reverse complementary strand it does not suffice to check the highest offset possible
					// 	if(uniSeqPos > sucUnitig->size - sucUnitig->getGraph()->getK()){
					// 		//Check the lowest position possible
					// 		offset = 0;
					// 	} else{
					// 		//Move k-1 positions to the left and calculate the corresponding offset on the other strand
					// 		offset = compOffset(uniSeqPos, sucUnitig->getGraph()->getK(), sucUnitig->size, false);
					// 	}
					// }

					// //Testing
					// // cout << "Do check for offset: " << offset << endl;

					// //Check whether quorum is still fulfilled
					// if((checkedPos = checkSearchCrit(*sucUnitig, quorum, offset, searchSet)) == 0){
						//Testing
						// cout << "2" << endl;
						// cout << "10 Option 2" << endl;
						// cout << "uniSeqPos: " << uniSeqPos << endl;

					break;

					// }

					//Testing
					// cout << "10 Option 1" << endl;
					// cout << "checkedPos: " << checkedPos << endl;
				}

				//Testing
				// cout << "2" << endl;

				//Check whether the score of our extension is positive
				if((tmpScore += compUScore(sucUniSeq[uniSeqPos], q[iniQoff + extLen + tmpSLen])) > 0){
					//Testing
					// cout << "11 Option 1" << endl;
					// cout << "12 Option ";

					// //Check weather we haven't created an extension seed for this unitig yet
					// if(extSeed == NULL){
					// 	//Testing
					// 	// cout << "not ";

					// 	//Create a new seed
					// 	extSeed = (struct seed*) malloc(sizeof(struct seed));
					// 	//Creating a new seed means we are right at the beginning of the unitig sequence
					// 	extSeed->offsetU = iniSeqPos;
					// 	//In q we start where the predecessive unitig's seed ended
					// 	extSeed->offsetQ = iniQoff + extLen + tmpSLen - (uniSeqPos - iniSeqPos);
					// 	//Seed's score is the score gained so far
					// 	extSeed->score = 0;
					// }

					// //Testing
					// // cout << "1" << endl;
					
					// //Update the seed info
					// extSeed->len = uniSeqPos - iniSeqPos + 1;//+1 because the index starts with 0
					// extSeed->score += tmpScore;

					//Update score
					score += tmpScore;
					//Update hit's length
					hitLen = extLen + tmpSLen + 1;
					//Reset the temporary score
					tmpScore = 0;
				} else{
					//Testing
					// cout << "11 Option 2" << endl;

					//Check if the current extension is already too bad
					if(tmpScore < -X){
						//Testing
						// cout << "13 Option 1" << endl;

						break;
					}

					//Testing
					// cout << "13 Option 2" << endl;
				}

				//Proceed with the next two positions
				++tmpSLen;
				++uniSeqPos;
				--checkedPos;
			} else{
				//Testing
				// cout << "8 Option 1" << endl;

				//Check if the current unitig has successors
				if(sucUnitig->getSuccessors().hasSuccessors()){
					//Testing
					// cout << "14 Option 1" << endl;

					//Calculate extension on successive unitigs
					score += extendAtNextUnitig_OnRevComp(sucUnitig->getSuccessors(), iniQoff, hitLen, extLen + tmpSLen, q, X, tmpScore, extPth, explCount, quorum, searchSet);
				}// else{
					//Testing
					// cout << "14 Option 2" << endl;
				//}

				break;
			}
		}
	}

	//Testing
	// cout << "1 Option 1" << endl;

	// //Check whether we had a successful extension within this unitig
	// if(extSeed != NULL){
	// 	//Testing
	// 	// cout << "15 Option 1" << endl;

	// 	//Check whether there was a successful extension while exploring successive seeds
	// 	if(extPtr == NULL){
	// 		//Testing
	// 		// cout << "16 Option 2" << endl;

	// 		//Create a new extension pointer
	// 		extPtr = (struct Ext_ptr*) malloc(sizeof(struct Ext_ptr));
	// 		//Link this unitig to the extension pointer
	// 		extPtr->origUnitig = *sucUnitig;
	// 		//Link the extension seed
	// 		extPtr->extSeed = extSeed;
	// 	} else{
	// 		//Testing
	// 		// cout << "16 Option 1" << endl;

	// 		//So far the score gain while exploring this unitig is not considered in the extension pointer seed
	// 		extPtr->extSeed->score += extSeed->score;
	// 		//Delete this extension seed
	// 		free(extSeed);
	// 	}
	// } else{
	// 	//Testing
	// 	// cout << "15 Option 2" << endl;
	// }

	// if(extPtr != NULL){
	// 	//Testing
	// 	// cout << "17 Option 1" << endl;

	// 	return extPtr->extSeed->score;
	// }

	//Testing
	// cout << "17 Option 2" << endl;

	return score;
}

// //This function performs an extension on all possible predecessors of a unitig considering a quorum and returns the maximum scoring one
// int32_t extendAtPrevUnitig(const BackwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> bwIter, uint32_t qPos, uint32_t &hitLen, const uint32_t &tmpExtLen, const string &q, const int16_t &X, const int32_t &lastExtSeedTmpScore, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum){
// 	uint16_t predID = 0;
// 	uint32_t tmpHitLen, bestHitLen;
// 	int32_t maxScore = 0, currScore;
// 	struct Ext_ptr *bExtSeed = NULL;
// 	//BackwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> bwIter = lastUni.getPredecessors();

// 	//Check whether we have reached the maximum recursion depth of an extension
// 	if(++explCount > MAXRECURSIONDEPTH){
// 		//Report this incident//TODO: If this happens the condition that an extended seed cannot be reached anymore is violated. Implement a procedure after the extension that merges such seeds with other seeds that would have reached them!
// 		//cerr << "Maximum recursion depth reached during extension" << endl;
// 		//Terminate this extension
// 		return 0;
// 	}

// 	// //Check whether we have predecessors
// 	// if(bwIter.hasPredecessors()){
		
// 	//Testing
// 	// if(report){
// 	// 	cerr << "extendAtPrevUnitig" << endl;
// 	// }
	
// 	//maxScore = 0;

// 	//Iterate over all predecessors
// 	for(neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> nI = bwIter.begin(); nI != bwIter.end(); ++nI){
// 		// //Temporary extension pointer to store the current extensions border seed and unitig
// 		// struct Ext_ptr *tmpExtSeed = NULL;

// 		//Temporary extension path
// 		list<uint16_t> tmpPth;
// 		//Note which predecessor we are on
// 		++predID;
// 		//Set tmpHitLen
// 		tmpHitLen = hitLen + tmpExtLen;

// 		//Calculate the score of an extension of a successor
// 		currScore = contLeftX_Drop(nI, qPos, tmpHitLen, q, X, lastExtSeedTmpScore, tmpPth, explCount, quorum);

// 		//Check whether the score of the current successors extension is the best one found so far
// 		if(currScore > maxScore){
// 			//Update maxScore
// 			maxScore = currScore;
// 			//Update hit length
// 			bestHitLen = tmpHitLen;
			
// 			// //Check whether we saved an extension list before already
// 			// if(bExtSeed != NULL){
// 			// 	//Delete old extension list
// 			// 	delExtPtr(bExtSeed);
// 			// }
// 			// //Save the new extension list
// 			// bExtSeed = tmpExtSeed;

// 			//Save the extension path
// 			extPth = tmpPth;
// 			//Save which predecessor we have chosen
// 			extPth.push_front(predID);
// 		}// else if(tmpExtSeed != NULL){
// 		// 	//Delete the new extension list
// 		// 	delExtPtr(tmpExtSeed);
// 		// }
// 	}

// 	// //Check whether we could find a good extension
// 	// if(maxScore > 0){
// 	// 	//Pass hit length of the best extension to the calling function
// 	// 	hitLen = bestHitLen;

// 	// 	return bExtSeed;
// 	// }
	
// 	//}

// 	//Nothing found
// 	return maxScore;
// }

// //This function performs an extension on all possible predecessors of a unitig considering a quorum and returns the maximum scoring one. This function is explicitly designed for seeds lying on the query's reverse complement (considering the overlap between unitigs in sequences' beginning)
// struct Ext_ptr* extendAtPrevUnitigOnRevComp(const BackwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> bwIter, uint32_t qPos, uint32_t &hitLen, const uint32_t &tmpExtLen, const string &q, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &explCount, const uint32_t &quorum){
// 	uint32_t tmpHitLen, bestHitLen;
// 	int32_t maxScore = 0, currScore;
// 	struct Ext_ptr *bExtSeed = NULL;

// 	//Check whether we have reached the maximum recursion depth of an extension
// 	if(++explCount > MAXRECURSIONDEPTH){
// 		//Report this incident//TODO: If this happens the condition that an extended seed cannot be reached anymore is violated. Implement a procedure after the extension that merges such seeds with other seeds that would have reached them!
// 		//cerr << "Maximum recursion depth reached during extension" << endl;
// 		//Terminate this extension
// 		return NULL;
// 	}

// 	//Iterate over all predecessors
// 	for(neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> nI = bwIter.begin(); nI != bwIter.end(); ++nI){
// 		//Temporary extension pointer to store the current extensions border seed and unitig
// 		struct Ext_ptr *tmpExtSeed = NULL;
// 		//Set tmpHitLen
// 		tmpHitLen = hitLen + tmpExtLen;
// 		//Calculate the score of an extension of a successor
// 		currScore = contLeftX_DropOnRevComp(nI, qPos, tmpHitLen, q, X, tmpExtSeed, lastExtSeedTmpScore, explCount, quorum);

// 		//Check whether the score of the current successors extension is the best one found so far
// 		if(currScore > maxScore){
// 			//Update maxScore
// 			maxScore = currScore;
// 			//Update hit length
// 			bestHitLen = tmpHitLen;
			
// 			//Check whether we saved an extension list before already
// 			if(bExtSeed != NULL){
// 				//Delete old extension list
// 				delExtPtr(bExtSeed);
// 			}
// 			//Save the new extension list
// 			bExtSeed = tmpExtSeed;
// 		} else if(tmpExtSeed != NULL){
// 			//Delete the new extension list
// 			delExtPtr(tmpExtSeed);
// 		}
// 	}

// 	//Check whether we could find a good extension
// 	if(maxScore > 0){
// 		//Pass hit length of the best extension to the calling function
// 		hitLen = bestHitLen;
// 		return bExtSeed;
// 	}

// 	//Nothing found
// 	return NULL;
// }

//This function performs an extension on all possible predecessors of a unitig considering a quorum and a search color set and returns the maximum scoring one
int32_t extendAtPrevUnitig(const BackwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> bwIter, uint32_t qPos, uint32_t &hitLen, const uint32_t &tmpExtLen, const string &q, const int16_t &X, const int32_t &lastExtSeedTmpScore, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	uint16_t predID = 0;
	uint32_t tmpHitLen, maxHitLen = hitLen;
	int32_t maxScore = 0, currScore;

	// struct Ext_ptr *bExtSeed = NULL;
	//BackwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> bwIter = lastUni.getPredecessors();

	//Check whether we have reached the maximum recursion depth of an extension
	if(++explCount > MAXRECURSIONDEPTH){
		//Report this incident//TODO: If this happens the condition that an extended seed cannot be reached anymore is violated. Implement a procedure after the extension that merges such seeds with other seeds that would have reached them!
		//cerr << "Maximum recursion depth reached during extension" << endl;
		//Terminate this extension
		return 0;
	}

	// //Check whether we have predecessors
	// if(bwIter.hasPredecessors()){
		// maxScore = 0;

	//Iterate over all predecessors
	for(neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> nI = bwIter.begin(); nI != bwIter.end(); ++nI){
		// //Temporary extension pointer to store the current extensions border seed and unitig
		// struct Ext_ptr *tmpExtSeed = NULL;

		//Temporary extension path
		list<uint16_t> tmpPth;
		//Note which predecessor we are on
		++predID;
		//Set tmpHitLen
		tmpHitLen = hitLen + tmpExtLen;

		//Testing
		// cout << "tmpHitLen: " << tmpHitLen << endl;

		//Calculate the score of an extension of a successor
		currScore = contLeftX_Drop(nI, qPos, tmpHitLen, q, X, lastExtSeedTmpScore, tmpPth, explCount, quorum, searchSet);

		//Testing
		// cout << "tmpHitLen: " << tmpHitLen << endl;

		//Check whether the score of the current successors extension is the best one found so far
		if(currScore > maxScore){
			//Update maxScore
			maxScore = currScore;
			//Update hit length
			maxHitLen = tmpHitLen;

			// //Check whether we saved an extension list before already
			// if(bExtSeed != NULL){
			// 	//Delete old extension list
			// 	delExtPtr(bExtSeed);
			// }
			// //Save the new extension list
			// bExtSeed = tmpExtSeed;

			//Save the extension path
			extPth = tmpPth;
			//Save which predecessor we have chosen
			extPth.push_front(predID);
		}// else if(tmpExtSeed != NULL){
		// 	//Delete the new extension list
		// 	delExtPtr(tmpExtSeed);
		// }
	}

	// //Check whether we could find a good extension
	// if(maxScore > 0){
	// 	//Pass hit length of the best extension to the calling function
	// 	hitLen = bestHitLen;

	// 	return bExtSeed;
	// }
	
	//}

	//Save hit length
	hitLen = maxHitLen;

	//Return score
	return maxScore;
}

//This function performs an extension on all possible predecessors of a unitig considering a quorum and a search color set and returns the maximum scoring one. This function is explicitly designed for seeds lying on the query's reverse complement (considering the overlap between unitigs in sequences' beginning)
int32_t extendAtPrevUnitigOnRevComp(const BackwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> bwIter, uint32_t qPos, uint32_t &hitLen, const uint32_t &tmpExtLen, const string &q, const int16_t &X, const int32_t &lastExtSeedTmpScore, list<uint16_t> &extPth, const uint32_t &lead, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	uint16_t predID = 0;
	uint32_t tmpHitLen, maxHitLen = hitLen;
	int32_t maxScore = 0, currScore;

	// struct Ext_ptr *bExtSeed = NULL;

	//Check whether we have reached the maximum recursion depth of an extension
	if(++explCount > MAXRECURSIONDEPTH){
		//Testing
		// cout << "1 Option 1" << endl;

		//Terminate this extension
		return 0;
	}

	//Testing
	// cout << "1 Option 2" << endl;

	//Iterate over all predecessors
	for(neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> nI = bwIter.begin(); nI != bwIter.end(); ++nI){
		// //Temporary extension pointer to store the current extensions border seed and unitig
		// struct Ext_ptr *tmpExtSeed = NULL;

		//Temporary extension path
		list<uint16_t> tmpPth;
		//Note which predecessor we are on
		++predID;
		//Set tmpHitLen
		tmpHitLen = hitLen + tmpExtLen;
		//Calculate the score of an extension of a successor
		currScore = contLeftX_DropOnRevComp(nI, qPos, tmpHitLen, q, X, lastExtSeedTmpScore, tmpPth, nI->size - lead, explCount, quorum, searchSet);

		//Check whether the score of the current successors extension is the best one found so far
		if(currScore > maxScore){
			//Testing
			// cout << "2 Option 1" << endl;

			//Update maxScore
			maxScore = currScore;
			//Update hit length
			maxHitLen = tmpHitLen;

			// //Check whether we saved an extension list before already
			// if(bExtSeed != NULL){
			// 	//Testing
			// 	// cout << "3 Option 1" << endl;

			// 	//Delete old extension list
			// 	delExtPtr(bExtSeed);
			// } else{
			// 	//Testing
			// 	// cout << "3 Option 2" << endl;
			// }

			// //Save the new extension list
			// bExtSeed = tmpExtSeed;

			//Save the extension path
			extPth = tmpPth;
			//Save which predecessor we have chosen
			extPth.push_front(predID);
		}// else if(tmpExtSeed != NULL){
		// 	//Testing
		// 	// cout << "4" << endl;

		// 	//Delete the new extension list
		// 	delExtPtr(tmpExtSeed);
		// }

		//Testing
		// if(currScore < maxScore) cout << "2 Option 2" << endl;
	}

	// //Check whether we could find a good extension
	// if(maxScore > 0){
	// 	//Testing
	// 	// cout << "5 Option 1" << endl;

	// 	//Pass hit length of the best extension to the calling function
	// 	hitLen = bestHitLen;
	// 	return bExtSeed;
	// }

	//Testing
	// cout << "5 Option 2" << endl;

	//Save hit length
	hitLen = maxHitLen;
	
	//Return score
	return maxScore;
}

// //This function starts the left extension for seeds lying on the query's reference strand considering a quorum
// void startLeftX_Drop(struct hit* hit, const string &q, const int16_t &X, const uint32_t &quorum){
// 	//Initialization of auxiliary variables
// 	int32_t tmpScore = 0;
// 	uint32_t tmpExtLen = 1, progress;
// 	//Counter to count tries to explore a further unitig
// 	uint32_t explCount;
// 	uint32_t checkedPos = 0;
// 	uint32_t k = hit->lUnitig.getGraph()->getK();
// 	size_t offset;
// 	string uSeq = hit->lUnitig.mappedSequenceToString();
// 	struct seed *nearestSeed, *prevSeed = NULL;
// 	struct Ext_ptr *extPtr = NULL;
// 	//Make a copy of our unitig
// 	UnitigColorMap<seedlist> cpyUnitig = hit->lUnitig;

// 	// } else{
// 	// 	//Change the unitig's strand
// 	// 	cpyUnitig.strand = !cpyUnitig.strand;
// 	// 	//Get the reverse complementary sequence
// 	// 	uSeq = cpyUnitig.referenceUnitigToString();
// 	// }

// 	//Find the nearest seed that we could reach
// 	nearestSeed = searchLeftNeighbor(cpyUnitig.getData()->getData(cpyUnitig)->getSeed(cpyUnitig.strand), hit->lSeedQoff, hit->lSeedUoff, prevSeed);

// 	//cout << "Neighbors searched" << endl;

// 	//Go through the query up to the beginning
// 	while(hit->lSeedQoff >= tmpExtLen){
// 		//Testing
// 		//cout << "We start the sequence comparision" << endl;
// 		// if(nearestSeed != NULL && nearestSeed->offsetU == 52 && nearestSeed->offsetQ == 5155 && nearestSeed->len == 17){
// 		// 	cout << "The seed has been reached before in startLeftX_Drop" << endl;
// 		// 	exit(0);
// 		// }

// 		//Check whether we have reached the next seed
// 		if(nearestSeed != NULL && hit->lSeedQoff - tmpExtLen < nearestSeed->offsetQ + nearestSeed->len){
// 			//Calculate the progress we have incorporating the reached seed
// 			progress = hit->lSeedQoff - tmpExtLen - nearestSeed->offsetQ + 1;//+1, because tmpExtLen is always 1 ahead
// 			//Temporary extension length
// 			tmpExtLen += progress;

// 			//Adjust number of positions checked for quorum fulfillment
// 			if(checkedPos > progress){
// 				//Testing
// 				//cout << "checkedPos:" << checkedPos << endl;

// 				checkedPos -= progress;
// 			} else{
// 				//Testing
// 				//cout << "1 Option 2" << endl;

// 				checkedPos = 0;
// 			}

// 			//Check if our temporary score becomes positive using this hit
// 			if((tmpScore += progress) > 0){//TODO Change this as soon as we use something else than a unit score!
// 				//Update hit length
// 				hit->length += tmpExtLen - 1;
// 				//Update hit's score
// 				hit->score += tmpScore;
// 				//Update hit's left border
// 				hit->lSeedUoff -= (tmpExtLen - 1);
// 				hit->lSeedQoff -= (tmpExtLen - 1);
// 				//Reset temporary hit length
// 				tmpExtLen = 1;
// 				//Reset temporary score
// 				tmpScore = 0;
// 			}

// 			//Check if the reached seed had a predecessor
// 			if(prevSeed != NULL){
// 				//Link predecessor and successor of the reached seed
// 				prevSeed->nextSeed = nearestSeed->nextSeed;
// 			} else{
// 				//Set the reached seed's successor as head of the seed list
// 				cpyUnitig.getData()->getData(cpyUnitig)->setSeed(nearestSeed->nextSeed, cpyUnitig.strand);
// 			}

// 			//Testing
// 			//cout << "Seed has been extracted" << endl;

// 			//Delete the reached seed
// 			free(nearestSeed);
// 			//Reset prevSeed
// 			prevSeed = NULL;
// 			//Search for the next neighbor
// 			nearestSeed = searchLeftNeighbor(cpyUnitig.getData()->getData(cpyUnitig)->getSeed(cpyUnitig.strand), hit->lSeedQoff,  hit->lSeedUoff, prevSeed);

// 			//Testing
// 			//cout << "New seed has been searched" << endl;

// 		} else if(hit->lSeedUoff >= tmpExtLen){//Check whether we have reached the unitig sequence's beginning
// 			//Testing
// 			//cout << "We want to compare the next bases" << endl;

// 			//Check whether quorum has to be checked
// 			if(checkedPos == 0){
// 				//Testing
// 				//cout << "We have to check whether the quorum is still fulfilled" << endl;
// 				//cout << "hit->lSeedUoff:" << hit->lSeedUoff << endl;

// 				//Calculate offset and make sure there is no overflow
// 				offset = (tmpExtLen + k - 1 > hit->lSeedUoff ? 0 : hit->lSeedUoff - (tmpExtLen + k - 1));

// 				//Check whether quorum is still fulfilled
// 				if((checkedPos = checkSearchCrit(cpyUnitig, quorum, offset)) == 0){
// 					//Testing
// 					//cout << "It's not fulfilled anymore" << endl;

// 					break;//TODO If we want to allow to check more than k positions at once we need a new function here!
// 				}

// 				//Testing
// 				// cout << "It is still fulfilled" << endl;
// 				// cout << "checkedPos:" << checkedPos << endl;
// 				// cout << "tmpExtLen + k - 1:" << (tmpExtLen + k - 1) << endl;
// 				// cout << "hit->lSeedUoff:" << hit->lSeedUoff << endl;
// 				// cout << "uSeq:" << uSeq << endl;
// 			}

// 			//Compare the next two bases and check whether this changes temporary score's sign
// 			if((tmpScore += compUScore(uSeq[hit->lSeedUoff - tmpExtLen], q[hit->lSeedQoff - tmpExtLen])) > 0){
// 				//Testing
// 				//cout << "We archieved a positive score" << endl;

// 				//Update hit length
// 				hit->length += tmpExtLen;
// 				//Update hit's score
// 				hit->score += tmpScore;
// 				//Update hit's left border
// 				hit->lSeedUoff -= (tmpExtLen);
// 				hit->lSeedQoff -= (tmpExtLen);
// 				//Reset temporary hit length
// 				tmpExtLen = 0;
// 				//Reset temporary score
// 				tmpScore = 0;
// 			} else if(tmpScore < -X){//Check if our score is already too bad
// 				//Testing
// 				//cout << "X dropped" << endl;

// 				break;
// 			}

// 			//Testing
// 			//cout << "Next 2 bases have been compared" << endl;

// 			++tmpExtLen;
// 			--checkedPos;
// 		} else{
// 			//Testing
// 			//cout << "We going to check the predecessive unitig" << endl;

// 			//Decrease temporary extension length
// 			--tmpExtLen;

// 			//Check if this unitig has a predecessor
// 			if(cpyUnitig.getPredecessors().hasPredecessors()){
// 				//Initialize explCount
// 				explCount = 0;
// 				//Continue the extension on the predecessive unitig
// 				extPtr = extendAtPrevUnitig(cpyUnitig.getPredecessors(), hit->lSeedQoff - tmpExtLen, hit->length, tmpExtLen, q, X, tmpScore, explCount, quorum);
// 			}

// 			//Testing
// 			//cout << "We have checked the predecessive unitig" << endl;

// 			break;
// 		}
// 	}

// 	//Check whether we could successfully extend our seed on the predecessive unitig
// 	if(extPtr != NULL){
// 		//Update hit score
// 		hit->score += extPtr->extSeed->score;
// 		//Update hit's left border info
// 		hit->lSeedUoff = extPtr->extSeed->offsetU;
// 		hit->lSeedQoff = extPtr->extSeed->offsetQ;
// 		hit->lUnitig = extPtr->origUnitig;

// 		// //Check whether we have to include the new border seed into a seed list
// 		// if(!extPtr->inSeedList){
// 		// 	//Incorporate seed from extension ptr into seed list
// 		// 	processExtPtr(extPtr);
// 		// }

// 		//Delete extension pointer
// 		delExtPtr(extPtr);

// 		//Testing
// 		//cout << "We could successfully extend on the predecessive unitig" << endl;

// 		//return true;
// 	}

// 	//Testing
// 	// if(record){
// 	// 	cerr << "Into the other direction the extension ended at u:" << hit->lSeedUoff << " q:" << hit->lSeedQoff << " with score " << hit->score << endl;
// 	// 	record = false;
// 	// }
// 	//cout << "Extension on predecessive unitig was not successful" << endl;

// 	//return false;
// }

// //This function starts the left extension for seeds lying on the query's reverse complement considering a quorum
// void startLeftX_Drop_OnRevComp(struct hit* hit, const string &q, const int16_t &X, const uint32_t &quorum){
// 	//Initialization of auxiliary variables
// 	int32_t tmpScore = 0;
// 	uint32_t tmpExtLen = 1, progress, endBuf = 0;
// 	//Counter to count tries to explore a further unitig
// 	uint32_t explCount;
// 	uint32_t checkedPos = 0;
// 	uint32_t k = hit->lUnitig.getGraph()->getK();
// 	size_t offset;
// 	string uSeq = hit->lUnitig.mappedSequenceToString();
// 	struct seed *nearestSeed, *prevSeed = NULL;
// 	struct Ext_ptr *extPtr = NULL;
// 	//Make a copy of our unitig
// 	UnitigColorMap<seedlist> cpyUnitig = hit->lUnitig;

// 	//Testing
// 	//bool seedReached = true;

// 	//Find the nearest seed that we could reach
// 	nearestSeed = searchLeftNeighbor(cpyUnitig.getData()->getData(cpyUnitig)->getSeed(cpyUnitig.strand), hit->lSeedQoff, hit->lSeedUoff, prevSeed);

// 	//Testing
// 	//if(nearestSeed != NULL) seedReached = false;

// 	//Check how far we should explore the unitig's sequence
// 	if(cpyUnitig.getPredecessors().hasPredecessors()){
// 		endBuf = k - 1;
// 	}

// 	//Go through the query up to the beginning
// 	while(hit->lSeedQoff >= tmpExtLen){
// 		//Check whether we have reached the next seed
// 		if(nearestSeed != NULL && hit->lSeedQoff - tmpExtLen < nearestSeed->offsetQ + nearestSeed->len){
// 			//Testing
// 			//seedReached = true;

// 			//Calculate the progress we have incorporating the reached seed
// 			progress = hit->lSeedQoff - tmpExtLen - nearestSeed->offsetQ + 1;//+1, because tmpExtLen is always 1 ahead
// 			//Temporary extension length
// 			tmpExtLen += progress;

// 			//Check if our temporary score becomes positive using this hit
// 			if((tmpScore += progress) > 0){//TODO Change this as soon as we use something else than a unit score!
// 				//Update hit length
// 				hit->length += tmpExtLen - 1;
// 				//Update hit's score
// 				hit->score += tmpScore;
// 				//Update hit's left border
// 				hit->lSeedUoff -= (tmpExtLen - 1);
// 				hit->lSeedQoff -= (tmpExtLen - 1);
// 				//Reset temporary hit length
// 				tmpExtLen = 1;
// 				//Reset temporary score
// 				tmpScore = 0;
// 			}

// 			//Check if the reached seed had a predecessor
// 			if(prevSeed != NULL){
// 				//Link predecessor and successor of the reached seed
// 				prevSeed->nextSeed = nearestSeed->nextSeed;
// 			} else{
// 				//Set the reached seed's successor as head of the seed list
// 				cpyUnitig.getData()->getData(cpyUnitig)->setSeed(nearestSeed->nextSeed, cpyUnitig.strand);
// 			}

// 			//Delete the reached seed
// 			free(nearestSeed);
// 			//Reset prevSeed
// 			prevSeed = NULL;
// 			//Search for the next neighbor
// 			nearestSeed = searchLeftNeighbor(cpyUnitig.getData()->getData(cpyUnitig)->getSeed(cpyUnitig.strand), hit->lSeedQoff,  hit->lSeedUoff, prevSeed);

// 			//Testing
// 			//if(nearestSeed != NULL) seedReached = false;
// 		} else if(hit->lSeedUoff >= tmpExtLen + endBuf){//Check whether we have reached the unitig sequence's beginning
// 			//Check whether quorum has to be checked
// 			if(checkedPos == 0){
// 				//Calculate offset and make sure there is no overflow
// 				offset = (tmpExtLen + k - 1 > hit->lSeedUoff ? 0 : hit->lSeedUoff - (tmpExtLen + k - 1));//TODO Theoretically we don't need this here!

// 				//Check whether quorum is still fulfilled
// 				if((checkedPos = checkSearchCrit(cpyUnitig, quorum, offset)) == 0){
// 					break;//TODO If we want to allow to check more than k positions at once we need a new function here!
// 				}
// 			}

// 			//Compare the next two bases and check whether this changes temporary score's sign
// 			if((tmpScore += compUScore(uSeq[hit->lSeedUoff - tmpExtLen], q[hit->lSeedQoff - tmpExtLen])) > 0){
// 				//Update hit length
// 				hit->length += tmpExtLen;
// 				//Update hit's score
// 				hit->score += tmpScore;
// 				//Update hit's left border
// 				hit->lSeedUoff -= (tmpExtLen);
// 				hit->lSeedQoff -= (tmpExtLen);
// 				//Reset temporary hit length
// 				tmpExtLen = 0;
// 				//Reset temporary score
// 				tmpScore = 0;
// 			} else if(tmpScore < -X){//Check if our score is already too bad
// 				break;
// 			}

// 			++tmpExtLen;
// 			--checkedPos;
// 		} else{
// 			//Decrease temporary extension length
// 			--tmpExtLen;

// 			//Check if this unitig has a predecessor
// 			if(endBuf != 0){
// 				//Testing
// 				// if(!seedReached){
// 				// 	cerr << "A reachable seed has not been reached during left extension on the unitig we started!" << endl;
// 				// 	cerr << "Seed is uoff:" << nearestSeed->offsetU << " qoff:" << nearestSeed->offsetQ << " len:" << nearestSeed->len << endl;
// 				// 	cerr << "lSeedQoff:" << hit->lSeedQoff << " tmpExtLen:" << tmpExtLen << endl;
// 				// 	//exit(0);
// 				// }

// 				//Initialize explCount
// 				explCount = 0;
// 				//Continue the extension on the predecessive unitig
// 				extPtr = extendAtPrevUnitigOnRevComp(cpyUnitig.getPredecessors(), hit->lSeedQoff - tmpExtLen, hit->length, tmpExtLen, q, X, tmpScore, explCount, quorum);
// 			}

// 			break;
// 		}
// 	}

// 	//Check whether we could successfully extend our seed on the predecessive unitig
// 	if(extPtr != NULL){
// 		//Update hit score
// 		hit->score += extPtr->extSeed->score;
// 		//Update hit's left border info
// 		hit->lSeedUoff = extPtr->extSeed->offsetU;
// 		hit->lSeedQoff = extPtr->extSeed->offsetQ;
// 		hit->lUnitig = extPtr->origUnitig;
// 		//Delete extension pointer
// 		delExtPtr(extPtr);
// 	}
// }

//This function starts the left extension for seeds lying on the query's reference strand considering a quorum and a search color set
void startLeftX_Drop(hit* hit, const string &q, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	//Initialization of auxiliary variables
	int32_t tmpScore = 0;
	uint32_t tmpExtLen = 1, progress, posU = hit->offU, posQ = hit->offQ;
	//Counter to count tries to explore a further unitig
	uint32_t explCount;
	//Get the number of covered positions coming from the sequence's end considering that we do not start in the very end
	int32_t checkedPos = getSrchCritCov(hit->origUni, quorum, searchSet, posU - tmpExtLen, false);// - (hit->origUni.size - posU);

	// uint32_t k = hit->origUni.getGraph()->getK();
	// size_t offset;

	//Testing
	// if(hit->origUni.mappedSequenceToString() == "CTCAGTTTACTCATACCCCATTCATCAGCAATTGAAAACAAA"){
	// 	cout << "We have found the interesting unitig and we have to move the starting offset" << endl;
	// }

	string uSeq = hit->origUni.mappedSequenceToString();
	list<uint16_t> extPth;
	struct seed *nearestSeed, *prevSeed = NULL;

	// struct Ext_ptr *extPtr = NULL;
	// //Make a copy of our unitig
	// UnitigColorMap<seedlist> cpyUnitig = hit->lUnitig;

	// //Check which sequence we have to look on
	// if(cpyUnitig.strand){
	// 	//Just use the sequence of our unitig
	// 	uSeq = cpyUnitig.referenceUnitigToString();
	// } else{
	// 	//Change the unitig's strand
	// 	cpyUnitig.strand = !cpyUnitig.strand;
	// 	//Get the reverse complementary sequence
	// 	uSeq = cpyUnitig.referenceUnitigToString();
	// }

	//Find the nearest seed that we could reach
	nearestSeed = searchLeftNeighbor(hit->origUni.getData()->getData(hit->origUni)->getSeed(hit->origUni.strand), hit->offQ, hit->offU, prevSeed);

	//Go through the query up to the beginning
	while(posQ >= tmpExtLen){
		//Check whether we have reached the next seed
		if(nearestSeed != NULL && posQ - tmpExtLen < nearestSeed->offsetQ + nearestSeed->len){
			//Calculate the progress we have incorporating the reached seed
			progress = posQ - tmpExtLen - nearestSeed->offsetQ + 1;//+1, because tmpExtLen is always 1 ahead
			//Temporary extension length
			tmpExtLen += progress;
			//Adjust number of remaining covered positions
			checkedPos -= progress;

			//Check if our temporary score becomes positive using this hit
			if((tmpScore += progress) > 0){//TODO Change this as soon as we use something else than a unit score!
				//Update hit length
				hit->length += tmpExtLen - 1;
				//Update hit's score
				hit->score += tmpScore;
				//Update hit's left border
				posU -= (tmpExtLen - 1);
				posQ -= (tmpExtLen - 1);
				//Reset temporary hit length
				tmpExtLen = 1;
				//Reset temporary score
				tmpScore = 0;
			}

			// //Delete the reached seed
			// extractUnprocSeed(currUnitig->getData()->getData(*currUnitig), nearestSeed, prevSeed);

			//Check if the reached seed had a predecessor
			if(prevSeed != NULL){
				//Testing
				// cout << "11 Option 1" << endl;

				//Link predecessor and successor of the reached seed
				prevSeed->nextSeed = nearestSeed->nextSeed;
			} else{
				//Testing
				// cout << "11 Option 2" << endl;

				//Set the reached seed's successor as head of the seed list
				hit->origUni.getData()->getData(hit->origUni)->setSeed(nearestSeed->nextSeed, hit->origUni.strand);
			}

			//Delete the reached seed
			free(nearestSeed);
			//Reset prevSeed
			prevSeed = NULL;
			//Search for the next neighbor
			nearestSeed = searchLeftNeighbor(hit->origUni.getData()->getData(hit->origUni)->getSeed(hit->origUni.strand), posQ,  posU, prevSeed);
		} else if(posU >= tmpExtLen){//Check whether we have reached the unitig sequence's beginning
			//Check whether quorum has to be checked
			if(checkedPos == 0){
				break;

				// //Check whether quorum is still fulfilled
				// if((checkedPos = checkSearchCrit(hit->origUni, quorum, offset, searchSet)) == 0) break;//TODO If we want to allow to check more than k positions at once we need a new function here!
			}

			//Compare the next two bases and check whether this changes temporary score's sign
			if((tmpScore += compUScore(uSeq[posU - tmpExtLen], q[posQ - tmpExtLen])) > 0){
				//Update hit length
				hit->length += tmpExtLen;
				//Update hit's score
				hit->score += tmpScore;
				//Update hit's left border
				posU -= (tmpExtLen);
				posQ -= (tmpExtLen);
				//Reset temporary hit length
				tmpExtLen = 0;
				//Reset temporary score
				tmpScore = 0;
			} else if(tmpScore < -X){//Check if our score is already too bad
				break;
			}
			++tmpExtLen;
			--checkedPos;
		} else{
			//Decrease temporary extension length
			--tmpExtLen;

			//Check if this unitig has a predecessor
			if(hit->origUni.getPredecessors().hasPredecessors()){
				//Testing
				// cout << "12 Option 1" << endl;

				//Initialize explCount
				explCount = 0;
				//Continue the extension on the predecessive unitig
				hit->score += extendAtPrevUnitig(hit->origUni.getPredecessors(), posQ - tmpExtLen, hit->length, tmpExtLen, q, X, tmpScore, extPth, explCount, quorum, searchSet);
			}// else{
			// 	//Testing
			// 	cout << "12 Option 2" << endl;
			// }

			break;
		}
	}

	// //Check whether we could successfully extend our seed on the predecessive unitig
	// if(extPtr != NULL){
	// 	//Update hit score
	// 	hit->score += extPtr->extSeed->score;
	// 	//Update hit's left border info
	// 	hit->lSeedUoff = extPtr->extSeed->offsetU;
	// 	hit->lSeedQoff = extPtr->extSeed->offsetQ;
	// 	hit->lUnitig = extPtr->origUnitig;

	// 	// //Check whether we have to include the new border seed into a seed list
	// 	// if(!extPtr->inSeedList){
	// 	// 	//Incorporate seed from extension ptr into seed list
	// 	// 	processExtPtr(extPtr);
	// 	// }

	// 	//Delete extension pointer
	// 	delExtPtr(extPtr);

	// 	//return true;
	// }

	hit->lExt = cmprExtPth(extPth);

	//return false;
}

//This function starts the left extension for seeds lying on the query's reverse complement considering a quorum and a search color set
void startLeftX_Drop_OnRevComp(hit* hit, const string &q, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	//Initialization of auxiliary variables
	int32_t tmpScore = 0;
	uint32_t tmpExtLen = 1, progress, overlap = 0;
	//Counter to count tries to explore a further unitig
	uint32_t explCount;
	int32_t checkedPos = getSrchCritCov(hit->origUni, quorum, searchSet, compOffset(hit->offU - tmpExtLen, 1, hit->origUni.size, false), true);
	uint32_t k = hit->origUni.getGraph()->getK();
	list<uint16_t> extPth;
	string uSeq = hit->origUni.mappedSequenceToString();
	struct seed *nearestSeed, *prevSeed = NULL;

	//Testing
	// cout << "Start startLeftX_Drop_OnRevComp" << endl;
	// if(hit->origUni.mappedSequenceToString() == "AACTTTTATATTTGTTTTCAATTGCTGATGAATGGGGTAT"){
	// 	report = true;
	// }

	// //Make a copy of our unitig
	// UnitigColorMap<seedlist> cpyUnitig = hit->lUnitig;

	//Find the nearest seed that we could reach
	nearestSeed = searchLeftNeighbor(hit->origUni.getData()->getData(hit->origUni)->getSeed(hit->origUni.strand), hit->offQ, hit->offU, prevSeed);

	//Testing
	// cout << "Neighbors searched" << endl;
	
	//Check how far we should explore the unitig's sequence
	if(hit->origUni.getPredecessors().hasPredecessors()){
		//Testing
		// cout << "1 Option 1" << endl;

		overlap = k - 1;
	}// else{
	// 	//Testing
	// 	// cout << "1 Option 2" << endl;
	// }

	//Go through the query up to the beginning
	while(hit->offQ >= tmpExtLen){
		//Testing
		// cout << "We are not at the query's beginning" << endl;

		//Check whether we have reached the next seed
		if(nearestSeed != NULL && hit->offQ - tmpExtLen < nearestSeed->offsetQ + nearestSeed->len){
			//Testing
			// cout << "2 Option 1" << endl;
			// cout << "Seed reached" << endl;

			//Calculate the progress we have incorporating the reached seed
			progress = hit->offQ - tmpExtLen - nearestSeed->offsetQ + 1;//+1, because tmpExtLen is always 1 ahead
			//Temporary extension length
			tmpExtLen += progress;

			// if(checkedPos > progress){

			//Adjust number of remaining covered positions
			checkedPos -= progress;

			// } else{
			// 	checkedPos = 0;
			// }

			//Check if our temporary score becomes positive using this hit
			if((tmpScore += progress) > 0){//TODO Change this as soon as we use something else than a unit score!
				//Testing
				// cout << "3 Option 1" << endl;

				//Update hit length
				hit->length += tmpExtLen - 1;
				//Update hit's score
				hit->score += tmpScore;
				//Update hit's left border
				hit->offU -= (tmpExtLen - 1);
				hit->offQ -= (tmpExtLen - 1);
				//Reset temporary hit length
				tmpExtLen = 1;
				//Reset temporary score
				tmpScore = 0;
			}// else{
			// 	//Testing
			// 	// cout << "3 Option 2" << endl;
			// }

			//Check if the reached seed had a predecessor
			if(prevSeed != NULL){
				//Testing
				// cout << "4 Option 1" << endl;

				//Link predecessor and successor of the reached seed
				prevSeed->nextSeed = nearestSeed->nextSeed;
			} else{
				//Testing
				// cout << "4 Option 2" << endl;

				//Set the reached seed's successor as head of the seed list
				hit->origUni.getData()->getData(hit->origUni)->setSeed(nearestSeed->nextSeed, hit->origUni.strand);
			}

			//Delete the reached seed
			free(nearestSeed);
			//Reset prevSeed
			prevSeed = NULL;
			//Search for the next neighbor
			nearestSeed = searchLeftNeighbor(hit->origUni.getData()->getData(hit->origUni)->getSeed(hit->origUni.strand), hit->offQ,  hit->offU, prevSeed);
		} else if(hit->offU >= overlap + tmpExtLen){//Check whether we have reached the unitig sequence's beginning
			//Testing
			// cout << "2 Option 2" << endl;
			// cout << "6 Option 2" << endl;
			// cout << "Let's compare bases" << endl;

			//Check whether quorum has to be checked
			if(checkedPos <= 0){
				//Testing
				// cout << "5 Option 1" << endl;

				// //Calculate offset and make sure there is no overflow
				// offset = (posU - tmpExtLen < k - 1 ? hit->origUni.size - (k - 1) : hit->origUni.size - (posU - tmpExtLen) - 1);

				// //Check whether quorum is still fulfilled
				// if((checkedPos = checkSearchCrit(hit->origUni, quorum, offset, searchSet)) == 0){
				// 	//Testing
				// 	// cout << "7 Option 2" << endl;

				break;

				// }// else{
				// 	//Testing
				// 	// cout << "7 Option 1" << endl;
				// }
			}// else{
			// 	//Testing
			// 	// cout << "5 Option 2" << endl;
			// }

			//Testing
			// if(hit->origUni.mappedSequenceToString() == "TCTTTACGGCGAAGTTCAGCGCCCTCATAGCC" && hit->offU == 20 && hit->offQ == 47){
			// 	cout << "Bases to be compared:" << endl;
			// 	cout << "u " << posU + tmpExtLen << " : " << uSeq[posU - tmpExtLen] << " " << q[posQ - tmpExtLen] << " : " << posQ - tmpExtLen << " q" << endl;
			// 	report = true;
			// }

			//Compare the next two bases and check whether this changes temporary score's sign
			if((tmpScore += compUScore(uSeq[hit->offU - tmpExtLen], q[hit->offQ - tmpExtLen])) > 0){
				//Testing
				// cout << "8 Option 1" << endl;

				//Update hit length
				hit->length += tmpExtLen;
				//Update hit's score
				hit->score += tmpScore;
				//Update hit's left border
				hit->offU -= tmpExtLen;
				hit->offQ -= tmpExtLen;
				//Reset temporary hit length
				tmpExtLen = 0;
				//Reset temporary score
				tmpScore = 0;
			} else if(tmpScore < -X){//Check if our score is already too bad
				//Testing
				// cout << "8 Option 3" << endl;

				break;
			}// else{
			// 	//Testing
			// 	// cout << "8 Option 2" << endl;
			// }

			++tmpExtLen;
			--checkedPos;
		} else{
			//Testing
			// cout << "2 Option 2" << endl;
			// cout << "6 Option 1" << endl;
			// cout << "Go to the next unitig" << endl;

			//Decrease temporary extension length
			--tmpExtLen;

			//Check if this unitig has a predecessor
			if(overlap != 0){
				//Testing
				// cout << "9 Option 1" << endl;
				// cout << "startLeftX_Drop_OnRevComp: Continue on predecessive unitig" << endl;

				//Initialize explCount
				explCount = 0;
				//Continue the extension on the predecessive unitig
				hit->score += extendAtPrevUnitigOnRevComp(hit->origUni.getPredecessors(), hit->offQ - tmpExtLen, hit->length, tmpExtLen, q, X, tmpScore, extPth, overlap - (hit->offU - tmpExtLen), explCount, quorum, searchSet);
			}// else{
			// 	//Testing
			// 	// cout << "9 Option 2" << endl;
			// }

			break;
		}
	}

	//Testing
	// cout << "Calculations finished" << endl;
	// if(hit->origUni.mappedSequenceToString() == "TTTGTTTTCAATTGCTGATGAATGGGGTATGAGTAAACTGAG"){
	// 	cout << "We have found the interesting unitig and we have to move the starting offset" << endl;
	// }

	// //Check whether we could successfully extend our seed on the predecessive unitig
	// if(extPtr != NULL){
	// 	//Testing
	// 	// cout << "10 Option 1" << endl;

	// 	//Update hit score
	// 	hit->score += extPtr->extSeed->score;
	// 	//Update hit's left border info
	// 	hit->lSeedUoff = extPtr->extSeed->offsetU;
	// 	hit->lSeedQoff = extPtr->extSeed->offsetQ;
	// 	hit->lUnitig = extPtr->origUnitig;
	// 	//Delete extension pointer
	// 	delExtPtr(extPtr);
	// } else{
	// 	//Testing
	// 	// cout << "10 Option 2" << endl;
	// }

	//Testing
	// cerr << "startLeftX_Drop_OnRevComp: origUni: " << hit->origUni.mappedSequenceToString() << endl;
	// if(hit->origUni.mappedSequenceToString() == "CAGTACGGTATCGGCCCCCAACGCAATCATGCGCACAACGTCCAGACCGTTACGGATCCCGCTGTCTGCCAGAATGGTGATGTCGCCTTTCACCGCATCGGCAATGGCGGG"){
	// 	cerr << "After left extension: length: " << hit->length << " score: " << hit->score << " left extension path: " << endl;
	// 	for(list<uint16_t>::const_iterator i = extPth.begin(); i != extPth.end(); ++i){
	// 		cerr << *i << " ";
	// 	}
	// 	cerr << endl;
	// 	report = true;
	// }

	//A hit's start offset must never be inside the overlap at a unitig sequence's end
	if(hit->offU > hit->origUni.size - k){
		//Testing
		// cout << "startLeftX_Drop_OnRevComp: Moving start" << endl;

		mvStartToValUni(hit, extPth);
	}

	//Testing
	// cout << "startLeftX_Drop_OnRevComp: Checked start offset" << endl;
	// if(hit->origUni.mappedSequenceToString() == "TTTCACCGCATCGGCAATGGCGGGCAGGGCG"){
	// 	cerr << "startLeftX_Drop_OnRevComp: We have found the interesting unitig" << endl;
	// 	exit(0);
	// }
	
	//If hit is invalid we do not need to compress its left extension path
	if(hit->score > 0){
		//Testing
		// cout << "startLeftX_Drop_OnRevComp: Compressing extension path" << endl;

		hit->lExt = cmprExtPth(extPth);
	}


	//Testing
	// cout << "Hit set" << endl;
	// cout << "startLeftX_Drop_OnRevComp: End of function" << endl;
}

// //This function continues a left extension on a predecessive unitig of a seed lying on the query's reference strand considering a quorum and returns the achieved score
// int32_t contLeftX_Drop(const neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> &prevUni, uint32_t qPos, uint32_t &hitLen, const string &q, const int16_t &X, struct Ext_ptr*& extPtr, const int32_t &lastSeedTmpScore, uint32_t &explCount, const uint32_t &quorum){
// 	//bool inSeedList;
// 	int32_t tmpScore;
// 	uint32_t tmpExtLen, uPos;
// 	uint32_t checkedPos = 0;
// 	uint32_t k = prevUni->getGraph()->getK();
// 	size_t offset;
// 	string prevUniSeq = prevUni->mappedSequenceToString();
// 	struct seed *nearestSeed, *prevSeed = NULL, *extSeed = NULL;

// 	//Calculate our offset position inside the unitig sequence (+1, because otherwise it is not possible to calculate the correct gain if reaching a seed)
// 	uPos = prevUni->size - prevUni->getGraph()->getK() + 1;

// 	// //Check whether we are working on the reference strand of the current unitig
// 	// if(prevUni->strand){

// 	//Testing
// 	// if(report){
// 	// 	cerr << "contLeftX_Drop" << endl;
// 	// }

// 	//Find the nearest seed that we might be able to reach during our extension
// 	nearestSeed = searchLeftNeighbor(prevUni->getData()->getData(*prevUni)->getSeed(prevUni->strand), qPos, uPos, prevSeed);

// 	// } else{
// 	// 	//So far we have not considered the reverse complementary strand for the seed detection //TODO: Change this!
// 	// 	nearestSeed = NULL;
// 	// }

// 	//Initialize our temporary score with the last temporary score of its successive unitig
// 	tmpScore = lastSeedTmpScore;
// 	//Initialize the temporary extension length
// 	tmpExtLen = 1;

// 	//We are done if we reach the beginning of the query
// 	while(qPos >= tmpExtLen){
// 		//Testing
// 		// if(nearestSeed != NULL && nearestSeed->offsetU == 52 && nearestSeed->offsetQ == 5155 && nearestSeed->len == 17){
// 		// 	cout << "The seed has been reached before in contLeftX_Drop" << endl;
// 		// 	exit(0);
// 		// }

// 		//Check whether we have reached a nearest seed
// 		if(nearestSeed != NULL && qPos - tmpExtLen <= nearestSeed->offsetQ + nearestSeed->len - 1){
// 			//Update the temporary score
// 			tmpScore += qPos - tmpExtLen - nearestSeed->offsetQ + 1;//TODO Change this as soon as we want to allow more than only unit score!
// 			//Calculate the gain we have by incorporating the seed into our extension and update the temporary extension length
// 			tmpExtLen = qPos - nearestSeed->offsetQ + 1;

// 			//Testing
// 			// if(prevUniSeq == "CGAGCGCCTGTACTGATCCACTCGACGCCTACGCGTCAGTGATCGCGCCTCAGCCCCGCCCCCGCC"){
// 			// 	cerr << "We have reached our interesting seed during left extension" << endl;
// 			// }

// 			//Check whether we get a score larger 0 by incorporating the reached seed
// 			if(tmpScore > 0){
// 				//Update hit's length
// 				hitLen += qPos - nearestSeed->offsetQ;

// 				//Check whether we already have an extension seed
// 				if(extSeed != NULL){
// 					//Update extension seed
// 					extSeed->offsetU = nearestSeed->offsetU;
// 					extSeed->offsetQ = nearestSeed->offsetQ;
// 					extSeed->score += tmpScore;
// 					extSeed->len += qPos - nearestSeed->offsetQ;
// 					//Update current position in q the unitig
// 					qPos = extSeed->offsetQ;
// 					uPos = extSeed->offsetU;

// 					//Check if the reached seed has a predecessor in its seed list
// 					if(prevSeed != NULL){
// 						//Link the reached seed's predecessor and successor
// 						prevSeed->nextSeed = nearestSeed->nextSeed;
// 					} else{
// 						//Set the reached seed's successor as the head of the seed list
// 						prevUni->getData()->getData(*prevUni)->setSeed(nearestSeed->nextSeed, prevUni->strand);
// 					}
// 					//Delete the reached seed
// 					free(nearestSeed);
					
// 					//extractUnprocSeed(prevUni->getData()->getData(*prevUni), nearestSeed, prevSeed);

// 					// //Reset prevSeed
// 					// prevSeed = NULL;
// 				} else{
// 					//Make the reached seed our new extension seed
// 					extSeed = nearestSeed;
// 					//Update seed's length and score (we have to substract 1 here because we initialize the temporary extension length with 1)
// 					extSeed->len = tmpExtLen - 1;
// 					extSeed->score = tmpScore;

// 					// //Mark that we recycle a seed instead of creating a new one for our extension
// 					// inSeedList = true;

// 					//TODO Check whether we can export this into a function!
// 					//Check whether the reached seed is the first in its list
// 					if(prevSeed != NULL){
// 						//Exclude the reached seed from its list
// 						prevSeed->nextSeed = extSeed->nextSeed;
// 					} else{
// 						//Exclude the reached seed from its list
// 						prevUni->getData()->getData(*prevUni)->setSeed(nearestSeed->nextSeed, prevUni->strand);
// 					}

// 					// //Incorporate seed into the list of processed seeds
// 					// extSeed->nextSeed = prevUni->getData()->getData(*prevUni)->getLastProcSeed();
// 					// prevUni->getData()->getData(*prevUni)->setLastProcSeed(extSeed);

// 					//Update current position in q the unitig
// 					qPos = extSeed->offsetQ;
// 					uPos = extSeed->offsetU;
// 				}

// 				tmpScore = 0;
// 			} else{
// 				//Check if the reached seed has a predecessor in its seed list
// 				if(prevSeed != NULL){
// 					//Link the reached seed's predecessor and successor
// 					prevSeed->nextSeed = nearestSeed->nextSeed;
// 				} else{
// 					//Set the reached seed's successor as the head of the seed list
// 					prevUni->getData()->getData(*prevUni)->setSeed(nearestSeed->nextSeed, prevUni->strand);
// 				}
// 				//Delete the reached seed
// 				free(nearestSeed);
				
// 				// //Delete the reached seed
// 				// extractUnprocSeed(prevUni->getData()->getData(*prevUni), nearestSeed, prevSeed);
// 			}

// 			//Reset prevSeed
// 			prevSeed = NULL;
// 			//Search for the next seed to reach
// 			nearestSeed = searchLeftNeighbor(prevUni->getData()->getData(*prevUni)->getSeed(prevUni->strand), qPos, uPos, prevSeed);
// 		} else if(uPos >= tmpExtLen){//Check whether we have already reached the beginning of the unitig sequence
// 			//Check whether quorum has to be checked
// 			if(checkedPos == 0){
// 				//Calculate offset and make sure there is no overflow
// 				offset = (tmpExtLen + k - 1 > uPos ? 0 : uPos - (tmpExtLen + k - 1));

// 				//Check whether quorum is still fulfilled
// 				if((checkedPos = checkSearchCrit(*prevUni, quorum, offset)) == 0) break;//TODO If we want to allow to check more than k positions at once we need a new function here!
// 			}

// 			//Compare the next two bases and check whether our temporary score is becoming > 0 by this
// 			if((tmpScore += compUScore(prevUniSeq[uPos - tmpExtLen], q[qPos - tmpExtLen])) > 0){
// 				//Check whether there is no extension seed yet
// 				if(extSeed == NULL){
// 					//Create a new extension seed
// 					extSeed = (struct seed*) malloc(sizeof(struct seed));
// 					extSeed->len = 0;
// 					extSeed->score = 0;

// 					// //Mark that this extension seed is not part of a seed list yet
// 					// inSeedList = false;
// 				}

// 				//Update positions in q and the unitig
// 				uPos -= tmpExtLen;
// 				qPos -= tmpExtLen;
// 				//Update extension seed's infos
// 				extSeed->offsetU = uPos;
// 				extSeed->offsetQ = qPos;
// 				extSeed->len += tmpExtLen;
// 				extSeed->score += tmpScore;
// 				//Update hit's length
// 				hitLen += tmpExtLen;
// 				//Reset temporary length and score
// 				tmpExtLen = 0;
// 				tmpScore = 0;
// 			} else if(tmpScore < -X){//Check whether our temporary score is already too negative
// 				break;
// 			}

// 			//Increment temporary extension length
// 			++tmpExtLen;
// 			--checkedPos;
// 		} else{
// 			//Decrease tmpExtLen 
// 			--tmpExtLen;

// 			//Check if this unitig has a predecessor
// 			if(prevUni->getPredecessors().hasPredecessors()){
// 				//Try to continue the extension on the predecessive unitigs
// 				extPtr = extendAtPrevUnitig(prevUni->getPredecessors(), qPos - tmpExtLen, hitLen, tmpExtLen, q, X, tmpScore, explCount, quorum);
// 			}

// 			break;
// 		}
// 	}

// 	//Check whether we have an extension seed
// 	if(extSeed != NULL){
// 		//Check whether there exists a successful extension on a predecessive unitig
// 		if(extPtr == NULL){
// 			//Create a new extension pointer
// 			extPtr = (struct Ext_ptr*) malloc(sizeof(struct Ext_ptr));
// 			//Link this unitig to the extension pointer
// 			extPtr->origUnitig = *prevUni;
// 			//Link the extension seed
// 			extPtr->extSeed = extSeed;

// 			// //Save whether the linked seed still needs to be incorporated into a seed list
// 			// extPtr->inSeedList = inSeedList;
// 		} else{
// 			//Add the score that we gained for this unitig to the score we gained extending its predecessors
// 			extPtr->extSeed->score += extSeed->score;
// 			//Delete extension seed for this unitig
// 			free(extSeed);

// 			//delExtPtrSeed(inSeedList, extSeed, prevUni->getData()->getData(*prevUni));
// 		}
// 	}

// 	//Check whether we have an extension pointer
// 	if(extPtr != NULL){ return extPtr->extSeed->score; }

// 	return 0;
// }

// //This function continues a left extension on a predecessive unitig of a seed lying on the query's reverse complement considering a quorum and returns the achieved score
// int32_t contLeftX_DropOnRevComp(const neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> &prevUni, uint32_t qPos, uint32_t &hitLen, const string &q, const int16_t &X, struct Ext_ptr*& extPtr, const int32_t &lastSeedTmpScore, uint32_t &explCount, const uint32_t &quorum){
// 	int32_t tmpScore;
// 	uint32_t tmpExtLen, uPos = prevUni->size, overlap = 0;
// 	uint32_t checkedPos = 0;
// 	uint32_t k = prevUni->getGraph()->getK();
// 	size_t offset;
// 	string prevUniSeq = prevUni->mappedSequenceToString();
// 	struct seed *nearestSeed, *prevSeed = NULL, *extSeed = NULL;

// 	//Testing
// 	//bool seedReached = true;

// 	//Find the nearest seed that we might be able to reach during our extension
// 	nearestSeed = searchLeftNeighbor(prevUni->getData()->getData(*prevUni)->getSeed(prevUni->strand), qPos, uPos, prevSeed);

// 	//Testing
// 	//if(nearestSeed != NULL) seedReached = false;

// 	//Check how far we should explore the unitig's sequence
// 	if(prevUni->getPredecessors().hasPredecessors()){
// 		overlap = k - 1;
// 	}

// 	//Initialize our temporary score with the last temporary score of its successive unitig
// 	tmpScore = lastSeedTmpScore;
// 	//Initialize the temporary extension length
// 	tmpExtLen = 1;

// 	//We are done if we reach the beginning of the query
// 	while(qPos >= tmpExtLen){
// 		//Check whether we have reached a nearest seed
// 		if(nearestSeed != NULL && qPos - tmpExtLen <= nearestSeed->offsetQ + nearestSeed->len - 1){
// 			//Testing
// 			//seedReached = true;

// 			//Update the temporary score
// 			tmpScore += qPos - tmpExtLen - nearestSeed->offsetQ + 1;//TODO Change this as soon as we want to allow more than only unit score!
// 			//Calculate the gain we have by incorporating the seed into our extension and update the temporary extension length
// 			tmpExtLen = qPos - nearestSeed->offsetQ + 1;

// 			//Check whether we get a score larger 0 by incorporating the reached seed
// 			if(tmpScore > 0){
// 				//Update hit's length
// 				hitLen += qPos - nearestSeed->offsetQ;

// 				//Check whether we already have an extension seed
// 				if(extSeed != NULL){
// 					//Update extension seed
// 					extSeed->offsetU = nearestSeed->offsetU;
// 					extSeed->offsetQ = nearestSeed->offsetQ;
// 					extSeed->score += tmpScore;
// 					extSeed->len += qPos - nearestSeed->offsetQ;
// 					//Update current position in q the unitig
// 					qPos = extSeed->offsetQ;
// 					uPos = extSeed->offsetU;

// 					//Check if the reached seed has a predecessor in its seed list
// 					if(prevSeed != NULL){
// 						//Link the reached seed's predecessor and successor
// 						prevSeed->nextSeed = nearestSeed->nextSeed;
// 					} else{
// 						//Set the reached seed's successor as the head of the seed list
// 						prevUni->getData()->getData(*prevUni)->setSeed(nearestSeed->nextSeed, prevUni->strand);
// 					}
// 					//Delete the reached seed
// 					free(nearestSeed);
// 				} else{
// 					//Make the reached seed our new extension seed
// 					extSeed = nearestSeed;
// 					//Update seed's length and score (we have to substract 1 here because we initialize the temporary extension length with 1)
// 					extSeed->len = tmpExtLen - 1;
// 					extSeed->score = tmpScore;

// 					//TODO Check whether we can export this into a function!
// 					//Check whether the reached seed is the first in its list
// 					if(prevSeed != NULL){
// 						//Exclude the reached seed from its list
// 						prevSeed->nextSeed = extSeed->nextSeed;
// 					} else{
// 						//Exclude the reached seed from its list
// 						prevUni->getData()->getData(*prevUni)->setSeed(nearestSeed->nextSeed, prevUni->strand);
// 					}

// 					//Update current position in q the unitig
// 					qPos = extSeed->offsetQ;
// 					uPos = extSeed->offsetU;
// 				}

// 				tmpScore = 0;
// 				tmpExtLen = 1;
// 			} else{
// 				//Check if the reached seed has a predecessor in its seed list
// 				if(prevSeed != NULL){
// 					//Link the reached seed's predecessor and successor
// 					prevSeed->nextSeed = nearestSeed->nextSeed;
// 				} else{
// 					//Set the reached seed's successor as the head of the seed list
// 					prevUni->getData()->getData(*prevUni)->setSeed(nearestSeed->nextSeed, prevUni->strand);
// 				}
// 				//Delete the reached seed
// 				free(nearestSeed);
// 			}
// 			//Testing
// 			//cout << "2" << endl;

// 			//Reset prevSeed
// 			prevSeed = NULL;
// 			//Search for the next seed to reach
// 			nearestSeed = searchLeftNeighbor(prevUni->getData()->getData(*prevUni)->getSeed(prevUni->strand), qPos, uPos, prevSeed);

// 			//Testing
// 			//if(nearestSeed != NULL) seedReached = false;
// 		} else if(uPos >= tmpExtLen + overlap){//Check whether we have already reached the beginning of the unitig sequence
// 			//Check whether quorum has to be checked
// 			if(checkedPos == 0){
// 				//Calculate offset and make sure there is no overflow
// 				offset = (tmpExtLen + k - 1 > uPos ? 0 : uPos - (tmpExtLen + k - 1));//TODO Theoretically, this is not necessary!

// 				//Check whether quorum is still fulfilled
// 				if((checkedPos = checkSearchCrit(*prevUni, quorum, offset)) == 0) break;//TODO If we want to allow to check more than k positions at once we need a new function here!
// 			}

// 			//Compare the next two bases and check whether our temporary score is becoming > 0 by this
// 			if((tmpScore += compUScore(prevUniSeq[uPos - tmpExtLen], q[qPos - tmpExtLen])) > 0){
// 				//Check whether there is no extension seed yet
// 				if(extSeed == NULL){
// 					//Create a new extension seed
// 					extSeed = (struct seed*) malloc(sizeof(struct seed));
// 					extSeed->len = 0;
// 					extSeed->score = 0;
// 				}

// 				//Update positions in q and the unitig
// 				uPos -= tmpExtLen;
// 				qPos -= tmpExtLen;
// 				//Update extension seed's infos
// 				extSeed->offsetU = uPos;
// 				extSeed->offsetQ = qPos;
// 				extSeed->len += tmpExtLen;
// 				extSeed->score += tmpScore;
// 				//Update hit's length
// 				hitLen += tmpExtLen;
// 				//Reset temporary length and score
// 				tmpExtLen = 0;
// 				tmpScore = 0;
// 			} else if(tmpScore < -X){//Check whether our temporary score is already too negative
// 				break;
// 			}

// 			//Increment temporary extension length
// 			++tmpExtLen;
// 			--checkedPos;
// 		} else{
// 			//Decrease tmpExtLen
// 			--tmpExtLen;

// 			//Testing
// 			// if(!seedReached){
// 			// 	cerr << "We have not reached a reachable seed during continued left extension!" << endl;
// 			// 	cerr << "Seed is uoff:" << nearestSeed->offsetU << " qoff:" << nearestSeed->offsetQ << " len:" << nearestSeed->len << endl;
// 			// 	cerr << "qPos:" << qPos << " tmpExtLen:" << tmpExtLen << endl << "We are on unitig " << prevUniSeq << endl << "We are " << (prevUni->strand ? "" : "not ") << "on the reference strand" << endl;
// 			// 	exit(0);
// 			// }

// 			//Check if this unitig has a predecessor
// 			if(overlap != 0){
// 				//Testing
// 				// if(!seedReached){
// 				// 	cerr << "The unitig has predecessors" << endl;
// 				// 	exit(0);
// 				// }

// 				//Try to continue the extension on the predecessive unitigs
// 				extPtr = extendAtPrevUnitigOnRevComp(prevUni->getPredecessors(), qPos - tmpExtLen, hitLen, tmpExtLen, q, X, tmpScore, explCount, quorum);
// 			}

// 			break;
// 		}
// 	}

// 	//Check whether we have an extension seed
// 	if(extSeed != NULL){
// 		//Check whether there exists a successful extension on a predecessive unitig
// 		if(extPtr == NULL){
// 			//Create a new extension pointer
// 			extPtr = (struct Ext_ptr*) malloc(sizeof(struct Ext_ptr));
// 			//Link this unitig to the extension pointer
// 			extPtr->origUnitig = *prevUni;
// 			//Link the extension seed
// 			extPtr->extSeed = extSeed;
// 		} else{
// 			//Add the score that we gained for this unitig to the score we gained extending its predecessors
// 			extPtr->extSeed->score += extSeed->score;
// 			//Delete extension seed for this unitig
// 			free(extSeed);
// 		}
// 	}
// 	//Testing
// 	//cout << "8 Option 2" << endl;

// 	//Check whether we have an extension pointer
// 	if(extPtr != NULL){ return extPtr->extSeed->score; }

// 	return 0;
// }

//This function continues a left extension on a predecessive unitig of a seed lying on the query's reference strand considering a quorum and a search color set and returns the achieved score
int32_t contLeftX_Drop(const neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> &prevUni, uint32_t qPos, uint32_t &hitLen, const string &q, const int16_t &X, const int32_t &lastSeedTmpScore, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	//bool inSeedList;
	int32_t tmpScore, score = 0;
	uint32_t tmpExtLen = 1;
	//Calculate our offset position inside the unitig sequence (+1, because otherwise it is not possible to calculate the correct gain if reaching a seed)
	uint32_t uPos= prevUni->size - prevUni->getGraph()->getK() + 1;
	int32_t checkedPos = getSrchCritCov(*prevUni, quorum, searchSet, compOffset(uPos - tmpExtLen, 1, prevUni->size, prevUni->strand), !prevUni->strand);// - (prevUni->size - (uPos + 1));

	//Testing
	// cout << "checkedPos: " << checkedPos << endl;
	// if(prevUni->mappedSequenceToString() == "CTCACCTTGGTTGAAAAATGCTAAAAGTGCCCCCTCAATACATTTATCAAGATATTGGTCT"){
	// 	cout << "contLeftX_DropOnRevComp: We are on unitig CTCACCTTGGTTGAAAAATGCTAAAAGTGCCCCCTCAATACATTTATCAAGATATTGGTCT" << endl;
	// 	cout << "checkedPos: " << checkedPos << endl;
	// 	return 0;
	// }

	// uint32_t k = prevUni->getGraph()->getK();
	// size_t offset;

	string prevUniSeq = prevUni->mappedSequenceToString();
	struct seed *nearestSeed, *prevSeed = NULL;//, *extSeed = NULL;

	// //Check whether we are working on the reference strand of the current unitig
	// if(prevUni->strand){

	//Find the nearest seed that we might be able to reach during our extension
	nearestSeed = searchLeftNeighbor(prevUni->getData()->getData(*prevUni)->getSeed(prevUni->strand), qPos, uPos, prevSeed);

	// } else{
	// 	//So far we have not considered the reverse complementary strand for the seed detection //TODO: Change this!
	// 	nearestSeed = NULL;
	// }

	//Initialize our temporary score with the last temporary score of its successive unitig
	tmpScore = lastSeedTmpScore;

	// //Initialize the temporary extension length
	// tmpExtLen = 1;

	//We are done if we reach the beginning of the query
	while(qPos >= tmpExtLen){
		//Check whether we have reached a nearest seed
		if(nearestSeed != NULL && qPos - tmpExtLen <= nearestSeed->offsetQ + nearestSeed->len - 1){
			//Update the temporary score
			tmpScore += qPos - tmpExtLen - nearestSeed->offsetQ + 1;//TODO Change this as soon as we want to allow more than only unit score!
			//Calculate the gain we have by incorporating the seed into our extension and update the temporary extension length
			tmpExtLen = qPos - nearestSeed->offsetQ + 1;
			//Adjust number of remaining covered positions
			checkedPos -= qPos - nearestSeed->offsetQ + 1;

			//Check whether we get a score larger 0 by incorporating the reached seed
			if(tmpScore > 0){
				//Update hit's length
				hitLen += qPos - nearestSeed->offsetQ;
				//Update score
				score += tmpScore;

				// //Check whether we already have an extension seed
				// if(extSeed != NULL){
				// 	//Update extension seed
				// 	extSeed->offsetU = nearestSeed->offsetU;
				// 	extSeed->offsetQ = nearestSeed->offsetQ;
				// 	extSeed->score += tmpScore;
				// 	extSeed->len += qPos - nearestSeed->offsetQ;

				//Update current position in q the unitig
				qPos = nearestSeed->offsetQ;
				uPos = nearestSeed->offsetU;
					
				// 	//Check if the reached seed has a predecessor in its seed list
				// 	if(prevSeed != NULL){
				// 		//Link the reached seed's predecessor and successor
				// 		prevSeed->nextSeed = nearestSeed->nextSeed;
				// 	} else{
				// 		//Set the reached seed's successor as the head of the seed list
				// 		prevUni->getData()->getData(*prevUni)->setSeed(nearestSeed->nextSeed, prevUni->strand);
				// 	}
				// 	//Delete the reached seed
				// 	free(nearestSeed);

				// 	// //Delete the reached seed
				// 	// extractUnprocSeed(prevUni->getData()->getData(*prevUni), nearestSeed, prevSeed);
				// 	// //Reset prevSeed
				// 	// prevSeed = NULL;
				// } else{
				// 	//Make the reached seed our new extension seed
				// 	extSeed = nearestSeed;
				// 	//Update seed's length and score (we have to substract 1 here because we initialize the temporary extension length with 1)
				// 	extSeed->len = tmpExtLen - 1;
				// 	extSeed->score = tmpScore;

				// 	// //Mark that we recycle a seed instead of creating a new one for our extension
				// 	// inSeedList = true;

				// 	//TODO Check whether we can export this into a function!
				// 	//Check whether the reached seed is the first in its list
				// 	if(prevSeed != NULL){
				// 		//Exclude the reached seed from its list
				// 		prevSeed->nextSeed = extSeed->nextSeed;
				// 	} else{
				// 		//Exclude the reached seed from its list
				// 		prevUni->getData()->getData(*prevUni)->setSeed(nearestSeed->nextSeed, prevUni->strand);
				// 	}

				// 	// //Incorporate seed into the list of processed seeds
				// 	// extSeed->nextSeed = prevUni->getData()->getData(*prevUni)->getLastProcSeed();
				// 	// prevUni->getData()->getData(*prevUni)->setLastProcSeed(extSeed);

				// 	//Update current position in q the unitig
				// 	qPos = extSeed->offsetQ;
				// 	uPos = extSeed->offsetU;
				// }

				//Reset temporary score and extension length
				tmpScore = 0;
				tmpExtLen = 1;
			}// else{

			//Check if the reached seed has a predecessor in its seed list
			if(prevSeed != NULL){
				//Link the reached seed's predecessor and successor
				prevSeed->nextSeed = nearestSeed->nextSeed;
			} else{
				//Set the reached seed's successor as the head of the seed list
				prevUni->getData()->getData(*prevUni)->setSeed(nearestSeed->nextSeed, prevUni->strand);
			}
			//Delete the reached seed
			free(nearestSeed);

				// //Delete the reached seed
				// extractUnprocSeed(prevUni->getData()->getData(*prevUni), nearestSeed, prevSeed);
				// //Reset prevSeed
				// prevSeed = NULL;
			// }

			//Reset prevSeed
			prevSeed = NULL;
			//Search for the next seed to reach
			nearestSeed = searchLeftNeighbor(prevUni->getData()->getData(*prevUni)->getSeed(prevUni->strand), qPos, uPos, prevSeed);
		} else if(uPos >= tmpExtLen){//Check whether we have already reached the beginning of the unitig sequence
			//Check whether quorum has to be checked
			if(checkedPos <= 0){
				break;

				// //Check on which strand we currently are
				// if(prevUni->strand){
				// 	//Calculate offset and make sure there is no overflow
				// 	if(uPos - tmpExtLen < k - 1){
				// 		//Testing
				// 		// cout << "11 Option 1" << endl;

				// 		offset = 0;
				// 	} else{
				// 		//Testing
				// 		// cout << "11 Option 2" << endl;

				// 		offset = uPos - tmpExtLen - (k - 1);
				// 	}
				// } else{
				// 	//Calculate offset and make sure there is no overflow
				// 	if(uPos - tmpExtLen < k - 1){
				// 		//Testing
				// 		// cout << "11 Option 3" << endl;

				// 		offset = prevUni->size - (k - 1);
				// 	} else{
				// 		//Testing
				// 		// cout << "11 Option 4" << endl;

				// 		offset = prevUni->size - (uPos - tmpExtLen) - 1;
				// 	}
				// }
				
				// // offset = (tmpExtLen + k - 1 > uPos ? 0 : uPos - (tmpExtLen + k - 1));

				// //Check whether quorum is still fulfilled
				// if((checkedPos = checkSearchCrit(*prevUni, quorum, offset, searchSet)) == 0) break;//TODO If we want to allow to check more than k positions at once we need a new function here!
			}

			//Testing
			//cout << "uPos - tmpExtLen: " << uPos - tmpExtLen << " qPos - tmpExtLen: " << qPos - tmpExtLen << endl;

			//Compare the next two bases and check whether our temporary score is becoming > 0 by this
			if((tmpScore += compUScore(prevUniSeq[uPos - tmpExtLen], q[qPos - tmpExtLen])) > 0){
				// //Check whether there is no extension seed yet
				// if(extSeed == NULL){
				// 	//Create a new extension seed
				// 	extSeed = (struct seed*) malloc(sizeof(struct seed));
				// 	extSeed->len = 0;
				// 	extSeed->score = 0;

				// 	// //Mark that this extension seed is not part of a seed list yet
				// 	// inSeedList = false;
				// }

				//Update positions in q and the unitig
				uPos -= tmpExtLen;
				qPos -= tmpExtLen;

				// //Update extension seed's infos
				// extSeed->offsetU = uPos;
				// extSeed->offsetQ = qPos;
				// extSeed->len += tmpExtLen;
				// extSeed->score += tmpScore;
				
				//Update score
				score += tmpScore;
				//Update hit's length
				hitLen += tmpExtLen;
				//Reset temporary length and score
				tmpExtLen = 0;
				tmpScore = 0;
			} else if(tmpScore < -X){//Check whether our temporary score is already too negative
				break;
			}

			//Increment temporary extension length
			++tmpExtLen;
			--checkedPos;
		} else{
			//Decrease temporary extension length
			--tmpExtLen;

			//Check if this unitig has a predecessor
			if(prevUni->getPredecessors().hasPredecessors()){
				//Try to continue the extension on the predecessive unitigs
				score += extendAtPrevUnitig(prevUni->getPredecessors(), qPos - tmpExtLen, hitLen, tmpExtLen, q, X, tmpScore, extPth, explCount, quorum, searchSet);
			}

			break;
		}
	}

	// //Check whether we have an extension seed
	// if(extSeed != NULL){
	// 	//Check whether there exists a successful extension on a predecessive unitig
	// 	if(extPtr == NULL){
	// 		//Create a new extension pointer
	// 		extPtr = (struct Ext_ptr*) malloc(sizeof(struct Ext_ptr));
	// 		//Link this unitig to the extension pointer
	// 		extPtr->origUnitig = *prevUni;
	// 		//Link the extension seed
	// 		extPtr->extSeed = extSeed;

	// 		// //Save whether the linked seed still needs to be incorporated into a seed list
	// 		// extPtr->inSeedList = inSeedList;
	// 	} else{
	// 		//Add the score that we gained for this unitig to the score we gained extending its predecessors
	// 		extPtr->extSeed->score += extSeed->score;
	// 		//Delete extension seed for this unitig
	// 		free(extSeed);

	// 		//delExtPtrSeed(inSeedList, extSeed, prevUni->getData()->getData(*prevUni));
	// 	}
	// }

	// //Check whether we have an extension pointer
	// if(extPtr != NULL){ return extPtr->extSeed->score; }

	return score;
}

//This function continues a left extension on a predecessive unitig of a seed lying on the query's reverse complement considering a quorum and a search color set and returns the achieved score
int32_t contLeftX_DropOnRevComp(const neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> &prevUni, uint32_t qPos, uint32_t &hitLen, const string &q, const int16_t &X, const int32_t &lastSeedTmpScore, list<uint16_t> &extPth, uint32_t uPos, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	int32_t tmpScore, score = 0;
	uint32_t tmpExtLen = 1, overlap = 0;

	int32_t checkedPos = getSrchCritCov(*prevUni, quorum, searchSet, compOffset(uPos - tmpExtLen, 1, prevUni->size, prevUni->strand), !prevUni->strand);// - (prevUni->size - uPos);
	uint32_t k = prevUni->getGraph()->getK();
	uint32_t progress;

	// size_t offset;

	string prevUniSeq = prevUni->mappedSequenceToString();
	struct seed *nearestSeed, *prevSeed = NULL;//, *extSeed = NULL;

	// //If we are on the reference strand we have to consider the overlap in the unitig sequence's end for checkedPos
	// if(prevUni->strand) checkedPos -= (k - 1);
	
	//Find the nearest seed that we might be able to reach during our extension
	nearestSeed = searchLeftNeighbor(prevUni->getData()->getData(*prevUni)->getSeed(prevUni->strand), qPos, uPos, prevSeed);

	//Testing
	// if(nearestSeed == NULL) cout << "2 Option 2" << endl;
	// cout << "uPos: " << uPos << endl;
	// if(prevUni->mappedSequenceToString() == "CTCACCTTGGTTGAAAAATGCTAAAAGTGCCCCCTCAATACATTTATCAAGATATTGGTCT"){
	// 	cout << "contLeftX_DropOnRevComp: We are on unitig CTCACCTTGGTTGAAAAATGCTAAAAGTGCCCCCTCAATACATTTATCAAGATATTGGTCT" << endl;
	// 	cout << "checkedPos: " << checkedPos << endl;
	// 	cout << "nearestSeed is " << (nearestSeed == NULL ? "" : "not ") << "NULL" << endl;
	// 	// exit(0);
	// }

	//Check how far we should explore the unitig's sequence
	if(prevUni->getPredecessors().hasPredecessors()){
		//Testing
		// cout << "1 Option 1" << endl;

		overlap = k - 1;
	}// else{
	// 	//Testing
	// 	// cout << "1 Option 2" << endl;
	// }

	//Initialize our temporary score with the last temporary score of its successive unitig
	tmpScore = lastSeedTmpScore;
	
	// //Initialize the temporary extension length
	// tmpExtLen = 1;

	//We are done if we reach the beginning of the query
	while(qPos >= tmpExtLen){
		//Check whether we have reached a nearest seed
		if(nearestSeed != NULL && qPos - tmpExtLen <= nearestSeed->offsetQ + nearestSeed->len - 1){
			//Calculate the progress with have by incorporating the seed
			progress = qPos - tmpExtLen - nearestSeed->offsetQ + 1;

			//Testing
			// cout << "tmpScore: " << tmpScore << endl;
			// cout << "progress: " << progress << endl;
			// if(prevUni->mappedSequenceToString() == "CTCACCTTGGTTGAAAAATGCTAAAAGTGCCCCCTCAATACATTTATCAAGATATTGGTCT") cout << "We should not go in here" << endl;

			//Update the temporary score
			tmpScore += progress;//TODO Change this as soon as we want to allow more than only unit score!
			//Calculate the gain we have by incorporating the seed into our extension and update the temporary extension length
			tmpExtLen = qPos - nearestSeed->offsetQ + 1;

			//Check whether we get a score larger 0 by incorporating the reached seed
			if(tmpScore > 0){
				//Testing
				//cout << "4 Option 1" << endl;
				// cout << "3 Option 1" << endl;

				//Update hit's length
				hitLen += qPos - nearestSeed->offsetQ;
				//Update score
				score += tmpScore;
				//Reset tmpScore
				tmpScore = 0;
				//Reset temporary extension length
				tmpExtLen = 1;

				//Testing
				// cout << "Score after seed incorporation: " << score << endl;

				//Update current position in q the unitig
				qPos = nearestSeed->offsetQ;
				uPos = nearestSeed->offsetU;
			}

			//Check if the reached seed has a predecessor in its seed list
			if(prevSeed != NULL){
					//Testing
					// cout << "7 Option 2" << endl;

				//Link the reached seed's predecessor and successor
				prevSeed->nextSeed = nearestSeed->nextSeed;
			} else{
					//Testing
					// cout << "7 Option 1" << endl;

				//Set the reached seed's successor as the head of the seed list
				prevUni->getData()->getData(*prevUni)->setSeed(nearestSeed->nextSeed, prevUni->strand);
			}
			//Delete the reached seed
			free(nearestSeed);

			// }

			//Testing
			// if(progress > checkedPos) cout << "18 Option 1" << endl;

			//Decrease number of remaining covered positions
			checkedPos -= progress; 
			//Reset prevSeed
			prevSeed = NULL;
			//Search for the next seed to reach
			nearestSeed = searchLeftNeighbor(prevUni->getData()->getData(*prevUni)->getSeed(prevUni->strand), qPos, uPos, prevSeed);
		} else if(uPos >= tmpExtLen + overlap){//Check whether we have already reached the beginning of the unitig sequence
			//Testing
			// cout << "8 Option 2" << endl;
			// if(prevUni->mappedSequenceToString() == "CTCACCTTGGTTGAAAAATGCTAAAAGTGCCCCCTCAATACATTTATCAAGATATTGGTCT") cout << "Let's compare the next bases" << endl;

			//Check whether quorum has to be checked
			if(checkedPos <= 0){
				//Testing
				// cout << "9 Option 1" << endl;
				// if(prevUni->mappedSequenceToString() == "CTCACCTTGGTTGAAAAATGCTAAAAGTGCCCCCTCAATACATTTATCAAGATATTGGTCT") cout << "checkedPos is there so we cannot do anything" << endl;

				// //Calculate offset and make sure there is no overflow
				// offset = (tmpExtLen + k - 1 > uPos ? 0 : uPos - (tmpExtLen + k - 1));//TODO Theoretically, this is not necessary!

				/*Calculate the correct offset for the check*/

				// //Make sure we do not calculate a negative offset
				// if((uPos - tmpExtLen) < (k - 1)){
				// 	//Testing
				// 	// cout << "16" << endl;

				// 	//Check the remaining sequence part
				// 	offset = 0;
				// } else{
				// 	if(prevUni->strand){
				// 		//Testing
				// 		// cout << "17 Option 1" << endl;

				// 		//On the reference strand we may just move k-1 positions to the left
				// 		offset = uPos - tmpExtLen - (k - 1);
				// 	} else{
				// 		//Testing
				// 		// cout << "17 Option 2" << endl;

				// 		//The offset on the reference strand corresponding to our position suffices already to check the next k positions
				// 		offset = prevUni->size - (uPos - tmpExtLen) - 1;
				// 	}
				// }

				// if(prevUni->strand){
				// 	//On the reference strand we may just move k-1 positions to the left
				// 	offset = uPos - tmpExtLen - (k - 1);
				// } else{
				// 	//The offset on the reference strand corresponding to our position suffices already to check the next k positions
				// 	offset = prevUni->size - (uPos - tmpExtLen) - 1;
				// }

				//Testing
				// cout << "We have checked the quorum for position " << offset << endl;
				// cout << "prevUni->size: " << prevUni->size << endl;
				// cout << "uPos: " << uPos << endl;
				// cout << "tmpExtLen: " << tmpExtLen << endl;

				//Check whether quorum is still fulfilled
				// if((checkedPos = checkSearchCrit(*prevUni, quorum, offset, searchSet)) == 0){
					//Testing
					//cout << "2" << endl;
					// cout << "10 Option 2" << endl;

				break;

				// }

				//Testing
				// cout << "10 Option 1" << endl;
			}// else{
			// 	//Testing
			// 	// cout << "9 Option 2" << endl;
			// }

			//Testing
			// cout << "Bases to compare:" << prevUniSeq[uPos - tmpExtLen] << " " << q[qPos - tmpExtLen] << endl;
			// cout << "uPos: " << uPos << endl;
			// cout << "tmpExtLen: " << tmpExtLen << endl;
			// cout << "checkedPos: " << checkedPos << endl;
			// cout << "Position u: " << uPos - tmpExtLen << " position q: " << qPos - tmpExtLen << endl;
			// if(report){
			// 	cout << "Bases to be compared:" << endl;
			// 	cout << "u " << uPos - tmpExtLen << " : " << prevUniSeq[uPos - tmpExtLen] << " " << q[qPos - tmpExtLen] << " : " << qPos - tmpExtLen << " q" << endl;
			// }

			//Compare the next two bases and check whether our temporary score is becoming > 0 by this
			if((tmpScore += compUScore(prevUniSeq[uPos - tmpExtLen], q[qPos - tmpExtLen])) > 0){
				//Testing
				// cout << "11 Option 1" << endl;

				// //Check whether there is no extension seed yet
				// if(extSeed == NULL){
				// 	//Create a new extension seed
				// 	extSeed = (struct seed*) malloc(sizeof(struct seed));
				// 	extSeed->len = 0;
				// 	extSeed->score = 0;
				// 	extSeed->nextSeed = NULL;
				// }

				//Update positions in q and the unitig
				uPos -= tmpExtLen;
				qPos -= tmpExtLen;

				// //Update extension seed's infos
				// extSeed->offsetU = uPos;
				// extSeed->offsetQ = qPos;
				// extSeed->len += tmpExtLen;
				// extSeed->score += tmpScore;

				//Update score
				score += tmpScore;
				//Update hit's length
				hitLen += tmpExtLen;
				//Reset temporary length and score
				tmpExtLen = 0;
				tmpScore = 0;
			} else if(tmpScore < -X){//Check whether our temporary score is already too negative
				//Testing
				// cout << "11 Option 3" << endl;

				break;
			}// else{
			// 	//Testing
			// 	// cout << "11 Option 2" << endl;
			// }
			//Testing
			//if(tmpScore < 0) cout << "8 Option 2" << endl;

			//Increment temporary extension length
			++tmpExtLen;
			--checkedPos;
		} else{
			//Testing
			// cout << "8 Option 1" << endl;

			//Decrease temporary extension length
			--tmpExtLen;

			//Check if this unitig has a predecessor
			if(overlap != 0){
				//Testing
				// cout << "12 Option 1" << endl;

				//Try to continue the extension on the predecessive unitigs
				score += extendAtPrevUnitigOnRevComp(prevUni->getPredecessors(), qPos - tmpExtLen, hitLen, tmpExtLen, q, X, tmpScore, extPth, overlap - (uPos - tmpExtLen), explCount, quorum, searchSet);
			}// else{
			// 	//Testing
			// 	// cout << "12 Option 2" << endl;
			// }

			break;
		}
	}

	//Testing
	// if(prevUni->mappedSequenceToString() == "CTCACCTTGGTTGAAAAATGCTAAAAGTGCCCCCTCAATACATTTATCAAGATATTGGTCT") exit(0);

	return score;

}

//This function checks if a hit's start offset for the gapped extension lies inside the overlap of a unitig sequence's end. If so it tries to move the start position on a unitig where it is not inside the overlap anymore. If this is not possible the hit is marked as invalid by setting its score to 0.
void mvStartToValUni(hit* h, list<uint16_t>& lExtPth){
	//Testing
	// cout << "mvStartToValUni: Start of function" << endl;

	//If we want to switch unitigs the right extension path must not be empty
	if(h->rExt.nbElem > 0){
		//Testing
		// cout << "1 Option 2" << endl;

		//Decompress right extension path
		list<uint16_t> rExtPth = decmprExtPth(h->rExt);

		//Switch unitigs
		switUni(h->offU, h->origUni, lExtPth, rExtPth);

		//Check if unitig switching was successful
		if(h->offU <= h->origUni.size - h->origUni.getGraph()->getK()){
			//Testing
			// cout << "2 Option 1" << endl;

			//Compress right extension path again
			h->rExt = cmprExtPth(rExtPth);
			return;
		}// else{
		// 	//Testing
		// 	cout << "2 Option 2" << endl;
		// }
	}// else{
	// 	//Testing
	// cout << "1 Option 1" << endl;
	// }

	//Mark hit as invalid
	h->score = 0;
}

//This function moves an offset from a given unitig to its successor if there is any while keeping left and right extension paths updated. If the offset at the successive unitig lies inside the overlap at the unitig sequence's end the function calls itself recursively rExtPth must not be empty
void switUni(uint32_t &offset, UnitigColorMap<seedlist> &currUni, list<uint16_t> &lExtPth, list<uint16_t> &rExtPth){
	uint16_t sucCount = 1, predCount = 1;

	//Testing
	// cout << "switUni: Start of function" << endl; 

	//Traverse successors of the current unitig
	for(neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> suc = currUni.getSuccessors().begin(); suc != currUni.getSuccessors().end(); ++suc, ++sucCount){
		//Check if we have found the required successor
		if(sucCount < rExtPth.front()){
			//Testing
			// cout << "1 Option 2" << endl;

			continue;
		}

		//Testing
		// cout << "1 Option 1" << endl;

		//Iterate over predecessors of found successor
		for(neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> pred = suc->getPredecessors().begin(); pred != suc->getPredecessors().end(); ++pred, ++predCount){
			//Check if we have found the unitig we came from
			if(*pred == currUni){
				//Testing
				// cout << "2 Option 1" << endl;

				break;
			}// else{
			// 	//Testing
			// 	// cout << "2 Option 2" << endl;
			// }
		}

		//Testing
		// cout << "switUni: Predecessor found" << endl;
		// cout << "switUni: rExtPth is " << (rExtPth.empty() ? "" : "not " ) << "empty" << endl;

		//Save count in left extension path
		lExtPth.push_front(predCount);
		//Delete successor from right extension path
		rExtPth.pop_front();

		//Testing
		// cout << "switUni: Extension paths modified" << endl;
		if(report) cerr << "switUni: Switch from unitig " << currUni.mappedSequenceToString() << " with len " << currUni.len << " offset was " << offset << endl;

		//Update offset
		offset -= currUni.len;
		//Update current unitig
		currUni = *suc;

		//Testing
		if(report) cerr << "switUni: Switch to unitig " << currUni.mappedSequenceToString() << " new offset: " << offset << endl;

		//End iteration
		break;

		//Testing
		// cout << "switUni: end of for loop" << endl;
	}

	//Check if it is necessary and possible to go on
	if(rExtPth.empty() || offset <= currUni.size - currUni.getGraph()->getK()){
		//Testing
		// if(rExtPth.empty()) cout << "3 Option 2" << endl;
		// if(offset <= currUni.size - currUni.getGraph()->getK()) cout << "3 Option 3" << endl;
		if(report){
			cerr << "switUni: Recursion terminated" << endl;
			report = false;
		}

		return;
	}

	//Testing
	// cout << "3 Option 1" << endl;

	//Move to next unitig
	switUni(offset, currUni, lExtPth, rExtPth);
}