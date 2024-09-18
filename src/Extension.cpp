#include <stack>
#include <deque>

#include "Extension.h"
#include "Hit.h"
#include "Search.h"
#include "UnitigInfo.cpp"

//Testing
// extern bool report = false;

//This function initiates the extension on all successors of a unitig and returns the best one considering a quorum and a search color set
int32_t extendAtNextUnitig(const ForwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &uniPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	uint16_t sucID;
	uint32_t tmpHitLen = 0;
	int32_t maxScore = 0, currScore;
	
	//Check whether we have reached the maximum recursion depth of an extension
	if(++explCount > EXPLORED_UNITIGS_MAX){
		//Report this incident
		//cerr << "Maximum recursion depth reached during extension. Position in q: " << iniQoff << endl;
		//Terminate this extension
		return 0;
	}

	maxScore = 0;
	sucID = 0;

	//Testing
	// if(report){
	// 	cout << "extendAtNextUnitig: Iterate over successors" << endl;
	// 	cout << "extendAtNextUnitig: sucIter.begin()->mappedSequenceToString(): " << sucIter.begin()->mappedSequenceToString() << endl;
	// }

	//Iterate over successors
	for(neighborIterator<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> nI = sucIter.begin(); nI != sucIter.end(); ++nI){
		//Temporary extention path
		list<uint16_t> tmpPth;
		//Note which successor we are on
		++sucID;

		//Testing
		// if(report) cout << "extendAtNextUnitig: Continue extension on unitig " << nI->mappedSequenceToString() << endl;

		//Calculate the score of an extension of a successor
		currScore = contRightX_Drop(nI, iniQoff, tmpHitLen, extLen, q, mscore, mmscore, X, lastExtSeedTmpScore, uniPos, tmpPth, explCount, quorum, searchSet, advIdx);

		//Testing
		// if(report) cout << "extendAtNextUnitig: Return from extension on unitig " << nI->mappedSequenceToString() << endl;

		//Check whether the score of the current successors extension is the best one found so far
		if(currScore > maxScore){
			//Update maxScore
			maxScore = currScore;
			//Update hit length (tmpHitLen does not have to be reseted, because it is only set to and never read from except here)
			hitLen = tmpHitLen;
			extPth = tmpPth;
			//Save which successor we have chosen
			extPth.push_front(sucID);
		}
	}

	//Nothing found
	return maxScore;
}

//This function  initiates the extension on all successors of a unitig and returns the best one considering a quorum and a search color set. This function is explicitly designed for seeds lying on the query's reverse complement (considering the overlap between unitigs in sequences' beginning)
int32_t extendAtNextUnitig_OnRevComp(const ForwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	uint16_t sucID;
	uint32_t tmpHitLen = 0;
	int32_t maxScore = 0, currScore;
	
	//Check whether we have reached the maximum recursion depth of an extension
	if(++explCount > EXPLORED_UNITIGS_MAX) return 0;//Terminate this extension

	maxScore = 0;
	sucID = 0;

	//Iterate over successors
	for(neighborIterator<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> nI = sucIter.begin(); nI != sucIter.end(); ++nI){
		//Temporary extention path
		list<uint16_t> tmpPth;
		//Note which successor we are on
		++sucID;
		//Calculate the score of an extension of a successor
		currScore = contRightX_Drop_OnRevComp(nI, iniQoff, tmpHitLen, extLen, q, mscore, mmscore, X, lastExtSeedTmpScore, tmpPth, explCount, quorum, searchSet, advIdx);

		//Check whether the score of the current successors extension is the best one found so far
		if(currScore > maxScore){
			//Update maxScore
			maxScore = currScore;
			//Update hit length
			hitLen = tmpHitLen;
			extPth = tmpPth;
			//Save which successor we have chosen
			extPth.push_front(sucID);
		}
	}

	//Nothing found
	return maxScore;
}

//The good old X-drop algorithm (extension to the right) for seeds matching the query's reference strand considering quorum and 
//search color set using an iterative approach
void perfRightX_Drop(Hit* hit, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const uint32_t 
	&quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	//The rank of a successor of a leading unitig (used to construct the extension path)
	uint16_t nr;
	//The last position in a unitig's sequence that needs to be compared with
	uint32_t lastSeqPos;
	//The number of times we have already switched unitigs during this extension
	uint32_t uniSwtchCnt = 0;
	//The k value are using
	int K = hit->origUni.getGraph()->getK();
	//The length of the currently explored extension on the leading unitig
	int32_t extLen;
	//A number of positions on the leading unitig fulfilling the search criteria
	int32_t nbffPos;
	//The progress we achieve due to reaching a seed
	int32_t progress;
	//The initial query position at the beginning of an exploration on the current leading unitig (might be different from 
	//currExt.offsetQ right after switching unitigs)
	int32_t posQ;
	//The leading unitig's sequence
	string uSeq;
	struct Seed *nearestSeed, *prevSeed;
	//An extension we are dealing with
	Ext currExt;
	//A vector of extensions to be continued on a successor of a leading unitig//TODO: This is only an artificial detour to allow 
	//the reproduction of results of the old recursive PLAST version for testing purposes. As soon as we are sure that the iterative
	// version of PLAST works this should be removed!
	stack<Ext> newExts;
	//Initialize DFS stack
	stack<Ext> exts;//TODO: We use a deque here as the container for the stack, because intuitively it sounds more efficient (e.g. compared to a vector). At some point we should try to verify this by a test!

	//Create the first extension right after the hit's end
	exts.push(Ext(hit->offU + hit->length, hit->offQ + hit->length, hit->offQ + hit->length, 0, 0, hit->origUni, list<uint16_t>()));

	while(!exts.empty()){
		//Get first extension on stack
		currExt = exts.top();
		//Remove first extension from stack
		exts.pop();

		//Try to extend current extension on the leading unitig as far as possible//

		//If the current unitig has successors we do not need to compare the k-1 overlap//TODO: Would not it actually be better to compare it on this unitig rather than first move to the next one?
		lastSeqPos = currExt.ldUni.size - (currExt.ldUni.getSuccessors().hasSuccessors() ? K : 1);
		//Set extension length on leading unitig
		extLen = 0;
		//Get the leading unitig's sequence
		uSeq = currExt.ldUni.mappedSequenceToString();
		//Compute the number of positions on which search criteria are fulfilled	
		nbffPos = getSrchCritCov(currExt.ldUni, quorum, searchSet, compOffset(currExt.offsetU, 1, currExt.ldUni.size, 
			currExt.ldUni.strand), currExt.ldUni.strand, advIdx);
		//Search for a seed we could reach during this extension (tmpQoff has to be used here to match with offsetU, because offsetQ
		// is only updated after an explored continuation of an extension is also accepted due to a positive score which might not 
		//be the case when switching to the next unitig whereas offsetU is always updated when switching)
		nearestSeed = searchRightNeighbor(currExt.ldUni.getData()->getData(currExt.ldUni)->getSeed(currExt.ldUni.strand), 
			(uint32_t) currExt.tmpQoff, (uint32_t) extLen, (uint32_t) currExt.offsetU, prevSeed);//TODO: As soon as the iterative extension is completely implemented, it might be possible to modify input parameter types of searchRightNeighbor() to avoid the necessity of these casts
		//The initial query position is always identical to currExt.tmpQoff in the beginning, but this might change during the ex-
		//tension on the leading unitig
		posQ = currExt.tmpQoff;

		//We are done if we have reached the end of the query or search criteria are no longer fulfilled
		while(nbffPos-- > 0 && posQ + extLen < q.length()){
			//Check if we have reached the next seed
			if(nearestSeed != NULL && posQ + extLen >= nearestSeed->offsetQ){
				//Calculate the number of positions on the unitig we do not need to compare
				progress = nearestSeed->offsetQ + nearestSeed->len - (posQ + extLen);
				//Update extension length
				extLen += progress;
				//Update nbffPos
				nbffPos -= progress;
				//Update extension with information about the reached seed (the score has to be increased by 1, because progress 
				//does not consider that we have just moved forward in the beginning of this iteration)
				updateExtension(currExt, posQ, (progress + 1) * mscore, extLen, true);//TODO: This function still needs to be tested!

				//Check if there is another seed to reach
				if(prevSeed != NULL){
					//Exclude the reached seed from its seed list
					prevSeed->nextSeed = nearestSeed->nextSeed;
					//Delete the reached seed
					free(nearestSeed);
					//Search for another neighbor
					nearestSeed = searchRightNeighbor(currExt.ldUni.getData()->getData(currExt.ldUni)->getSeed(currExt.ldUni.strand)
						, (uint32_t) posQ, (uint32_t) extLen, (uint32_t) currExt.offsetU, prevSeed)
					;
				} else{
					//Make seed's successor the head of the seed list
					currExt.ldUni.getData()->getData(currExt.ldUni)->setSeed(nearestSeed->nextSeed, currExt.ldUni.strand);
					//Delete the reached seed
					free(nearestSeed);
					//Reset nearestSeed
					nearestSeed = NULL;
				}

				//Possibly this has to be done after the nearest neighbor search
				++extLen;
			//Check if there is still sequence left to compare on the leading unitig
			} else if(currExt.offsetU + extLen <= lastSeqPos){
				//Update extension with information about the next two compared positions
				updateExtension(currExt, posQ, compUScore(uSeq[currExt.offsetU + extLen], q[posQ + extLen], mscore, mmscore), 
					extLen, true);

				//Terminate exploration if the score has become too bad
				if(currExt.tmpScore < -X) break;

				//Move on
				++extLen;
			} else{
				//Check if leading unitig has successors (the way it is done is probably cheaper than to call hasSuccessors()) and 
				//if we have not yet reached the maximum number of unitig switches during this extension
				if(lastSeqPos ==  currExt.ldUni.size - K && uniSwtchCnt++ < EXPLORED_UNITIGS_MAX){
					//Reset successor's rank
					nr = 0;

					//Iterate over successors
					for(auto n = currExt.ldUni.getSuccessors().begin(); n != currExt.ldUni.getSuccessors().end(); ++n, ++nr){
						//Copy current extension
						Ext newExt = Ext(currExt);
						//Calculate sequence position at which we start on the next unitig
						newExt.offsetU = currExt.offsetU + extLen - lastSeqPos;
						//Set tmpQoff correctly
						newExt.tmpQoff = posQ + extLen;
						//Update leading unitig
						newExt.ldUni = *n;
						//Extend extension path by new leading unitig
						newExt.pth.push_back(nr);
						//Add the new extension to the auxiliar vector to preserve the correct order
						newExts.push(newExt);
					}

					//Push the new extensions to the stack in the correct order
					while(!newExts.empty()){
						exts.push(newExts.top());
						newExts.pop();
					}
				}

				//Terminate current extension
				break;
			}
		}

		//Update info about the best extension found so far//

		//Check if the score of the current extension is better than what we found so far
		if(currExt.score > hit->score){
			//Update max score found
			hit->score = currExt.score;
			//Recalculate hit's length
			hit->length = (uint32_t) currExt.offsetQ - hit->offQ + 1;//TODO: We should think about whether it would not also make sense to change the type of offets inside the Hit Class
			//Clear extension path if it already exists
			decmprExtPth(hit->rExt);
			//Copy extension path
			hit->rExt = cmprExtPth(currExt.pth);
		}
	}
}

//The good old X-drop algorithm (extension to the right) for seeds matching the query's reference strand considering quorum and 
//search color set using recursive function calls
void startRightX_Drop(Hit* hit, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	//Initialization of auxiliary variables
	int32_t tmpScore = 0;
	uint32_t tmpSeedLen = hit->length;
	uint32_t overlap = 0;
	uint32_t iniUniPos;
	int32_t checkedPos = getSrchCritCov(hit->origUni, quorum, searchSet, hit->length, true, advIdx);
	//Counter to count tries to explore a further unitig
	uint32_t explCount;
	string uSeq = hit->origUni.mappedSequenceToString();
	list<uint16_t> extPth;

	//Testing
	// if(report) cout << "startRightX_Drop: Start extension to the right" << endl;

	//Calculate hit's initial score
	hit->score = hit->length * mscore;

	//Check how far we should explore the unitig's sequence
	if(hit->origUni.getSuccessors().hasSuccessors()) overlap = hit->origUni.getGraph()->getK() - 1;

	//We are done if we have reached the end of the query
	while(hit->offQ + tmpSeedLen < q.length()){
		//Check whether we have reached the end of the unitig's sequence
		if(hit->offU + tmpSeedLen < hit->origUni.size - overlap){
			//Ensure that search criteria are still fulfilled
			if(checkedPos == 0) break;

			//Check whether the score of our extension is positive
			if((tmpScore += compUScore(uSeq[hit->offU + tmpSeedLen], q[hit->offQ + tmpSeedLen], mscore, mmscore)) > 0){
				//Update the seed info
				hit->score += tmpScore;
				hit->length = tmpSeedLen + 1;//+1 because we haven't increased tmpSLen yet
				//Reset the temporary score
				tmpScore = 0;
			} else{
				//Check if the current extension is already too bad
				if(tmpScore < -X) break;
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

				//Testing
				// if(report) cout << "startRightX_Drop: Continue extension on successors" << endl;

				//Explore unitig's successors
				hit->score += extendAtNextUnitig(hit->origUni.getSuccessors(), hit->offQ, hit->length, tmpSeedLen, q, mscore, mmscore, X, tmpScore, iniUniPos, extPth, explCount, quorum, searchSet, advIdx);
			}

			break;
		}
	}

	//Testing
	// if(report) cout << "startRightX_Drop: Terminating right extension" << endl;

	//Compress extension path
	hit->rExt = cmprExtPth(extPth);
}

//The good old X-drop algorithm (extension to the right) for seeds matching the query's reverse complement considering quorum and search color set. Returns an extension pointer storing the extension path through the graph
void startRightX_Drop_OnRevComp(Hit* hit, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	//Initialization of auxiliary variables
	int32_t tmpScore = 0;
	uint32_t tmpSeedLen = hit->length;
	//Counter to count tries to explore a further unitig
	uint32_t explCount;
	string uSeq = hit->origUni.mappedSequenceToString();
	list<uint16_t> extPth;

	//Get the number of checked positions considering that we do not start in the sequence's very end
	int32_t checkedPos = getSrchCritCov(hit->origUni, quorum, searchSet, compOffset(hit->offU + hit->length, 1, hit->origUni.size, false), false, advIdx);
	//Calculate hit's initial score
	hit->score = hit->length * mscore;

	//We are done if we have reached the end of the query
	while(hit->offQ + tmpSeedLen < q.length()){
		//Check whether we have reached the end of the unitig's sequence
		if(hit->offU + tmpSeedLen < hit->origUni.size){
			//Check whether the current position still fulfills the search criteria
			if(checkedPos == 0) break;

			//Check whether the score of our extension is positive
			if((tmpScore += compUScore(uSeq[hit->offU + tmpSeedLen], q[hit->offQ + tmpSeedLen], mscore, mmscore)) > 0){
				//Update the seed info
				hit->score += tmpScore;
				hit->length = tmpSeedLen + 1;//+1 because we haven't increased tmpSLen yet
				//Reset the temporary score
				tmpScore = 0;
			} else{
				//Check if the current extension is already too bad
				if(tmpScore < -X) break;
			}
			//Proceed with the next two positions
			++tmpSeedLen;
			--checkedPos;
		} else{
			//Check if the current unitig has successors
			if(hit->origUni.getSuccessors().hasSuccessors()){
				//Initialize explCount
				explCount = 0;
				//Explore unitig's successors
				hit->score += extendAtNextUnitig_OnRevComp(hit->origUni.getSuccessors(), hit->offQ, hit->length, tmpSeedLen, q, mscore, mmscore, X, tmpScore, extPth, explCount, quorum, searchSet, advIdx);
			}

			break;
		}
	}

	//Compress extension path
	hit->rExt = cmprExtPth(extPth);
}

//This function continues an extension to the right on a successive unitig of a seed lying on the query's reference strand considering a quorum and a search color set. Returns the maximum score reached.
int32_t contRightX_Drop(const neighborIterator<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> &sucUnitig, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t &extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastSeedTmpScore, uint32_t uniSeqPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	int32_t tmpScore, progress, score = 0;
	int32_t overlap = sucUnitig->getGraph()->getK() - 1;
	uint32_t iniSeqPos;
	uint32_t tmpSLen;
	//Even if we are on the reverse complementary strand and no position is covered on this unitig, because we have checked the first k - 1 position on the last unitig already
	int32_t checkedPos = getSrchCritCov(*sucUnitig, quorum, searchSet, compOffset(uniSeqPos, 1, sucUnitig->size, sucUnitig->strand), sucUnitig->strand, advIdx);
	string sucUniSeq = sucUnitig->mappedSequenceToString();
	struct Seed *nearestSeed, *prevSeed;

	//Testing
	// if(report){
	// 	cout << "contRightX_Drop: Continue extension on unitig " << sucUniSeq << endl;
	// 	cout << "contRightX_Drop: iniQoff: " << iniQoff << " hitLen: " << hitLen << " extLen: " << extLen << " lastSeedTmpScore: " << lastSeedTmpScore;
	// 	cout << " uniSeqPos: " << uniSeqPos << " explCount: " << explCount << endl;
	// 	UnitigColorMap<UnitigInfo> probUni = sucUnitig->getGraph()->find(Kmer("GCTATTGGCACACCAATCTATTCACCAGCAGATGG"));
	// 	cout << "contRightX_Drop: probUni.getSuccessors().hasSuccessors(): " << probUni.getSuccessors().hasSuccessors() << endl;
	// 	cout << "contRightX_Drop: sucUnitig->getSuccessors().hasSuccessors(): " << sucUnitig->getSuccessors().hasSuccessors() << endl;
	// 	// report = true;
	// }

	//Save the initial offset in the current unitig which we need for all nearest neighbor calculations
	iniSeqPos = uniSeqPos;
	//Find the nearest seed that we might be able to reach during our extension
	nearestSeed = searchRightNeighbor(sucUnitig->getData()->getData(*sucUnitig)->getSeed(sucUnitig->strand), iniQoff, extLen, iniSeqPos, prevSeed);
	//Perform the X-drop algorithm on the successive unitig
	tmpScore = lastSeedTmpScore;
	tmpSLen = 0;

	//Testing
	// if(report){
	// 	cout << "contRightX_Drop: Nearest neighbor search done" << endl;
	// }

	//We are done if we have reached the end of the query
	while(iniQoff + extLen + tmpSLen < q.length()){
		//Testing
		// if(report) cout << "contRightX_Drop: Current position in query is " << iniQoff + extLen + tmpSLen << endl;

		//Check whether we have reached the next seed
		if(nearestSeed != NULL && iniQoff + extLen + tmpSLen >= nearestSeed->offsetQ){
			//Testing
			// if(report){
			// 	cout << "contRightX_Drop: We have reached a seed" << endl;
			// 	cout << "contRightX_Drop: nearestSeed->offsetQ: " << nearestSeed->offsetQ << " iniQoff + extLen + tmpSLen: ";
			// 	cout << iniQoff + extLen + tmpSLen << endl;
			// }

			//Calculate the gain we get by incorporating the reached seed
			progress = nearestSeed->offsetQ + nearestSeed->len - (iniQoff + extLen + tmpSLen);
			//Update temporary seed length
			tmpSLen += progress;
			//Update current position in the unitig sequence
			uniSeqPos += progress;
			//Adjust number of remaining covered positions
			checkedPos -= progress;

			//If we have reached a seed check if it suffices to get a positive tempScore
			if((tmpScore += progress * mscore) > 0){
				//Update hit's length
				hitLen = extLen + tmpSLen;
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
			//Testing
			// if(report){
			// 	cout << "contRightX_Drop: We have not reached a seed" << endl;
			// 	cout << "contRightX_Drop: sucUnitig->getSuccessors().hasSuccessors(): " << sucUnitig->getSuccessors().hasSuccessors() << endl;
			// 	cout << "contRightX_Drop: sucUnitig->mappedSequenceToString(): " << sucUnitig->mappedSequenceToString() << endl;
			// }

			//Check up to which point we have to compare the unitig sequence
			if(!sucUnitig->getSuccessors().hasSuccessors()){
				overlap = 0;
			}

			//Testing
			// if(report) cout << "contRightX_Drop: Sequence overlap check done" << endl;

			//Check whether we have reached the end of the unitig's sequence
			if(uniSeqPos < sucUnitig->size - overlap){
				//Testing
				// if(report) cout << "contRightX_Drop: We have not yet reached the unitig's sequence end" << endl;

				//Are search criteria still fulfilled?
				if(checkedPos <= 0) break;

				//Check whether the score of our extension is positive
				if((tmpScore += compUScore(sucUniSeq[uniSeqPos], q[iniQoff + extLen + tmpSLen], mscore, mmscore)) > 0){
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
				//Testing
				// if(report) cout << "contRightX_Drop: We have reached the unitig's sequence end" << endl;

				//Check if the current unitig has successors
				if(overlap != 0){
					//Calculate the position in the next unitig's sequence we have to start with
					uniSeqPos = uniSeqPos - sucUnitig->size + overlap;

					//Testing
					// if(report) cout << "contRightX_Drop: Continue extension on successors" << endl;

					//Check out next unitig
					score += extendAtNextUnitig(sucUnitig->getSuccessors(), iniQoff, hitLen, extLen + tmpSLen, q, mscore, mmscore, X, tmpScore, uniSeqPos, extPth, explCount, quorum, searchSet, advIdx);
				}

				break;
			}
		}
	}

	//Testing
	// if(report) cout << "contRightX_Drop: Extension done" << endl;

	return score;
}

//This function continues an extension to the right on a successive unitig of a seed lying on the query's reverse complement considering a quorum and a search color set. Returns the maximum score reached.
int32_t contRightX_Drop_OnRevComp(const neighborIterator<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> &sucUnitig, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t &extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastSeedTmpScore, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	int32_t tmpScore, progress, score = 0;
	//Since we consider the overlap in the sequence's beginning for seed's on the reverse complementary strand the initial unitig position is always the same
	uint32_t iniSeqPos = sucUnitig->getGraph()->getK() - 1;
	uint32_t uniSeqPos = iniSeqPos;
	uint32_t tmpSLen;
	//Get the number of covered positions
	int32_t checkedPos = getSrchCritCov(*sucUnitig, quorum, searchSet, compOffset(uniSeqPos, 1, sucUnitig->size, sucUnitig->strand), sucUnitig->strand, advIdx);
	string sucUniSeq = sucUnitig->mappedSequenceToString();
	struct Seed *nearestSeed, *prevSeed;

	//Find the nearest seed that we might be able to reach during our extension
	nearestSeed = searchRightNeighbor(sucUnitig->getData()->getData(*sucUnitig)->getSeed(sucUnitig->strand), iniQoff, extLen, iniSeqPos, prevSeed);
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
			if((tmpScore += progress * mscore) > 0){
				//Update hit's length
				hitLen = extLen + tmpSLen;

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
			//Check whether we have reached the end of the unitig's sequence
			if(uniSeqPos < sucUnitig->size){
				//Check if the current position still fulfills the search criteria
				if(checkedPos <= 0) break;

				//Check whether the score of our extension is positive
				if((tmpScore += compUScore(sucUniSeq[uniSeqPos], q[iniQoff + extLen + tmpSLen], mscore, mmscore)) > 0){
					//Update score
					score += tmpScore;
					//Update hit's length
					hitLen = extLen + tmpSLen + 1;
					//Reset the temporary score
					tmpScore = 0;
				} else{
					//Check if the current extension is already too bad
					if(tmpScore < -X) break;
				}

				//Proceed with the next two positions
				++tmpSLen;
				++uniSeqPos;
				--checkedPos;
			} else{
				//Check if the current unitig has successors
				if(sucUnitig->getSuccessors().hasSuccessors()){
					//Calculate extension on successive unitigs
					score += extendAtNextUnitig_OnRevComp(sucUnitig->getSuccessors(), iniQoff, hitLen, extLen + tmpSLen, q, mscore, mmscore, X, tmpScore, extPth, explCount, quorum, searchSet, advIdx);
				}

				break;
			}
		}
	}

	return score;
}

//This function performs an extension on all possible predecessors of a unitig considering a quorum and a search color set and returns the maximum scoring one
int32_t extendAtPrevUnitig(const BackwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> bwIter, uint32_t qPos, uint32_t &hitLen, const uint32_t &tmpExtLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	uint16_t predID = 0;
	uint32_t tmpHitLen, maxHitLen = hitLen;
	int32_t maxScore = 0, currScore;

	//Check whether we have reached the maximum recursion depth of an extension
	if(++explCount > EXPLORED_UNITIGS_MAX){
		//Report this incident//TODO: If this happens the condition that an extended seed cannot be reached anymore is violated. Implement a procedure after the extension that merges such seeds with other seeds that would have reached them!
		//cerr << "Maximum recursion depth reached during extension" << endl;
		//Terminate this extension
		return 0;
	}

	//Iterate over all predecessors
	for(neighborIterator<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> nI = bwIter.begin(); nI != bwIter.end(); ++nI){
		//Temporary extension path
		list<uint16_t> tmpPth;
		//Note which predecessor we are on
		++predID;
		//Set tmpHitLen
		tmpHitLen = hitLen + tmpExtLen;

		//Calculate the score of an extension of a successor
		currScore = contLeftX_Drop(nI, qPos, tmpHitLen, q, mscore, mmscore, X, lastExtSeedTmpScore, tmpPth, explCount, quorum, searchSet, advIdx);

		//Check whether the score of the current successors extension is the best one found so far
		if(currScore > maxScore){
			//Update maxScore
			maxScore = currScore;
			//Update hit length
			maxHitLen = tmpHitLen;
			//Save the extension path
			extPth = tmpPth;
			//Save which predecessor we have chosen
			extPth.push_front(predID);
		}
	}

	//Save hit length
	hitLen = maxHitLen;

	//Return score
	return maxScore;
}

//This function performs an extension on all possible predecessors of a unitig considering a quorum and a search color set and returns the maximum scoring one. This function is explicitly designed for seeds lying on the query's reverse complement (considering the overlap between unitigs in sequences' beginning)
int32_t extendAtPrevUnitigOnRevComp(const BackwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> bwIter, uint32_t qPos, uint32_t &hitLen, const uint32_t &tmpExtLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, list<uint16_t> &extPth, const uint32_t &lead, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	uint16_t predID = 0;
	uint32_t tmpHitLen, maxHitLen = hitLen;
	int32_t maxScore = 0, currScore;

	//Check whether we have reached the maximum recursion depth of an extension
	if(++explCount > EXPLORED_UNITIGS_MAX){
		//Terminate this extension
		return 0;
	}

	//Iterate over all predecessors
	for(neighborIterator<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> nI = bwIter.begin(); nI != bwIter.end(); ++nI){
		//Temporary extension path
		list<uint16_t> tmpPth;
		//Note which predecessor we are on
		++predID;
		//Set tmpHitLen
		tmpHitLen = hitLen + tmpExtLen;
		//Calculate the score of an extension of a successor
		currScore = contLeftX_DropOnRevComp(nI, qPos, tmpHitLen, q, mscore, mmscore, X, lastExtSeedTmpScore, tmpPth, nI->size - lead, explCount, quorum, searchSet, advIdx);

		//Check whether the score of the current successors extension is the best one found so far
		if(currScore > maxScore){
			//Update maxScore
			maxScore = currScore;
			//Update hit length
			maxHitLen = tmpHitLen;
			//Save the extension path
			extPth = tmpPth;
			//Save which predecessor we have chosen
			extPth.push_front(predID);
		}
	}

	//Save hit length
	hitLen = maxHitLen;	
	//Return score
	return maxScore;
}

//This function starts the left extension for seeds lying on the query's reference strand considering a quorum and a search color 
//set using an iterative approach
void perfLeftX_Drop(Hit* hit, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const uint32_t 
	&quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	//The rank of a predecessor of a leading unitig (used to construct the extension path)
	uint16_t nr;
	//The length of the currently explored extension on the leading unitig
	int32_t extLen;
	//The number of times we have already switched unitigs during this extension
	uint32_t uniSwtchCnt = 0;
	//The k value are using
	int K = hit->origUni.getGraph()->getK();
	//A number of positions on the leading unitig fulfilling the search criteria
	int32_t nbffPos;
	//The progress we achieve due to reaching a seed
	int32_t progress;
	//The initial query position at the beginning of an exploration on the current leading unitig (might be different from 
	//currExt.offsetQ right after switching unitigs)
	int32_t posQ;
	//The leading unitig's sequence
	string uSeq;
	struct Seed *nearestSeed, *prevSeed;
	//An extension we are dealing with
	Ext currExt;
	//A vector of extensions to be continued on a predecessor of a leading unitig//TODO: This is only an artificial detour to allow 
	//the reproduction of results of the old recursive PLAST version for testing purposes. As soon as we are sure that the iterative
	// version of PLAST works this should be removed!
	stack<Ext> newExts;
	//Initialize DFS stack
	stack<Ext> exts;//TODO: We use a deque here as the container for the stack, because intuitively it sounds more efficient (e.g. compared to a vector). At some point we should try to verify this by a test!

	//Create the first extension right in front of the hit's start
	exts.push(Ext(-1 + hit->offU, -1 + hit->offQ, -1 + hit->offQ, hit->score, 0, hit->origUni, list<uint16_t>()));

	//Perform DFS through the graph while extending the seed
	while(!exts.empty()){
		//Get first extension on stack
		currExt = exts.top();
		//Remove first extension from stack
		exts.pop();

		//Try to extend current extension on the leading unitig as far as possible//

		//Set extension length on leading unitig
		extLen = 0;
		//Get the leading unitig's sequence
		uSeq = currExt.ldUni.mappedSequenceToString();
		//Compute the number of positions on which search criteria are fulfilled	
		nbffPos = getSrchCritCov(currExt.ldUni, quorum, searchSet, compOffset(currExt.offsetU, 1, currExt.ldUni.size, 
			currExt.ldUni.strand), !currExt.ldUni.strand, advIdx);
		//Search for a seed we could reach during this extension (tmpQoff has to be used here to match with offsetU, because offsetQ
		// is only updated after an explored continuation of an extension is also accepted due to a positive score which might not 
		//be the case when switching to the next unitig whereas offsetU is always updated when switching)
		nearestSeed = searchLeftNeighbor(currExt.ldUni.getData()->getData(currExt.ldUni)->getSeed(currExt.ldUni.strand), (uint32_t) 
			currExt.tmpQoff, (uint32_t) currExt.offsetU, prevSeed);//TODO: As soon as the iterative extension is completely implemented, it might be possible to modify input parameter types of searchLeftNeighbor() to avoid the necessity of these casts
		//The initial query position is always identical to currExt.tmpQoff in the beginning, but this might change during the ex-
		//tension on the leading unitig
		posQ = currExt.tmpQoff;

		//Go through the query up to the beginning or until search criteria are no longer fulfilled
		while(nbffPos-- > 0 && posQ >= extLen){
			//Check if we have reached the next seed
			if(nearestSeed != NULL && posQ - extLen < nearestSeed->offsetQ + nearestSeed->len){
				//Calculate the number of positions on the unitig we do not need to compare
				progress = posQ - extLen - nearestSeed->offsetQ;
				//Update extension length
				extLen += progress;
				//Update nbffPos
				nbffPos -= progress;
				//Update extension with information about the reached seed (the score has to be increased by 1, because progress 
				//does not consider that we have just moved forward in the beginning of this iteration)
				updateExtension(currExt, posQ, (progress + 1) * mscore, extLen, false);

				//Exclude reached seed from seed list
				if(prevSeed != NULL){
					//Link predecessor and successor of the reached seed
					prevSeed->nextSeed = nearestSeed->nextSeed;
				} else{
					//Set the reached seed's successor as head of the seed list
					currExt.ldUni.getData()->getData(currExt.ldUni)->setSeed(nearestSeed->nextSeed, currExt.ldUni.strand);
				}

				//Delete the reached seed
				free(nearestSeed);
				//Reset prevSeed
				prevSeed = NULL;
				//Search for next seed
				nearestSeed = searchLeftNeighbor(currExt.ldUni.getData()->getData(currExt.ldUni)->getSeed(currExt.ldUni.strand), (uint32_t) posQ, (uint32_t) currExt.offsetU, prevSeed);
				//Extend further
				++extLen;
			//Make sure there is still unitig sequence left to compare
			} else if(currExt.offsetU >= extLen){
				//Update extension with information about the next two compared positions
				updateExtension(currExt, posQ, compUScore(uSeq[currExt.offsetU - extLen], q[posQ - extLen], mscore, mmscore), extLen, false);

				//Terminate exploration if the score has become too bad
				if(currExt.tmpScore < -X) break;

				//Move on
				++extLen;
			} else{
				//Check if the leading unitig has predecessors and make sure we do not exceed the maximum number of unitig switches 
				//during this extension
				if(currExt.ldUni.getPredecessors().hasPredecessors() && uniSwtchCnt++ < EXPLORED_UNITIGS_MAX){
					//Reset predecessor's rank
					nr = 0;

					//Iterate over predecessors
					for(auto n = currExt.ldUni.getPredecessors().begin(); n != currExt.ldUni.getPredecessors().end(); ++n, ++nr){
						//Copy current extension
						Ext newExt = Ext(currExt);
						//Calculate sequence position at which we start on the next unitig
						newExt.offsetU = n->size - K;
						//Set tmpQoff correctly
						newExt.tmpQoff = posQ - extLen;
						//Update leading unitig
						newExt.ldUni = *n;
						//Extend extension path by new leading unitig
						newExt.pth.push_back(nr);
						//Add the new extension to the auxiliar vector to preserve the correct order
						newExts.push(newExt);
					}

					//Push the new extensions to the stack in the correct order
					while(!newExts.empty()){
						exts.push(newExts.top());
						newExts.pop();
					}
				}

				//Terminate current extension
				break;
			}
		}

		//Update info about the best extension found so far//

		//Check if the score of the current extension is better than what we found so far
		if(currExt.score > hit->score){
			//Update max score found
			hit->score = currExt.score;
			//Recalculate hit's length
			hit->length += hit->offQ - currExt.offsetQ;//TODO: We should think about whether it would not also make sense to change the type of offets inside the Hit Class
			//Clear extension path if it already exists
			decmprExtPth(hit->lExt);
			//Copy extension path
			hit->lExt = cmprExtPth(currExt.pth);
		}
	}
}

//This function starts the left extension for seeds lying on the query's reference strand considering a quorum and a search color set
void startLeftX_Drop(Hit* hit, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	//Initialization of auxiliary variables
	int32_t tmpScore = 0;
	uint32_t tmpExtLen = 1, progress, posU = hit->offU, posQ = hit->offQ;
	//Counter to count tries to explore a further unitig
	uint32_t explCount;
	//Get the number of covered positions coming from the sequence's end considering that we do not start in the very end
	int32_t checkedPos = getSrchCritCov(hit->origUni, quorum, searchSet, posU - tmpExtLen, false, advIdx);
	string uSeq = hit->origUni.mappedSequenceToString();
	list<uint16_t> extPth;
	struct Seed *nearestSeed, *prevSeed = NULL;

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
			if((tmpScore += progress * mscore) > 0){
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

			//Check if the reached seed had a predecessor
			if(prevSeed != NULL){
				//Link predecessor and successor of the reached seed
				prevSeed->nextSeed = nearestSeed->nextSeed;
			} else{
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
			if(checkedPos == 0) break;

			//Compare the next two bases and check whether this changes temporary score's sign
			if((tmpScore += compUScore(uSeq[posU - tmpExtLen], q[posQ - tmpExtLen], mscore, mmscore)) > 0){
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
				//Initialize explCount
				explCount = 0;
				//Continue the extension on the predecessive unitig
				hit->score += extendAtPrevUnitig(hit->origUni.getPredecessors(), posQ - tmpExtLen, hit->length, tmpExtLen, q, mscore, mmscore, X, tmpScore, extPth, explCount, quorum, searchSet, advIdx);
			}

			break;
		}
	}

	//Compress extension path
	hit->lExt = cmprExtPth(extPth);
}

//This function starts the left extension for seeds lying on the query's reverse complement considering a quorum and a search color set
void startLeftX_Drop_OnRevComp(Hit* hit, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	//Initialization of auxiliary variables
	int32_t tmpScore = 0;
	uint32_t tmpExtLen = 1, progress, overlap = 0;
	//Counter to count tries to explore a further unitig
	uint32_t explCount;
	int32_t checkedPos = getSrchCritCov(hit->origUni, quorum, searchSet, compOffset(hit->offU - tmpExtLen, 1, hit->origUni.size, false), true, advIdx);
	uint32_t k = hit->origUni.getGraph()->getK();
	list<uint16_t> extPth;
	string uSeq = hit->origUni.mappedSequenceToString();
	struct Seed *nearestSeed, *prevSeed = NULL;

	//Find the nearest seed that we could reach
	nearestSeed = searchLeftNeighbor(hit->origUni.getData()->getData(hit->origUni)->getSeed(hit->origUni.strand), hit->offQ, hit->offU, prevSeed);

	//Check how far we should explore the unitig's sequence
	if(hit->origUni.getPredecessors().hasPredecessors()) overlap = k - 1;

	//Go through the query up to the beginning
	while(hit->offQ >= tmpExtLen){
		//Check whether we have reached the next seed
		if(nearestSeed != NULL && hit->offQ - tmpExtLen < nearestSeed->offsetQ + nearestSeed->len){
			//Calculate the progress we have incorporating the reached seed
			progress = hit->offQ - tmpExtLen - nearestSeed->offsetQ + 1;
			//Temporary extension length
			tmpExtLen += progress;
			//Adjust number of remaining covered positions
			checkedPos -= progress;

			//Check if our temporary score becomes positive using this hit
			if((tmpScore += progress * mscore) > 0){
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
			}

			//Check if the reached seed had a predecessor
			if(prevSeed != NULL){
				//Link predecessor and successor of the reached seed
				prevSeed->nextSeed = nearestSeed->nextSeed;
			} else{
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
			//Check whether quorum has to be checked
			if(checkedPos <= 0) break;

			//Compare the next two bases and check whether this changes temporary score's sign
			if((tmpScore += compUScore(uSeq[hit->offU - tmpExtLen], q[hit->offQ - tmpExtLen], mscore, mmscore)) > 0){
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
				break;
			}

			++tmpExtLen;
			--checkedPos;
		} else{
			//Decrease temporary extension length
			--tmpExtLen;

			//Check if this unitig has a predecessor
			if(overlap != 0){
				//Initialize explCount
				explCount = 0;
				//Continue the extension on the predecessive unitig
				hit->score += extendAtPrevUnitigOnRevComp(hit->origUni.getPredecessors(), hit->offQ - tmpExtLen, hit->length, tmpExtLen, q, mscore, mmscore, X, tmpScore, extPth, overlap - (hit->offU - tmpExtLen), explCount, quorum, searchSet, advIdx);
			}

			break;
		}
	}

	//A hit's start offset must never be inside the overlap at a unitig sequence's end
	if(hit->offU > hit->origUni.size - k) mvStartToValUni(hit, extPth);

	//If hit is invalid we do not need to compress its left extension path
	if(hit->score > 0) hit->lExt = cmprExtPth(extPth);
}

//This function continues a left extension on a predecessive unitig of a seed lying on the query's reference strand considering a quorum and a search color set and returns the achieved score
int32_t contLeftX_Drop(const neighborIterator<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> &prevUni, uint32_t qPos, uint32_t &hitLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastSeedTmpScore, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	//bool inSeedList;
	int32_t tmpScore, score = 0;
	uint32_t tmpExtLen = 1;
	//Calculate our offset position inside the unitig sequence (+1, because otherwise it is not possible to calculate the correct gain if reaching a seed)
	uint32_t uPos= prevUni->size - prevUni->getGraph()->getK() + 1;
	int32_t checkedPos = getSrchCritCov(*prevUni, quorum, searchSet, compOffset(uPos - tmpExtLen, 1, prevUni->size, prevUni->strand), !prevUni->strand, advIdx);
	string prevUniSeq = prevUni->mappedSequenceToString();
	struct Seed *nearestSeed, *prevSeed = NULL;

	//Find the nearest seed that we might be able to reach during our extension
	nearestSeed = searchLeftNeighbor(prevUni->getData()->getData(*prevUni)->getSeed(prevUni->strand), qPos, uPos, prevSeed);
	//Initialize our temporary score with the last temporary score of its successive unitig
	tmpScore = lastSeedTmpScore;

	//We are done if we reach the beginning of the query
	while(qPos >= tmpExtLen){
		//Check whether we have reached a nearest seed
		if(nearestSeed != NULL && qPos - tmpExtLen <= nearestSeed->offsetQ + nearestSeed->len - 1){
			//Update the temporary score
			tmpScore += (qPos - tmpExtLen - nearestSeed->offsetQ + 1) * mscore;
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
				//Update current position in q the unitig
				qPos = nearestSeed->offsetQ;
				uPos = nearestSeed->offsetU;	
				//Reset temporary score and extension length
				tmpScore = 0;
				tmpExtLen = 1;
			}

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
			//Reset prevSeed
			prevSeed = NULL;
			//Search for the next seed to reach
			nearestSeed = searchLeftNeighbor(prevUni->getData()->getData(*prevUni)->getSeed(prevUni->strand), qPos, uPos, prevSeed);
		} else if(uPos >= tmpExtLen){//Check whether we have already reached the beginning of the unitig sequence
			//Check whether quorum has to be checked
			if(checkedPos <= 0) break;

			//Compare the next two bases and check whether our temporary score is becoming > 0 by this
			if((tmpScore += compUScore(prevUniSeq[uPos - tmpExtLen], q[qPos - tmpExtLen], mscore, mmscore)) > 0){
				//Update positions in q and the unitig
				uPos -= tmpExtLen;
				qPos -= tmpExtLen;
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
				score += extendAtPrevUnitig(prevUni->getPredecessors(), qPos - tmpExtLen, hitLen, tmpExtLen, q, mscore, mmscore, X, tmpScore, extPth, explCount, quorum, searchSet, advIdx);
			}

			break;
		}
	}

	return score;
}

//This function continues a left extension on a predecessive unitig of a seed lying on the query's reverse complement considering a quorum and a search color set and returns the achieved score
int32_t contLeftX_DropOnRevComp(const neighborIterator<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> &prevUni, uint32_t qPos, uint32_t &hitLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastSeedTmpScore, list<uint16_t> &extPth, uint32_t uPos, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	int32_t tmpScore, score = 0;
	uint32_t tmpExtLen = 1, overlap = 0;

	int32_t checkedPos = getSrchCritCov(*prevUni, quorum, searchSet, compOffset(uPos - tmpExtLen, 1, prevUni->size, prevUni->strand), !prevUni->strand, advIdx);
	uint32_t k = prevUni->getGraph()->getK();
	uint32_t progress;
	string prevUniSeq = prevUni->mappedSequenceToString();
	struct Seed *nearestSeed, *prevSeed = NULL;

	//Find the nearest seed that we might be able to reach during our extension
	nearestSeed = searchLeftNeighbor(prevUni->getData()->getData(*prevUni)->getSeed(prevUni->strand), qPos, uPos, prevSeed);

	//Check how far we should explore the unitig's sequence
	if(prevUni->getPredecessors().hasPredecessors()) overlap = k - 1;

	//Initialize our temporary score with the last temporary score of its successive unitig
	tmpScore = lastSeedTmpScore;
	
	//We are done if we reach the beginning of the query
	while(qPos >= tmpExtLen){
		//Check whether we have reached a nearest seed
		if(nearestSeed != NULL && qPos - tmpExtLen <= nearestSeed->offsetQ + nearestSeed->len - 1){
			//Calculate the progress with have by incorporating the seed
			progress = qPos - tmpExtLen - nearestSeed->offsetQ + 1;
			//Update the temporary score
			tmpScore += progress * mscore;
			//Calculate the gain we have by incorporating the seed into our extension and update the temporary extension length
			tmpExtLen = qPos - nearestSeed->offsetQ + 1;

			//Check whether we get a score larger 0 by incorporating the reached seed
			if(tmpScore > 0){
				//Update hit's length
				hitLen += qPos - nearestSeed->offsetQ;
				//Update score
				score += tmpScore;
				//Reset tmpScore
				tmpScore = 0;
				//Reset temporary extension length
				tmpExtLen = 1;
				//Update current position in q the unitig
				qPos = nearestSeed->offsetQ;
				uPos = nearestSeed->offsetU;
			}

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
			//Decrease number of remaining covered positions
			checkedPos -= progress; 
			//Reset prevSeed
			prevSeed = NULL;
			//Search for the next seed to reach
			nearestSeed = searchLeftNeighbor(prevUni->getData()->getData(*prevUni)->getSeed(prevUni->strand), qPos, uPos, prevSeed);
		} else if(uPos >= tmpExtLen + overlap){//Check whether we have already reached the beginning of the unitig sequence
			//Check whether quorum has to be checked
			if(checkedPos <= 0) break;

			//Compare the next two bases and check whether our temporary score is becoming > 0 by this
			if((tmpScore += compUScore(prevUniSeq[uPos - tmpExtLen], q[qPos - tmpExtLen], mscore, mmscore)) > 0){
				//Update positions in q and the unitig
				uPos -= tmpExtLen;
				qPos -= tmpExtLen;
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
			if(overlap != 0){
				//Try to continue the extension on the predecessive unitigs
				score += extendAtPrevUnitigOnRevComp(prevUni->getPredecessors(), qPos - tmpExtLen, hitLen, tmpExtLen, q, mscore, mmscore, X, tmpScore, extPth, overlap - (uPos - tmpExtLen), explCount, quorum, searchSet, advIdx);
			}

			break;
		}
	}

	return score;
}

//This function checks if a hit's start offset for the gapped extension lies inside the overlap of a unitig sequence's end. If so it moves the start position to a unitig where it is not inside the overlap anymore.
void mvStartToValUni(Hit* h, list<uint16_t>& lExtPth){
	list<uint16_t> rExtPth;

	//Check if there is exists a right extension path to follow when switching unitigs and decompress it
	if(h->rExt.nbElem > 0) rExtPth = decmprExtPth(h->rExt);

	//Switch unitigs
	switUni(h->offU, h->origUni, lExtPth, rExtPth);

	//If there is still a right extension path it needs to be compressed again
	if(h->rExt.nbElem > 0) h->rExt = cmprExtPth(rExtPth);
}

//This function moves an offset from a given unitig to its successor while keeping left and right extension paths updated. If the offset at the successive unitig lies inside the overlap at the unitig sequence's end the function calls itself recursively
void switUni(uint32_t &offset, UnitigColorMap<UnitigInfo> &currUni, list<uint16_t> &lExtPth, list<uint16_t> &rExtPth){
	uint16_t sucCount = 1, predCount = 1;

	//Traverse successors of the current unitig
	for(neighborIterator<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> suc = currUni.getSuccessors().begin(); suc != currUni.getSuccessors().end(); ++suc, ++sucCount){
		//Check if we have found the required successor (if such exists)
		if(!rExtPth.empty() && sucCount < rExtPth.front()) continue;

		//Iterate over predecessors of found successor
		for(neighborIterator<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> pred = suc->getPredecessors().begin(); pred != suc->getPredecessors().end(); ++pred, ++predCount){
			//Check if we have found the unitig we came from
			if(*pred == currUni) break;
		}

		//Save count in left extension path
		lExtPth.push_front(predCount);
		//Delete successor from right extension path if such exists
		if(!rExtPth.empty()) rExtPth.pop_front();

		//Update offset
		offset -= currUni.len;
		//Update current unitig
		currUni = *suc;
		//End iteration
		break;
	}

	//Check if it is not necessary to go on anymore (which is the case if either 1. the current unitig has no successor or 2. the offset is not inside the overlap anymore)
	if(!currUni.getSuccessors().hasSuccessors() || offset <= currUni.size - currUni.getGraph()->getK()) return;

	//Move to next unitig
	switUni(offset, currUni, lExtPth, rExtPth);
}