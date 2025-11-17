#include "Extension.h"
#include "Hit.h"
#include "Search.h"
#include "UnitigInfo.cpp"
#include <bifrost/CompactedDBG.hpp>
#include <bifrost/UnitigMap.hpp>
#include <queue>
#include <tuple>
#include <algorithm>
#include <vector>

//This function initiates the extension on all successors of a unitig and returns the best one considering a quorum and a search color set
int32_t extendAtNextUnitig(const ForwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &uniPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	uint16_t sucID;
	uint32_t tmpHitLen = 0;
	int32_t maxScore = 0, currScore;
	
	//Check whether we have reached the maximum recursion depth of an extension
	if(++explCount > MAXRECURSIONDEPTH){
		//Report this incident
		//cerr << "Maximum recursion depth reached during extension. Position in q: " << iniQoff << endl;
		//Terminate this extension
		return 0;
	}

	maxScore = 0;
	sucID = 0;

	
	//Iterate over successors
	for(neighborIterator<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> nI = sucIter.begin(); nI != sucIter.end(); ++nI){
		//Temporary extention path
		list<uint16_t> tmpPth;
		//Note which successor we are on
		++sucID;
		//Calculate the score of an extension of a successor
		//cout << "" << endl;
		//cout << "unitig right now: " << nI->mappedSequenceToString() << endl;
		currScore = contRightX_Drop(nI, iniQoff, tmpHitLen, extLen, q, mscore, mmscore, X, lastExtSeedTmpScore, uniPos, tmpPth, explCount, quorum, searchSet, advIdx);

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


int32_t extendAtNextUnitig_BFS(const ForwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &uniPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	uint16_t sucID = 0;
	int32_t maxScore = 0;
	bool check = true;

	//cout << "  " << endl;
	//cout << "  " << endl;
	//cout << "  " << endl;

	//cout << "first unitig succesor" << endl;


	list<uint16_t> tempPath;


	queue<tuple<shorterTemp, uint32_t, int32_t, pathList, uint32_t, uint32_t, uint32_t>> queueTest;
	list<uint16_t> bestPath;


	for(shorterTemp nI = sucIter.begin(); nI != sucIter.end(); ++nI){
		//cout << "successor: " << nI->mappedSequenceToString() << endl;
		++sucID;
		uint32_t tempextLen = extLen;
		uint32_t tmpHitLen = hitLen;
		uint32_t tmpuniPos = uniPos;
		int32_t tmpScore = lastExtSeedTmpScore;
		//cout << "uniPos: " << uniPos << endl;
        tempPath.clear();
		int32_t startScore = contRightX_Drop_BFS(nI, iniQoff, tmpHitLen, tempextLen, q, mscore, mmscore, X, tmpScore, tmpuniPos, tempPath, explCount, quorum, searchSet, advIdx, check);
		tempPath.push_back(sucID);
		if(!check){
			//cout << "check failed" << endl;
		}
		if(check){
			//cout << "Add from startRight_XDrop" << endl;
			//cout << sucID << endl;
			queueTest.push(std::make_tuple(nI,startScore,tmpScore,tempPath,tmpHitLen,tempextLen, tmpuniPos));
		}
		if (startScore > maxScore) {
            maxScore = startScore;
            //bestPath = tempPath;
			hitLen = tmpHitLen;
        }
	}


	//cout << get<1>(queueTest.front()) << endl;


	uint16_t stop = 0;
	int counter = 0;
	while (!queueTest.empty()) {

		auto tempFront = queueTest.front();
		queueTest.pop();

		//cout << " " << endl;
		//cout << " " << endl;
		//cout << " " << endl;
		//cout << "accessing unitig from queue" << endl;
		//cout << "std::get<0>(tempFront)->mappedSequenceToString(): " << std::get<0>(tempFront)->mappedSequenceToString() << endl;

		shorterTemp currUnitig = std::get<0>(tempFront);
		uint32_t currScore = std::get<1>(tempFront);
		int32_t currtmpScore = std::get<2>(tempFront);
		pathList currPath = std::get<3>(tempFront);
		uint32_t currHitLen = std::get<4>(tempFront);
		uint32_t currextLen = std::get<5>(tempFront);
		uint32_t curruniPos = std::get<6>(tempFront);

		//cout << "iniQoff: " << iniQoff << endl;
		//cout << "HitLen: " << currHitLen << endl;
		//cout << "queue length: " << q.length() << endl;
		//cout << "currextLen: " << currextLen << endl;

		
		auto& tempcurrUnitig = *currUnitig;
		shorterContainer sucIter2 = tempcurrUnitig.getSuccessors();

		/*
		if((currextLen + iniQoff >= q.length())) {
			cout << " " << endl;
			cout << "removed from queue longer than queue" << endl;
		}
		*/
		if((currextLen + iniQoff < q.length())) {
			sucID = 0;
			for(shorterTemp nI = sucIter2.begin(); nI != sucIter2.end(); ++nI){
				check = true;
				++sucID;
				//cout << " " << endl;
				//cout << "going trough successor: " << nI->mappedSequenceToString() << endl;
				uint32_t tmpHitLen = currHitLen;
				int32_t tmpScore = currtmpScore;
                tempPath = currPath;
				uint32_t tmpExtLen = currextLen;
				uint32_t nextUniPos = curruniPos; // + (tempcurrUnitig.getGraph()->getK() - 1);
				//cout << "extLen_before: " << tmpExtLen << endl;
				//cout << "HitLen_before: " << tmpHitLen << endl;
				//cout << "nextUniPos: " << nextUniPos << endl;
				int32_t tempScore = contRightX_Drop_BFS(nI, iniQoff, tmpHitLen, tmpExtLen, q, mscore, mmscore, X, tmpScore, nextUniPos, tempPath, explCount, quorum, searchSet, advIdx, check);
				//cout << tempScore << endl;
				//tempScore += currScore;
				//std::cout << "nextUniPos after call = " << nextUniPos << std::endl;
				//cout << "extLen_after " << tmpExtLen << endl;
				//cout << "HitLen_after: " << tmpHitLen << endl;
				tempPath.push_back(sucID);
				//cout << boolalpha;
				//cout << check << endl;
				if(check == true){
					//cout << "Added in while loop " << counter << endl;
					//cout << sucID << endl;
					queueTest.push(std::make_tuple(nI, (currScore+tempScore), tmpScore, tempPath, tmpHitLen, tmpExtLen, nextUniPos));
                }
				int32_t scoreNow = currScore+tempScore;
				if (scoreNow > maxScore) {
                    maxScore = currScore+tempScore;
                    //bestPath = tempPath;
                    //hitLen = tmpHitLen;
					//uniPos = nextUniPos;
				}
				//cout << "maxScore rn: " << maxScore << endl;
			}
		}
		++stop;
		++counter;
	}

	//cout << " " << endl;

	//cout << "final maxScore: " << maxScore << endl;
	//cout << "Final maxScore: " << maxScore << ", hitLen: " << hitLen << ", uniPos: " << uniPos << endl;
	extPth = bestPath;
	return maxScore;
}

int32_t extendAtNextUnitig_BFS_SMART1(const ForwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &uniPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	uint16_t sucID = 0;
	int32_t maxScore = 0;
	bool check = true;
	uint numOfUnitig = 5;

	shorterPrioQueue bestUnitigsQueue(prioLongest);


	shorterVector bestUnitigs;
	list<uint16_t> tempPath;
	list<uint16_t> bestPath;


	shorterVector tempBestUnitigs;
	for(shorterTemp nI = sucIter.begin(); nI != sucIter.end(); ++nI){
		++sucID;
		uint32_t tempextLen = extLen;
		uint32_t tmpHitLen = hitLen;
		uint32_t tmpuniPos = uniPos;
		int32_t tmpScore = lastExtSeedTmpScore;
		tempPath.clear();
		int32_t tempScore = contRightX_Drop_BFS(nI, iniQoff, tmpHitLen, tempextLen, q, mscore, mmscore, X, tmpScore, tmpuniPos, tempPath, explCount, quorum, searchSet, advIdx, check);
		tempPath.push_back(sucID);
		if(check){
			//cout << "Added from startRight_XDrop" << endl;
			//cout << sucID << endl;
			bestUnitigsQueue.push(make_tuple(nI,tempScore,tmpScore,tempPath,tmpHitLen,tempextLen, tmpuniPos));
		}
		if(tempScore > maxScore){
			maxScore = tempScore;
            bestPath = tempPath;
            hitLen = tmpHitLen;
		}
	}

	for(uint i = 0; i < numOfUnitig; i++) {

		if(bestUnitigsQueue.empty()){
			break;
		} else {
			bestUnitigs.push_back(bestUnitigsQueue.top());
			bestUnitigsQueue.pop(); 
		}

	}

	//cout << get<1>(bestUnitigs[0]) << endl;

	int counter = 0;
	while(!(bestUnitigs.empty())){
		shorterVector allSucUnitigs;
		shorterPrioQueue allSucUnitigsQueue(prioLongest);

		for(uint i = 0; i < bestUnitigs.size(); i++) {
			shorterTuple temptuple = bestUnitigs[i];
			shorterTemp currUnitig = std::get<0>(temptuple);
			uint32_t currScore = std::get<1>(temptuple);
			int32_t currtmpScore = std::get<2>(temptuple);
			pathList currPath = std::get<3>(temptuple);
			uint32_t currHitLen = std::get<4>(temptuple);
			uint32_t currextLen = std::get<5>(temptuple);
			uint32_t currUniPos = std::get<6>(temptuple);
			auto& tempcurrUnitig = *currUnitig;
			shorterContainer sucIter2 = tempcurrUnitig.getSuccessors();
			if((currextLen + iniQoff < q.length())) {
				sucID = 0;
				for(shorterTemp nI = sucIter2.begin(); nI != sucIter2.end(); ++nI){
					++sucID;
					check = true;
					uint32_t tmpHitLen = currHitLen;
					int32_t tmpScore = currtmpScore;
                	tempPath = currPath;
					uint32_t tmpExtLen = currextLen;
					uint32_t nextUniPos = currUniPos;
					int32_t tempScore = contRightX_Drop_BFS(nI, iniQoff, tmpHitLen, tmpExtLen, q, mscore, mmscore, X, tmpScore, nextUniPos, tempPath, explCount, quorum, searchSet, advIdx, check);
					//cout << tempScore << endl;
					tempPath.push_back(sucID);
					if(check){
						//cout << "Added in while loop " << counter << endl;
						//cout << sucID << endl;
						allSucUnitigsQueue.push(make_tuple(nI, (currScore+tempScore), tmpScore, tempPath, tmpHitLen, tmpExtLen, nextUniPos));
					}
					int32_t scoreNow = (currScore+tempScore);
					if(scoreNow > maxScore){
						maxScore = currScore+tempScore;
                    	bestPath = tempPath;
                    	hitLen = tmpHitLen;
						uniPos = nextUniPos;
					}
				}
			}
		}
		bestUnitigs.clear();
		for(uint i = 0; i < numOfUnitig; i++) {
			if(allSucUnitigsQueue.empty()){
				break;
			} else{
				bestUnitigs.push_back(allSucUnitigsQueue.top());
				allSucUnitigsQueue.pop();
			}
		}
	++counter;
	}
	extPth = bestPath;
	return maxScore;
}



int32_t extendAtNextUnitig_BFS_SMART2(const ForwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &uniPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	uint16_t sucID = 0;
	int32_t maxScore = 0;
	bool check = true;
	uint numOfUnitig = 1000;

	shorterVector bestUnitigs;
	vector<int> bestScores;
	list<uint16_t> tempPath;
	list<uint16_t> bestPath;


	for(shorterTemp nI = sucIter.begin(); nI != sucIter.end(); ++nI){
		++sucID;
		uint32_t tempextLen = extLen;
		uint32_t tmpHitLen = hitLen;
		uint32_t tmpuniPos = uniPos;
		int32_t tmpScore = lastExtSeedTmpScore;
		tempPath.clear();
		int32_t tempScore = contRightX_Drop_BFS(nI, iniQoff, tmpHitLen, tempextLen, q, mscore, mmscore, X, tmpScore, tmpuniPos, tempPath, explCount, quorum, searchSet, advIdx, check);
		tempPath.push_back(sucID);
		int temptempScore = tempScore;
		if(check && bestUnitigs.size() <= numOfUnitig){
			bestUnitigs.push_back(std::make_tuple(nI,tempScore,tmpScore,tempPath,tmpHitLen,tempextLen, tmpuniPos));
			bestScores.push_back(tempScore);
		} else {
			if(check) {
				int tempMinElement = *min_element(bestScores.begin(), bestScores.end());
				if(tempMinElement < temptempScore){
					for(uint i = 0; i < bestScores.size(); i++){
						if(bestScores[i] == tempMinElement){
							bestScores[i] = temptempScore;
							bestUnitigs[i] = make_tuple(nI,tempScore,tmpScore,tempPath,tmpHitLen,tempextLen, tmpuniPos);
							break;
						}	
					}
				}
			}	
		}
		if(tempScore > maxScore){
			maxScore = tempScore;
            bestPath = tempPath;
            hitLen = tmpHitLen;
		}
	}

	for(list<uint16_t>::iterator n = bestPath.begin(); n != bestPath.end(); ++n){
		//cout << *n << endl;
	}
	//cout << (bestUnitigs.empty()?"Liste empty":"Liste full") << endl;
	//cout << (bestUnitigs.size()) << endl;

	while(!(bestUnitigs.empty())){
		for(uint i = 0; i < bestUnitigs.size(); i++) {
			shorterTuple temptuple = bestUnitigs[0];
			shorterTemp currUnitig = std::get<0>(temptuple);
			uint32_t currScore = std::get<1>(temptuple);
			int32_t currtmpScore = std::get<2>(temptuple);
			pathList currPath = std::get<3>(temptuple);
			uint32_t currHitLen = std::get<4>(temptuple);
			uint32_t currextLen = std::get<5>(temptuple);
			uint32_t currUniPos = std::get<6>(temptuple);
			auto& tempcurrUnitig = *currUnitig;
			shorterContainer sucIter2 = tempcurrUnitig.getSuccessors();
			bestUnitigs.erase(bestUnitigs.begin() + 0);
			bestScores.erase(bestScores.begin() + 0);
			if((currextLen + iniQoff < q.length())) {
				sucID = 1;
				for(shorterTemp nI = sucIter2.begin(); nI != sucIter2.end(); ++nI){
					check = true;
					uint32_t tmpHitLen = currHitLen;
					int32_t tmpScore = currtmpScore;
                	tempPath = currPath;
					uint32_t tmpExtLen = currextLen;
					uint32_t nextUniPos = currUniPos;
					int32_t tempScore = contRightX_Drop_BFS(nI, iniQoff, tmpHitLen, tmpExtLen, q, mscore, mmscore, X, tmpScore, nextUniPos, tempPath, explCount, quorum, searchSet, advIdx, check);
					//cout << "sucID: " << sucID << endl;
					//cout << tempScore << endl;
					tempPath.push_back(sucID);
					int temptempScore = currScore+tempScore;
					if(check &&  bestUnitigs.size() <= numOfUnitig){
						bestUnitigs.push_back(make_tuple(nI, (currScore+tempScore), tmpScore, tempPath, tmpHitLen, tmpExtLen, nextUniPos));
						bestScores.push_back((currScore+tempScore));
					} else {
						if(check) {
							int tempMinElement = *min_element(bestScores.begin(), bestScores.end());
							if(tempMinElement < temptempScore){
								for(uint i = 0; i < bestScores.size(); i++){
									if(bestScores[i] == tempMinElement){
										bestScores[i] = temptempScore;
										bestUnitigs[i] = make_tuple(nI, (currScore+tempScore), tmpScore, tempPath, tmpHitLen, tmpExtLen, nextUniPos);
										break;
									}	
								}
							}
						}	
					}
					++sucID;
					if(temptempScore > maxScore){
						maxScore = temptempScore;
						bestPath = currPath;
                    	hitLen = currHitLen;
						uniPos = currUniPos;
					}
				}
			}
		}
	}
	extPth = bestPath;
	return maxScore;
}



//  							input Queue       iniQof    q      mscore   mmscore  X      extPath        explCount  quorum   searchSet            advIdx maxScore numOfBases hitLen
//using inputTypes = tuple<queue<shorterTuple2>,uint32_t,string,uint16_t,int16_t,int16_t,list<uint16_t>, uint32_t,uint32_t,list<pair<string, size_t>>,bool,int32_t,int32_t,uint32_t>;

//	  						PrioQueue		maxScore  hitLen
//using outputTypes = tuple<shorterPrioQueue2,int32_t, uint32_t>;

outputTypes calcUnitigsMitBasenberechnung(inputTypes extensionCalcInputs) {
	queue<shorterTuple2> extensionQueue 	= get<0>(extensionCalcInputs);
	uint32_t iniQoff 						= get<1>(extensionCalcInputs);
	string q 								= get<2>(extensionCalcInputs);
	uint16_t mscore 						= get<3>(extensionCalcInputs);
	int16_t mmscore 						= get<4>(extensionCalcInputs);
	int16_t X 								= get<5>(extensionCalcInputs);
	list<uint16_t> extPath 					= get<6>(extensionCalcInputs);
	uint32_t explCount 						= get<7>(extensionCalcInputs);
	uint32_t quorum 						= get<8>(extensionCalcInputs);
	list<pair<string, size_t>> searchSet 	= get<9>(extensionCalcInputs);
	bool advIdx 							= get<10>(extensionCalcInputs);
	int32_t maxScore 						= get<11>(extensionCalcInputs);
	int32_t numOfBases 						= get<12>(extensionCalcInputs);
	uint32_t hitLen 						= get<13>(extensionCalcInputs);

	shorterPrioQueue2 bestUnitigsPrioQueue(prioLongest2);

	while(!extensionQueue.empty()){
		shorterTuple2 currExtension = extensionQueue.front();
		extensionQueue.pop();

		shorterTemp currUnitig 		= get<0>(currExtension);
		uint32_t currScore 			= get<1>(currExtension);
		int32_t currtmpScore 		= get<2>(currExtension);
		pathList currPath 			= get<3>(currExtension);
		uint32_t currHitLen 		= get<4>(currExtension);
		uint32_t currextLen 		= get<5>(currExtension);
		uint32_t curruniPos 		= get<6>(currExtension);

		auto& tempcurrUnitig = *currUnitig;
		shorterContainer sucIter2 = tempcurrUnitig.getSuccessors();

		if((currextLen + iniQoff < q.length())){
			uint16_t sucID = 0;
			for(shorterTemp nI = sucIter2.begin(); nI != sucIter2.end(); ++nI){
				++sucID;
				bool check = false;
				uint32_t nextHitLen 	= currHitLen;
				int32_t nextScore 		= currtmpScore;
                list<uint16_t>tempPath 	= currPath;
				uint32_t nextExtLen 	= currextLen;
				uint32_t nextUniPos 	= curruniPos;
				int32_t tempNumOfBases 	= numOfBases;

				int32_t addScore = contRightX_Drop_BFS_2(nI, iniQoff, nextHitLen, nextExtLen, q, mscore, mmscore, X, nextScore, nextUniPos, tempPath, explCount, quorum, searchSet, advIdx, check, tempNumOfBases);
				tempPath.push_back(sucID);
				int32_t fullScore = currScore+addScore;

				if((tempNumOfBases == 0) && !check){
					//cout << "(tempNumOfBases == 0) && !check" << endl;
					bestUnitigsPrioQueue.push(make_tuple(nI,fullScore,nextScore,tempPath,nextHitLen,nextExtLen, nextUniPos, 1));
				} else if(check && tempNumOfBases != 0){
					//cout << "check && tempNumOfBases != 0" << endl;
					extensionQueue.push(make_tuple(nI,fullScore,nextScore,tempPath,nextHitLen,nextExtLen,nextUniPos,tempNumOfBases));
				}
				if (fullScore > maxScore) {
                    maxScore = fullScore;
                    hitLen = nextHitLen;
				}
			}
		}
	}

	outputTypes returnValues = make_tuple(bestUnitigsPrioQueue,maxScore,hitLen);
	return returnValues;
}



 tuple<shorterVector2> getBestUnitigs(shorterPrioQueue2 bestUnitigsPrioQueue, uint numOfAcceptedExtensions){
	shorterVector2 bestUnitigs;
	
	for(uint i = 0; i < numOfAcceptedExtensions; i++) {
		if(bestUnitigsPrioQueue.empty()){
			break;
		} else {
			bestUnitigs.push_back(bestUnitigsPrioQueue.top());
			bestUnitigsPrioQueue.pop(); 
		}
	}

	return bestUnitigs;
}




int32_t extendAtNextUnitig_BFS_SMART3(const ForwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &uniPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	uint16_t sucID = 0;
	int32_t maxScore = 0;
	uint numOfUnitig = 1000000;
	int32_t numOfBases = 30;
	bool check = false;


	shorterPrioQueue2 bestUnitigsQueue(prioLongest2);

	queue<tuple<shorterTemp, uint32_t, int32_t, pathList, uint32_t, uint32_t, uint32_t, int>> finishBases;

	shorterVector2 bestUnitigs;
	list<uint16_t> tempPath;
	list<uint16_t> bestPath;


	for(shorterTemp nI = sucIter.begin(); nI != sucIter.end(); ++nI){
		++sucID;
		uint32_t tempextLen = extLen;
		uint32_t tmpHitLen = hitLen;
		uint32_t tmpuniPos = uniPos;
		int32_t tmpScore = lastExtSeedTmpScore;
		tempPath.clear();
		int32_t tmpNumOfBases = numOfBases;
		int32_t tempScore = contRightX_Drop_BFS_2(nI, iniQoff, tmpHitLen, tempextLen, q, mscore, mmscore, X, tmpScore, tmpuniPos, tempPath, explCount, quorum, searchSet, advIdx, check, tmpNumOfBases);
		tempPath.push_back(sucID);
		if(check && tmpNumOfBases == 0){
			//cout << "check && tmpNumOfBases == 0" << endl;
			bestUnitigsQueue.push(make_tuple(nI,tempScore,tmpScore,tempPath,tmpHitLen,tempextLen, tmpuniPos, 0));
		} else if((tmpNumOfBases == 0) && !check){
			//cout << "(tmpNumOfBases == 0) && !check" << endl;
			bestUnitigsQueue.push(make_tuple(nI,tempScore,tmpScore,tempPath,tmpHitLen,tempextLen, tmpuniPos, 1));
		} else if(check && tmpNumOfBases != 0){
			//cout << "check && tmpNumOfBases != 0" << endl;
			finishBases.push(make_tuple(nI,tempScore,tmpScore,tempPath,tmpHitLen,tempextLen,tmpuniPos,tmpNumOfBases));
		}
		if(tempScore > maxScore){
			maxScore = tempScore;
            bestPath = tempPath;
            hitLen = tmpHitLen;
		}
	}

	//cout << "zwischen" << endl;

	while(!(finishBases.empty())){
		auto tempFront = finishBases.front();
		finishBases.pop();

		shorterTemp currUnitig = get<0>(tempFront);
		uint32_t currScore = get<1>(tempFront);
		int32_t currtmpScore = get<2>(tempFront);
		pathList currPath = get<3>(tempFront);
		uint32_t currHitLen = get<4>(tempFront);
		uint32_t currextLen = get<5>(tempFront);
		uint32_t curruniPos = get<6>(tempFront);

		auto& tempcurrUnitig = *currUnitig;
		shorterContainer sucIter2 = tempcurrUnitig.getSuccessors();

		if((currextLen + iniQoff < q.length())){
			sucID = 0;
			for(shorterTemp nI = sucIter2.begin(); nI != sucIter2.end(); ++nI){
				++sucID;
				uint32_t tmpHitLen = currHitLen;
				int32_t tmpScore = currtmpScore;
                tempPath = currPath;
				uint32_t tmpExtLen = currextLen;
				uint32_t nextUniPos = curruniPos;
				int tmpNumOfBases = numOfBases;
				int32_t tempScore = contRightX_Drop_BFS_2(nI, iniQoff, tmpHitLen, tmpExtLen, q, mscore, mmscore, X, tmpScore, nextUniPos, tempPath, explCount, quorum, searchSet, advIdx, check, tmpNumOfBases);
				tempPath.push_back(sucID);
				if(check && tmpNumOfBases == 0){
					//cout << "check && tmpNumOfBases == 0" << endl;
					bestUnitigsQueue.push(make_tuple(nI,(currScore+tempScore),tmpScore,tempPath,tmpHitLen,tmpExtLen, nextUniPos, 0));
				} else if((tmpNumOfBases == 0) && !check){
					//cout << "(tmpNumOfBases == 0) && !check" << endl;
					bestUnitigsQueue.push(make_tuple(nI,tempScore,tmpScore,tempPath,tmpHitLen,tmpExtLen, nextUniPos, 1));
				} else if(check && tmpNumOfBases != 0){
					//cout << "check && tmpNumOfBases != 0" << endl;
					finishBases.push(make_tuple(nI,(currScore+tempScore),tmpScore,tempPath,tmpHitLen,tmpExtLen,nextUniPos,tmpNumOfBases));
				}
				int32_t scoreNow = currScore+tempScore;
				if (scoreNow > maxScore) {
                    	maxScore = scoreNow;
                    	bestPath = tempPath;
                    	hitLen = tmpHitLen;
						uniPos = nextUniPos;
				}
			}
		}
	}


	for(uint i = 0; i < numOfUnitig; i++) {

		if(bestUnitigsQueue.empty()){
			break;
		} else {
			bestUnitigs.push_back(bestUnitigsQueue.top());
			bestUnitigsQueue.pop(); 
		}

	}

	//cout << "zwische2" << endl;

	while(!(bestUnitigs.empty())){
		shorterVector2 allSucUnitigs;
		shorterPrioQueue2 allSucUnitigsQueue(prioLongest2);

		for(uint i = 0; i < bestUnitigs.size(); i++) {
			shorterTuple2 temptuple = bestUnitigs[i];
			shorterTemp currUnitig = std::get<0>(temptuple);
			uint32_t currScore = std::get<1>(temptuple);
			int32_t currtmpScore = std::get<2>(temptuple);
			pathList currPath = std::get<3>(temptuple);
			uint32_t currHitLen = std::get<4>(temptuple);
			uint32_t currextLen = std::get<5>(temptuple);
			uint32_t currUniPos = std::get<6>(temptuple);
			int currnumOfBases = get<7>(temptuple);

			if(currnumOfBases == 1){
				//cout << "Erweiterung markiert" << endl;
				int tmpNumOfBases = numOfBases;
				int32_t tempScore = contRightX_Drop_BFS_2(currUnitig, iniQoff, currHitLen, currextLen, q, mscore, mmscore, X, currtmpScore, currUniPos, currPath, explCount, quorum, searchSet, advIdx, check, tmpNumOfBases);
				if(check && tmpNumOfBases == 0){
					//cout << "check && tmpNumOfBases == 0" << endl;
					bestUnitigsQueue.push(make_tuple(currUnitig,(currScore+tempScore),currtmpScore,currPath,currHitLen,currextLen, currUniPos, 0));
				} else if((tmpNumOfBases == 0) && !check){
					//cout << "(tmpNumOfBases == 0) && !check" << endl;
					bestUnitigsQueue.push(make_tuple(currUnitig,tempScore,currtmpScore,currPath,currHitLen,currextLen, currUniPos, 1));
				} else if(check && tmpNumOfBases != 0){
					//cout << "check && tmpNumOfBases != 0" << endl;
					finishBases.push(make_tuple(currUnitig,(currScore+tempScore),currtmpScore,currPath,currHitLen,currextLen,currUniPos,tmpNumOfBases));
				}
				int32_t scoreNow = currScore+tempScore;
				if (scoreNow > maxScore) {
               		maxScore = scoreNow;
                	bestPath = currPath;
                	hitLen = currHitLen;
					uniPos = currUniPos;
				}
			} else {
				//cout << "Erweiterung nicht markiert" << endl;
				auto& tempcurrUnitig = *currUnitig;
				shorterContainer sucIter2 = tempcurrUnitig.getSuccessors();
				if((currextLen + iniQoff < q.length())) {
					sucID = 0;
					for(shorterTemp nI = sucIter2.begin(); nI != sucIter2.end(); ++nI){
						++sucID;
						uint32_t tmpHitLen = currHitLen;
						int32_t tmpScore = currtmpScore;
                		tempPath = currPath;
						uint32_t tmpExtLen = currextLen;
						uint32_t nextUniPos = currUniPos;
						int32_t tmpNumOfBases = numOfBases;
						int32_t tempScore = contRightX_Drop_BFS_2(nI, iniQoff, tmpHitLen, tmpExtLen, q, mscore, mmscore, X, tmpScore, nextUniPos, tempPath, explCount, quorum, searchSet, advIdx, check, tmpNumOfBases);
						tempPath.push_back(sucID);
						if(check && tmpNumOfBases == 0){
							//cout << "check && tmpNumOfBases == 0" << endl;
							bestUnitigsQueue.push(make_tuple(nI,(currScore+tempScore),tmpScore,tempPath,tmpHitLen,tmpExtLen, nextUniPos, 0));
						} else if((tmpNumOfBases == 0) && !check){
							//cout << "(tmpNumOfBases == 0) && !check" << endl;
							bestUnitigsQueue.push(make_tuple(nI,tempScore,tmpScore,tempPath,tmpHitLen,tmpExtLen, nextUniPos, 1));
						} else if(check && tmpNumOfBases != 0){
							//cout << "check && tmpNumOfBases != 0" << endl;
							finishBases.push(make_tuple(nI,(currScore+tempScore),tmpScore,tempPath,tmpHitLen,tmpExtLen,nextUniPos,tmpNumOfBases));
						}
						int32_t scoreNow = (currScore+tempScore);
						if(scoreNow > maxScore){
							maxScore = scoreNow;
                    		bestPath = tempPath;
                    		hitLen = tmpHitLen;
							uniPos = nextUniPos;
						}
					}
				}
			}
		}

		//cout << "zwischen3" << endl;

		while(!(finishBases.empty())){
			auto tempFront = finishBases.front();
			finishBases.pop();

			shorterTemp currUnitig = get<0>(tempFront);
			uint32_t currScore = get<1>(tempFront);
			int32_t currtmpScore = get<2>(tempFront);
			pathList currPath = get<3>(tempFront);
			uint32_t currHitLen = get<4>(tempFront);
			uint32_t currextLen = get<5>(tempFront);
			uint32_t curruniPos = get<6>(tempFront);

			auto& tempcurrUnitig = *currUnitig;
			shorterContainer sucIter2 = tempcurrUnitig.getSuccessors();

			if((currextLen + iniQoff < q.length())){
				sucID = 0;
				for(shorterTemp nI = sucIter2.begin(); nI != sucIter2.end(); ++nI){
					++sucID;
					uint32_t tmpHitLen = currHitLen;
					int32_t tmpScore = currtmpScore;
                	tempPath = currPath;
					uint32_t tmpExtLen = currextLen;
					uint32_t nextUniPos = curruniPos;
					int tmpNumOfBases = numOfBases;
					int32_t tempScore = contRightX_Drop_BFS_2(nI, iniQoff, tmpHitLen, tmpExtLen, q, mscore, mmscore, X, tmpScore, nextUniPos, tempPath, explCount, quorum, searchSet, advIdx, check, tmpNumOfBases);
						//cout << tmpNumOfBases << endl;
					tempPath.push_back(sucID);
					if(check && tmpNumOfBases == 0){
						//cout << "check && tmpNumOfBases == 0" << endl;
						bestUnitigsQueue.push(make_tuple(nI,(currScore+tempScore),tmpScore,tempPath,tmpHitLen,tmpExtLen, nextUniPos, 0));
					} else if((tmpNumOfBases == 0) && !check){
						//cout << "(tmpNumOfBases == 0) && !check" << endl;
						bestUnitigsQueue.push(make_tuple(nI,tempScore,tmpScore,tempPath,tmpHitLen,tmpExtLen, nextUniPos, 0));
					} else if(check && tmpNumOfBases != 0){
						//cout << "check && tmpNumOfBases != 0" << endl;
						finishBases.push(make_tuple(nI,(currScore+tempScore),tmpScore,tempPath,tmpHitLen,tmpExtLen,nextUniPos,tmpNumOfBases));
					}
					int32_t scoreNow = currScore+tempScore;
					if (scoreNow > maxScore) {
                    	maxScore = scoreNow;
                    	bestPath = tempPath;
                    	hitLen = tmpHitLen;
						uniPos = nextUniPos;
					}
				}
			}

		}

		bestUnitigs.clear();
		for(uint i = 0; i < numOfUnitig; i++) {
			if(allSucUnitigsQueue.empty()){
				break;
			} else{
				bestUnitigs.push_back(allSucUnitigsQueue.top());
				allSucUnitigsQueue.pop();
			}
		}
	}
	//cout << "while 3" << endl;
	extPth = bestPath;
	return maxScore;
}



//This function  initiates the extension on all successors of a unitig and returns the best one considering a quorum and a search color set. This function is explicitly designed for seeds lying on the query's reverse complement (considering the overlap between unitigs in sequences' beginning)
int32_t extendAtNextUnitig_OnRevComp(const ForwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	uint16_t sucID;
	uint32_t tmpHitLen = 0;
	int32_t maxScore = 0, currScore;
	
	//Check whether we have reached the maximum recursion depth of an extension
	if(++explCount > MAXRECURSIONDEPTH) return 0;//Terminate this extension

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

//The good old X-drop algorithm (extension to the right) for seeds matching the query's reference strand considering quorum and search color set. Returns an extension pointer storing the extension path through the graph
void startRightX_Drop(Hit* hit, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx, const int16_t extend_modus){
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

	//cout << "startUnitig :" << uSeq << endl;

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
			//cout << "check succesor" << endl;
			//Check if the current unitig has successors
			//cout << "startScore: " << hit->score << endl;
			if(overlap != 0){
				//Calculate the unitig sequence position we have to start with in the successive unitig
				iniUniPos = hit->offU + tmpSeedLen - hit->origUni.size + overlap;
				//Initialize explCount
				explCount = 0;
				//Explore unitig's successors
				switch(extend_modus) {
					case 0:
						hit->score += extendAtNextUnitig(hit->origUni.getSuccessors(), hit->offQ, hit->length, tmpSeedLen, q, mscore, mmscore, X, tmpScore, iniUniPos, extPth, explCount, quorum, searchSet, advIdx);
						break;
					case 1:
						hit->score += extendAtNextUnitig_BFS(hit->origUni.getSuccessors(), hit->offQ, hit->length, tmpSeedLen, q, mscore, mmscore, X, tmpScore, iniUniPos, extPth, explCount, quorum, searchSet, advIdx);
						break;
					case 2:
						hit->score += extendAtNextUnitig_BFS_SMART1(hit->origUni.getSuccessors(), hit->offQ, hit->length, tmpSeedLen, q, mscore, mmscore, X, tmpScore, iniUniPos, extPth, explCount, quorum, searchSet, advIdx);
						break;
					case 3:
						hit->score += extendAtNextUnitig_BFS_SMART2(hit->origUni.getSuccessors(), hit->offQ, hit->length, tmpSeedLen, q, mscore, mmscore, X, tmpScore, iniUniPos, extPth, explCount, quorum, searchSet, advIdx);
						break;
					case 4:
						hit->score += extendAtNextUnitig_BFS_SMART3(hit->origUni.getSuccessors(), hit->offQ, hit->length, tmpSeedLen, q, mscore, mmscore, X, tmpScore, iniUniPos, extPth, explCount, quorum, searchSet, advIdx);
						break;
				}
				
			break;
			}
		}
	}

	//cout << "startRightX_drop score after" << endl;
	//cout << hit->score << endl;
	for(list<uint16_t>::iterator i = extPth.begin(); i != extPth.end(); i++){
		//cout << *i << endl;
	}
	//cout << endl;

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

	//Save the initial offset in the current unitig which we need for all nearest neighbor calculations
	iniSeqPos = uniSeqPos;
	//Find the nearest seed that we might be able to reach during our extension
	nearestSeed = searchRightNeighbor(sucUnitig->getData()->getData(*sucUnitig)->getSeed(sucUnitig->strand), iniQoff, extLen, iniSeqPos, prevSeed);
	//Perform the X-drop algorithm on the successive unitig
	tmpScore = lastSeedTmpScore;
	tmpSLen = 0;

	//We are done if we have reached the end of the query
	while(iniQoff + extLen + tmpSLen < q.length()){
		//Check whether we have reached the next seed
		//Testing
		//cout << "Not at query's end" << endl;
		//cout << "nearestSeed is " << (nearestSeed == NULL ? "NULL" : "not NULL") << endl;
		if(nearestSeed != NULL && iniQoff + extLen + tmpSLen >= nearestSeed->offsetQ){
			//Calculate the gain we get by incorporating the reached seed

			//Testing
			//cout << "Found seed" << endl;

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
			//Check up to which point we have to compare the unitig sequence
			if(!sucUnitig->getSuccessors().hasSuccessors()){
				overlap = 0;
			}

			//Check whether we have reached the end of the unitig's sequence
			if(uniSeqPos < sucUnitig->size - overlap){
				//Are search criteria still fulfilled?
				if(checkedPos <= 0) break;

				//TEsting
				//cout << "Compare bases" << endl;

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
				//cout << "tmpScore: " << tmpScore << endl;

				//cout << "sucUniSeq[uniSeqPos]: " << sucUniSeq[uniSeqPos] << "  q[iniQoff + extLen + tmpSLen]: " << q[iniQoff + extLen + tmpSLen] << "  score: " << score << endl;

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
					//cout << " " << endl;
					//cout << "next unitig" << endl;
					int32_t tempScore = extendAtNextUnitig(sucUnitig->getSuccessors(), iniQoff, hitLen, extLen + tmpSLen, q, mscore, mmscore, X, tmpScore, uniSeqPos, extPth, explCount, quorum, searchSet, advIdx);
					score += tempScore;
				}

				break;
			}
		}
	}
	
	return score;
}


int32_t contRightX_Drop_BFS(const neighborIterator<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> &sucUnitig, const uint32_t &iniQoff, uint32_t &hitLen, uint32_t &extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, int32_t &tmpScore, uint32_t &uniSeqPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx, bool &check){
	int32_t progress, score = 0;
	int32_t overlap = sucUnitig->getGraph()->getK() - 1;
	uint32_t iniSeqPos;
	uint32_t tmpSLen;
	//Even if we are on the reverse complementary strand and no position is covered on this unitig, because we have checked the first k - 1 position on the last unitig already
	int32_t checkedPos = getSrchCritCov(*sucUnitig, quorum, searchSet, compOffset(uniSeqPos, 1, sucUnitig->size, sucUnitig->strand), sucUnitig->strand, advIdx);
	string sucUniSeq = sucUnitig->mappedSequenceToString();
	struct Seed *nearestSeed, *prevSeed;

	//cout << "starttmpScore: " << tmpScore << endl;

	//Save the initial offset in the current unitig which we need for all nearest neighbor calculations
	iniSeqPos = uniSeqPos;
	//Find the nearest seed that we might be able to reach during our extension
	nearestSeed = searchRightNeighbor(sucUnitig->getData()->getData(*sucUnitig)->getSeed(sucUnitig->strand), iniQoff, extLen, iniSeqPos, prevSeed);
	//Perform the X-drop algorithm on the successive unitig
	tmpSLen = 0;
	check = false;

	//We are done if we have reached the end of the query
	while(iniQoff + extLen + tmpSLen < q.length()){
		//Testing
		//cout << "Not at query's end" << endl;
		//cout << "nearestSeed is " << (nearestSeed == NULL ? "NULL" : "not NULL") << endl;

		//Check whether we have reached the next seed
		if(nearestSeed != NULL && iniQoff + extLen + tmpSLen >= nearestSeed->offsetQ){
			//Testing
			//cout << "Found seed" << endl;

			/*
			if(nearestSeed->offsetQ == 5017 && nearestSeed->offsetU == 31){
				cout << "treffer hat score 1 " << endl;
				cout << sucUnitig->mappedSequenceToString() << endl;
				cout << "iniQoff: " << iniQoff << endl;
				//exit(0);
			}
			*/

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
			//Check up to which point we have to compare the unitig sequence
			if(!sucUnitig->getSuccessors().hasSuccessors()){
				overlap = 0;
			}

			//Check whether we have reached the end of the unitig's sequence
			if(uniSeqPos < sucUnitig->size - overlap){
				//Are search criteria still fulfilled?
				if(checkedPos <= 0){
					break;
				}

				//TEsting
				//cout << "Compare bases" << endl;

				//cout << "sucUniSeq[uniSeqPos]: " << sucUniSeq[uniSeqPos] << "  q[iniQoff + extLen + tmpSLen]: " << q[iniQoff + extLen + tmpSLen] << endl;

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
					if(tmpScore < -X){
						break;
					}
				}
				//cout << "tmpScore: " << tmpScore << endl;
				//cout << "sucUniSeq[uniSeqPos]: " << sucUniSeq[uniSeqPos] << "  q[iniQoff + extLen + tmpSLen]: " << q[iniQoff + extLen + tmpSLen] << "  score: " << score << endl;
				//cout << "iniQoff + extLen + tmpSLen: " << iniQoff + extLen + tmpSLen << ", uniSeqPos: " << uniSeqPos << endl;

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
					extLen = extLen + tmpSLen;
					//cout << "tmpSLen: " << tmpSLen << endl;
					check = true;
				}

				break;
			}
		}
	}

	//cout << "iniQoff: " << iniQoff << ", uniPos: " << uniSeqPos << ", extLen: " << extLen << ", tmpSLen: " << tmpSLen << ", iniQoff + extLen + tmpSLen: " << iniQoff + extLen + tmpSLen << endl;

	//std::cout << "contRightX_Drop_BFS Final score: " << score << std::endl;
	//cout << "Unitig: " << sucUnitig->mappedSequenceToString() << ", hitLen: " << hitLen << ", score: " << score << ", uniPos: " << uniSeqPos << endl;
	//cout << "final uniSeqPos = " << uniSeqPos << std::endl;
	//cout << "finaltmpScore: " << tmpScore << endl;
	return score;
}

int32_t contRightX_Drop_BFS_2(const neighborIterator<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> &sucUnitig, const uint32_t &iniQoff, uint32_t &hitLen, uint32_t &extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, int32_t &tmpScore, uint32_t &uniSeqPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx, bool &check, int &numOfBases){
	int32_t progress, score = 0;
	int32_t overlap = sucUnitig->getGraph()->getK() - 1;
	uint32_t iniSeqPos;
	uint32_t tmpSLen;
	//Even if we are on the reverse complementary strand and no position is covered on this unitig, because we have checked the first k - 1 position on the last unitig already
	int32_t checkedPos = getSrchCritCov(*sucUnitig, quorum, searchSet, compOffset(uniSeqPos, 1, sucUnitig->size, sucUnitig->strand), sucUnitig->strand, advIdx);
	string sucUniSeq = sucUnitig->mappedSequenceToString();
	struct Seed *nearestSeed, *prevSeed;

	//cout << "starttmpScore: " << tmpScore << endl;

	//Save the initial offset in the current unitig which we need for all nearest neighbor calculations
	iniSeqPos = uniSeqPos;
	//Find the nearest seed that we might be able to reach during our extension
	nearestSeed = searchRightNeighbor(sucUnitig->getData()->getData(*sucUnitig)->getSeed(sucUnitig->strand), iniQoff, extLen, iniSeqPos, prevSeed);
	//Perform the X-drop algorithm on the successive unitig
	tmpSLen = 0;
	check = false;

	//We are done if we have reached the end of the query
	while(iniQoff + extLen + tmpSLen < q.length()){
		//Testing
		//cout << "Not at query's end" << endl;
		//cout << "nearestSeed is " << (nearestSeed == NULL ? "NULL" : "not NULL") << endl;

		if(numOfBases == 0){
			break;
		}


		//Check whether we have reached the next seed
		if(nearestSeed != NULL && iniQoff + extLen + tmpSLen >= nearestSeed->offsetQ){
			//Testing
			//cout << "Found seed" << endl;

			//Calculate the gain we get by incorporating the reached seed
			progress = nearestSeed->offsetQ + nearestSeed->len - (iniQoff + extLen + tmpSLen);

			//cout << "before Seed calc: " << numOfBases << endl;
			if(progress >= numOfBases){
				progress = numOfBases;
				numOfBases = 0;
			}else{
				numOfBases -= progress;
			}
			//cout << "after Seed calc: " << numOfBases << endl;

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
			//Check up to which point we have to compare the unitig sequence
			if(!sucUnitig->getSuccessors().hasSuccessors()){
				overlap = 0;
			}

			//Check whether we have reached the end of the unitig's sequence
			if(uniSeqPos < sucUnitig->size - overlap){
				//Are search criteria still fulfilled?
				if(checkedPos <= 0){
					break;
				}

				//TEsting
				//cout << "Compare bases" << endl;

				//cout << "sucUniSeq[uniSeqPos]: " << sucUniSeq[uniSeqPos] << "  q[iniQoff + extLen + tmpSLen]: " << q[iniQoff + extLen + tmpSLen] << endl;

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
					if(tmpScore < -X){
						break;
					}
				}
				//cout << "tmpScore: " << tmpScore << endl;
				//cout << "sucUniSeq[uniSeqPos]: " << sucUniSeq[uniSeqPos] << "  q[iniQoff + extLen + tmpSLen]: " << q[iniQoff + extLen + tmpSLen] << "  score: " << score << endl;
				//cout << "iniQoff + extLen + tmpSLen: " << iniQoff + extLen + tmpSLen << ", uniSeqPos: " << uniSeqPos << endl;

				//Proceed with the next two positions
				++tmpSLen;
				++uniSeqPos;
				--checkedPos;
				--numOfBases;
				//cout << "after Base calc: " << numOfBases << endl;
			} else{
				//Check if the current unitig has successors
				if(overlap != 0){
					//Calculate the position in the next unitig's sequence we have to start with
					uniSeqPos = uniSeqPos - sucUnitig->size + overlap;
					//Check out next unitig
					extLen = extLen + tmpSLen;
					//cout << "tmpSLen: " << tmpSLen << endl;
					check = true;
				}

				break;
			}
		}
	}

	//cout << "iniQoff: " << iniQoff << ", uniPos: " << uniSeqPos << ", extLen: " << extLen << ", tmpSLen: " << tmpSLen << ", iniQoff + extLen + tmpSLen: " << iniQoff + extLen + tmpSLen << endl;

	//std::cout << "contRightX_Drop_BFS Final score: " << score << std::endl;
	//cout << "Unitig: " << sucUnitig->mappedSequenceToString() << ", hitLen: " << hitLen << ", score: " << score << ", uniPos: " << uniSeqPos << endl;
	//cout << "final uniSeqPos = " << uniSeqPos << std::endl;
	//cout << "finaltmpScore: " << tmpScore << endl;
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
	if(++explCount > MAXRECURSIONDEPTH){
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
	if(++explCount > MAXRECURSIONDEPTH){
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