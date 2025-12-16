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

outputTypes2 stepUnitig(inputTypes2 inputStructure){
	queue<shorterTuple> extensionQueue		= inputStructure.extensionQueue;
	uint32_t iniQoff						= inputStructure.iniQoff;
	string q								= inputStructure.q;
	uint16_t mscore							= inputStructure.mscore;
	int16_t mmscore							= inputStructure.mmscore;
	int16_t X								= inputStructure.X;
	uint32_t explCount						= inputStructure.explCount;
	uint32_t quorum							= inputStructure.quorum;
	list<pair<string, size_t>> searchSet	= inputStructure.searchSet;
	bool advIdx								= inputStructure.advIdx;
	int32_t maxScore						= inputStructure.maxScore;
	int32_t numOfBases 						= inputStructure.numOfBases;
	uint32_t hitLen							= inputStructure.hitLen;
	list<uint16_t> bestPath					= inputStructure.bestPath;
	bool compareBases						= inputStructure.compareBases;

	shorterPrioQueue unitigsPrioQueue(prioLongest);

	unitigsPrioQueue						= inputStructure.priorityQueue;



	outputTypes2 p;

	p.unitigsPrioQueue = unitigsPrioQueue;
	p.extensionQueue = extensionQueue;
	p.maxScore = maxScore;
	p.hitLen = hitLen;
	p.bestPath = bestPath;

	return p;
}


//	does a BFS-step for all saved unitigs
//	and returns them in a priority queue
outputTypes calcUnitigs(inputTypes inputStructure) {

	//	get all information contained and given in the input structure
	queue<shorterTuple> extensionQueue		= inputStructure.extensionQueue;
	uint32_t iniQoff						= inputStructure.iniQoff;
	string q								= inputStructure.q;
	uint16_t mscore							= inputStructure.mscore;
	int16_t mmscore							= inputStructure.mmscore;
	int16_t X								= inputStructure.X;
	uint32_t explCount						= inputStructure.explCount;
	uint32_t quorum							= inputStructure.quorum;
	list<pair<string, size_t>> searchSet	= inputStructure.searchSet;
	bool advIdx								= inputStructure.advIdx;
	int32_t maxScore						= inputStructure.maxScore;
	int32_t numOfBases 						= inputStructure.numOfBases;
	uint32_t hitLen							= inputStructure.hitLen;
	list<uint16_t> bestPath					= inputStructure.bestPath;
	bool compareBases						= inputStructure.compareBases;

	//	create a priority queue sorted decreasingly by score for all finished extensions
	shorterPrioQueue unitigsPrioQueue(prioLongest);

	//	process all extensions
	while(!(extensionQueue.empty())){

		//	get next extension
		shorterTuple currExtension = extensionQueue.front();
		extensionQueue.pop();
		
		//	extract all relevant information of the extension
		UnitigColorMap<UnitigInfo> currUnitig 		= get<0>(currExtension);
		uint32_t currScore 			= get<1>(currExtension);
		int32_t currtmpScore 		= get<2>(currExtension);
		list<uint16_t> currPath		= get<3>(currExtension);
		uint32_t currHitLen 		= get<4>(currExtension);
		uint32_t currextLen 		= get<5>(currExtension);
		uint32_t curruniPos 		= get<6>(currExtension);
		int currnumOfBases 			= get<7>(currExtension);
		int32_t tempNumOfBases 		= numOfBases;

		//	get the successors of the unitig
		shorterContainer sucIter2 = currUnitig.getSuccessors();

		//	stop if we reached the end of the query
		if((currextLen + iniQoff < q.length())){
			uint16_t sucID = 0;

			if(!compareBases){
				//	iterate over all successors
				for(shorterTemp nI = sucIter2.begin(); nI != sucIter2.end(); ++nI){
					++sucID;

					//	copy current extension for multiple iterations
					//	because the variables will be overwritten when calculation the extension of the first successor
					bool check = false;
					uint32_t nextHitLen 	= currHitLen;
					int32_t nextScore 		= currtmpScore;
        			list<uint16_t>tempPath 	= currPath;
					uint32_t nextExtLen 	= currextLen;
					uint32_t nextUniPos 	= curruniPos;

					//	calculate the score of an extension of a successor and add it to its current score
					int32_t addScore = contRightX_Drop_BFS(*nI, iniQoff, nextHitLen, nextExtLen, q, mscore, mmscore, X, nextScore, nextUniPos, tempPath, explCount, quorum, searchSet, advIdx, check, tempNumOfBases, false);
					tempPath.push_back(sucID);
					int32_t fullScore = currScore+addScore;

					//	if we reached the end of the unitig add extension to the priority queue
					if(check){
						unitigsPrioQueue.push(make_tuple(*nI,fullScore,nextScore,tempPath,nextHitLen,nextExtLen, nextUniPos, 0));
					}

					//	update the best scoring extension
					if (fullScore > maxScore) {
          				maxScore = fullScore;
          				hitLen = nextHitLen;
						bestPath = tempPath;
					}
				}
			} else {
				//	if the extension is marked we don't lock at successors
				if(currnumOfBases == -1){
					bool check = false;
					list<uint16_t>tempPath 	= currPath;

					//	continue extension on the current unitig and calculate the new score
					int32_t addScore = contRightX_Drop_BFS(currUnitig, iniQoff, currHitLen, currextLen, q, mscore, mmscore, X, currtmpScore, curruniPos, tempPath, explCount, quorum, searchSet, advIdx, check, tempNumOfBases, true);
					int32_t fullScore = currScore+addScore;

					//	if having compared enough bases while still not reaching the end of the unitig
					//	add marked extension to the priority queue
					if((tempNumOfBases == 0) && !check){
						unitigsPrioQueue.push(make_tuple(currUnitig,fullScore,currtmpScore,tempPath,currHitLen,currextLen, curruniPos, -1));

					//	if not having compared enough bases and we reached the end of the unitig
					//	the extension is added back to the queue unmarked
					} else if(check && tempNumOfBases != 0){
						extensionQueue.push(make_tuple(currUnitig,fullScore,currtmpScore,tempPath,currHitLen,currextLen,curruniPos,tempNumOfBases));
					} else {
						//cout << "removed" << endl;
					}

					//	update the best scoring extension
					if (fullScore > maxScore) {
                    	maxScore = fullScore;
                    	hitLen = currHitLen;
						bestPath = tempPath;
					}

				//	if the extension is not marked go through the succesors of the unitig	
				} else {
					uint16_t sucID = 0;

					//	iterate over all successors
					for(shorterTemp nI = sucIter2.begin(); nI != sucIter2.end(); ++nI){
						++sucID;

						//	copy current extension for multiple iterations
						//	because the variables will be overwritten when calculation the extension of the first successor
						bool check = false;
						uint32_t nextHitLen 	= currHitLen;
						int32_t nextScore 		= currtmpScore;
               			list<uint16_t>tempPath 	= currPath;
						uint32_t nextExtLen 	= currextLen;
						uint32_t nextUniPos 	= curruniPos;
						int32_t tempNumOfBases 	= numOfBases;

						//	calculate the score of an extension of a successor and add it to its current score
						int32_t addScore = contRightX_Drop_BFS(*nI, iniQoff, nextHitLen, nextExtLen, q, mscore, mmscore, X, nextScore, nextUniPos, tempPath, explCount, quorum, searchSet, advIdx, check, tempNumOfBases, true);
						tempPath.push_back(sucID);
						int32_t fullScore = currScore+addScore;

						//	if having compared enough bases while still not reaching the end of the unitig
						//	add marked extension to the priority queue
						if((tempNumOfBases == 0) && !check){
							unitigsPrioQueue.push(make_tuple(*nI,fullScore,nextScore,tempPath,nextHitLen,nextExtLen, nextUniPos, -1));

						//	if not having compared enough bases and we reached the end of the unitig
						//	the extension is added back to the queue unmarked
						} else if(check && tempNumOfBases != 0){
							extensionQueue.push(make_tuple(*nI,fullScore,nextScore,tempPath,nextHitLen,nextExtLen,nextUniPos,tempNumOfBases));
						} else {
							//cout << "removed" << endl;
						}

						//	update the best scoring extension
						if (fullScore > maxScore) {
                    		maxScore = fullScore;
                    		hitLen = nextHitLen;
							bestPath = tempPath;
						}
					}
				}
			}
		}
	}

	//	return the best extension and the priority queue
	outputTypes p;

	p.bestUnitigsPrioQueue = unitigsPrioQueue;
	p.maxScore = maxScore;
	p.hitLen = hitLen;
	p.bestPath = bestPath;

	return p;
}

//	this function initiates the extension on all successors of a unitig and returns the best one considering a quorum and a search color set
//	using a exhaustive BFS algorithm
int32_t extendAtNextUnitig_BFS(const UnitigColorMap<UnitigInfo> startUnitig, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &uniPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	
	//	best score
	int32_t maxScore = 0;

	//	temporary and best global path
	list<uint16_t> tempPath;
	list<uint16_t> bestPath;

	//	priority queue ordered decreasingly by score
	shorterPrioQueue unitigsPrioQueue(prioLongest);

	//	queue for next iterations extensions
	queue<shorterTuple> unitigsQueue;

	//	add starting extension to the queue
	shorterTuple startTuple = make_tuple(startUnitig, 0, lastExtSeedTmpScore, tempPath, hitLen, extLen, uniPos, 0);
	unitigsQueue.push(startTuple);

	//	prepare the parameters given to the BFS extension function
	inputTypes inputStruct;
	inputStruct.iniQoff     	= iniQoff;
	inputStruct.q           	= q;
	inputStruct.mscore      	= mscore;
	inputStruct.mmscore     	= mmscore;
	inputStruct.X           	= X;
	inputStruct.explCount   	= explCount;
	inputStruct.quorum      	= quorum;
	inputStruct.searchSet   	= searchSet;
	inputStruct.advIdx      	= advIdx;
	inputStruct.numOfBases  	= 0;
	inputStruct.compareBases	= false;

	outputTypes outputStruct;

	//	repeat BFS extension until no more extensions are possible
	while(!(unitigsQueue.empty())){
		//	update extension queue and best result found so far
		inputStruct.extensionQueue = unitigsQueue;
		inputStruct.maxScore = maxScore;
		inputStruct.hitLen = hitLen;
		inputStruct.bestPath = bestPath;

		//	extend and score
		outputStruct = calcUnitigs(inputStruct);

		//	update with new results
		unitigsPrioQueue = outputStruct.bestUnitigsPrioQueue;
		maxScore = outputStruct.maxScore;
		hitLen = outputStruct.hitLen;
		bestPath = outputStruct.bestPath;

		//	empty the queue
		//while (!unitigsQueue.empty()) {
      		//unitigsQueue.pop();
    	//}

		unitigsQueue = queue<shorterTuple>();

		//	transfer element from priority queue to normal queue
		while (!unitigsPrioQueue.empty()) {
    		unitigsQueue.push(unitigsPrioQueue.top());
    		unitigsPrioQueue.pop();
		}
	}

	//	update path of best extension and return best score
	extPth = bestPath;
	return maxScore;
}

	
/*
	//	best extension overall
	list<uint16_t> bestPath;

	//	initial expansion from the given succesors
	for(shorterTemp nI = sucIter.begin(); nI != sucIter.end(); ++nI){
		++sucID;

		//	copy for the extension important parameters
		uint32_t tempextLen = extLen;
		uint32_t tmpHitLen 	= hitLen;
		uint32_t tmpuniPos 	= uniPos;
		int32_t tmpScore 	= lastExtSeedTmpScore;
        tempPath.clear();

		//	calculate the score of an extension of a successor
		int numOfBases = 0;
		int32_t startScore = contRightX_Drop_BFS(*nI, iniQoff, tmpHitLen, tempextLen, q, mscore, mmscore, X, tmpScore, tmpuniPos, tempPath, explCount, quorum, searchSet, advIdx, check, numOfBases, false);
		tempPath.push_back(sucID);

		//	if extensions is valid it is added for further BFS processing
		if(check){
			queueTest.push(std::make_tuple(nI,startScore,tmpScore,tempPath,tmpHitLen,tempextLen, tmpuniPos));
		}

		//	tracks best extension
		if (startScore > maxScore) {
            maxScore = startScore;
            bestPath = tempPath;
			hitLen = tmpHitLen;
        }
	}

	//	process all extensions
	while (!queueTest.empty()) {

		//	get next extension
		auto tempFront = queueTest.front();
		queueTest.pop();

		//	get current extension parameters
		shorterTemp currUnitig 	= get<0>(tempFront);
		uint32_t currScore 		= get<1>(tempFront);
		int32_t currtmpScore 	= get<2>(tempFront);
		pathList currPath 		= get<3>(tempFront);
		uint32_t currHitLen 	= get<4>(tempFront);
		uint32_t currextLen 	= get<5>(tempFront);
		uint32_t curruniPos 	= get<6>(tempFront);

		//	get the succesors of the current unitig
		UnitigColorMap<UnitigInfo> tempcurrUnitig = *currUnitig;
		shorterContainer sucIter2 = tempcurrUnitig.getSuccessors();

		//	continue extension only if it does not exceed query length
		if((currextLen + iniQoff < q.length())) {
			sucID = 0;

			//	try extending every succesor
			for(shorterTemp nI = sucIter2.begin(); nI != sucIter2.end(); ++nI){
				check = true;
				++sucID;

				//	copy current extension information
				uint32_t tmpHitLen 	= currHitLen;
				int32_t tmpScore 	= currtmpScore;
                tempPath 			= currPath;
				uint32_t tmpExtLen 	= currextLen;
				uint32_t nextUniPos = curruniPos;

				//	calculate the score of an extension of a successor
				int numOfBases = 0;
				int32_t tempScore = contRightX_Drop_BFS(*nI, iniQoff, tmpHitLen, tmpExtLen, q, mscore, mmscore, X, tmpScore, nextUniPos, tempPath, explCount, quorum, searchSet, advIdx, check, numOfBases, false);
				int32_t scoreNow = currScore+tempScore;
				tempPath.push_back(sucID);

				//	if extensions is valid it is added for further BFS processing
				if(check == true){
					queueTest.push(std::make_tuple(nI, (scoreNow), tmpScore, tempPath, tmpHitLen, tmpExtLen, nextUniPos));
                }

				//	tracks best extension
				if (scoreNow > maxScore) {
                    maxScore = scoreNow;
                    bestPath = tempPath;
                    hitLen = tmpHitLen;
					uniPos = nextUniPos;
				}
			}
		}
	}

	// update the best path and return the max score
	extPth = bestPath;
	return maxScore;
}

*/

//	takes up to numOfAcceptedExtensions elements from the given priority queue
//	and returns them as a queue
queue<shorterTuple> getBestUnitigs(shorterPrioQueue bestUnitigsPrioQueue, uint numOfAcceptedExtensions){
	queue<shorterTuple> bestUnitigs;
	
	for(uint i = 0; i < numOfAcceptedExtensions; i++) {
		if(bestUnitigsPrioQueue.empty()){
			break;
		} else {
			bestUnitigs.push(bestUnitigsPrioQueue.top());
			bestUnitigsPrioQueue.pop(); 
		}
	}

	return bestUnitigs;
}


//	this function initiates the extension on all successors of a unitig and returns the best one considering a quorum and a search color set
//	using a BFS heuristic where each iteration a number of best scoring extensions are pushed to the next iteration 
int32_t extendAtNextUnitig_BFS_SMART1(const UnitigColorMap<UnitigInfo> startUnitig, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &uniPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	
	//	best score
	int32_t maxScore = 0;

	//	max number of unitigs pushed to the next iteration
	uint numOfUnitig = 1000000;

	//	priority queue ordered decreasingly by score
	shorterPrioQueue unitigsPrioQueue(prioLongest);

	//	queue for next iterations extensions
	queue<shorterTuple> bestUnitigs;

	//	temporary and best global path
	list<uint16_t> tempPath;
	list<uint16_t> bestPath;

	//	add starting extension to the queue
	shorterTuple startTuple = make_tuple(startUnitig, 0, lastExtSeedTmpScore, tempPath, hitLen, extLen, uniPos, 0);
	bestUnitigs.push(startTuple);

	//	prepare the parameters given to the BFS extension function
	inputTypes inputStruct;
	inputStruct.iniQoff 		= iniQoff;
	inputStruct.q 				= q;
	inputStruct.mscore 			= mscore;
	inputStruct.mmscore 		= mmscore;
	inputStruct.X 				= X;
	inputStruct.explCount 		= explCount;
	inputStruct.quorum 			= quorum;
	inputStruct.searchSet 		= searchSet;
	inputStruct.advIdx 			= advIdx;
	inputStruct.numOfBases 		= 0;
	inputStruct.compareBases	= false;

	outputTypes outputStruct;

	//	repeat BFS extension until no more extensions are possible
	while(!(bestUnitigs.empty())){

		//	update extension queue and best result found so far
		inputStruct.extensionQueue = bestUnitigs;
		inputStruct.maxScore = maxScore;
		inputStruct.hitLen = hitLen;
		inputStruct.bestPath = bestPath;

		//	extend and score
		outputStruct = calcUnitigs(inputStruct);

		//	update with new results
		unitigsPrioQueue = outputStruct.bestUnitigsPrioQueue;
		maxScore = outputStruct.maxScore;
		hitLen = outputStruct.hitLen;
		bestPath = outputStruct.bestPath;

		//	choose extensions for the next iteration
		bestUnitigs = getBestUnitigs(unitigsPrioQueue,numOfUnitig);
	}

	//	update path of best extension and return best score
	extPth = bestPath;
	return maxScore;
}

/*
	//	how many extensions are kept
	uint numOfUnitig = 1000000;

	//	priority queue sorted by score (highes-lowest)
	shorterPrioQueue bestUnitigsQueue(prioLongest);

	//	vector of the chosen extensions for the enxt iteration
	shorterVector bestUnitigs;

	//	temporary path and best global path
	list<uint16_t> tempPath;
	list<uint16_t> bestPath;

	//	initial expansion from the given succesors
	for(shorterTemp nI = sucIter.begin(); nI != sucIter.end(); ++nI){
		++sucID;

		//	copy for the extension important parameters
		uint32_t tempextLen = extLen;
		uint32_t tmpHitLen 	= hitLen;
		uint32_t tmpuniPos 	= uniPos;
		int32_t tmpScore 	= lastExtSeedTmpScore;
		tempPath.clear();

		//	calculate the score of an extension of a successor
		int numOfBases = 0;
		int32_t tempScore = contRightX_Drop_BFS(*nI, iniQoff, tmpHitLen, tempextLen, q, mscore, mmscore, X, tmpScore, tmpuniPos, tempPath, explCount, quorum, searchSet, advIdx, check, numOfBases, false);
		tempPath.push_back(sucID);

		//	if extensions is valid it is added for further BFS processing
		if(check){
			bestUnitigsQueue.push(make_tuple(nI,tempScore,tmpScore,tempPath,tmpHitLen,tempextLen, tmpuniPos));
		}

		//	tracks best extension
		if(tempScore > maxScore){
			maxScore = tempScore;
            bestPath = tempPath;
            hitLen = tmpHitLen;
		}
	}

	//	takes up to numOfUnitig amount of best scoring unitigs
	for(uint i = 0; i < numOfUnitig; i++) {

		if(bestUnitigsQueue.empty()){
			break;
		} else {
			bestUnitigs.push_back(bestUnitigsQueue.top());
			bestUnitigsQueue.pop(); 
		}

	}

	//	process the best extensions
	while(!(bestUnitigs.empty())){
		shorterPrioQueue allSucUnitigsQueue(prioLongest);

		//	go through all extensions
		for(uint i = 0; i < bestUnitigs.size(); i++) {

			//	get current extension parameters
			shorterTuple temptuple 	= bestUnitigs[i];
			shorterTemp currUnitig 	= get<0>(temptuple);
			uint32_t currScore 		= get<1>(temptuple);
			int32_t currtmpScore 	= get<2>(temptuple);
			pathList currPath 		= get<3>(temptuple);
			uint32_t currHitLen 	= get<4>(temptuple);
			uint32_t currextLen 	= get<5>(temptuple);
			uint32_t currUniPos 	= get<6>(temptuple);

			//	get the succesors of the current unitig
			UnitigColorMap<UnitigInfo> tempcurrUnitig = *currUnitig;
			shorterContainer sucIter2 = tempcurrUnitig.getSuccessors();

			//	continue extension only if it does not exceed query length
			if((currextLen + iniQoff < q.length())) {
				sucID = 0;

				//	try extending every succesor
				for(shorterTemp nI = sucIter2.begin(); nI != sucIter2.end(); ++nI){
					++sucID;

					//	copy current extension information
					uint32_t tmpHitLen 	= currHitLen;
					int32_t tmpScore 	= currtmpScore;
                	tempPath 			= currPath;
					uint32_t tmpExtLen 	= currextLen;
					uint32_t nextUniPos = currUniPos;

					//	calculate the score of an extension of a successor
					int numOfBases = 0;
					int32_t tempScore = contRightX_Drop_BFS(*nI, iniQoff, tmpHitLen, tmpExtLen, q, mscore, mmscore, X, tmpScore, nextUniPos, tempPath, explCount, quorum, searchSet, advIdx, check, numOfBases, false);
					int32_t scoreNow = currScore+tempScore;
					tempPath.push_back(sucID);

					//	if extensions is valid it is added for further BFS processing
					if(check){
						allSucUnitigsQueue.push(make_tuple(nI, (scoreNow), tmpScore, tempPath, tmpHitLen, tmpExtLen, nextUniPos));
					}

					//	tracks best extension
					if(scoreNow > maxScore){
						maxScore = currScore+tempScore;
                    	bestPath = tempPath;
                    	hitLen = tmpHitLen;
						uniPos = nextUniPos;
					}
				}
			}
		}
		
		//	takes up to numOfUnitig amount of best scoring unitigs
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

	// update the best path and return the max score
	extPth = bestPath;
	return maxScore;
}
*/

//	this function initiates the extension on all successors of a unitig and returns the best one considering a quorum and a search color set
//	using a BFS heuristic where at all times only a certain number of the best scoring extensions are kept
int32_t extendAtNextUnitig_BFS_SMART2(const UnitigColorMap<UnitigInfo> startUnitig, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &uniPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	
	//	best score
	int32_t maxScore = 0;

	//	max number of unitigs pushed to the next iteration
	uint numOfUnitig = 1000000;

	//	temporary and best global path
	list<uint16_t> tempPath;
	list<uint16_t> bestPath;

	//	list of the unitigs
	shorterVector bestUnitigs;

	//	list of the scores of the unitigs
	vector<int> bestScores;

	//	add starting extension to the starting list
	shorterTuple startTuple = make_tuple(startUnitig, 0, lastExtSeedTmpScore, tempPath, hitLen, extLen, uniPos, 0);
	bestUnitigs.push_back(startTuple);
	bestScores.push_back(0);

	//	repeat BFS extension until no more extensions are possible
	while(!(bestUnitigs.empty())){

		//	go through every extension of the list
		for(uint i = 0; i < bestUnitigs.size(); i++) {

			//	extract all relevant information of the extension
			shorterTuple currExtension 	= bestUnitigs[0];
			UnitigColorMap<UnitigInfo> currUnitig 	= get<0>(currExtension);
			uint32_t currScore 		= get<1>(currExtension);
			int32_t currtmpScore 	= get<2>(currExtension);
			pathList currPath 		= get<3>(currExtension);
			uint32_t currHitLen 	= get<4>(currExtension);
			uint32_t currextLen 	= get<5>(currExtension);
			uint32_t curruniPos 	= get<6>(currExtension);

			//	get the successors of the unitig
			shorterContainer sucIter2 = currUnitig.getSuccessors();

			//	remove the unitig from the list
			bestUnitigs.erase(bestUnitigs.begin() + 0);
			bestScores.erase(bestScores.begin() + 0);
			
			//	stop if we reached the end of the query
			if((currextLen + iniQoff < q.length())){
				uint16_t sucID = 0;

				//	iterate over all successors
				for(shorterTemp nI = sucIter2.begin(); nI != sucIter2.end(); ++nI){
					++sucID;

					//	copy current extension for multiple iterations
					//	because the variables will be overwritten when calculation the extension of the first successor
					bool check = false;
					uint32_t nextHitLen 	= currHitLen;
					int32_t nextScore 		= currtmpScore;
               		list<uint16_t>tempPath 	= currPath;
					uint32_t nextExtLen 	= currextLen;
					uint32_t nextUniPos 	= curruniPos;
					int32_t tempNumOfBases 	= 0;

					//	calculate the score of an extension of a successor and add it to its current score
					int32_t addScore = contRightX_Drop_BFS(*nI, iniQoff, nextHitLen, nextExtLen, q, mscore, mmscore, X, nextScore, nextUniPos, tempPath, explCount, quorum, searchSet, advIdx, check, tempNumOfBases, false);
					tempPath.push_back(sucID);
					int32_t fullScore = currScore+addScore;

					//	if we reached the end of the unitig and still place in the list
					//	add extension to the list
					if(check && bestUnitigs.size() <= numOfUnitig){
						bestUnitigs.push_back(make_tuple(*nI,fullScore,nextScore,tempPath,nextHitLen,nextExtLen, nextUniPos, 0));
						bestScores.push_back(fullScore);

					//	if we reached the end of the unitig and there is no place in the list
					} else if(check){

						//	get minimum score from the score list
						int minElement = *min_element(bestScores.begin(), bestScores.end());

						//	if the calculated score is higher than the minimum score then replace it
						if(minElement < fullScore){
							for(uint i = 0; i < bestScores.size(); i++){
								if(bestScores[i] == minElement){
									bestScores[i] = fullScore;
									bestUnitigs[i] = make_tuple(*nI,fullScore,nextScore,tempPath,nextHitLen,nextExtLen, nextUniPos, 0);
									break;
								}
							}
						}	
					}

					//	update the best scoring extension
					if (fullScore > maxScore) {
                    	maxScore = fullScore;
                    	hitLen = nextHitLen;
						bestPath = tempPath;

					}
				}
			}
		}
	}

	//	update path of best extension and return best score
	extPth = bestPath;
	return maxScore;
}

/*

	shorterVector bestUnitigs;
	vector<int> bestScores;
	for(shorterTemp nI = sucIter.begin(); nI != sucIter.end(); ++nI){
		++sucID;

		uint32_t tempextLen = extLen;
		uint32_t tmpHitLen = hitLen;
		uint32_t tmpuniPos = uniPos;
		int32_t tmpScore = lastExtSeedTmpScore;
		tempPath.clear();

		int numOfBases = 0;
		int32_t tempScore = contRightX_Drop_BFS(*nI, iniQoff, tmpHitLen, tempextLen, q, mscore, mmscore, X, tmpScore, tmpuniPos, tempPath, explCount, quorum, searchSet, advIdx, check, numOfBases, false);
		tempPath.push_back(sucID);
		int fullScore = tempScore;

		if(check && bestUnitigs.size() < numOfUnitig){
			bestUnitigs.push_back(make_tuple(nI,fullScore,tmpScore,tempPath,tmpHitLen,tempextLen, tmpuniPos));
			bestScores.push_back(fullScore);

		} else {
			if(check) {
				int tempMinElement = *min_element(bestScores.begin(), bestScores.end());
				if(tempMinElement < fullScore){
					for(uint i = 0; i < bestScores.size(); i++){
						if(bestScores[i] == tempMinElement){
							bestScores[i] = fullScore;
							bestUnitigs[i] = make_tuple(nI,fullScore,tmpScore,tempPath,tmpHitLen,tempextLen, tmpuniPos);
							break;
						}	
					}
				}
			}	
		}

		if(fullScore > maxScore){
			maxScore = fullScore;
            bestPath = tempPath;
            hitLen = tmpHitLen;
		}
	}

	while(!(bestUnitigs.empty())){

		for(uint i = 0; i < bestUnitigs.size(); i++) {

			shorterTuple temptuple 	= bestUnitigs[0];
			shorterTemp currUnitig 	= std::get<0>(temptuple);
			uint32_t currScore 		= std::get<1>(temptuple);
			int32_t currtmpScore 	= std::get<2>(temptuple);
			pathList currPath 		= std::get<3>(temptuple);
			uint32_t currHitLen 	= std::get<4>(temptuple);
			uint32_t currextLen 	= std::get<5>(temptuple);
			uint32_t currUniPos 	= std::get<6>(temptuple);

			UnitigColorMap<UnitigInfo> tempcurrUnitig = *currUnitig;
			shorterContainer sucIter2 = tempcurrUnitig.getSuccessors();

			bestUnitigs.erase(bestUnitigs.begin() + 0);
			bestScores.erase(bestScores.begin() + 0);

			if((currextLen + iniQoff < q.length())) {
				sucID = 0;
				for(shorterTemp nI = sucIter2.begin(); nI != sucIter2.end(); ++nI){
					check = false;
					++sucID;

					uint32_t tmpHitLen = currHitLen;
					int32_t tmpScore = currtmpScore;
                	tempPath = currPath;
					uint32_t tmpExtLen = currextLen;
					uint32_t nextUniPos = currUniPos;

					int numOfBases = 0;
					int32_t tempScore = contRightX_Drop_BFS(*nI, iniQoff, tmpHitLen, tmpExtLen, q, mscore, mmscore, X, tmpScore, nextUniPos, tempPath, explCount, quorum, searchSet, advIdx, check, numOfBases, false);
					tempPath.push_back(sucID);
					int fullScore = currScore+tempScore;

					if(check && bestUnitigs.size() < numOfUnitig){
						bestUnitigs.push_back(make_tuple(nI, fullScore, tmpScore, tempPath, tmpHitLen, tmpExtLen, nextUniPos));
						bestScores.push_back(fullScore);
					} else {
						if(check) {
							int tempMinElement = *min_element(bestScores.begin(), bestScores.end());
							if(tempMinElement < fullScore){
								for(uint i = 0; i < bestScores.size(); i++){
									if(bestScores[i] == tempMinElement){
										bestScores[i] = fullScore;
										bestUnitigs[i] = make_tuple(nI, fullScore, tmpScore, tempPath, tmpHitLen, tmpExtLen, nextUniPos);
										break;
									}	
								}
							}
						}	
					}

					if(fullScore > maxScore){
						maxScore = fullScore;
						bestPath = tempPath;
                    	hitLen = tmpHitLen;
						uniPos = nextUniPos;
					}
				}
			}
		}
	}

	extPth = bestPath;
	return maxScore;
}
*/



//	does a BFS-step of a given number of bases for all saved unitigs
//	and returns them in a priority queue

/*
outputTypes calcUnitigsMitBasenberechnung(inputTypes inputStructure) {

	//	get all information contained and given in the input structure
	queue<shorterTuple> extensionQueue 		= inputStructure.extensionQueue;
	uint32_t iniQoff 						= inputStructure.iniQoff;
	string q 								= inputStructure.q;
	uint16_t mscore 						= inputStructure.mscore;
	int16_t mmscore 						= inputStructure.mmscore;
	int16_t X 								= inputStructure.X;
	uint32_t explCount 						= inputStructure.explCount;
	uint32_t quorum 						= inputStructure.quorum;
	list<pair<string, size_t>> searchSet 	= inputStructure.searchSet;
	bool advIdx 							= inputStructure.advIdx;
	int32_t maxScore 						= inputStructure.maxScore;
	int32_t numOfBases 						= inputStructure.numOfBases;
	uint32_t hitLen 						= inputStructure.hitLen;
	list<uint16_t> bestPath					= inputStructure.bestPath;

	//	create a priority queue sorted decreasingly by score for all finished extensions
	shorterPrioQueue bestUnitigsPrioQueue(prioLongest);

	//	process all extensions
	while(!extensionQueue.empty()){

		//	get next extension
		shorterTuple currExtension = extensionQueue.front();
		extensionQueue.pop();

		//	extract all relevant information of the extension
		UnitigColorMap<UnitigInfo> currUnitig 		= get<0>(currExtension);
		uint32_t currScore 			= get<1>(currExtension);
		int32_t currtmpScore 		= get<2>(currExtension);
		list<uint16_t> currPath		= get<3>(currExtension);
		uint32_t currHitLen 		= get<4>(currExtension);
		uint32_t currextLen 		= get<5>(currExtension);
		uint32_t curruniPos 		= get<6>(currExtension);
		int currnumOfBases 			= get<7>(currExtension);
		int32_t tempNumOfBases 		= numOfBases;

		//	get the successors of the unitig
		shorterContainer sucIter2 = currUnitig.getSuccessors();

		//	stop if we reached the end of the query
		if((currextLen + iniQoff < q.length())){

			//	if the extension is marked we don't lock at successors
			if(currnumOfBases == -1){
				bool check = false;
				list<uint16_t>tempPath 	= currPath;

				//	continue extension on the current unitig and calculate the new score
				int32_t addScore = contRightX_Drop_BFS(currUnitig, iniQoff, currHitLen, currextLen, q, mscore, mmscore, X, currtmpScore, curruniPos, tempPath, explCount, quorum, searchSet, advIdx, check, tempNumOfBases, true);
				int32_t fullScore = currScore+addScore;

				//	if having compared enough bases while still not reaching the end of the unitig
				//	add marked extension to the priority queue
				if((tempNumOfBases == 0) && !check){
					bestUnitigsPrioQueue.push(make_tuple(currUnitig,fullScore,currtmpScore,tempPath,currHitLen,currextLen, curruniPos, -1));

				//	if not having compared enough bases and we reached the end of the unitig
				//	the extension is added back to the queue unmarked
				} else if(check && tempNumOfBases != 0){
					extensionQueue.push(make_tuple(currUnitig,fullScore,currtmpScore,tempPath,currHitLen,currextLen,curruniPos,tempNumOfBases));
				} else {
					//cout << "removed" << endl;
				}

				//	update the best scoring extension
				if (fullScore > maxScore) {
                    maxScore = fullScore;
                    hitLen = currHitLen;
					bestPath = tempPath;
				}

			//	if the extension is not marked go through the succesors of the unitig	
			} else {
				uint16_t sucID = 0;

				//	iterate over all successors
				for(shorterTemp nI = sucIter2.begin(); nI != sucIter2.end(); ++nI){
					++sucID;

					//	copy current extension for multiple iterations
					//	because the variables will be overwritten when calculation the extension of the first successor
					bool check = false;
					uint32_t nextHitLen 	= currHitLen;
					int32_t nextScore 		= currtmpScore;
               		list<uint16_t>tempPath 	= currPath;
					uint32_t nextExtLen 	= currextLen;
					uint32_t nextUniPos 	= curruniPos;
					int32_t tempNumOfBases 	= numOfBases;

					//	calculate the score of an extension of a successor and add it to its current score
					int32_t addScore = contRightX_Drop_BFS(*nI, iniQoff, nextHitLen, nextExtLen, q, mscore, mmscore, X, nextScore, nextUniPos, tempPath, explCount, quorum, searchSet, advIdx, check, tempNumOfBases, true);
					tempPath.push_back(sucID);
					int32_t fullScore = currScore+addScore;

					//	if having compared enough bases while still not reaching the end of the unitig
					//	add marked extension to the priority queue
					if((tempNumOfBases == 0) && !check){
						bestUnitigsPrioQueue.push(make_tuple(*nI,fullScore,nextScore,tempPath,nextHitLen,nextExtLen, nextUniPos, -1));

					//	if not having compared enough bases and we reached the end of the unitig
					//	the extension is added back to the queue unmarked
					} else if(check && tempNumOfBases != 0){
						extensionQueue.push(make_tuple(*nI,fullScore,nextScore,tempPath,nextHitLen,nextExtLen,nextUniPos,tempNumOfBases));
					} else {
						//cout << "removed" << endl;
					}

					//	update the best scoring extension
					if (fullScore > maxScore) {
                    	maxScore = fullScore;
                    	hitLen = nextHitLen;
						bestPath = tempPath;
					}
				}
			}
		}
	}

	//	return the best extension and the priority queue
	outputTypes p;

	p.bestUnitigsPrioQueue = bestUnitigsPrioQueue;
	p.maxScore = maxScore;
	p.hitLen = hitLen;
	p.bestPath = bestPath;

	return p;
}


//	takes up to numOfAcceptedExtensions elements from the given priority queue
//	and returns them as a queue
queue<shorterTuple> getBestUnitigs(shorterPrioQueue bestUnitigsPrioQueue, uint numOfAcceptedExtensions){
	queue<shorterTuple> bestUnitigs;
	
	for(uint i = 0; i < numOfAcceptedExtensions; i++) {
		if(bestUnitigsPrioQueue.empty()){
			break;
		} else {
			bestUnitigs.push(bestUnitigsPrioQueue.top());
			bestUnitigsPrioQueue.pop(); 
		}
	}

	return bestUnitigs;
}
*/

//	This function initiates the extension on all successors of a unitig and returns the best one considering a quorum and a search color set
//	using a BFS heuristic where x bases are compared and the best scoring extensions are followed upon
int32_t extendAtNextUnitig_BFS_SMART3(const UnitigColorMap<UnitigInfo> startUnitig, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &uniPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx){
	
	//	best score
	int32_t maxScore = 0;

	//	max number of unitigs pushed to the next iteration
	uint numOfUnitig = 1000000;

	//	number of bases considered in each iteration
	int32_t numOfBases = 30;

	//	priority queue ordered decreasingly by score
	shorterPrioQueue bestUnitigsPrioQueue(prioLongest);

	//	queue for next iterations extensions
	queue<shorterTuple> bestUnitigs;

	//	temporary and best global path
	list<uint16_t> tempPath;
	list<uint16_t> bestPath;

	//	queue for starting extensions
	//queue<shorterTuple> startQueue;

	//	initialize start extension parameters
	//int32_t tempScore = 0;
	//int32_t tmpScore = lastExtSeedTmpScore;
	//uint32_t tmpHitLen = hitLen;
	//uint32_t tempextLen = extLen;
	//uint32_t tmpuniPos = uniPos;

	//	add starting extension to the starting queue
	shorterTuple startTuple = make_tuple(startUnitig, 0, lastExtSeedTmpScore, tempPath, hitLen, extLen, uniPos, 0);
	bestUnitigs.push(startTuple);

	//	prepare the parameters given to the BFS extension function
	inputTypes inputStruct;
	//inputStruct.extensionQueue = startQueue;
	inputStruct.iniQoff 		= iniQoff;
	inputStruct.q 				= q;
	inputStruct.mscore 			= mscore;
	inputStruct.mmscore 		= mmscore;
	inputStruct.X 				= X;
	inputStruct.explCount 		= explCount;
	inputStruct.quorum 			= quorum;
	inputStruct.searchSet 		= searchSet;
	inputStruct.advIdx 			= advIdx;
	//inputStruct.maxScore = maxScore;
	inputStruct.numOfBases 		= numOfBases;
	//inputStruct.hitLen = hitLen;
	//inputStruct.bestPath = bestPath;
	inputStruct.compareBases	= true;

	outputTypes outputStruct;

	//	first BFS extension
	//outputStruct = calcUnitigsMitBasenberechnung(inputStruct);

	//	extract updated information
	//bestUnitigsPrioQueue = outputStruct.bestUnitigsPrioQueue;
	//maxScore = outputStruct.maxScore;
	//hitLen = outputStruct.hitLen;
	//bestPath = outputStruct.bestPath;

	//	get the best extensions
	//bestUnitigs = getBestUnitigs(bestUnitigsPrioQueue,numOfUnitig);

	//	repeat BFS extension until no more extensions are possible
	while(!(bestUnitigs.empty())){

		//	update extension queue and best result found so far
		inputStruct.extensionQueue = bestUnitigs;
		inputStruct.maxScore = maxScore;
		inputStruct.hitLen = hitLen;
		inputStruct.bestPath = bestPath;

		//	extend and score again
		outputStruct = calcUnitigs(inputStruct);

		//	update with new results
		bestUnitigsPrioQueue = outputStruct.bestUnitigsPrioQueue;
		maxScore = outputStruct.maxScore;
		hitLen = outputStruct.hitLen;
		bestPath = outputStruct.bestPath;

		//	choose extensions for the next iteration
		bestUnitigs = getBestUnitigs(bestUnitigsPrioQueue,numOfUnitig);
	}

	//	update path of best extension and return best score
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
				//cout << "hit->score before:" << hit->score << endl;
				switch(extend_modus) {
					case 0:
						hit->score += extendAtNextUnitig(hit->origUni.getSuccessors(), hit->offQ, hit->length, tmpSeedLen, q, mscore, mmscore, X, tmpScore, iniUniPos, extPth, explCount, quorum, searchSet, advIdx);
						break;
					case 1:
						hit->score += extendAtNextUnitig_BFS(hit->origUni, hit->offQ, hit->length, tmpSeedLen, q, mscore, mmscore, X, tmpScore, iniUniPos, extPth, explCount, quorum, searchSet, advIdx);
						break;
					case 2:
						hit->score += extendAtNextUnitig_BFS_SMART1(hit->origUni, hit->offQ, hit->length, tmpSeedLen, q, mscore, mmscore, X, tmpScore, iniUniPos, extPth, explCount, quorum, searchSet, advIdx);
						break;
					case 3:
						hit->score += extendAtNextUnitig_BFS_SMART2(hit->origUni, hit->offQ, hit->length, tmpSeedLen, q, mscore, mmscore, X, tmpScore, iniUniPos, extPth, explCount, quorum, searchSet, advIdx);
						break;
					case 4:
						hit->score += extendAtNextUnitig_BFS_SMART3(hit->origUni, hit->offQ, hit->length, tmpSeedLen, q, mscore, mmscore, X, tmpScore, iniUniPos, extPth, explCount, quorum, searchSet, advIdx);
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

//	This function continues an extension to the right on a successive unitig of a seed lying on the query's reference strand considering a quorum and a search color set. Returns the maximum score reached.
//	using an iterative BFS algorithm and can compare a certain number of bases
int32_t contRightX_Drop_BFS(const UnitigColorMap<UnitigInfo> &sucUnitig, const uint32_t &iniQoff, uint32_t &hitLen, uint32_t &extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, int32_t &tmpScore, uint32_t &uniSeqPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx, bool &check, int &numOfBases, const bool &compBases){
	int32_t progress, score = 0;
	int32_t overlap = sucUnitig.getGraph()->getK() - 1;
	uint32_t iniSeqPos;
	uint32_t tmpSLen;
	//Even if we are on the reverse complementary strand and no position is covered on this unitig, because we have checked the first k - 1 position on the last unitig already
	int32_t checkedPos = getSrchCritCov(sucUnitig, quorum, searchSet, compOffset(uniSeqPos, 1, sucUnitig.size, sucUnitig.strand), sucUnitig.strand, advIdx);
	string sucUniSeq = sucUnitig.mappedSequenceToString();
	struct Seed *nearestSeed, *prevSeed;

	//cout << "starttmpScore: " << tmpScore << endl;

	//Save the initial offset in the current unitig which we need for all nearest neighbor calculations
	iniSeqPos = uniSeqPos;
	//Find the nearest seed that we might be able to reach during our extension
	nearestSeed = searchRightNeighbor(sucUnitig.getData()->getData(sucUnitig)->getSeed(sucUnitig.strand), iniQoff, extLen, iniSeqPos, prevSeed);
	//Perform the X-drop algorithm on the successive unitig
	tmpSLen = 0;
	check = false;

	//We are done if we have reached the end of the query
	while(iniQoff + extLen + tmpSLen < q.length()){
		//Testing
		//cout << "Not at query's end" << endl;
		//cout << "nearestSeed is " << (nearestSeed == NULL ? "NULL" : "not NULL") << endl;

		//	stop the extension if the number of considered bases is reached
		
		if(numOfBases == 0 && compBases){
			extLen = extLen + tmpSLen;
			break;
		}
		

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
					nearestSeed = searchRightNeighbor(sucUnitig.getData()->getData(sucUnitig)->getSeed(sucUnitig.strand), iniQoff, extLen, iniSeqPos, prevSeed);
				} else{
					//Set seed's successor as head of the seed list
					sucUnitig.getData()->getData(sucUnitig)->setSeed(nearestSeed->nextSeed, sucUnitig.strand);
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
					nearestSeed = searchRightNeighbor(sucUnitig.getData()->getData(sucUnitig)->getSeed(sucUnitig.strand), iniQoff, extLen, iniSeqPos, prevSeed);
				} else{
					//Set seed's successor as head of the seed list
					sucUnitig.getData()->getData(sucUnitig)->setSeed(nearestSeed->nextSeed, sucUnitig.strand);
					//Delete the reached seed
					free(nearestSeed);
					//Reset nearestSeed
					nearestSeed = NULL;
				}
			}
		} else{
			//Check up to which point we have to compare the unitig sequence
			if(!sucUnitig.getSuccessors().hasSuccessors()){
				overlap = 0;
			}

			//Check whether we have reached the end of the unitig's sequence
			if(uniSeqPos < sucUnitig.size - overlap){
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
				if(compBases){
					--numOfBases;
				}
			} else{
				//Check if the current unitig has successors
				if(overlap != 0){
					//Calculate the position in the next unitig's sequence we have to start with
					uniSeqPos = uniSeqPos - sucUnitig.size + overlap;
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

	//cout << "contRightX_Drop_BFS Final score: " << score << std::endl;
	//cout << "Unitig: " << sucUnitig->mappedSequenceToString() << ", hitLen: " << hitLen << ", score: " << score << ", uniPos: " << uniSeqPos << endl;
	//cout << "final uniSeqPos = " << uniSeqPos << std::endl;
	//cout << "finaltmpScore: " << tmpScore << endl;
	return score;
}

//	This function continues an extension to the right on a successive unitig of a seed lying on the query's reference strand considering a quorum and a search color set. Returns the maximum score reached.
//	using an iterative BFS algorithm and comparing a certain number of bases
int32_t contRightX_Drop_BFS_2(const UnitigColorMap<UnitigInfo> &sucUnitig, const uint32_t &iniQoff, uint32_t &hitLen, uint32_t &extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, int32_t &tmpScore, uint32_t &uniSeqPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx, bool &check, int &numOfBases){
	int32_t progress, score = 0;
	int32_t overlap = sucUnitig.getGraph()->getK() - 1;
	uint32_t iniSeqPos;
	uint32_t tmpSLen;
	//Even if we are on the reverse complementary strand and no position is covered on this unitig, because we have checked the first k - 1 position on the last unitig already
	int32_t checkedPos = getSrchCritCov(sucUnitig, quorum, searchSet, compOffset(uniSeqPos, 1, sucUnitig.size, sucUnitig.strand), sucUnitig.strand, advIdx);
	string sucUniSeq = sucUnitig.mappedSequenceToString();
	struct Seed *nearestSeed, *prevSeed;

	//cout << "starttmpScore: " << tmpScore << endl;

	//Save the initial offset in the current unitig which we need for all nearest neighbor calculations
	iniSeqPos = uniSeqPos;
	//Find the nearest seed that we might be able to reach during our extension
	nearestSeed = searchRightNeighbor(sucUnitig.getData()->getData(sucUnitig)->getSeed(sucUnitig.strand), iniQoff, extLen, iniSeqPos, prevSeed);
	//Perform the X-drop algorithm on the successive unitig
	tmpSLen = 0;
	check = false;

	//We are done if we have reached the end of the query
	while(iniQoff + extLen + tmpSLen < q.length()){
		//Testing
		//cout << "Not at query's end" << endl;
		//cout << "nearestSeed is " << (nearestSeed == NULL ? "NULL" : "not NULL") << endl;

		//	stop the extension if the number of considered bases is reached
		if(numOfBases == 0){
			extLen = extLen + tmpSLen;
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
					nearestSeed = searchRightNeighbor(sucUnitig.getData()->getData(sucUnitig)->getSeed(sucUnitig.strand), iniQoff, extLen, iniSeqPos, prevSeed);
				} else{
					//Set seed's successor as head of the seed list
					sucUnitig.getData()->getData(sucUnitig)->setSeed(nearestSeed->nextSeed, sucUnitig.strand);
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
					nearestSeed = searchRightNeighbor(sucUnitig.getData()->getData(sucUnitig)->getSeed(sucUnitig.strand), iniQoff, extLen, iniSeqPos, prevSeed);
				} else{
					//Set seed's successor as head of the seed list
					sucUnitig.getData()->getData(sucUnitig)->setSeed(nearestSeed->nextSeed, sucUnitig.strand);
					//Delete the reached seed
					free(nearestSeed);
					//Reset nearestSeed
					nearestSeed = NULL;
				}
			}
		} else{
			//Check up to which point we have to compare the unitig sequence
			if(!sucUnitig.getSuccessors().hasSuccessors()){
				overlap = 0;
			}

			//Check whether we have reached the end of the unitig's sequence
			if(uniSeqPos < sucUnitig.size - overlap){
				//Are search criteria still fulfilled?
				if(checkedPos <= 0){
					break;
				}

				//TEsting
				//cout << "tmpScore: " << tmpScore << endl;
				//cout << "compUScore(sucUniSeq[uniSeqPos], q[iniQoff + extLen + tmpSLen], mscore, mmscore): " << compUScore(sucUniSeq[uniSeqPos], q[iniQoff + extLen + tmpSLen], mscore, mmscore) << endl;
				//cout << "score: " << score << endl;
				//cout << "extLen: " << extLen << endl;
				//cout << "tmpSLen: " << tmpSLen << endl;
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
					uniSeqPos = uniSeqPos - sucUnitig.size + overlap;
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
	//cout << "Unitig: " << sucUnitig.mappedSequenceToString() << ", hitLen: " << hitLen << ", score: " << score << ", uniPos: " << uniSeqPos << endl;
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