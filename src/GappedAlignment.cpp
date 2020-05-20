#include "GappedAlignment.h"
#include "Hit.h"
#include "Search.h"

//This function calculates a semi-global alignment of a unitig and the query sequence considering a quorum and a search color set. Returns true if calculations on unitig could be finished (i.e. the unitig sequence was not too long to be stored inside the edit matrix).
bool calcSemiGlobAlignment(const UnitigColorMap<seedlist> &uni, const string &q, uint32_t& posU, uint32_t& posQ, const uint16_t &X, const uint32_t &maxGaps, uint32_t& maxPosQ, uint32_t& maxPosU, struct Algn &maxAlgn, struct Algn &brdAlgn, int32_t& maxScore, int32_t& eMax, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool srchCritCheck){
	bool termCalcs, matShrunk = false;
	uint32_t matBrth, matHgth, i, j, lCalcBorder, rCalcBorder, endBuf = 0, uCmpSeqLen, qCmpSeqLen, maxMatPosI = 0, maxMatPosJ = 0, eMaxPosU = posU, eMaxPosQ = posQ, edgeMaxPosI = 0, edgeMaxPosJ = 0;
	int32_t covPos;
	string uSeq;

	//Initialize the current max position
	eMax = -X;
	//Get unitig's sequence
	uSeq = uni.mappedSequenceToString();

	//Check whether current unitig has a successor
	if(uni.getSuccessors().hasSuccessors()){
		endBuf = uni.getGraph()->getK() - 1;
	}

	//Calculate length of unitig sequence to compare
	uCmpSeqLen = uni.size - (posU + endBuf);

	//Do a search criteria check only if necessary
	if(srchCritCheck){
		//Get the number of positions fulfilling our search crit on this unitig reduced by our starting offset
		covPos = getSrchCritCov(uni, quorum, searchSet, posU, true);

		//Check if there is at least one position fulfilling the search crit
		if(covPos <= 0){
			//We don't need to proceed
			uCmpSeqLen = 0;
		} else if((uint32_t) covPos < uCmpSeqLen){
			//We only need to calculate until our search crits are not fulfilled anymore
			uCmpSeqLen = covPos;
		}	
	}

	//Calculate length of query sequence to compare
	qCmpSeqLen = q.length() - posQ;

	//Make sure that the matrix is not becoming too big
	if(uCmpSeqLen > MAX_MATRIX_SIZE){
		//If we reach the lower matrix border we definitely do not miss any sequence of the unitig
		matBrth = min((uint32_t) MAX_MATRIX_SIZE, qCmpSeqLen + maxGaps) + 1;
		//Do not matter which matrix border we reach and the matrix should not become too big
		matHgth = min(qCmpSeqLen, (uint32_t) MAX_MATRIX_SIZE) + 1;
		//If we still have uncompared query sequence mark that we are not done
		matShrunk = true;
	} else{
		//Consider as much unitig sequence as necessary but not more
		matHgth = min(uCmpSeqLen + maxGaps, qCmpSeqLen) + 1;
		//Make sure we cannot reach the lower matrix border except if our calculations are done
		matBrth = min(qCmpSeqLen + maxGaps, uCmpSeqLen) + 1;

	}

	//Initialize dynamic programming matrix
	int32_t **mat = (int32_t**) malloc(matHgth * sizeof(int32_t*));
	for(i = 0; i < matHgth; ++i) mat[i] = (int32_t*) malloc(matBrth * sizeof(int32_t));

	//Initialize upper left cell
	mat[0][0] = maxScore;

	//Initialize first column
	for(i = 1; i < matHgth; ++i){
		if(i <= maxGaps){
			//GAP_SCORE is supposed to be negative
			mat[i][0] = maxScore + i * GAP_SCORE;//TODO: This has to be changed as soon as we are not using unit score anymore!
		} else{
			mat[i][0] = INT32_MIN - GAP_SCORE;//TODO: This has to be changed as soon as we are not using unit score anymore!
			//We want to save time by not treating the whole matrix
			break;
		}
	}

	//Initialize first row
	for(j = 1; j < matBrth; ++j){
		if(j <= maxGaps){
			//GAP_SCORE is supposed to be negative
			mat[0][j] = maxScore + j * GAP_SCORE;//TODO: This has to be changed as soon as we are not using unit score anymore!
		} else{
			mat[0][j] = INT32_MIN - GAP_SCORE;//TODO: This has to be changed as soon as we are not using unit score anymore!
			//We want to save time by not treating the whole matrix
			break;
		}
	}

	//Set borders
	lCalcBorder = 1;

	//Fill the matrix
	for(i = 1; i < matHgth; ++i){
		//Reset abort var
		termCalcs = true;

		//Check if we can adjust left alignment calculation border
		if(maxGaps < i){
			lCalcBorder = i - maxGaps;
			//Initialize cell to the left of the first one to calculate as we will need it there
			mat[i][lCalcBorder - 1] = INT32_MIN - GAP_SCORE;//TODO: This has to be changed as soon as we are not using unit score anymore!
		}

		//Adjust right border
		rCalcBorder = min(i + maxGaps, matBrth - 1);
		//If existing, initialize cell to the right of the last calculated one as we might need it for calculations in the next iteration
		if(rCalcBorder < matBrth - 1) mat[i][rCalcBorder + 1] = INT32_MIN - GAP_SCORE;//TODO: This has to be changed as soon as we are not using unit score anymore!

		for(j = lCalcBorder; j <= rCalcBorder; ++j){
			//Calculate current cells value and check we don't drop
			if((mat[i][j] = max(mat[i - 1][j - 1] + compUScore(q[posQ + i - 1], uSeq[posU + j - 1]), max(mat[i][j - 1] + GAP_SCORE, mat[i - 1][j] + GAP_SCORE))) > -X){//TODO: This has to be changed as soon as we are not using unit score anymore!
				//Make sure we do not terminate our calculations too soon
				termCalcs = false;

				//Check whether we have a new max score
				if(maxScore < mat[i][j]){
					//Update max score
					maxScore = mat[i][j];
					//Update maximum position
					maxPosQ = posQ + i - 1;
					maxPosU = posU + j - 1;
					//Save maximum position in matrix
					maxMatPosI = i;
					maxMatPosJ = j;
				}

				//Check whether we are at the right or lower matrix edge and whether we have found a new maximum
				if((i == matHgth - 1 || j == matBrth - 1) && eMax <= mat[i][j]){//TODO Using <= instead of < here let's the edge maximum be as close to the diagonal as possible. This minimizes the number of gaps in our alignment, but also prevents us from finding the best possible score sometimes (like in test case 5 of the startRightGappedAlignment test). Maybe it is worth to change this, but could also be that we have to consider this when we want do implement hybrid matrices' alignment
					//Update edge maximum
					eMax = mat[i][j];
					eMaxPosQ = posQ + i - 1;
					eMaxPosU = posU + j - 1;
					//Save edge maximum in matrix
					edgeMaxPosI = i;
					edgeMaxPosJ = j;
				}
			}
		}

		//Check whether calculations may be terminated
		if(termCalcs) break;
	}

	//Check if we had to terminate prematurely due to a lack of search criteria coverage
	if(uni.size - (posU + endBuf) > uCmpSeqLen){
		//Reset edge maximum to prevent an extension on a successive unitig
		eMax = -X;
		//Set edge maximum matrix positions to overall maximum postions to prevent a trace back of the border alignment
		edgeMaxPosI = maxMatPosI;
		edgeMaxPosJ = maxMatPosJ;
	}

	//Check if maximum score is at the matrix border
	if(maxMatPosI == edgeMaxPosI && maxMatPosJ == edgeMaxPosJ){
		//Do backtracing to get the alignment maximum scoring alignment
		traceBack(maxMatPosI, maxPosQ, q, maxMatPosJ, maxPosU, uSeq, mat, maxAlgn);
		//Copy maximum scoring alignment to save time
		brdAlgn = maxAlgn;
	} else{
		//Do backtracing to get the alignment maximum scoring alignment
		traceBack(maxMatPosI, maxPosQ, q, maxMatPosJ, maxPosU, uSeq, mat, maxAlgn);
		//Backtrace alignment to go on with on next unitig
		traceBack(edgeMaxPosI, eMaxPosQ, q, edgeMaxPosJ, eMaxPosU, uSeq, mat, brdAlgn);
	}

	//Free matrix
	for(i = 0; i < matHgth; ++i) free(mat[i]);
	free(mat);

	//Update positions in unitig und q
	if(eMax > -X){
		posU = eMaxPosU;
		posQ = eMaxPosQ;
	} else{
		posU = maxPosU;
		posQ = maxPosQ;
	}

	if(eMax > -X && matShrunk && posQ < q.length() - 1) return false;
	
	return true;
}

//This function calculates a banded alignment of a unitig and the query sequence to the left considering a quorum and a search color set. Returns true if calculations on unitig could be finished (i.e. the unitig sequence was not to long to be stored inside the edit matrix).
bool calcLeftGlobAlignment(const UnitigColorMap<seedlist> &uni, const string &q, uint32_t& posU, uint32_t& posQ, const uint16_t &X, const uint32_t &maxGaps, uint32_t& maxPosQ, uint32_t& maxPosU, struct Algn &maxAlgn, struct Algn &brdAlgn, int32_t& maxScore, int32_t& eMax, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool srchCritCheck){
	bool termCalcs, matShrunk = false;
	uint32_t matBrth, matHgth, i, j, lCalcBorder, rCalcBorder, maxMatPosI = 0, maxMatPosJ = 0, eMaxPosU = posU, eMaxPosQ = posQ, edgeMaxPosI = 0, edgeMaxPosJ = 0, uCmpSeqLen = posU;
	int32_t covPos = 0;
	string uSeq;

	//Get unitig's sequence
	uSeq = uni.mappedSequenceToString();

	//Do we have to check search criteria on this unitig?
	if(srchCritCheck){
		covPos = getSrchCritCov(uni, quorum, searchSet, posU, false);

		//Check if there is at least one position fulfilling the search crit
		if(covPos > 0 && (uint32_t) covPos <= uCmpSeqLen){
			//We only need to calculate until our search crits are not fulfilled anymore (-1, because offsets start at 0)
			uCmpSeqLen = (uint32_t) covPos - 1;
		}
	}

	//Make sure that the matrix is not becoming too big
	if(uCmpSeqLen > MAX_MATRIX_SIZE){
		//If we reach the lower matrix border we definitely do not miss any sequence of the unitig
		matBrth = min((uint32_t) MAX_MATRIX_SIZE, posQ + maxGaps) + 2;//+2 since we need 1 column and row for the empty string and offset indexes start with 0
		//Do not matter which matrix border we reach and the matrix should not become too big
		matHgth = min(posQ, (uint32_t) MAX_MATRIX_SIZE) + 2;
		//Mark that we do not consider the whole unitig sequence
		matShrunk = true;
	} else{
		//Make sure we cannot reach the lower matrix border except if our calculations are done
		matBrth = min(uCmpSeqLen, posQ + maxGaps) + 2;
		//Consider as much unitig sequence as necessary but not more
		matHgth = min(posQ, uCmpSeqLen + maxGaps) + 2;
	}

	//Borders are not set correctly if search crit is not fulfilled
	if(srchCritCheck && covPos <= 0){
		//If the search criteria is not fulfilled we are not going to do anything!
		matHgth = 1;
		matBrth = 1;
	}

	//Initialize dynamic programming matrix
	int32_t **mat = (int32_t**) malloc(matHgth * sizeof(int32_t*));
	for(i = 0; i < matHgth; ++i) mat[i] = (int32_t*) malloc(matBrth * sizeof(int32_t));

	//Initialize upper left cell
	mat[0][0] = maxScore;
	//Initialize edge maximum
	eMax = -X;

	//Initialize first column
	for(i = 1; i < matHgth; ++i){
		if(i <= maxGaps){
			//GAP_SCORE is supposed to be negative
			mat[i][0] = maxScore + i * GAP_SCORE;//TODO: This has to be changed as soon as we are not using unit score anymore!
		} else{
			mat[i][0] = INT32_MIN - GAP_SCORE;//TODO: This has to be changed as soon as we are not using unit score anymore!
		}
	}

	//Initialize first row
	for(j = 1; j < matBrth; ++j){
		if(j <= maxGaps){
			//GAP_SCORE is supposed to be negative
			mat[0][j] = maxScore + j * GAP_SCORE;//TODO: This has to be changed as soon as we are not using unit score anymore!
		} else{
			mat[0][j] = INT32_MIN - GAP_SCORE;//TODO: This has to be changed as soon as we are not using unit score anymore!
		}
	}

	//Set left border
	lCalcBorder = 1;

	//Fill the matrix
	for(i = 1; i < matHgth; ++i){
		//Reset abort var
		termCalcs = true;

		//Check if we can adjust left alignment calculation border
		if(maxGaps < i){
			lCalcBorder = i - maxGaps;
			//Initialize cell to the left of the first one to calculate as we will need it there
			mat[i][lCalcBorder - 1] = INT32_MIN - GAP_SCORE;//TODO: This has to be changed as soon as we are not using unit score anymore!
		}

		//Adjust right border
		rCalcBorder = min(i + maxGaps, matBrth - 1);
		//If existing, initialize cell to the right of the last calculated one as we might need it for calculations in the next iteration
		if(rCalcBorder < matBrth - 1) mat[i][rCalcBorder + 1] = INT32_MIN - GAP_SCORE;//TODO: This has to be changed as soon as we are not using unit score anymore!

		for(j = lCalcBorder; j <= rCalcBorder; ++j){
			//Calculate current cells value and check we don't drop
			if((mat[i][j] = max(mat[i - 1][j - 1] + compUScore(q[posQ - i + 1], uSeq[posU - j + 1]), max(mat[i][j - 1] + GAP_SCORE, mat[i - 1][j] + GAP_SCORE))) > -X){//TODO: This has to be changed as soon as we are not using unit score anymore!
				//Make sure we do not terminate our calculations too soon
				termCalcs = false;

				//Check whether we have a new max score (NOTE: We could check only for better scores here as well. Do it this way will in doubt lead to larger alignments.)//TODO: Maybe we don't what this because this might lead to less significant alignments!
				if(maxScore < mat[i][j]){
					//Update max score
					maxScore = mat[i][j];
					//Update maximum position
					maxPosQ = posQ - i + 1;
					maxPosU = posU - j + 1;
					maxMatPosI = i;
					maxMatPosJ = j;
				}

				//Check whether we are at the right or lower matrix edge and whether we have found a new maximum
				if((i == matHgth - 1 || j == matBrth - 1) && eMax <= mat[i][j]){
					//Update edge maximum
					eMax = mat[i][j];
					eMaxPosQ = posQ - i + 1;
					eMaxPosU = posU - j + 1;
					edgeMaxPosI = i;
					edgeMaxPosJ = j;
				}
			}
		}

		//Check whether calculations may be terminated
		if(termCalcs) break;
	}

	//Check if we had to terminate prematurely due to a lack of search criteria coverage
	if(posU > uCmpSeqLen){
		//Reset edge maximum to prevent an extension on a predecessive unitig
		eMax = -X;
		//Set edge maximum matrix positions to overall maximum postions to prevent a trace back of the border alignment
		edgeMaxPosI = maxMatPosI;
		edgeMaxPosJ = maxMatPosJ;
	}

	//Check if maximum score is at the matrix border
	if(maxMatPosI == edgeMaxPosI && maxMatPosJ == edgeMaxPosJ){
		//Do backtracing to get the alignment maximum scoring alignment
		traceForth(maxMatPosI, maxPosQ, q, maxMatPosJ, maxPosU, uSeq, mat, maxAlgn);
		//Copy maximum scoring alignment to save time
		brdAlgn = maxAlgn;
	} else{
		//Do backtracing to get the alignment maximum scoring alignment
		traceForth(maxMatPosI, maxPosQ, q, maxMatPosJ, maxPosU, uSeq, mat, maxAlgn);
		//Backtrace alignment to go on with on next unitig
		traceForth(edgeMaxPosI, eMaxPosQ, q, edgeMaxPosJ, eMaxPosU, uSeq, mat, brdAlgn);
	}

	//Free matrix
	for(i = 0; i < matHgth; ++i) free(mat[i]);
	free(mat);

	//Update positions in unitig und q
	if(eMax > -X){
		posU = eMaxPosU;
		posQ = eMaxPosQ;
	} else{
		posU = maxPosU;
		posQ = maxPosQ;
	}

	if(eMax > -X && matShrunk && posQ > 0) return false;
	
	return true;
}

//This function calculates the continuation of a gapped alignment on a successive unitig or the next peace of the query considering a quorum and a search color set
void contRightGappedAlignment(UnitigColorMap<seedlist> &uni, list<uint16_t> &extPth, const string &q, uint32_t &posQ, uint32_t &posU, const uint16_t &X, const uint32_t &maxGaps, struct Algn &algn, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	bool explSuc = false, unitigDone;
	uint32_t maxPosQ = posQ, maxPosU = posU;
	int32_t maxBorderScore = -X;
	struct Algn globAlgn;
	UnitigColorMap<seedlist> maxUni = uni;

	//Calculate gapped alignment
	unitigDone = calcSemiGlobAlignment(uni, q, posU, posQ, X, maxGaps, maxPosQ, maxPosU, algn, globAlgn, score, maxBorderScore, quorum, searchSet, extPth.empty());

	//Check outcome of alignment calculation
	if(!unitigDone){
		explSuc = contGappedOnSameUni(uni, extPth, q, posQ, posU, X, maxGaps, globAlgn, maxBorderScore, explCount, quorum, searchSet);
	} else if(maxBorderScore > -X && posQ < q.length() - 1){
		//Continue calculations on the next unitig
		explSuc = contGappedOnSuccUni(uni, extPth, q, posQ, posU, X, maxGaps, globAlgn, maxBorderScore, explCount, quorum, searchSet);
	}

	//Check whether exploration of next unitig was successful (this is only the case if a new maximum score has been found as well)
	if(explSuc && score < maxBorderScore){
		//Score has to be updated if we have found a better one on a successive unitig
		score = maxBorderScore;
		//Overwrite previous best alignment
		algn = globAlgn;
	} else{
		//Save new offsets
		posU = maxPosU;
		posQ = maxPosQ;
		//Reset the unitig on which the maximum has been found
		uni = maxUni;
	}
}

//This function calculates the continuation of a gapped alignment on a predecessive unitig or the next peace of the query considering a quorum and a search color set
void contLeftGappedAlignment(UnitigColorMap<seedlist> &uni, list<uint16_t> &extPth, const string &q, uint32_t &posQ, uint32_t &posU, const uint16_t &X, const uint32_t &maxGaps, struct Algn &maxAlgn, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	bool unitigDone, explSuc = false;
	uint32_t maxPosQ = posQ, maxPosU = posU;
	int32_t maxBorderScore;
	struct Algn brdAlgn;
	UnitigColorMap<seedlist> maxUni = uni;

	//Calculate gapped alignment
	unitigDone = calcLeftGlobAlignment(uni, q, posU, posQ, X, maxGaps, maxPosQ, maxPosU, maxAlgn, brdAlgn, score, maxBorderScore, quorum, searchSet, extPth.empty());

	//Check outcome of alignment calculation
	if(!unitigDone){
		explSuc = contLeftOnSameUni(uni, extPth, q, posQ, posU, X, maxGaps, brdAlgn, maxBorderScore, explCount, quorum, searchSet);
	} else if(maxBorderScore > -X && posQ > 0){
		//Continue calculations on the next unitig
		explSuc = contGappedOnPredUni(uni, extPth, q, posQ, posU, X, maxGaps, brdAlgn, maxBorderScore, explCount, quorum, searchSet);
	}

	//Check whether exploration of next unitig was successful (this is only the case if a new maximum score has been found as well)
	if(explSuc && score < maxBorderScore){
		//Score has to be updated if we have found a better one on a successive unitig
		score = maxBorderScore;
		//Save new maximum scoring alignment
		maxAlgn = brdAlgn;
	} else{
		//Save new offsets in hit
		posU = maxPosU;
		posQ = maxPosQ;
		//Reset the unitig on which the maximum has been found
		uni = maxUni;
	}
}

//Continues a gapped alignment on the same unitig as before considering a quorum and a search color set. Returns true if a new maximal score has been found.
bool contGappedOnSameUni(UnitigColorMap<seedlist> &uni, list<uint16_t> &extPth, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, struct Algn &algn, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	uint32_t tmpQOff = qOff + 1, tmpUOff = uOff + 1;
	int32_t tmpScore = score;
	struct Algn tmpAlgn;

	contRightGappedAlignment(uni, extPth, q, tmpQOff, tmpUOff, X, maxGaps, tmpAlgn, tmpScore, explCount, quorum, searchSet);

	//Check whether we have found a new maximum
	if(score < tmpScore){
		//Update maximum
		score = tmpScore;
		qOff = tmpQOff;
		uOff = tmpUOff;
		algn.aSeqG += tmpAlgn.aSeqG;
		algn.aSeqQ += tmpAlgn.aSeqQ;
		return true;
	}

	return false;
}

//This function continues a gapped alignment calculation to the left on the same unitig considering a quorum and a search color set. Returns true if a new maximal score has been found.
bool contLeftOnSameUni(UnitigColorMap<seedlist> &uni, list<uint16_t> &extPth, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, struct Algn &algn, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	uint32_t tmpQOff = qOff - 1, tmpUOff = uOff - 1;
	int32_t tmpScore = score;
	struct Algn tmpAlgn;

	contLeftGappedAlignment(uni, extPth, q, tmpQOff, tmpUOff, X, maxGaps, tmpAlgn, tmpScore, explCount, quorum, searchSet);

	//Check whether we have found a new maximum
	if(score < tmpScore){
		//Update maximum
		score = tmpScore;
		qOff = tmpQOff;
		uOff = tmpUOff;
		algn.aSeqG = tmpAlgn.aSeqG + algn.aSeqG;
		algn.aSeqQ = tmpAlgn.aSeqQ + algn.aSeqQ;
		return true;
	}

	return false;
}

//Initiates a gapped alignment calculation on all successive unitigs considering a quorum and a search color set and returns the alignment's end position in the unitig
bool contGappedOnSuccUni(UnitigColorMap<seedlist> &uni, list<uint16_t> &extPth, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, struct Algn &algn, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	bool suc = false;
	uint16_t sucCount, nextSuc;
	uint32_t tmpPosU, tmpPosQ, maxPosQ, maxPosU;
	int32_t maxScore, tmpScore;
	struct Algn maxAlgn;
	UnitigColorMap<seedlist> sucUni;
	ForwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> sucIter = uni.getSuccessors();

	//Check whether we have reached the maximum recursion depth of an extension
	if(++explCount > MAXRECURSIONDEPTH){
		//Report this incident
		//cerr << "Maximum recursion depth reached during extension. Position in q: " << qOff << endl;
		//Terminate this extension
		return suc;
	}

	//Check whether there are successors to continue
	if(sucIter.hasSuccessors()){
		maxScore = score;
		//Set counter
		sucCount = 1;

		//Check if there still exists an extension path to follow
		if(!extPth.empty()){
			//Get next successor
			nextSuc = extPth.front();
			//Delete that successor from path
			extPth.pop_front();
			//We do not need to count how many successors we have explored yet
			explCount = 0;
		} else{
			//We have no successor given to go on with
			nextSuc = 0;
		}

		//Explore each neighbor
		for(neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> nI = sucIter.begin(); nI != sucIter.end(); ++nI, ++sucCount){
			struct Algn tmpAlgn;

			//Check if we still have to wait for the correct successor to go on with
			if(sucCount < nextSuc) continue;

			//Get next unitig
			sucUni = *nI;
			tmpPosQ = qOff + 1;
			tmpPosU = 0;
			tmpScore = score;
			//Calculate gapped alignment on the next unitig
			contRightGappedAlignment(sucUni, extPth, q, tmpPosQ, tmpPosU, X, maxGaps, tmpAlgn, tmpScore, explCount, quorum, searchSet);

			//Check for new maximum
			if(maxScore < tmpScore){
				//Update maximum
				maxScore = tmpScore;
				maxPosU = tmpPosU;
				maxPosQ = tmpPosQ;
				uni = sucUni;
				maxAlgn = tmpAlgn;
				suc = true;
			}

			//If we reach this point and had a specific successor to go on with we are done
			if(nextSuc > 0) break;
		}

		//Only update score and qOff if exploration was successful
		if(suc){
			//Give back best result found
			score = maxScore;
			qOff = maxPosQ;
			uOff = maxPosU;
			algn.aSeqG += maxAlgn.aSeqG;
			algn.aSeqQ += maxAlgn.aSeqQ;
		}
	}

	return suc;
}

//This function initiates a gapped alignment calculation on all predecessive unitigs considering a quorum and a search color set. Returns true if calculations have been successful
bool contGappedOnPredUni(UnitigColorMap<seedlist> &uni, list<uint16_t> &extPth, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, struct Algn &algn, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	bool suc = false;
	uint16_t predCount, nextPred;
	uint32_t tmpPosU, tmpPosQ, maxPosQ, maxPosU;
	int32_t maxScore, tmpScore;
	struct Algn maxAlgn;
	UnitigColorMap<seedlist> predUni;
	BackwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> predIter = uni.getPredecessors();

	//Check whether we have reached the maximum recursion depth of an extension
	if(++explCount > MAXRECURSIONDEPTH){
		//Report this incident
		//cerr << "Maximum recursion depth reached during extension. Position in q: " << qOff << endl;
		//Terminate this extension
		return suc;
	}

	//Check whether there are predecessors to continue
	if(predIter.hasPredecessors()){
		maxScore = score;
		//Set counter
		predCount = 1;

		//Check if there still exists an extension path to follow
		if(!extPth.empty()){
			//Get next predecessor
			nextPred = extPth.front();
			//Delete that predecessor from path
			extPth.pop_front();
			//We do not need to count how many successors we have explored yet
			explCount = 0;
		} else{
			//We have no predecessor given to go on with
			nextPred = 0;
		}

		//Explore each neighbor
		for(neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> nI = predIter.begin(); nI != predIter.end(); ++nI, ++predCount){
			//Check if we still have to wait for the correct predecessor to go on with
			if(predCount < nextPred) continue;

			struct Algn tmpAlgn;
			//Get next unitig
			predUni = *nI;
			tmpPosQ = qOff - 1;
			tmpPosU = predUni.size - predUni.getGraph()->getK();
			tmpScore = score;
			//Calculate gapped alignment on the predecessive unitig
			contLeftGappedAlignment(predUni, extPth, q, tmpPosQ, tmpPosU, X, maxGaps, tmpAlgn, tmpScore, explCount, quorum, searchSet);

			//Check for new maximum
			if(maxScore < tmpScore){
				//Update maximum
				maxScore = tmpScore;
				maxPosU = tmpPosU;
				maxPosQ = tmpPosQ;
				maxAlgn = tmpAlgn;
				uni = predUni;
				suc = true;
			}

			//If we reach this point and had a specific predecessor to go on with we are done
			if(nextPred > 0) break;
		}

		//Only update score and qOff if exploration was successful
		if(suc){
			//Give back best result found
			score = maxScore;
			qOff = maxPosQ;
			uOff = maxPosU;
			algn.aSeqG = maxAlgn.aSeqG + algn.aSeqG;
			algn.aSeqQ = maxAlgn.aSeqQ + algn.aSeqQ;
		}
	}

	return suc;
}

//This function calculates a gapped alignment to the right side of the starting position considering a quorum and a search color set. ATTENTION: Hit's length attribute will be deprecated after function call!
void startRightGappedAlignment(Hit *h, const string &q, const uint16_t &X, const uint32_t maxGaps, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	bool explSuc = false, unitigDone = true;
	uint32_t posU = h->offU, posQ = h->offQ, maxPosU = h->offU, maxPosQ = h->offQ, explCount = 0;
	int32_t maxScore = 0, maxBorderScore;
	struct Algn globAlgn;
	list<uint16_t> extPth = decmprExtPth(h->rExt);
	UnitigColorMap<seedlist> currUni = h->origUni;

	//Calculate gapped alignment
	unitigDone = calcSemiGlobAlignment(currUni, q, posU, posQ, X, maxGaps, maxPosQ, maxPosU, h->gAlgn, globAlgn, maxScore, maxBorderScore, quorum, searchSet, extPth.empty());

	//Check whether we have reached the end of the current unitig
	if(!unitigDone){
		explSuc = contGappedOnSameUni(currUni, extPth, q, posQ, posU, X, maxGaps, globAlgn, maxBorderScore, explCount, quorum, searchSet);
	} else if(maxBorderScore > -X && posQ < q.length() - 1){
		//Continue calculations on the next unitig
		explSuc = contGappedOnSuccUni(currUni, extPth, q, posQ, posU, X, maxGaps, globAlgn, maxBorderScore, explCount, quorum, searchSet);
	}

	//Check whether exploration of next unitig was successful
	if(explSuc && maxScore < maxBorderScore){
		//Save the achieved score
		h->score = maxBorderScore;
		//Save alignment
		h->gAlgn = globAlgn;
	} else{
		//Save the achieved score
		h->score = maxScore;
	}
}

//This function calculates a gapped alignment to the left side of the starting position considering a quorum and a search color set. ATTENTION: Hit's length attribute will be deprecated after function call!
void startLeftGappedAlignment(Hit *h, const string &q, const uint16_t &X, const uint32_t maxGaps, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	bool explSuc = false, unitigDone = true;
	uint32_t posU = h->offU, posQ = h->offQ, maxPosU = h->offU, maxPosQ = h->offQ, explCount = 0;
	int32_t maxScore = 0, maxBorderScore;
	list<uint16_t> extPth = decmprExtPth(h->lExt);
	struct Algn globAlgn = h->gAlgn;
	UnitigColorMap<seedlist> currUni = h->origUni;

	//Check whether we can start the calculation right on this unitig or whether we are already at the begin of unitig or query sequence
	if(posQ > 0 && posU > 0){
		//Set alignment start positions
		--posU;
		--posQ;
		//Calculate gapped alignment
		unitigDone = calcLeftGlobAlignment(currUni, q, posU, posQ, X, maxGaps, maxPosQ, maxPosU, h->gAlgn, globAlgn, maxScore, maxBorderScore, quorum, searchSet, extPth.empty());
	} else{
		//Set variable for exploration of predecessive unitigs
		maxBorderScore = 0;
	}

	//Check whether we have reached the beginning of the current unitig
	if(!unitigDone){
		explSuc = contLeftOnSameUni(currUni, extPth, q, posQ, posU, X, maxGaps, globAlgn, maxBorderScore, explCount, quorum, searchSet);
	} else if(maxBorderScore > -X && posQ > 0){
		//Continue calculations on the next unitig
		explSuc = contGappedOnPredUni(currUni, extPth, q, posQ, posU, X, maxGaps, globAlgn, maxBorderScore, explCount, quorum, searchSet);
	}

	//Check whether exploration of predecessive unitig was successful
	if(explSuc && maxScore < maxBorderScore){
		//Save new offsets in hit
		h->offU = posU;
		h->offQ = posQ;
		//Save unitig we end at
		h->origUni = currUni;
		//Save the achieved score
		h->score += maxBorderScore;
		//Save the better alignment
		h->gAlgn.aSeqG = globAlgn.aSeqG;
		h->gAlgn.aSeqQ = globAlgn.aSeqQ;
	} else{
		//Save new offsets in hit
		h->offU = maxPosU;
		h->offQ = maxPosQ;
		//Save the achieved score
		h->score += maxScore;
	}
}

//This function traces back the optimal path from a given start position within an edit matrix to get the corresponding alignment
void traceBack(uint32_t &matPosI, uint32_t posQ, const string &query, uint32_t &matPosJ, uint32_t posU, const string &uniSeq, int32_t **mat, struct Algn &algn){
	//We walk through the matrix until we have reached the upper, left corner
	while(matPosI != 0 || matPosJ != 0){
		//Check if we have reached the left or upper matrix border already
		if(matPosI == 0){
			//We can only walk on with a gap in the query sequence
			algn.aSeqG = uniSeq[posU] + algn.aSeqG;
			algn.aSeqQ = GAP_SYMB + algn.aSeqQ;
			//Move one cell up
			--matPosJ;
			--posU;
		} else if(matPosJ == 0){
			//We can only walk on with a gap in the unitig sequence
			algn.aSeqG = GAP_SYMB + algn.aSeqG;
			algn.aSeqQ = query[posQ] + algn.aSeqQ;
			//Move one cell to the left
			--matPosI;
			--posQ;
		} else{
			//Find the cell with the maximum score
			if(mat[matPosI - 1][matPosJ - 1] + compUScore(uniSeq[posU], query[posQ]) == mat[matPosI][matPosJ]){
				//We can add a base from both sequences
				algn.aSeqG = uniSeq[posU] + algn.aSeqG;
				algn.aSeqQ = query[posQ] + algn.aSeqQ;
				//Move to upper left cell
				--matPosI;
				--matPosJ;
				--posQ;
				--posU;
			} else if(mat[matPosI - 1][matPosJ] + GAP_SCORE == mat[matPosI][matPosJ]){
				//We can only walk on with a gap in the unitig sequence
				algn.aSeqG = GAP_SYMB + algn.aSeqG;
				algn.aSeqQ = query[posQ] + algn.aSeqQ;
				//Move one cell to the left
				--matPosI;
				--posQ;
			} else{
				//We can only walk on with a gap in the query sequence
				algn.aSeqG = uniSeq[posU] + algn.aSeqG;
				algn.aSeqQ = GAP_SYMB + algn.aSeqQ;
				//Move one cell up
				--matPosJ;
				--posU;
			}
		}
	}
}

//This function traces back the optimal path of a left gapped extension from a given start position within an edit matrix to get the corresponding alignment
void traceForth(uint32_t &matPosI, uint32_t posQ, const string &query, uint32_t &matPosJ, uint32_t posU, const string &uniSeq, int32_t **mat, struct Algn &algn){
	// int32_t mScore;
	string gSeq, qSeq;

	//We walk through the matrix until we have reached the upper, left corner
	while(matPosI != 0 || matPosJ != 0){
		//Check if we have reached the left or upper matrix border already
		if(matPosI == 0){
			//We can only walk on with a gap in the query sequence
			gSeq = gSeq + uniSeq[posU];
			qSeq = qSeq + GAP_SYMB;
			//Move one cell up
			--matPosJ;
			++posU;
		} else if(matPosJ == 0){
			//We can only walk on with a gap in the unitig sequence
			gSeq = gSeq + GAP_SYMB;
			qSeq = qSeq + query[posQ];
			//Move one cell to the left
			--matPosI;
			++posQ;
		} else{
			//Find the cell with the maximum score
			if(mat[matPosI - 1][matPosJ - 1]  + compUScore(uniSeq[posU], query[posQ]) == mat[matPosI][matPosJ]){
				//We can add a base from both sequences
				gSeq = gSeq + uniSeq[posU];
				qSeq = qSeq + query[posQ];
				//Move to upper left cell
				--matPosI;
				--matPosJ;
				++posQ;
				++posU;
			} else if(mat[matPosI - 1][matPosJ] + GAP_SCORE == mat[matPosI][matPosJ]){
				//We can only walk on with a gap in the unitig sequence
				gSeq = gSeq + GAP_SYMB;
				qSeq = qSeq + query[posQ];
				//Move one cell to the left
				--matPosI;
				++posQ;
			} else{
				//We can only walk on with a gap in the query sequence
				gSeq = gSeq + uniSeq[posU];
				qSeq = qSeq + GAP_SYMB;
				//Move one cell up
				--matPosJ;
				++posU;
			}
		}
	}

	algn.aSeqG = gSeq + algn.aSeqG;
	algn.aSeqQ = qSeq + algn.aSeqQ;
}