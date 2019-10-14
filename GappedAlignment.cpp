#include "GappedAlignment.h"
#include "Hit.h"
#include "Search.h"

// bool report = false;

// //This function calculates a semi-global (actually its rather half-global) alignment of a unitig and the query sequence considering a quorum. Returns true if calculations on unitig could be finished (i.e. the unitig sequence was not too long to be stored inside the edit matrix).
// bool calcSemiGlobAlignment(const UnitigColorMap<seedlist> &uni, const string &q, uint32_t& posU, uint32_t& posQ, const uint16_t &X, const uint32_t &maxGaps, uint32_t& maxPosQ, uint32_t& maxPosU, struct Algn &maxAlgn, struct Algn &brdAlgn, int32_t& maxScore, int32_t& eMax, const uint32_t &quorum){//TODO: So far we do not adjust the maximal number of gaps. This has to be changed as soon as we calculate real alignments!//TODO: This function does more than it should!
// 	bool termCalcs, matShriked = false;
// 	uint32_t matBrth, matHgth, i, j, lCalcBorder, rCalcBorder, endBuf = 0, uCmpSeqLen, qCmpSeqLen, maxMatPosI = 0, maxMatPosJ = 0, eMaxPosU = posU, eMaxPosQ = posQ, edgeMaxPosI = 0, edgeMaxPosJ = 0;
// 	string uSeq;

// 	//Initialize the current max position
// 	eMax = -X;
// 	//Get unitig's sequence
// 	uSeq = uni.mappedSequenceToString();

// 	//Check whether current unitig has a successor
// 	if(uni.getSuccessors().hasSuccessors()){
// 		endBuf = uni.getGraph()->getK() - 1;
// 	}

// 	//Testing
// 	// cout << "maxGaps:" << maxGaps << endl;
// 	// cout << "uSeq:" << uSeq << endl;

// 	//Calculate length of unitig sequence to compare
// 	uCmpSeqLen = uni.size - (posU + endBuf);

// 	//Testing
// 	// cout << "uCmpSeqLen:" << uCmpSeqLen << endl;

// 	//TODO This is just preliminar and has to be changed a.s.a.p.!
// 	//Check whether the unitig sequence which has to be compared fulfills the quorum
// 	if(!quorumFulfilled(uni, posU, uni.len, quorum)){
// 		//If the quorum is not fulfilled we are not going to do anything!
// 		uCmpSeqLen = 0;
// 	}

// 	//Calculate length of query sequence to compare
// 	qCmpSeqLen = q.length() - posQ;

// 	//Testing
// 	// cout << "qCmpSeqLen:" << qCmpSeqLen << endl;
// 	// cout << "posQ:" << posQ << endl;

// 	//Make sure that the matrix is not becoming too big
// 	if(uCmpSeqLen > MAX_MATRIX_SIZE){
// 		//Testing
// 		//cout << "Do we get here?" << endl;

// 		//If we reach the lower matrix border we definitely do not miss any sequence of the unitig
// 		matBrth = min((uint32_t) MAX_MATRIX_SIZE, qCmpSeqLen + maxGaps) + 1;
// 		//Do not matter which matrix border we reach and the matrix should not become too big
// 		matHgth = min(qCmpSeqLen, (uint32_t) MAX_MATRIX_SIZE) + 1;
// 		//If we still have uncompared query sequence mark that we are not done
// 		if(matBrth == (uint32_t) MAX_MATRIX_SIZE) matShriked = true;
// 	} else{
// 		//Consider as much unitig sequence as necessary but not more
// 		matHgth = min(uCmpSeqLen + maxGaps, qCmpSeqLen) + 1;
// 		//Make sure we cannot reach the lower matrix border except if our calculations are done
// 		matBrth = min(qCmpSeqLen + maxGaps, uCmpSeqLen) + 1;

// 	}

// 	//Testing
// 	// cout << "matBrth:" << matBrth << endl;

// 	//Initialize dynamic programming matrix

// 	// int32_t mat[dim][dim];

// 	int32_t **mat = (int32_t**) malloc(matHgth * sizeof(int32_t*));//TOOD: Is this the most efficient way to handle this?
// 	for(i = 0; i < matHgth; ++i) mat[i] = (int32_t*) malloc(matBrth * sizeof(int32_t));

// 	//Initialize upper left cell
// 	mat[0][0] = maxScore;

// 	//If one of the to compare sequences is the empty sequence we need to set eMax correctly because we are not entering the main for loop at all
// 	if(matBrth == 1) eMax = mat[0][0];

// 	//Testing
// 	// cout << "First cell of matrix is initialized" << endl;

// 	//Initialize first column
// 	for(i = 1; i < matHgth; ++i){
// 		if(i <= maxGaps){
// 			mat[i][0] = maxScore + i * GAP_SCORE;//GAP_SCORE is supposed to be negative//TODO: This has to be changed as soon as we are not using unit score anymore!
// 		} else{
// 			mat[i][0] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
// 			//We want to save time by not treating the whole matrix
// 			break;
// 		}

// 		//Testing
// 		// cout << "We initialized position " << i << endl;
// 	}

// 	//Testing
// 	// cout << "First matrix column is initialized" << endl;

// 	//Initialize first row
// 	for(j = 1; j < matBrth; ++j){
// 		if(j <= maxGaps){
// 			mat[0][j] = maxScore + j * GAP_SCORE;//GAP_SCORE is supposed to be negative//TODO: This has to be changed as soon as we are not using unit score anymore!
// 		} else{
// 			mat[0][j] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
// 			//We want to save time by not treating the whole matrix
// 			break;
// 		}
// 	}

// 	//Set left border
// 	lCalcBorder = 1;

// 	// rCalcBorder = min(matBrth, maxGaps);

// 	//Testing
// 	// cout << "STarte Matrixberechnung" << endl;
// 	// cout << "posQ: " << posQ << " posU: " << posU << endl; 
// 	// cout << "Last characters in this matrix are: q:" << q[posQ + dim - 2] << " u:" << uSeq[posU + dim - 2] << endl;

// 	//Fill the matrix
// 	for(i = 1; i < matHgth; ++i){
// 		//Reset abort var
// 		termCalcs = true;

// 		//Check if we can adjust left alignment calculation border
// 		if(maxGaps < i){
// 			lCalcBorder = i - maxGaps;
// 			//Initialize cell to the left of the first one to calculate as we will need it there
// 			mat[i][lCalcBorder - 1] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
// 		}

// 		//Adjust right border
// 		rCalcBorder = min(i + maxGaps, matBrth - 1);
// 		//If existing, initialize cell to the right of the last calculated one as we might need it for calculations in the next iteration
// 		if(rCalcBorder < matBrth - 1) mat[i][rCalcBorder + 1] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!

// 		for(j = lCalcBorder; j <= rCalcBorder; ++j){
// 			//Testing
// 			// cout << "lCalcBorder:" << lCalcBorder << " rCalcBorder:" << rCalcBorder << endl;

// 			// //Make sure we calculate a bended alignment if matrix is too big
// 			// if(j < lCalcBorder || j > rCalcBorder){
// 			// 	//Initialize cell to the right of the last calculated one as we might need it for calculations in the next iteration
// 			// 	mat[i][j] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
// 			// 	//We want to save time by not treating the whole matrix
// 			// 	break;
// 			// } else{

// 			//Testing
// 			// cout << "We compare mat[i - 1][j - 1] + compUScore(q[posQ + i - 1], uSeq[posU + j - 1]):" << mat[i - 1][j - 1] + compUScore(q[posQ + i - 1], uSeq[posU + j - 1]) << " mat[i][j - 1] - 1:" << mat[i][j - 1] - 1 << " mat[i - 1][j] - 1:" << mat[i - 1][j] - 1 << endl;

// 			//Calculate current cells value and check we don't drop
// 			if((mat[i][j] = max(mat[i - 1][j - 1] + compUScore(q[posQ + i - 1], uSeq[posU + j - 1]), max(mat[i][j - 1] + GAP_SCORE, mat[i - 1][j] + GAP_SCORE))) > -X){//TODO: This has to be changed as soon as we are not using unit score anymore!
// 				//Testing
// 				// cout << "mat[i][j]:" << mat[i][j] << endl;

// 				//Make sure we do not terminate our calculations too soon
// 				termCalcs = false;

// 				//Check whether we have a new max score (NOTE: We could check only for better scores here as well. Do it this way will in doubt lead to larger alignments. <- In favour for alignments with smaller E-values, I changed this)
// 				if(maxScore < mat[i][j]){
// 					//Update max score
// 					maxScore = mat[i][j];
// 					//Update maximum position
// 					maxPosQ = posQ + i - 1;
// 					maxPosU = posU + j - 1;
// 					//Save maximum position in matrix
// 					maxMatPosI = i;
// 					maxMatPosJ = j;

// 					//Testing
// 					//cout << "maxPosQ inside calcSemiGlobAlignment: " << maxPosQ << endl;
// 				}

// 				//Testing
// 				//cout << mat[i][j] << endl;

// 				//Check whether we are at the right or lower matrix edge and whether we have found a new maximum
// 				if((i == matHgth - 1 || j == matBrth - 1) && eMax <= mat[i][j]){
// 					//Update edge maximum
// 					eMax = mat[i][j];
// 					eMaxPosQ = posQ + i - 1;
// 					eMaxPosU = posU + j - 1;
// 					//Save edge maximum in matrix
// 					edgeMaxPosI = i;
// 					edgeMaxPosJ = j;
// 				}
// 			}
// 		}

// 		//Check whether calculations may be terminated
// 		if(termCalcs) break;
// 	}

// 	//Testing
// 	// cout << "Alignment matrix:" << endl;
// 	// for(uint32_t k = 0; k < matHgth; ++k){
// 	// 	for(uint32_t l = 0; l < matBrth; ++l){
// 	// 		cout << mat[k][l] << " ";
// 	// 	}
// 	// 	cout << endl;
// 	// }
// 	// cout << "i:" << i << " posQ:" << posQ << endl;
// 	// if(counter == 1) exit(0);
// 	// ++counter;
// 	// cout << "maxPosQ:" << maxPosQ << "maxPosU:" << maxPosU << endl;

// 	//Do backtracing to get the alignment maximum scoring alignment
// 	traceBack(maxMatPosI, maxPosQ, q, maxMatPosJ, maxPosU, uSeq, mat, maxAlgn);//TODO Implement this function!

// 	//Testing
// 	// cout << "Backtracing done" << endl;

// 	//Check if maximum score is at the matrix border
// 	if(maxPosU == eMaxPosU && maxPosQ == eMaxPosQ){
// 		//Copy maximum scoring alignment to save time
// 		brdAlgn = maxAlgn;
// 	} else{
// 		//Backtrace alignment to go on with on next unitig
// 		traceBack(edgeMaxPosI, eMaxPosQ, q, edgeMaxPosJ, eMaxPosU, uSeq, mat, brdAlgn);
// 	}

// 	//Testing
// 	// cout << "are we already here?" << endl;

// 	//Free matrix
// 	for(i = 0; i < matHgth; ++i) free(mat[i]);
// 	free(mat);

// 	//Testing
// 	// cout << "Backtracing done" << endl;

// 	//Update positions in unitig und q
// 	if(eMax > -X){
// 		posU = eMaxPosU;
// 		posQ = eMaxPosQ;
// 	} else{
// 		posU = maxPosU;
// 		posQ = maxPosQ;

// 		//Testing
// 		// cout << "posQ is newly set here and maxPosQ is " << maxPosQ << endl;
// 	}

// 	//Testing
// 	// cout << "posQ in the end of calcSemiGlobAlignment:" << posQ << endl;
// 	// cout << ":" << maxAlgn.aSeqG.length() << endl;

// 	if(eMax > -X && matShriked) return false;
	
// 	return true;
// }

//This function calculates a semi-global alignment of a unitig and the query sequence considering a quorum and a search color set. Returns true if calculations on unitig could be finished (i.e. the unitig sequence was not too long to be stored inside the edit matrix).
bool calcSemiGlobAlignment(const UnitigColorMap<seedlist> &uni, const string &q, uint32_t& posU, uint32_t& posQ, const uint16_t &X, const uint32_t &maxGaps, uint32_t& maxPosQ, uint32_t& maxPosU, struct Algn &maxAlgn, struct Algn &brdAlgn, int32_t& maxScore, int32_t& eMax, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool srchCritCheck){//TODO: So far we do not adjust the maximal number of gaps. This has to be changed as soon as we calculate real alignments!//TODO: This function does more than it should!
	bool termCalcs, matShrunk = false;
	uint32_t matBrth, matHgth, i, j, lCalcBorder, rCalcBorder, endBuf = 0, uCmpSeqLen, qCmpSeqLen, maxMatPosI = 0, maxMatPosJ = 0, eMaxPosU = posU, eMaxPosQ = posQ, edgeMaxPosI = 0, edgeMaxPosJ = 0;
	int32_t covPos;
	string uSeq;

	//Testing
	// cout << "calcSemiGlobAlignment: Start of function" << endl;

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

	//Testing
	//cout << "8 Option ";
	// cout << "Do the actual calculations..." << endl;

	//Do a search criteria check only if necessary
	if(srchCritCheck){
		//Testing
		// cout << "9 Option 1" << endl;
		// cout << "quorum: " << quorum << endl;
		// cout << "searchSet is " << (searchSet.empty() ? "" : "not ") << "empty" << endl;
		// cout << "uni: " << uni.mappedSequenceToString() << endl;
		// cout << "posU: " << posU << endl;

		//Get the number of positions fulfilling our search crit on this unitig reduced by our starting offset (-1, because we do not have to compare the offset position itself anymore)
		covPos = getSrchCritCov(uni, quorum, searchSet, posU, true);// - posU - 1;

		//Testing
		// cout << "covPos: " << covPos << endl;
		// cout << "Got number of covered positions" << endl;

		//Check if there is at least one position fulfilling the search crit
		if(covPos <= 0){
			//Testing
			// cout << "10 Option 2" << endl;

			//We don't need to proceed
			uCmpSeqLen = 0;
		} else if((uint32_t) covPos < uCmpSeqLen){
			//Testing
			// cout << "10 Option 1" << endl;

			//We only need to calculate until our search crits are not fulfilled anymore
			uCmpSeqLen = covPos;
		}	
	}// else{
	// 	//Testing
	// 	 cout << "9 Option 2" << endl;
	// }

	// //Check whether the unitig sequence which has to be compared fulfills the quorum
	// if(!quorumFulfilled(uni, posU, uni.len, quorum, searchSet)){
	// 	//Testing
	// 	//cout << "not ";

	// 	//If the quorum is not fulfilled we are not going to do anything!
	// 	uCmpSeqLen = 0;
	// }

	//Calculate length of query sequence to compare
	qCmpSeqLen = q.length() - posQ;

	// //Calculate dynamic programming matrix borders
	// dim = min(uCmpSeqLen, qCmpSeqLen) + 1;

	//Testing
	//cout << "dim: " << dim << endl;
	//cout << "1" << endl;
	// cout << "Search criteria done" << endl;

	//Make sure that the matrix is not becoming too big
	if(uCmpSeqLen > MAX_MATRIX_SIZE){
		//Testing
		//cout << "Do we get here?" << endl;

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

	//Testing
	// if(report){
	// 	cout << "Inside calcSemiGlobAlignment:  posU: " << posU << " endBuf: " << endBuf << " unitig size: " << uni.size << " uCmpSeqLen: " << uCmpSeqLen << endl;
	// }
	// cout << "qCmpSeqLen + maxGaps: " << qCmpSeqLen + maxGaps << endl;

	//Initialize dynamic programming matrix
	int32_t **mat = (int32_t**) malloc(matHgth * sizeof(int32_t*));//TOOD: Is this the most efficient way to handle this?
	for(i = 0; i < matHgth; ++i) mat[i] = (int32_t*) malloc(matBrth * sizeof(int32_t));

	// int32_t mat[dim][dim];

	//Initialize upper left cell
	mat[0][0] = maxScore;

	// //If one of the two compared sequences is the empty sequence we need to set eMax correctly because we are not entering the main for loop at all
	// if(matBrth == 1){
	// 	//Testing
	// 	// cout << "11 Option 1" << endl;

	// 	eMax = mat[0][0];
	// }// else{
	// 	//Testing
	// 	cout << "11 Option 2" << endl;
	// }

	//Initialize first column
	for(i = 1; i < matHgth; ++i){
		if(i <= maxGaps){
			//GAP_SCORE is supposed to be negative
			mat[i][0] = maxScore + i * GAP_SCORE;//TODO: This has to be changed as soon as we are not using unit score anymore!
		} else{
			mat[i][0] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
			//We want to save time by not treating the whole matrix
			break;
		}
	}

	//Initialize first row
	for(j = 1; j < matBrth; ++j){
		if(j <= maxGaps){
			mat[0][j] = maxScore + j * GAP_SCORE;//GAP_SCORE is supposed to be negative//TODO: This has to be changed as soon as we are not using unit score anymore!
		} else{
			mat[0][j] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
			//We want to save time by not treating the whole matrix
			break;
		}
	}

	//Set borders
	lCalcBorder = 1;

	// rCalcBorder = min(dim, maxGaps);

	//Testing
	// cout << "Starte Matrixberechnung" << endl;

	//Fill the matrix
	for(i = 1; i < matHgth; ++i){
		//Reset abort var
		termCalcs = true;

		//Check if we can adjust left alignment calculation border
		if(maxGaps < i){
			lCalcBorder = i - maxGaps;
			//Initialize cell to the left of the first one to calculate as we will need it there
			mat[i][lCalcBorder - 1] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
		}

		//Adjust right border
		rCalcBorder = min(i + maxGaps, matBrth - 1);
		//If existing, initialize cell to the right of the last calculated one as we might need it for calculations in the next iteration
		if(rCalcBorder < matBrth - 1) mat[i][rCalcBorder + 1] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!

		for(j = lCalcBorder; j <= rCalcBorder; ++j){

			// //Make sure we calculate a bended alignment if matrix is too big
			// if(j < lCalcBorder || j > rCalcBorder){
			// 	mat[i][j] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
			// } else{

			//Calculate current cells value and check we don't drop
			if((mat[i][j] = max(mat[i - 1][j - 1] + compUScore(q[posQ + i - 1], uSeq[posU + j - 1]), max(mat[i][j - 1] + GAP_SCORE, mat[i - 1][j] + GAP_SCORE))) > -X){//TODO: This has to be changed as soon as we are not using unit score anymore!
				//Make sure we do not terminate our calculations too soon
				termCalcs = false;

				//Check whether we have a new max score (NOTE: We could check only for better scores here as well. Do it this way will in doubt lead to larger alignments. <- In favour for alignments with smaller E-values, I changed this)
				if(maxScore < mat[i][j]){
					//Update max score
					maxScore = mat[i][j];
					//Update maximum position
					maxPosQ = posQ + i - 1;
					maxPosU = posU + j - 1;
					//Save maximum position in matrix
					maxMatPosI = i;
					maxMatPosJ = j;

					//Testing
					//cout << "maxPosQ inside calcSemiGlobAlignment: " << maxPosQ << endl;
				}

				//Testing
				//cout << mat[i][j] << endl;

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

	//Testing
	// if(report){
		// cout << "Inside calcSemiGlobAlignment: Matrix filled" << endl;
	// 	cout << "Alignment matrix:" << endl;
	// 	for(i = 0; i < matHgth; ++i){
	// 		for(j = 0; j < matBrth; ++j){
	// 			cout << mat[i][j] << " ";
	// 		}
	// 		cout << endl;
	// 	}
	// // 	cout << "maxMatPosI: " << maxMatPosI << " maxMatPosJ: " << maxMatPosJ << " maxPosU: " << maxPosU << " maxPosQ: " << maxPosQ << endl;
	// // 	cout << "edgeMaxPosI: " << edgeMaxPosI << " edgeMaxPosJ: " << edgeMaxPosJ << " eMaxPosU: " << eMaxPosU << " eMaxPosQ: " << eMaxPosQ << endl;
 // 	}

	//Do backtracing to get the alignment maximum scoring alignment
	traceBack(maxMatPosI, maxPosQ, q, maxMatPosJ, maxPosU, uSeq, mat, maxAlgn);

	//Check if maximum score is at the matrix border
	if(maxMatPosI == edgeMaxPosI && maxMatPosJ == edgeMaxPosJ){
		//Copy maximum scoring alignment to save time
		brdAlgn = maxAlgn;
	} else{
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

//This function calculates a global, banded alignment of the reverse of two sequences. Calculations end as soon as the end of one of the sequences has been reached or the score drops below a specified value X. Returns true if a calculation was successful.
/*bool calcLeftGlobAlignment(const string hseq, const string vseq, const uint16_t& X, const uint32_t& bandRadius, uint32_t& maxPos_hseq, uint32_t& maxPos_vseq, int32_t& maxScore, uint32_t& ePos_hseq, uint32_t& ePos_vseq, int32_t& eMaxScore){
	bool 
}*/

//This function calculates a global, bended alignment of the reverse of a unitig and the query sequence. Returns true if calculations on this unitig could be finished (i.e. the unitig sequence was not to long to be stored inside the edit matrix).
/*bool calcLeftGlobAlignment(const UnitigMap<seedlist> &uni, const string &q, uint32_t& startPosU, uint32_t& startPosQ, const uint16_t &X, const uint32_t &maxGaps, uint32_t& maxPosQ, uint32_t& maxPosU, int32_t& maxScore, int32_t& eMax){
	bool termCalcs, matShriked = false;
	uint32_t dim, i, j, lCalcBorder, rCalcBorder, uCmpSeqLen, qCmpSeqLen, eMaxPosU = 0, eMaxPosQ = 0;
	string uSeq;

	//Get unitig's sequence
	uSeq = uni.toString();

	//Check whether current unitig has a predecessor
	if(uni.{
		endBuf = uni.getGraph()->getK() - 1;
	}
}*/

//This function calculates a banded alignment of a unitig and the query sequence to the left considering a quorum. Returns true if calculations on unitig could be finished (i.e. the unitig sequence was not to long to be stored inside the edit matrix).
bool calcLeftGlobAlignment(const UnitigColorMap<seedlist> &uni, const string &q, uint32_t& posU, uint32_t& posQ, const uint16_t &X, const uint32_t &maxGaps, uint32_t& maxPosQ, uint32_t& maxPosU, struct Algn &maxAlgn, struct Algn &brdAlgn, int32_t& maxScore, int32_t& eMax, const uint32_t &quorum){//TODO: So far we do not adjust the maximal number of gaps. This has to be changed as soon as we calculate real alignments!//TODO: This function does more than it should!
	bool termCalcs, matShriked = false;
	uint32_t matBrth, matHgth, i, j, lCalcBorder, rCalcBorder, maxMatPosI = 0, maxMatPosJ = 0, eMaxPosU = posU, eMaxPosQ = posQ, edgeMaxPosI = 0, edgeMaxPosJ = 0;
	string uSeq;

	//Get unitig's sequence
	uSeq = uni.mappedSequenceToString();

	//Testing
	// cout << "Unitig sequence to compare:" << uSeq << " from pos " << posU << endl;
	// cout << "Query sequence to compare: " << q << " from pos " << posQ << endl;

	//Calculate dynamic programming matrix borders
	// matHgth = min(posU + maxGaps, posQ) + 2;//+2 since we need 1 column and row for the empty string and offset indexes start with 0
	// matBrth = min(pos)<-
	//dim = min(posU, posQ) + 2;//+2 since we need 1 column and row for the empty string and offset indexes start with 0

	//Make sure that the matrix is not becoming too big
	if(posU > MAX_MATRIX_SIZE){
		//If we reach the lower matrix border we definitely do not miss any sequence of the unitig
		matBrth = min((uint32_t) MAX_MATRIX_SIZE, posQ + maxGaps) + 2;//+2 since we need 1 column and row for the empty string and offset indexes start with 0
		//Do not matter which matrix border we reach and the matrix should not become too big
		matHgth = min(posQ, (uint32_t) MAX_MATRIX_SIZE) + 2;
		matShriked = true;
	} else{
		//Make sure we cannot reach the lower matrix border except if our calculations are done
		matBrth = min(posU, posQ + maxGaps) + 2;
		//Consider as much unitig sequence as necessary but not more
		matHgth = min(posQ, posU + maxGaps) + 2;
	}

	//TODO This is just preliminar and has to be changed a.s.a.p.!
	//Check whether the unitig sequence which has to be compared fulfills the quorum
	if(!quorumFulfilled(uni, 0, posU + 1, quorum)){
		//If the quorum is not fulfilled we are not going to do anything!
		matHgth = 1;
		matBrth = 1;
	}

	//Initialize dynamic programming matrix
	int32_t **mat = (int32_t**) malloc(matHgth * sizeof(int32_t*));//TOOD: Is this the most efficient way to handle this?
	for(i = 0; i < matHgth; ++i) mat[i] = (int32_t*) malloc(matBrth * sizeof(int32_t));

	//Testing
	// cout << "We are getting here" << endl;
	// cout << "matHgth:" << matHgth << " matBrth:" << matBrth << endl;

	//Initialize upper left cell
	mat[0][0] = maxScore;

	//Testing
	// cout << "mat[0][0]:" << mat[0][0] << endl;
	
	//If one of the to compare sequences is the empty sequence we need to set eMax correctly because we are not entering the main for loop at all
	if(matBrth == 1){
		eMax = mat[0][0];
	} else{
		//Initialize edge maximum
		eMax = -X;
	}

	//Initialize first column
	for(i = 1; i < matHgth; ++i){
		if(i <= maxGaps){
			mat[i][0] = maxScore + i * GAP_SCORE;//GAP_SCORE is supposed to be negative//TODO: This has to be changed as soon as we are not using unit score anymore!
		} else{
			mat[i][0] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
			//We want to save time by not treating the whole matrix
			break;
		}
	}

	//Testing
	// cout << "But do we get here?" << endl;

	//Initialize first row
	for(j = 1; j < matBrth; ++j){
		if(j <= maxGaps){
			mat[0][j] = maxScore + j * GAP_SCORE;//GAP_SCORE is supposed to be negative//TODO: This has to be changed as soon as we are not using unit score anymore!
		} else{
			mat[0][j] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
		}
	}

	//Set left border
	lCalcBorder = 1;
	//rCalcBorder = min(dim, maxGaps);

	//Fill the matrix
	for(i = 1; i < matHgth; ++i){
		//Reset abort var
		termCalcs = true;

		//Check if we can adjust left alignment calculation border
		if(maxGaps < i){
			lCalcBorder = i - maxGaps;
			//Initialize cell to the left of the first one to calculate as we will need it there
			mat[i][lCalcBorder - 1] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
		}

		//Adjust right border
		rCalcBorder = min(i + maxGaps, matBrth - 1);
		//If existing, initialize cell to the right of the last calculated one as we might need it for calculations in the next iteration
		if(rCalcBorder < matBrth - 1) mat[i][rCalcBorder + 1] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!

		for(j = lCalcBorder; j <= rCalcBorder; ++j){
			// //Make sure we calculate a banded alignment if matrix is too big
			// if(j < lCalcBorder || j > rCalcBorder){
			// 	mat[i][j] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
			// } else{

			//Calculate current cells value and check we don't drop
			if((mat[i][j] = max(mat[i - 1][j - 1] + compUScore(q[posQ - i + 1], uSeq[posU - j + 1]), max(mat[i][j - 1] + GAP_SCORE, mat[i - 1][j] + GAP_SCORE))) > -X){//TODO: This has to be changed as soon as we are not using unit score anymore!
				//Make sure we do not terminate our calculations too soon
				termCalcs = false;

				//Testing
				//cout << "Inside calcLeftGlobAlignment: Characters that have been compared were " << q[posQ - i + 1] << " in the query and " << uSeq[posU - j + 1] << " in the unitig" << endl;

				//Check whether we have a new max score (NOTE: We could check only for better scores here as well. Do it this way will in doubt lead to larger alignments.<- In favour for alignments with smaller E-values, I changed this)//TODO: Maybe we don't what this because this might lead to less significant alignments!
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

			//Testing
			// cout << "mat[" << i << "][" << j << "]:" << mat[i][j] << " ";
		}

		//Testing
		// cout << endl;

		//Check whether calculations may be terminated
		if(termCalcs) break;
	}

	//Testing
	// cout << "Filled matrix:" << endl;
	// for(i = 0; i < matHgth; ++i){
	// 	for(j = 0; j < matBrth; ++j){
	// 		cout << mat[i][j] << " ";
	// 	}
	// 	cout << endl;
	// }
	// cout << "Do forthtracing" << endl;

	//Do backtracing to get the alignment maximum scoring alignment
	traceForth(maxMatPosI, maxPosQ, q, maxMatPosJ, maxPosU, uSeq, mat, maxAlgn);

	//Check if maximum score is at the matrix border
	if(maxPosU == eMaxPosU && maxPosQ == eMaxPosQ){
		//Copy maximum scoring alignment to save time
		brdAlgn = maxAlgn;
	} else{
		//Backtrace alignment to go on with on next unitig
		traceForth(edgeMaxPosI, eMaxPosQ, q, edgeMaxPosJ, eMaxPosU, uSeq, mat, brdAlgn);
	}

	//Testing
	// cout << "Free the matrix" << endl;

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

	if(eMax > -X && matShriked) return false;
	
	return true;
}

//This function calculates a banded alignment of a unitig and the query sequence to the left considering a quorum and a search color set. Returns true if calculations on unitig could be finished (i.e. the unitig sequence was not to long to be stored inside the edit matrix).
bool calcLeftGlobAlignment(const UnitigColorMap<seedlist> &uni, const string &q, uint32_t& posU, uint32_t& posQ, const uint16_t &X, const uint32_t &maxGaps, uint32_t& maxPosQ, uint32_t& maxPosU, struct Algn &maxAlgn, struct Algn &brdAlgn, int32_t& maxScore, int32_t& eMax, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool srchCritCheck){//TODO: So far we do not adjust the maximal number of gaps. This has to be changed as soon as we calculate real alignments!//TODO: This function does more than it should!
	bool termCalcs, matShrunk = false;
	uint32_t matBrth, matHgth, i, j, lCalcBorder, rCalcBorder, maxMatPosI = 0, maxMatPosJ = 0, eMaxPosU = posU, eMaxPosQ = posQ, edgeMaxPosI = 0, edgeMaxPosJ = 0, uCmpSeqLen = posU;
	int32_t covPos = 0;
	string uSeq;

	//Get unitig's sequence
	uSeq = uni.mappedSequenceToString();

	//Calculate dynamic programming matrix borders
	// dim = min(posU, posQ) + 2;//+2 since we need 1 column and row for the empty string and offset indexes start with 0

	//Testing
	// cout << "1 Option";
	// if(report){
		// cout << "Inside calcLeftGlobAlignment: posU: " << posU << " posQ: " << posQ << " unitig is " << uni.mappedSequenceToString() << endl;
	// 	UnitigColorMap<seedlist> bla = uni.getGraph()->find(Kmer(uni.mappedSequenceToString().c_str()));
	// 	cout << "This unitig can " << (bla.isEmpty ? "not " : "") << "be found in the graph" << endl;
	// 	// exit(0);
	// }
	
	//Do we have to check search criteria on this unitig?
	if(srchCritCheck){
		//Testing
		// cout << "2 Option 1" << endl;

		covPos = getSrchCritCov(uni, quorum, searchSet, posU, false);// - (uni.size - posU - 1);

		//Check if there is at least one position fulfilling the search crit
		if(covPos > 0 && (uint32_t) covPos <= uCmpSeqLen){
			//Testing
			// cout << "3 Option 1" << endl;
			// cout << "covPos: " << covPos << endl;
			// if(uSeq == "CTCACCTTGGTTGAAAAATGCTAAAAGTGCCCCCTCAATACATTTATCAAGATATTGGTCT"){
			// 	cout << "covPos: " << covPos << endl;
			// }

			//We only need to calculate until our search crits are not fulfilled anymore (-1, because offsets start at 0)
			uCmpSeqLen = (uint32_t) covPos - 1;
		}

		// 	//Testing
		// cout << "3 Option 2" << endl;
		// }
	}// else{
	// 	//Testing
	// 	cout << "2 Option 2" << endl;
	// }

	//Testing
	// cout << "uCmpSeqLen: " << uCmpSeqLen << endl;
	// if(uSeq == "CTCACCTTGGTTGAAAAATGCTAAAAGTGCCCCCTCAATACATTTATCAAGATATTGGTCT"){
	// 	cout << "maxAlgn before calculations on this unitig: " << maxAlgn.aSeqG << endl;
	// 	// exit(0);	
	// }

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

	// //Check whether the unitig sequence which has to be compared fulfills the quorum
	// if(!quorumFulfilled(uni, 0, posU + 1, quorum, searchSet)){
	// 	//Testing
	// 	// cout << "not ";

	// 	//If the search criteria is not fulfilled we are not going to do anything!
	// 	matHgth = 1;
	// 	matBrth = 1;
	// }

	//Testing
	// cout << "1" << endl;

	// //Make sure that the matrix is not becoming too big
	// if(dim > MAX_MATRIX_SIZE){
	// 	dim = MAX_MATRIX_SIZE;
	// 	matShrunk = true;
	// }

	//Initialize dynamic programming matrix
	int32_t **mat = (int32_t**) malloc(matHgth * sizeof(int32_t*));//TOOD: Is this the most efficient way to handle this?
	for(i = 0; i < matHgth; ++i) mat[i] = (int32_t*) malloc(matBrth * sizeof(int32_t));

	// int32_t mat[dim][dim];

	//Initialize upper left cell
	mat[0][0] = maxScore;

	// //If one of the to compare sequences is the empty sequence we need to set eMax correctly because we are not entering the main for loop at all
	// if(matBrth == 1){
	// 	eMax = mat[0][0];
	// } else{

	//Initialize edge maximum
	eMax = -X;

	// }

	//Initialize first column
	for(i = 1; i < matHgth; ++i){
		if(i <= maxGaps){
			mat[i][0] = maxScore + i * GAP_SCORE;//GAP_SCORE is supposed to be negative//TODO: This has to be changed as soon as we are not using unit score anymore!
		} else{
			mat[i][0] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
		}
	}

	//Initialize first row
	for(j = 1; j < matBrth; ++j){
		if(j <= maxGaps){
			mat[0][j] = maxScore + j * GAP_SCORE;//GAP_SCORE is supposed to be negative//TODO: This has to be changed as soon as we are not using unit score anymore!
		} else{
			mat[0][j] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
		}
	}

	//Set left border
	lCalcBorder = 1;

	//rCalcBorder = min(dim, maxGaps);

	//Testing
	// cout << "Start main calculations" << endl;

	//Fill the matrix
	for(i = 1; i < matHgth; ++i){
		//Reset abort var
		termCalcs = true;

		//Check if we can adjust left alignment calculation border
		if(maxGaps < i){
			lCalcBorder = i - maxGaps;
			//Initialize cell to the left of the first one to calculate as we will need it there
			mat[i][lCalcBorder - 1] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
		}

		//Adjust right border
		rCalcBorder = min(i + maxGaps, matBrth - 1);
		//If existing, initialize cell to the right of the last calculated one as we might need it for calculations in the next iteration
		if(rCalcBorder < matBrth - 1) mat[i][rCalcBorder + 1] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!

		for(j = lCalcBorder; j <= rCalcBorder; ++j){
			// //Make sure we calculate a banded alignment if matrix is too big
			// if(j < lCalcBorder || j > rCalcBorder){
			// 	mat[i][j] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
			// } else{

			//Testing
			// cout << "Do we get here?" << endl;
			// cout << "i: " << i << " j: " << j << endl;
			// cout << "posQ: " << posQ << " posU: " << posU << endl;
			// if(i == j) cout << "Bases to compare unitig: " << uSeq[posU - j + 1] << " query: " << q[posQ - i + 1] << endl;

			//Calculate current cells value and check we don't drop
			if((mat[i][j] = max(mat[i - 1][j - 1] + compUScore(q[posQ - i + 1], uSeq[posU - j + 1]), max(mat[i][j - 1] + GAP_SCORE, mat[i - 1][j] + GAP_SCORE))) > -X){//TODO: This has to be changed as soon as we are not using unit score anymore!
				//Make sure we do not terminate our calculations too soon
				termCalcs = false;

				//Testing
				// cout << "Inside calcLeftGlobAlignment: Characters that have been compared were " << q[posQ - i + 1] << " in the query and " << uSeq[posU - j + 1] << " in the unitig" << endl;

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

			//Testing
			// cout << "mat[i][j]: " << mat[i][j] << " ";
		}

		//Testing
		// cout << endl;

		//Check whether calculations may be terminated
		if(termCalcs) break;
	}

	//Testing
	// cout << "Filled matrix:" << endl;
	// for(i = 0; i < matHgth; ++i){
	// 	for(j = 0; j < matBrth; ++j){
	// 		cout << mat[i][j] << " ";
	// 	}
	// 	cout << endl;
	// }

	//Do backtracing to get the alignment maximum scoring alignment
	traceForth(maxMatPosI, maxPosQ, q, maxMatPosJ, maxPosU, uSeq, mat, maxAlgn);

	//Check if maximum score is at the matrix border
	if(maxMatPosI == edgeMaxPosI && maxMatPosJ == edgeMaxPosJ){
		//Testing
		// if(uSeq == "CTCACCTTGGTTGAAAAATGCTAAAAGTGCCCCCTCAATACATTTATCAAGATATTGGTCT"){
		// 	cout << "We can copy the maximum alignment" << endl;
		// }

		//Copy maximum scoring alignment to save time
		brdAlgn = maxAlgn;
	} else{
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
		//Testing
		// cout << "Do we get into this branch?" << endl;
		
		posU = maxPosU;
		posQ = maxPosQ;
	}

	//Testing
	// if(uSeq == "CTCACCTTGGTTGAAAAATGCTAAAAGTGCCCCCTCAATACATTTATCAAGATATTGGTCT"){
	// 	cout << "maxAlgn: " << maxAlgn.aSeqG << endl;
	// 	exit(0);
	// }
	// cout << "calcLeftGlobAlignment: posQ before leaving: " << posQ << endl;
	// cout << "eMax: " << eMax << " matShrunk? " << (matShrunk ? "Yes" : "No") << endl;

	if(eMax > -X && matShrunk && posQ > 0) return false;
	
	return true;
}

// //This function calculates the continuation of a gapped alignment on a successive unitig or the next peace of the query considering a quorum
// void contRightGappedAlignment(UnitigColorMap<seedlist> &uni, list<uint16_t> &extPth, const string &q, uint32_t &posQ, uint32_t &posU, const uint16_t &X, const uint32_t &maxGaps, struct Algn &algn, int32_t &score, uint32_t &explCount, const uint32_t &quorum){
// 	bool explSuc = false, unitigDone;
// 	uint32_t maxPosQ = posQ, maxPosU = posU;
// 	int32_t maxBorderScore = -X;
// 	struct Algn globAlgn;
// 	UnitigColorMap<seedlist> maxUni = uni;

// 	//
// 	//++posQ;

// 	//Calculate gapped alignment
// 	unitigDone = calcSemiGlobAlignment(uni, q, posU, posQ, X, maxGaps, maxPosQ, maxPosU, algn, globAlgn, score, maxBorderScore, quorum);

// 	//Testing
// 	//cout << "posU: " << posU << "\nposQ: " << posQ << endl;
// 	// cout << "Inside contRightGappedAlignment after calculation of gapped alignment on the current unitig:\nunitigDone: " << unitigDone << "\nmaxBorderScore: " << maxBorderScore << "\nposU: " << posU << "\nposQ: " << posQ << "\nmaxPosQ: " << maxPosQ << "\nmaxPosU: " << maxPosU << "\nscore: " << score << "\nQuery length: " << q.length() << endl;

// 	//Check outcome of alignment calculation
// 	if(!unitigDone){
// 		explSuc = contGappedOnSameUni(uni, extPth, q, posQ, posU, X, maxGaps, globAlgn, maxBorderScore, explCount, quorum);
// 	} else if(maxBorderScore > -X && posQ < q.length() - 1){
// 		//Testing
// 		// cout << "We are on unitig " << uni.mappedSequenceToString() << " and continue on its successors" << endl;

// 		//Continue calculations on the next unitig
// 		explSuc = contGappedOnSuccUni(uni, extPth, q, posQ, posU, X, maxGaps, globAlgn, maxBorderScore, explCount, quorum);
// 	}

// 	//Check whether exploration of next unitig was successful (this is only the case if a new maximum score has been found as well)
// 	if(explSuc && score < maxBorderScore){
// 		//Score has to be updated if we have found a better one on a successive unitig
// 		score = maxBorderScore;
// 		//Overwrite previous best alignment
// 		algn = globAlgn;
// 	} else{
// 		//Save new offsets in hit
// 		posU = maxPosU;
// 		posQ = maxPosQ;
// 		//Reset the unitig on which the maximum has been found
// 		uni = maxUni;
// 	}
// }

//This function calculates the continuation of a gapped alignment on a successive unitig or the next peace of the query considering a quorum and a search color set
void contRightGappedAlignment(UnitigColorMap<seedlist> &uni, list<uint16_t> &extPth, const string &q, uint32_t &posQ, uint32_t &posU, const uint16_t &X, const uint32_t &maxGaps, struct Algn &algn, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	bool explSuc = false, unitigDone;
	uint32_t maxPosQ = posQ, maxPosU = posU;
	int32_t maxBorderScore = -X;
	struct Algn globAlgn;
	UnitigColorMap<seedlist> maxUni = uni;

	//Testing
	// cout << "Start calculations on unitig " << uni.mappedSequenceToString() << endl;

	//
	//++posQ;
	//Calculate gapped alignment
	unitigDone = calcSemiGlobAlignment(uni, q, posU, posQ, X, maxGaps, maxPosQ, maxPosU, algn, globAlgn, score, maxBorderScore, quorum, searchSet, extPth.empty());

	//Testing
	// cout << "Finished calculations on this unitig" << endl;
	// if(report){
	// 	cout << "Inside contRightGappedAlignment after calculations on this unitig: algn:" << "aSeqG: " << algn.aSeqG << endl << "aSeqQ: " << algn.aSeqQ << endl << "globAlgn:" << "aSeqG: " << globAlgn.aSeqG << endl << "aSeqQ: " << globAlgn.aSeqQ << endl;
	// 	// cout << "maxBorderScore: " << maxBorderScore << " unitigDone: " << unitigDone << endl;
	// }

	//Testing
	// cout << "posU: " << posU << " posQ: " << posQ << endl;
	//cout << "Inside contRightGappedAlignment after calculation of gapped alignment on the current unitig:\nunitigDone: " << unitigDone << "\nmaxBorderScore: " << maxBorderScore << "\nposU: " << posU << "\nposQ: " << posQ << "\nmaxPosQ: " << maxPosQ << "\nmaxPosU: " << maxPosU << "\nscore: " << score << "\nQuery length: " << q.length() << endl;

	//Check if query sequence is left
	// if(posQ < q.length() - 1){

	//Check outcome of alignment calculation
	if(!unitigDone){
		explSuc = contGappedOnSameUni(uni, extPth, q, posQ, posU, X, maxGaps, globAlgn, maxBorderScore, explCount, quorum, searchSet);
	} else if(maxBorderScore > -X && posQ < q.length() - 1){
		//Testing
		// cout << "Check out next unitig" << endl;

		//Continue calculations on the next unitig
		explSuc = contGappedOnSuccUni(uni, extPth, q, posQ, posU, X, maxGaps, globAlgn, maxBorderScore, explCount, quorum, searchSet);
	}

	// }

	//Testing
	// cout << "Did calculation continuation stuff" << endl;

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
		uni = maxUni;//TODO I guess this is not needed anymore because we do not store it in the end anymore anyways!
	}
}

//This function calculates the continuation of a gapped alignment on a predecessive unitig or the next peace of the query considering a quorum
void contLeftGappedAlignment(UnitigColorMap<seedlist> &uni, list<uint16_t> &extPth, const string &q, uint32_t &posQ, uint32_t &posU, const uint16_t &X, const uint32_t &maxGaps, struct Algn &maxAlgn, int32_t &score, uint32_t &explCount, const uint32_t &quorum){
	bool unitigDone, explSuc = false;
	uint32_t maxPosQ = posQ, maxPosU = posU;
	int32_t maxBorderScore;
	struct Algn brdAlgn;
	UnitigColorMap<seedlist> maxUni = uni;

	//Calculate gapped alignment
	unitigDone = calcLeftGlobAlignment(uni, q, posU, posQ, X, maxGaps, maxPosQ, maxPosU, maxAlgn, brdAlgn, score, maxBorderScore, quorum);

	//Testing
	//cout << "Exploration counter: " << explCount << " posQ: " << posQ << endl;

	//Check outcome of alignment calculation
	if(!unitigDone){
		explSuc = contLeftOnSameUni(uni, extPth, q, posQ, posU, X, maxGaps, brdAlgn, maxBorderScore, explCount, quorum);
	} else if(maxBorderScore > -X && posQ > 0){
		//Testing
		//cout << "Call of contGappedOnPredUni" << endl;

		//Continue calculations on the next unitig
		explSuc = contGappedOnPredUni(uni, extPth, q, posQ, posU, X, maxGaps, brdAlgn, maxBorderScore, explCount, quorum);
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

//This function calculates the continuation of a gapped alignment on a predecessive unitig or the next peace of the query considering a quorum and a search color set
void contLeftGappedAlignment(UnitigColorMap<seedlist> &uni, list<uint16_t> &extPth, const string &q, uint32_t &posQ, uint32_t &posU, const uint16_t &X, const uint32_t &maxGaps, struct Algn &maxAlgn, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	bool unitigDone, explSuc = false;
	uint32_t maxPosQ = posQ, maxPosU = posU;
	int32_t maxBorderScore;
	struct Algn brdAlgn;
	UnitigColorMap<seedlist> maxUni = uni;

	//Testing
	// cout << "Starting contLeftGappedAlignment" << endl;
	// if(uni.mappedSequenceToString() == "CTCACCTTGGTTGAAAAATGCTAAAAGTGCCCCCTCAATACATTTATCAAGATATTGGTCT"){
	// 	cout << "We have found the interesting unitig" << endl;
	// 	cout << "Extension path is " << (extPth.empty() ? "" : "not ") << "empty" << endl;
	// }

	//Calculate gapped alignment
	unitigDone = calcLeftGlobAlignment(uni, q, posU, posQ, X, maxGaps, maxPosQ, maxPosU, maxAlgn, brdAlgn, score, maxBorderScore, quorum, searchSet, extPth.empty());

	//Testing
	// cout << "Exploration counter: " << explCount << " posQ: " << posQ << endl;
	// cout << "Unitig is " << (unitigDone ? "" : "not ") << "done" << endl;

	// //If we are at the queries beginning we are done
	// if(posQ > 0){

	//Check outcome of alignment calculation
	if(!unitigDone){
		//Testing
		// cout << "Further calculations on this unitig" << endl;

		explSuc = contLeftOnSameUni(uni, extPth, q, posQ, posU, X, maxGaps, brdAlgn, maxBorderScore, explCount, quorum, searchSet);
	} else if(maxBorderScore > -X && posQ > 0){
		//Testing
		// cout << "Call of contGappedOnPredUni" << endl;

		//Continue calculations on the next unitig
		explSuc = contGappedOnPredUni(uni, extPth, q, posQ, posU, X, maxGaps, brdAlgn, maxBorderScore, explCount, quorum, searchSet);
	}

	// }

	//Testing
	// cout << "Setting score alignment and borders" << endl;

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

// //Continues a gapped alignment on the same unitig as before considering a quorum. Returns true if a new maximal score has been found.
// bool contGappedOnSameUni(UnitigColorMap<seedlist> &uni, list<uint16_t> &extPth, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, struct Algn &algn, int32_t &score, uint32_t &explCount, const uint32_t &quorum){
// 	uint32_t tmpQOff = qOff + 1, tmpUOff = uOff + 1;//TODO Wouldn't it be better not to add 1 here? When the actual matrix is calculated positions in q and u are always lowered by one
// 	int32_t tmpScore = score;
// 	struct Algn tmpAlgn;

// 	//++qOff;
// 	//++uOff;

// 	contRightGappedAlignment(uni, extPth, q, tmpQOff, tmpUOff, X, maxGaps, tmpAlgn, tmpScore, explCount, quorum);

// 	//Check whether we have found a new maximum
// 	if(score < tmpScore){
// 		//Update maximum
// 		score = tmpScore;
// 		qOff = tmpQOff;
// 		uOff = tmpUOff;
// 		algn.aSeqG += tmpAlgn.aSeqG;
// 		algn.aSeqQ += tmpAlgn.aSeqQ;
// 		return true;
// 	}

// 	return false;
// }

//Continues a gapped alignment on the same unitig as before considering a quorum and a search color set. Returns true if a new maximal score has been found.
bool contGappedOnSameUni(UnitigColorMap<seedlist> &uni, list<uint16_t> &extPth, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, struct Algn &algn, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	uint32_t tmpQOff = qOff + 1, tmpUOff = uOff + 1;
	int32_t tmpScore = score;
	struct Algn tmpAlgn;

	//Testing
	// cout << "Exploration on same unitig" << endl;
	// cout << "tmpQOff: " << tmpQOff << endl;
	// cout << "qOff: " << qOff << endl;

	//++qOff;
	//++uOff;

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

//This function continues a gapped alignment calculation to the left on the same unitig considering a quorum. Returns true if a new maximal score has been found.
bool contLeftOnSameUni(UnitigColorMap<seedlist> &uni, list<uint16_t> &extPth, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, struct Algn &algn, int32_t &score, uint32_t &explCount, const uint32_t &quorum){
	uint32_t tmpQOff = qOff - 1, tmpUOff = uOff - 1;
	int32_t tmpScore = score;
	struct Algn tmpAlgn;

	contLeftGappedAlignment(uni, extPth, q, tmpQOff, tmpUOff, X, maxGaps, tmpAlgn, tmpScore, explCount, quorum);

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

// //Initiates a gapped alignment calculation on all successive unitigs considering a quorum and returns the alignment's end position in the unitig
// bool contGappedOnSuccUni(UnitigColorMap<seedlist> &uni, list<uint16_t> &extPth, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, struct Algn &algn, int32_t &score, uint32_t &explCount, const uint32_t &quorum){
// 	bool suc = false;
// 	uint16_t sucCount, nextSuc;
// 	uint32_t tmpPosU, tmpPosQ, maxPosQ, maxPosU;
// 	int32_t maxScore, tmpScore;
// 	struct Algn maxAlgn;
// 	UnitigColorMap<seedlist> sucUni;
// 	ForwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> sucIter = uni.getSuccessors();

// 	//Testing
// 	//cout << "We are coming from unitig " << uni.mappedSequenceToString() << endl;

// 	//Check whether we have reached the maximum recursion depth of an extension
// 	if(++explCount > MAXRECURSIONDEPTH){
// 		//Report this incident//TODO: If this happens the condition that an extended seed cannot be reached anymore is violated. Implement a procedure after the extension that merges such seeds with other seeds that would have reached them!
// 		//cerr << "Maximum recursion depth reached during extension. Position in q: " << qOff << endl;
// 		//Terminate this extension
// 		return suc;
// 	}

// 	//Check whether there are successors to continue
// 	if(sucIter.hasSuccessors()){
// 		maxScore = score;
// 		//Set counter
// 		sucCount = 1;

// 		//Check if there still exists an extension path to follow
// 		if(!extPth.empty()){
// 			//Get next successor
// 			nextSuc = extPth.front();
// 			//Delete that successor from path
// 			extPth.pop_front();
// 			//We do not need to count how many successors we have explored yet
// 			explCount = 0;

// 			//Testing
// 			// cout << "We continue gapped on successor " << nextSuc << endl;
// 		} else{
// 			//We have no successor given to go on with
// 			nextSuc = 0;
// 		}

// 		//Explore each neighbor
// 		for(neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> nI = sucIter.begin(); nI != sucIter.end(); ++nI, ++sucCount){
// 			struct Algn tmpAlgn;

// 			//Testing
// 			// cout << "Starting next iteration" << endl;

// 			//Check if we still have to wait for the correct successor to go on with
// 			if(sucCount < nextSuc){
// 				//Testing
// 				// cout << "We skip successor " << nI->mappedSequenceToString() << endl;
// 				//++nI;
// 				// cout << "Next successor would be " << nI->mappedSequenceToString() << endl;
// 				// cout << "(before continue) sucCount:" << sucCount << endl;

// 				continue;	
// 			}

// 			//Testing
// 			// cout << "We have found the correct successor: " << nI->mappedSequenceToString() << endl;
// 			// cout << "sucCount:" << sucCount << endl;

// 			//Get next unitig
// 			sucUni = *nI;
// 			tmpPosQ = qOff + 1;
// 			tmpPosU = 0;
// 			tmpScore = score;
// 			//Calculate gapped alignment on the next unitig
// 			contRightGappedAlignment(sucUni, extPth, q, tmpPosQ, tmpPosU, X, maxGaps, tmpAlgn, tmpScore, explCount, quorum);

// 			//Testing
// 			//cout << "We return from further calculating right extension" << endl << "sucUni: " << sucUni.mappedSequenceToString() << " maxScore: " << maxScore << " tmpScore: " << tmpScore << endl;

// 			//Check for new maximum
// 			if(maxScore < tmpScore){
// 				//Testing
// 				// cout << "Found best alignment on " << sucUni.mappedSequenceToString() << " score is " << tmpScore << endl;

// 				//Update maximum
// 				maxScore = tmpScore;
// 				maxPosU = tmpPosU;
// 				maxPosQ = tmpPosQ;
// 				uni = sucUni;
// 				maxAlgn = tmpAlgn;
// 				suc = true;
// 			}

// 			//If we reach this point and had a specific successor to go on with we are done
// 			if(nextSuc > 0) break;
// 		}

// 		//Testing
// 		// cout << "Do we get here?" << endl;

// 		//Only update score and qOff if exploration was successful
// 		if(suc){
// 			//Give back best result found
// 			score = maxScore;
// 			qOff = maxPosQ;
// 			uOff = maxPosU;
// 			algn.aSeqG += maxAlgn.aSeqG;
// 			algn.aSeqQ += maxAlgn.aSeqQ;
// 		}
// 	}

// 	return suc;
// }

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
		//Testing
		// cout << "4 Option 1" << endl;

		//Report this incident//TODO: If this happens the condition that an extended seed cannot be reached anymore is violated. Implement a procedure after the extension that merges such seeds with other seed that would have reached them!
		//cerr << "Maximum recursion depth reached during extension. Position in q: " << qOff << endl;
		//Terminate this extension
		return suc;
	}

	//Testing
	// cout << "4 Option 2" << endl;
	// if(report){
	// 	cout << "Inside contGappedOnSuccUni: Maximum recursion depth has not yet been reached" << endl;
	// }

	//Check whether there are successors to continue
	if(sucIter.hasSuccessors()){
		//Testing
		// if(report){
		// 	cout << "Inside contGappedOnSuccUni: There are successors to visit" << endl;
		// 	cout << "Inside contGappedOnSuccUni: Extension path is " << (extPth.empty() ? "" : "not ") << "empty" << endl;
		// }

		maxScore = score;
		//Set counter
		sucCount = 1;

		//Check if there still exists an extension path to follow
		if(!extPth.empty()){
			//Testing
			// cout << "5 Option 1" << endl;

			//Get next successor
			nextSuc = extPth.front();
			//Delete that successor from path
			extPth.pop_front();
			//We do not need to count how many successors we have explored yet
			explCount = 0;
		} else{
			//Testing
			// cout << "5 Option 2" << endl;

			//We have no successor given to go on with
			nextSuc = 0;
		}

		//Explore each neighbor
		for(neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> nI = sucIter.begin(); nI != sucIter.end(); ++nI, ++sucCount){
			struct Algn tmpAlgn;

			//Testing
			// if(report){
			// 	cout << "Inside contGappedOnSuccUni: nextSuc: " << nextSuc << " sucCount: " << sucCount << endl;
			// }

			//Check if we still have to wait for the correct successor to go on with
			if(sucCount < nextSuc){
				//Testing
				// cout << "6" << endl;

				continue;
			}

			//Get next unitig
			sucUni = *nI;
			tmpPosQ = qOff + 1;
			tmpPosU = 0;
			tmpScore = score;

			//Testing
			// cout << "Do we get here?" << endl;
			// if(report){
			// 	cout << "Last unitig was " << uni.mappedSequenceToString() << endl;
			// 	cout << "Let's try unitig " << nI->mappedSequenceToString() << endl;
			// }

			//Calculate gapped alignment on the next unitig
			contRightGappedAlignment(sucUni, extPth, q, tmpPosQ, tmpPosU, X, maxGaps, tmpAlgn, tmpScore, explCount, quorum, searchSet);

			//Testing
			// cout << "Do we get out there as well?" << endl;
			// if(report){
			// 	cout << "We came back." << endl; //" Alignment is:" << "aSeqG: " << tmpAlgn.aSeqG << endl << "aSeqQ: " << tmpAlgn.aSeqQ << endl;
			// }

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
			if(nextSuc > 0){
				//Testing
				// cout << "7 Option 1" << endl;

				break;
			}

			//Testing
			// cout << "7 Option 2" << endl;
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

//This function initiates a gapped alignment calculation on all predecessive unitigs considering a quorum. Returns true if calculations have been successful
bool contGappedOnPredUni(UnitigColorMap<seedlist> &uni, list<uint16_t> &extPth, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, struct Algn &algn, int32_t &score, uint32_t &explCount, const uint32_t &quorum){
	bool suc = false;
	uint16_t predCount, nextPred;
	uint32_t tmpPosU, tmpPosQ, maxPosQ, maxPosU;
	int32_t maxScore, tmpScore;
	struct Algn maxAlgn;
	UnitigColorMap<seedlist> predUni;
	BackwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> predIter = uni.getPredecessors();

	//Testing
	//cout << "Last unitig was " << uni.toString() << endl;
	// cout << "We start on the predecessor directly" << endl;

	//Check whether we have reached the maximum recursion depth of an extension
	if(++explCount > MAXRECURSIONDEPTH){
		//Report this incident//TODO: If this happens the condition that an extended seed cannot be reached anymore is violated. Implement a procedure after the extension that merges such seeds with other seed that would have reached them!
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
			//Testing
			// cout << "Extension path is not empty" << endl << "Next comes " << extPth.front() << endl;

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

		//Testing
		// cout << "We are about to switch on the predecessive unitig" << endl;
		// cout << "nextPred:" << nextPred << endl;

		//Explore each neighbor
		for(neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> nI = predIter.begin(); nI != predIter.end(); ++nI, ++predCount){
			//Check if we still have to wait for the correct predecessor to go on with
			if(predCount < nextPred){
				//Testing
				// cout << "We skip unitig " << nI->mappedSequenceToString() << endl;

				continue;	
			}

			struct Algn tmpAlgn;
			//Get next unitig
			predUni = *nI;
			tmpPosQ = qOff - 1;
			tmpPosU = predUni.size - predUni.getGraph()->getK();
			tmpScore = score;

			//Testing
			// cout << "Next unitig is " << predUni.mappedSequenceToString() << endl;
			// cout << "tmpPosU:" << tmpPosU << "tmpPosQ:" << tmpPosQ << "tmpScore:" << tmpScore << endl;

			//Calculate gapped alignment on the predecessive unitig
			contLeftGappedAlignment(predUni, extPth, q, tmpPosQ, tmpPosU, X, maxGaps, tmpAlgn, tmpScore, explCount, quorum);

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

//This function initiates a gapped alignment calculation on all predecessive unitigs considering a quorum and a search color set. Returns true if calculations have been successful
bool contGappedOnPredUni(UnitigColorMap<seedlist> &uni, list<uint16_t> &extPth, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, struct Algn &algn, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	bool suc = false;
	uint16_t predCount, nextPred;
	uint32_t tmpPosU, tmpPosQ, maxPosQ, maxPosU;
	int32_t maxScore, tmpScore;
	struct Algn maxAlgn;
	UnitigColorMap<seedlist> predUni;
	BackwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> predIter = uni.getPredecessors();

	//Testing
	// cout << "Last unitig was " << uni.mappedSequenceToString() << endl;

	//Check whether we have reached the maximum recursion depth of an extension
	if(++explCount > MAXRECURSIONDEPTH){
		//Testing
		// cout << "1 Option 1" << endl;

		//Report this incident//TODO: If this happens the condition that an extended seed cannot be reached anymore is violated. Implement a procedure after the extension that merges such seeds with other seed that would have reached them!
		//cerr << "Maximum recursion depth reached during extension. Position in q: " << qOff << endl;
		//Terminate this extension
		return suc;
	}

	//Testing
	// cout << "1 Option 2" << endl;

	//Check whether there are predecessors to continue
	if(predIter.hasPredecessors()){
		maxScore = score;
		//Set counter
		predCount = 1;

		//Check if there still exists an extension path to follow
		if(!extPth.empty()){
			//Testing
			// cout << "2 Option 1" << endl;

			//Get next predecessor
			nextPred = extPth.front();
			//Delete that predecessor from path
			extPth.pop_front();
			//We do not need to count how many successors we have explored yet
			explCount = 0;
		} else{
			//Testing
			// cout << "2 Option 2" << endl;

			//We have no predecessor given to go on with
			nextPred = 0;
		}

		//Explore each neighbor
		for(neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> nI = predIter.begin(); nI != predIter.end(); ++nI, ++predCount){
			//Check if we still have to wait for the correct predecessor to go on with
			if(predCount < nextPred){
				//Testing
				// cout << "3 Option 1" << endl;

				continue;
			}

			struct Algn tmpAlgn;
			//Get next unitig
			predUni = *nI;
			tmpPosQ = qOff - 1;
			tmpPosU = predUni.size - predUni.getGraph()->getK();
			tmpScore = score;

			//Testing
			// if(report){
			// cout << "Next unitig is " << predUni.mappedSequenceToString() << endl;
			// }

			//Calculate gapped alignment on the predecessive unitig
			contLeftGappedAlignment(predUni, extPth, q, tmpPosQ, tmpPosU, X, maxGaps, tmpAlgn, tmpScore, explCount, quorum, searchSet);

			//Testing
			// cout << "tmpScore: " << tmpScore << endl;

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
			if(nextPred > 0){
				//Testing
				// cout << "3 Option 2" << endl;

				break;
			}
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

// //This function calculates a gapped alignment to the right side of the starting position considering a quorum. ATTENTION: Hit's length attribute will be deprecated after function call!
// void startRightGappedAlignment(hit *h, const string &q, const uint16_t &X, const uint32_t maxGaps, const uint32_t &quorum){
// 	bool explSuc = false, unitigDone = true;
// 	uint32_t posU = h->offU, posQ = h->offQ, maxPosU = h->offU, maxPosQ = h->offQ, explCount = 0;//, endBuf = 0
// 	int32_t maxScore = 0, maxBorderScore;
// 	struct Algn globAlgn;
// 	list<uint16_t> extPth = decmprExtPth(h->rExt);
// 	UnitigColorMap<seedlist> currUni = h->origUni;

// 	//Testing
// 	// cout << "We want to do gapped extension now" << endl << "Extension path is" << endl;
// 	// for(list<uint16_t>::const_iterator i = extPth.begin(); i != extPth.end(); ++i) cout << *i << endl;
// 	// cout << "(Beginning of right gapped extension) We are on unitig " << currUni.mappedSequenceToString() << endl;

// 	// //Check whether we can start the calculation right on this unitig or whether we are already at the end of unitig or query sequence
// 	// if(posQ < q.length() && posU < currUni.referenceUnitigToString().length() - endBuf){<-Actually this has to be the case always
// 		// //Set alignment start positions (+1, because we want to start right behind the hit)
// 		// ++posU;
// 		// ++posQ;

// 	//Calculate gapped alignment
// 	unitigDone = calcSemiGlobAlignment(currUni, q, posU, posQ, X, maxGaps, maxPosQ, maxPosU, h->gAlgn, globAlgn, maxScore, maxBorderScore, quorum);
	
// 	// } else{
// 	// 	//Set variable for exploration of successive unitigs
// 	// 	maxBorderScore = 0;
// 	// }

// 	//Check whether we have reached the end of the current unitig
// 	if(!unitigDone){
// 		//Testing
// 		//cout << "Do we go in here?" << endl;

// 		//Continue calculatiions on the same unitig
// 		explSuc = contGappedOnSameUni(currUni, extPth, q, posQ, posU, X, maxGaps, globAlgn, maxBorderScore, explCount, quorum);
// 	} else if(maxBorderScore > -X && posQ < q.length() - 1){
// 		//Testing
// 		//cout << "Gehe zum nächsten unitig" << endl;

// 		//Continue calculations on the next unitig
// 		explSuc = contGappedOnSuccUni(currUni, extPth, q, posQ, posU, X, maxGaps, globAlgn, maxBorderScore, explCount, quorum);

// 		//Testing
// 		// cout << "Weiter zum Rest" << endl;
// 		// cout << "currUni: " << currUni.mappedSequenceToString() << endl;
// 		// cout << "explSuc: " << explSuc << " maxScore: " << maxScore << " maxBorderScore: " << maxBorderScore << endl;
// 	}

// 	//Check whether exploration of next unitig was successful
// 	if(explSuc && maxScore < maxBorderScore){
// 		// //Save new offsets in hit
// 		// h->rSeedUoff = posU;
// 		// h->rSeedQoff = posQ;

// 		//Save the achieved score
// 		h->score = maxBorderScore;

// 		// //Reset right border unitig//TODO This was in the else case so far which does not make any sense to me (we should assimilate a new unitig if the extension on a successive unitig was successful right?). However nothing happens in this case if we cannot become better using gaps anyways. Maybe this case never appeared when this function was tested the last time. The test of this function should be revised to check if this is actually a fix or not!  <-Not necessary rUnitig will be removed soon anyways
// 		// h->rUnitig = currUni;

// 		//Save alignment
// 		h->gAlgn = globAlgn;
// 	} else{
// 		//Save the achieved score
// 		h->score = maxScore;

// 		// //Reset right border unitig
// 		// h->rUnitig = currUni;
// 	}

// 	//Testing
// 	// cout << "(End of right gapped extension) We are on unitig " << h->origUni.mappedSequenceToString() << endl;
// 	// cout << "h->rUnitig: " << h->rUnitig.mappedSequenceToString() << endl;
// }

//This function calculates a gapped alignment to the right side of the starting position considering a quorum and a search color set. ATTENTION: Hit's length attribute will be deprecated after function call!
void startRightGappedAlignment(hit *h, const string &q, const uint16_t &X, const uint32_t maxGaps, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	bool explSuc = false, unitigDone = true;
	uint32_t posU = h->offU, posQ = h->offQ, maxPosU = h->offU, maxPosQ = h->offQ, explCount = 0;
	int32_t maxScore = 0, maxBorderScore;
	struct Algn globAlgn;
	list<uint16_t> extPth = decmprExtPth(h->rExt);
	UnitigColorMap<seedlist> currUni = h->origUni;

	//Testing
	// cout << "2 Option ";

	// //Check whether current unitig has a successor
	// if(currUni.getSuccessors().hasSuccessors()){
	// 	//Testing
	// 	// cout << "not ";

	// 	endBuf = currUni.getGraph()->getK();
	// }

	//Testing
	// cout << "2" << endl;
	// if(h->origUni.mappedSequenceToString() == "TCTTTACGGCGAAGTTCAGCGCCCTCATAGCC" && h->offU == 20 && h->offQ == 47){
	// 	// cout << "Seed of the unitig:" << endl << "offsetU: " << currSeed->offsetU << " offsetQ: " << currSeed->offsetQ << " score: " << currSeed->score << " length: " << currSeed->len << endl;
	// 	cout << "Inside startRightGappedAlignment before calculations: posU: " << posU << "posQ: " << posQ << endl << "aSeqG: " << (*it)->gAlgn.aSeqG << endl << "aSeqQ: " << (*it)->gAlgn.aSeqQ << endl;
	// }
	// if(h->origUni.mappedSequenceToString() == "AATTTGTCGCATATTATTTTTTGTATTTTTGGC"){
	// 	// cout << "Seed of the unitig:" << endl << "offsetU: " << currSeed->offsetU << " offsetQ: " << currSeed->offsetQ << " score: " << currSeed->score << " length: " << currSeed->len << endl;
	// 	cout << "Inside startRightGappedAlignment before calculations on this unitig: posU: " << posU << " posQ: " << posQ << endl;
	// 	report = true;
	// }

	// //Check whether we can start the calculation right on this unitig or whether we are already at the end of unitig or query sequence
	// if(posQ < q.length() - 1 && posU < currUni.referenceUnitigToString().length() - endBuf){
	// 	//Set alignment start positions (+1, because we want to start right behind the hit)
	// 	++posU;
	// 	++posQ;

		//Testing
		// cout << "1 Option 1" << endl;

	//Calculate gapped alignment
	unitigDone = calcSemiGlobAlignment(currUni, q, posU, posQ, X, maxGaps, maxPosQ, maxPosU, h->gAlgn, globAlgn, maxScore, maxBorderScore, quorum, searchSet, extPth.empty());

	// } else{
	// 	//Testing
	// 	// cout << "1 Option 2" << endl;

	// 	//Set variable for exploration of successive unitigs
	// 	maxBorderScore = 0;
	// }

	//Testing
	// cout << "Inside startRightGappedAlignment after calculation of gapped alignment on the first unitig:\nunitigDone: " << unitigDone << "\nmaxBorderScore: " << maxBorderScore << "\nposU: " << posU << "\nposQ: " << posQ << "\nmaxPosQ: " << maxPosQ << "\nmaxPosU: " << maxPosU << "\nmaxScore: " << maxScore << endl;
	// if(h->origUni.mappedSequenceToString() == "GCAAGAGCCGCTGTTTCTTGAACAATATCTCG"){
	// 	// cout << "Seed of the unitig:" << endl << "offsetU: " << currSeed->offsetU << " offsetQ: " << currSeed->offsetQ << " score: " << currSeed->score << " length: " << currSeed->len << endl;
	// 	cout << "Inside startRightGappedAlignment after calculations on this unitig: gAlgn:" << endl << "aSeqG: " << h->gAlgn.aSeqG << endl << "aSeqQ: " << h->gAlgn.aSeqQ << endl << "globAlgn: aSeqG: " << globAlgn.aSeqG << endl << "globAlgn: aSeqQ: " << globAlgn.aSeqQ << endl;
	// 	cout << "Inside startRightGappedAlignment after calculations on this unitig: Extension path:" << endl;
	// 	for(list<uint16_t>::const_iterator i = extPth.begin(); i != extPth.end(); ++i) cout << *i << " ";
	// 	cout << endl;
	// }

	// //Check if query sequence is left
	// if(){

	//Check whether we have reached the end of the current unitig
	if(!unitigDone){
		explSuc = contGappedOnSameUni(currUni, extPth, q, posQ, posU, X, maxGaps, globAlgn, maxBorderScore, explCount, quorum, searchSet);
	} else if(maxBorderScore > -X && posQ < q.length() - 1){
		//Testing
		//cout << "Gehe zum nächsten unitig" << endl;

		//Continue calculations on the next unitig
		explSuc = contGappedOnSuccUni(currUni, extPth, q, posQ, posU, X, maxGaps, globAlgn, maxBorderScore, explCount, quorum, searchSet);

		//Testing
		//cout << "Weiter zum Rest" << endl;
	}

	// }

	//Testing
	// if(h->origUni.mappedSequenceToString() == "AATTTGTCGCATATTATTTTTTGTATTTTTGGC"){
	// 	cout << "Inside startRightGappedAlignment after all calculations: gAlgn:" << endl << "aSeqG: " << h->gAlgn.aSeqG << endl << "aSeqQ: " << h->gAlgn.aSeqQ << endl << "globAlgn: aSeqG: " << globAlgn.aSeqG << endl << "globAlgn: aSeqQ: " << globAlgn.aSeqQ << endl;
	// 	// exit(0);
	// }
	cout << "maxBorderScore: " << maxBorderScore << endl;
	cout << "maxScore: " << maxScore << endl;

	//Check whether exploration of next unitig was successful
	if(explSuc && maxScore < maxBorderScore){
		// //Save new offsets in hit
		// h->offU = posU;
		// h->offQ = posQ;

		//Save the achieved score
		h->score = maxBorderScore;
		//Save alignment
		h->gAlgn = globAlgn;
	} else{
		// //Save new offsets in hit
		// h->offU = maxPosU;
		// h->offQ = maxPosQ;

		//Save the achieved score
		h->score = maxScore;

		// //Reset right border unitig
		// h->origUni = currUni;
	}
}

//This function calculates a gapped alignment to the left side of the starting position considering a quorum. ATTENTION: Hit's length attribute will be deprecated after function call!
void startLeftGappedAlignment(hit *h, const string &q, const uint16_t &X, const uint32_t maxGaps, const uint32_t &quorum){
	bool explSuc = false, unitigDone = true;
	uint32_t posU = h->offU, posQ = h->offQ, maxPosU = h->offU, maxPosQ = h->offQ, explCount = 0;
	int32_t maxScore = 0, maxBorderScore;
	list<uint16_t> extPth = decmprExtPth(h->lExt);
	struct Algn globAlgn;
	UnitigColorMap<seedlist> currUni = h->origUni;

	//Testing
	// cout << "OrigUni:" << currUni.mappedSequenceToString() << endl;
	// cout << "Start of left gapped extension" << endl << "Initial posU: " << posU << " initial posQ: " << posQ << endl;
	// cout << "Hit score:" << h->score << endl;

	//Check whether we can start the calculation right on this unitig or whether we are already at the begin of unitig or query sequence
	if(posQ > 0 && posU > 0){
		//Testing
		// cout << "Do we get here?" << endl;

		//Set alignment start positions (-1, because we want to start right in front of the hit)
		--posU;
		--posQ;
		//Calculate gapped alignment
		unitigDone = calcLeftGlobAlignment(currUni, q, posU, posQ, X, maxGaps, maxPosQ, maxPosU, h->gAlgn, globAlgn, maxScore, maxBorderScore, quorum);
	} else{
		//Set variable for exploration of predecessive unitigs
		maxBorderScore = 0;
	}

	//Testing
	// cout << "Left extension on the same unitig has finished" << endl << "posU:" << posU << "posQ:" << posQ << "maxScore:" << maxScore << endl;

	//Check whether we have reached the beginning of the current unitig
	if(!unitigDone){
		explSuc = contLeftOnSameUni(currUni, extPth, q, posQ, posU, X, maxGaps, globAlgn, maxBorderScore, explCount, quorum);
	} else if(maxBorderScore > -X && posQ > 0){
		//Continue calculations on the next unitig
		explSuc = contGappedOnPredUni(currUni, extPth, q, posQ, posU, X, maxGaps, globAlgn, maxBorderScore, explCount, quorum);
	}

	//Check whether exploration of predecessive unitig was successful
	if(explSuc && maxScore < maxBorderScore){
		//Save new offsets in hit
		h->offU = posU;
		h->offQ = posQ;

		//Save the achieved score
		h->score += maxBorderScore;
		//Save the better alignment
		h->gAlgn.aSeqG = globAlgn.aSeqG + h->gAlgn.aSeqG;
		h->gAlgn.aSeqQ = globAlgn.aSeqQ + h->gAlgn.aSeqQ;
	} else{
		//Save new offsets in hit
		h->offU = maxPosU;
		h->offQ = maxPosQ;

		//Save the achieved score
		h->score += maxScore;

		// //Reset left border unitig
		// h->lUnitig = currUni;
	}
}

//This function calculates a gapped alignment to the left side of the starting position considering a quorum and a search color set. ATTENTION: Hit's length attribute will be deprecated after function call!
void startLeftGappedAlignment(hit *h, const string &q, const uint16_t &X, const uint32_t maxGaps, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	bool explSuc = false, unitigDone = true;
	uint32_t posU = h->offU, posQ = h->offQ, maxPosU = h->offU, maxPosQ = h->offQ, explCount = 0;
	int32_t maxScore = 0, maxBorderScore;
	list<uint16_t> extPth = decmprExtPth(h->lExt);//TODO After the has been read out we should free the path!
	struct Algn globAlgn = h->gAlgn;
	UnitigColorMap<seedlist> currUni = h->origUni;

	//Testing
	// cout << "Inside startLeftGappedAlignment: Start calculations on this unitig" << endl;

	//Check whether we can start the calculation right on this unitig or whether we are already at the begin of unitig or query sequence
	if(posQ > 0 && posU > 0){
		//Testing
		//cout << "Do we get here?" << endl;

		//Set alignment start positions (-1, because we want to start right in front of the hit)
		--posU;
		--posQ;
		//Calculate gapped alignment
		unitigDone = calcLeftGlobAlignment(currUni, q, posU, posQ, X, maxGaps, maxPosQ, maxPosU, h->gAlgn, globAlgn, maxScore, maxBorderScore, quorum, searchSet, extPth.empty());
	} else{
		//Set variable for exploration of predecessive unitigs
		maxBorderScore = 0;
	}

	//Testing                                  TCTTTACGGCGAAGTTCAGCGCCCTCATAGCC
	// if(h->origUni.mappedSequenceToString() == "CTCACCTTGGTTGAAAAATGCTAAAAGTGCCCCCTCAATACATTTATCAAGATATTGGTCT"){
	// 	// cout << "Inside startLeftGappedAlignment after calculations on same unitig: posU: " << posU << " posQ: " << posQ << " maxPosU: " << maxPosU << " maxPosQ: " << maxPosQ << endl << "aSeqG: " << h->gAlgn.aSeqG << endl << "aSeqQ: " << h->gAlgn.aSeqQ << endl;
	// 	cout << "We have found the interesting unitig" << endl;
	// // 	report = true;
	// 	exit(0);
	// }
	// cout << "Inside startLeftGappedAlignment: Calculations on this unitig done" << endl;

	// //Check if beginning of query is reached
	// if(posQ > 0){

	//Check whether we have reached the beginning of the current unitig
	if(!unitigDone){
		//Testing
		// cout << "Inside startLeftGappedAlignment: Continuing calculations on this unitig" << endl;

		explSuc = contLeftOnSameUni(currUni, extPth, q, posQ, posU, X, maxGaps, globAlgn, maxBorderScore, explCount, quorum, searchSet);
	} else if(maxBorderScore > -X && posQ > 0){
		//Testing
		// cout << "Inside startLeftGappedAlignment: Go on with the next unitig" << endl;

		//Continue calculations on the next unitig
		explSuc = contGappedOnPredUni(currUni, extPth, q, posQ, posU, X, maxGaps, globAlgn, maxBorderScore, explCount, quorum, searchSet);
	}

	// }

	//Testing
	// cout << "Inside startLeftGappedAlignment: All calculations done" << endl;

	//Check whether exploration of predecessive unitig was successful
	if(explSuc && maxScore < maxBorderScore){
		//Save new offsets in hit
		h->offU = posU;
		h->offQ = posQ;
		//Save the achieved score
		h->score += maxBorderScore;

		//Testing                                  TCTTTACGGCGAAGTTCAGCGCCCTCATAGCC
		// if(h->origUni.mappedSequenceToString() == "TCTTTACGGCGAAGTTCAGCGCCCTCATAGCC" && h->offU == posU && h->offQ == posQ){
		// 	cout << "Inside startLeftGappedAlignment before updating the alignment: aSeqG: " << h->gAlgn.aSeqG << endl << "aSeqQ: " << h->gAlgn.aSeqQ << endl;
		// 	// report = true;
		// }	

		//Save the better alignment
		h->gAlgn.aSeqG = globAlgn.aSeqG;// + h->gAlgn.aSeqG;
		h->gAlgn.aSeqQ = globAlgn.aSeqQ;// + h->gAlgn.aSeqQ;

		//Testing                                  TCTTTACGGCGAAGTTCAGCGCCCTCATAGCC
		// if(h->origUni.mappedSequenceToString() == "TCTTTACGGCGAAGTTCAGCGCCCTCATAGCC" && h->offU == posU && h->offQ == posQ){
		// 	cout << "Inside startLeftGappedAlignment after updating the alignment: aSeqG: " << h->gAlgn.aSeqG << endl << "aSeqQ: " << h->gAlgn.aSeqQ << endl;
		// 	//exit(0);
		// }
	} else{
		//Save new offsets in hit
		h->offU = maxPosU;
		h->offQ = maxPosQ;
		//Save the achieved score
		h->score += maxScore;
		
		// //Reset left border unitig
		// h->origUni = currUni;
	}
}

//This function traces back the optimal path from a given start position within an edit matrix to get the corresponding alignment
void traceBack(uint32_t &matPosI, uint32_t posQ, const string &query, uint32_t &matPosJ, uint32_t posU, const string &uniSeq, int32_t **mat, struct Algn &algn){
	// int32_t mScore;

	//We walk through the matrix until we have reached the upper, left corner
	while(matPosI != 0 || matPosJ != 0){
		//Testing
		// cout << "We haven't reached the matrix corner yet" << endl;

		//Check if we have reached the left or upper matrix border already
		if(matPosI == 0){
			//Testing
			// cout << "1 Option 1" << endl;

			//We can only walk on with a gap in the query sequence
			algn.aSeqG = uniSeq[posU] + algn.aSeqG;
			algn.aSeqQ = GAP_SYMB + algn.aSeqQ;
			//Move one cell up
			--matPosJ;
			--posU;
		} else if(matPosJ == 0){
			//Testing
			// cout << "1 Option 2" << endl;

			//We can only walk on with a gap in the unitig sequence
			algn.aSeqG = GAP_SYMB + algn.aSeqG;
			algn.aSeqQ = query[posQ] + algn.aSeqQ;
			//Move one cell to the left
			--matPosI;
			--posQ;
		} else{
			//Testing
			// cout << "compUScore(uniSeq[posU], query[posQ]):" << compUScore(uniSeq[posU], query[posQ]) << endl;
			// cout << "matPosI:" << matPosI << " matPosJ:" << matPosJ << endl;
			// cout << "mat[matPosI - 1][matPosJ - 1]: " << mat[matPosI - 1][matPosJ - 1] << endl;

			// //Find the maximum score from all cells we can reach
			// mScore = max(mat[matPosI - 1][matPosJ], max(mat[matPosI - 1][matPosJ - 1], mat[matPosI][matPosJ - 1]));

			//Testing
			// cout << "Maximum score is " << mScore << endl;

			//Find the cell with the maximum score
			if(mat[matPosI - 1][matPosJ - 1] + compUScore(uniSeq[posU], query[posQ]) == mat[matPosI][matPosJ]){
				//Testing
				// cout << "Let's do a diagonal step" << endl;
				// if(mat[matPosI - 1][matPosJ - 1] < mat[matPosI][matPosJ]){
				// 	// cout << "2 Option 1" << endl;
				// } else{
				// 	cout << "2 Option 2" << endl;
				// }

				//We can add a base from both sequences
				algn.aSeqG = uniSeq[posU] + algn.aSeqG;
				algn.aSeqQ = query[posQ] + algn.aSeqQ;
				//Move to upper left cell
				--matPosI;
				--matPosJ;
				--posQ;
				--posU;
			} else if(mat[matPosI - 1][matPosJ] + GAP_SCORE == mat[matPosI][matPosJ]){
				//Testing
				// cout << "We shouldn't get here" << endl;
				// cout << "2 Option 4" << endl;

				//We can only walk on with a gap in the unitig sequence
				algn.aSeqG = GAP_SYMB + algn.aSeqG;
				algn.aSeqQ = query[posQ] + algn.aSeqQ;
				//Move one cell to the left
				--matPosI;
				--posQ;
			} else{
				//Testing
				// cout << "We shouldn't get here" << endl;
				// cout << "2 Option 3" << endl;

				//We can only walk on with a gap in the query sequence
				algn.aSeqG = uniSeq[posU] + algn.aSeqG;
				algn.aSeqQ = GAP_SYMB + algn.aSeqQ;
				//Move one cell up
				--matPosJ;
				--posU;
			}
		}
		//Testing
		// cout << "We got on cell further" << endl;
		// cout << "matPosI:" << matPosI << " matPosJ:" << matPosJ << endl;
	}
}

//This function traces back the optimal path of a left gapped extension from a given start position within an edit matrix to get the corresponding alignment
void traceForth(uint32_t &matPosI, uint32_t posQ, const string &query, uint32_t &matPosJ, uint32_t posU, const string &uniSeq, int32_t **mat, struct Algn &algn){
	// int32_t mScore;
	string gSeq, qSeq;

	//We walk through the matrix until we have reached the upper, left corner
	while(matPosI != 0 || matPosJ != 0){
		//Testing
		// cout << "We haven't reached the matrix corner yet" << endl;

		//Check if we have reached the left or upper matrix border already
		if(matPosI == 0){
			//Testing
			// cout << "1 Option 1" << endl;

			//We can only walk on with a gap in the query sequence
			gSeq = gSeq + uniSeq[posU];
			qSeq = qSeq + GAP_SYMB;
			//Move one cell up
			--matPosJ;
			++posU;
		} else if(matPosJ == 0){
			//Testing
			// cout << "1 Option 2" << endl;

			//We can only walk on with a gap in the unitig sequence
			gSeq = gSeq + GAP_SYMB;
			qSeq = qSeq + query[posQ];
			//Move one cell to the left
			--matPosI;
			++posQ;
		} else{
			//Testing
			// cout << "matPosI:" << matPosI << " matPosJ:" << matPosJ << endl;
			// cout << "mat[matPosI - 1][matPosJ - 1]: " << mat[matPosI - 1][matPosJ - 1] << endl;

			// //Find the maximum score from all cells we can reach
			// mScore = max(mat[matPosI - 1][matPosJ], max(mat[matPosI - 1][matPosJ - 1], mat[matPosI][matPosJ - 1]));

			//Testing
			// cout << "Maximum score is " << mScore << endl;

			//Find the cell with the maximum score
			if(mat[matPosI - 1][matPosJ - 1]  + compUScore(uniSeq[posU], query[posQ]) == mat[matPosI][matPosJ]){
				//Testing
				// cout << "Let's do a diagonal step" << endl;
				// if(mat[matPosI - 1][matPosJ - 1] < mat[matPosI][matPosJ]){
				// 	// cout << "2 Option 1" << endl;
				// } else{
				// 	cout << "2 Option 2" << endl;
				// }

				//We can add a base from both sequences
				gSeq = gSeq + uniSeq[posU];
				qSeq = qSeq + query[posQ];

				// algn.aSeqG = uniSeq[posU] + algn.aSeqG;
				// algn.aSeqQ = query[posQ] + algn.aSeqQ;

				//Move to upper left cell
				--matPosI;
				--matPosJ;
				++posQ;
				++posU;
			} else if(mat[matPosI - 1][matPosJ] + GAP_SCORE == mat[matPosI][matPosJ]){
				//Testing
				// cout << "We shouldn't get here" << endl;
				// cout << "2 Option 4" << endl;

				//We can only walk on with a gap in the unitig sequence
				gSeq = gSeq + GAP_SYMB;
				qSeq = qSeq + query[posQ];
				//Move one cell to the left
				--matPosI;
				++posQ;
			} else{
				//Testing
				// cout << "We shouldn't get here" << endl;
				// cout << "2 Option 3" << endl;

				//We can only walk on with a gap in the query sequence
				gSeq = gSeq + uniSeq[posU];
				qSeq = qSeq + GAP_SYMB;
				//Move one cell up
				--matPosJ;
				++posU;
			}
		}
		//Testing
		// cout << "We got on cell further" << endl;
	}

	algn.aSeqG = gSeq + algn.aSeqG;
	algn.aSeqQ = qSeq + algn.aSeqQ;
}