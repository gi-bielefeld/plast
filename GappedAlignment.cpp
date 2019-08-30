#include "GappedAlignment.h"
#include "Hit.h"
#include "Search.h"

// #include "Extension.h"

//This function calculates a semi-global alignment of a unitig and the query sequence considering a quorum. Returns true if calculations on unitig could be finished (i.e. the unitig sequence was not too long to be stored inside the edit matrix).
bool calcSemiGlobAlignment(const UnitigColorMap<seedlist> &uni, const string &q, uint32_t& posU, uint32_t& posQ, const uint16_t &X, const uint32_t &maxGaps, uint32_t& maxPosQ, uint32_t& maxPosU, int32_t& maxScore, int32_t& eMax, const uint32_t &quorum){//TODO: So far we do not adjust the maximal number of gaps. This has to be changed as soon as we calculate real alignments!//TODO: This function does more than it should!
	bool termCalcs, matShriked = false;
	uint32_t dim, i, j, lCalcBorder, rCalcBorder, endBuf = 0, uCmpSeqLen, qCmpSeqLen, eMaxPosU = posU, eMaxPosQ = posQ;
	string uSeq;

	//Initialize the current max position
	eMax = -X;
	//Get unitig's sequence
	uSeq = uni.mappedSequenceToString();

	//Testing
	//cout << "uSeq:" << uSeq << "\nStart pos u: " << posU << "\nStart pos q: " << posQ << endl;

	//Check whether current unitig has a successor
	if(uni.getSuccessors().hasSuccessors()){
		endBuf = uni.getGraph()->getK() - 1;
	}

	//Testing
	//cout << "Initialisiere den Kram" << endl;

	//Calculate length of unitig sequence to compare
	uCmpSeqLen = uni.size - (posU + endBuf);

	//TODO This is just preliminar and has to be changed a.s.a.p.!
	//Check whether the unitig sequence which has to be compared fulfills the quorum
	if(!quorumFulfilled(uni, posU, uni.len, quorum)){
		//If the quorum is not fulfilled we are not going to do anything!
		uCmpSeqLen = 0;
	}

	//Calculate length of query sequence to compare
	qCmpSeqLen = q.length() - posQ;
	//Calculate dynamic programming matrix borders
	dim = min(uCmpSeqLen, qCmpSeqLen) + 1;

	//Testing
	//cout << "dim: " << dim << endl;

	//Make sure that the matrix is not becoming too big
	if(dim > MAX_MATRIX_SIZE){
		//Testing
		//cout << "Do we get here?" << endl;

		dim = MAX_MATRIX_SIZE;
		matShriked = true;
	}

	//Testing
	//cout << "Martixdimension:" << dim << endl;

	//Initialize dynamic programming matrix
	int32_t mat[dim][dim];

	//Initialize upper left cell
	mat[0][0] = maxScore;

	//If one of the to compare sequences is the empty sequence we need to set eMax correctly because we are not entering the main for loop at all
	if(dim == 1) eMax = mat[0][0];

	//Initialize first column
	for(i = 1; i < dim; ++i){
		if(i <= maxGaps){
			mat[i][0] = maxScore + i * GAP_SCORE;//GAP_SCORE is supposed to be negative//TODO: This has to be changed as soon as we are not using unit score anymore!
		} else{
			mat[i][0] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
		}
	}

	//Initialize first row
	for(j = 1; j < dim; ++j){
		if(j <= maxGaps){
			mat[0][j] = maxScore + j * GAP_SCORE;//GAP_SCORE is supposed to be negative//TODO: This has to be changed as soon as we are not using unit score anymore!
		} else{
			mat[0][j] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
		}
	}

	//Set borders
	lCalcBorder = 1;
	rCalcBorder = min(dim, maxGaps);

	//Testing
	//cout << "STarte Matrixberechnung" << endl;

	//Fill the matrix
	for(i = 1; i < dim; ++i){
		//Reset abort var
		termCalcs = true;

		//Check if we can adjust left alignment calculation border
		if(maxGaps < i){
			lCalcBorder = i - maxGaps;
		}

		//Adjust right border
		rCalcBorder = min(i + maxGaps, dim);

		for(j = 1; j < dim; ++j){
			//Make sure we calculate a bended alignment if matrix is too big
			if(j < lCalcBorder || j > rCalcBorder){
				mat[i][j] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
			} else{
				//Calculate current cells value and check we don't drop
				if((mat[i][j] = max(mat[i - 1][j - 1] + compUScore(q[posQ + i - 1], uSeq[posU + j - 1]), max(mat[i][j - 1] - 1, mat[i - 1][j] - 1))) > -X){//TODO: This has to be changed as soon as we are not using unit score anymore!
					//Make sure we do not terminate our calculations too soon
					termCalcs = false;

					//Check whether we have a new max score (NOTE: We could check only for better scores here as well. Do it this way will in doubt lead to larger alignments.)
					if(maxScore <= mat[i][j]){
						//Update max score
						maxScore = mat[i][j];
						//Update maximum position
						maxPosQ = posQ + i - 1;
						maxPosU = posU + j - 1;

						//Testing
						//cout << "maxPosQ inside calcSemiGlobAlignment: " << maxPosQ << endl;
					}

					//Testing
					//cout << mat[i][j] << endl;

					//Check whether we are at the right or lower matrix edge and whether we have found a new maximum
					if((i == dim - 1 || j == dim - 1) && eMax <= mat[i][j]){
						//Update edge maximum
						eMax = mat[i][j];
						eMaxPosQ = posQ + i - 1;
						eMaxPosU = posU + j - 1;
					}
				}
			}
		}

		//Check whether calculations may be terminated
		if(termCalcs) break;
	}

	//Testing
	/*cout << "Alignment matrix:" << endl;
	for(i = 0; i < dim; ++i){
		for(j = 0; j < dim; ++j){
			cout << mat[i][j] << " ";
		}
		cout << endl;
	}*/

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

//This function calculates a semi-global alignment of a unitig and the query sequence considering a quorum and a search color set. Returns true if calculations on unitig could be finished (i.e. the unitig sequence was not too long to be stored inside the edit matrix).
bool calcSemiGlobAlignment(const UnitigColorMap<seedlist> &uni, const string &q, uint32_t& posU, uint32_t& posQ, const uint16_t &X, const uint32_t &maxGaps, uint32_t& maxPosQ, uint32_t& maxPosU, int32_t& maxScore, int32_t& eMax, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){//TODO: So far we do not adjust the maximal number of gaps. This has to be changed as soon as we calculate real alignments!//TODO: This function does more than it should!
	bool termCalcs, matShriked = false;
	uint32_t dim, i, j, lCalcBorder, rCalcBorder, endBuf = 0, uCmpSeqLen, qCmpSeqLen, eMaxPosU = posU, eMaxPosQ = posQ;
	string uSeq;

	//Initialize the current max position
	eMax = -X;
	//Get unitig's sequence
	uSeq = uni.mappedSequenceToString();

	//Check whether current unitig has a successor
	if(uni.getSuccessors().hasSuccessors()){
		endBuf = uni.getGraph()->getK() - 1;
	}

	//Testing
	//cout << "Initialisiere den Kram" << endl;

	//Calculate length of unitig sequence to compare
	uCmpSeqLen = uni.size - (posU + endBuf);

	//Testing
	//cout << "8 Option ";

	//TODO This is just preliminar and has to be changed a.s.a.p.! Especially because it does not work for hits on the reverse complementary strand!
	//Check whether the unitig sequence which has to be compared fulfills the quorum
	if(!quorumFulfilled(uni, posU, uni.len, quorum, searchSet)){
		//Testing
		//cout << "not ";

		//If the quorum is not fulfilled we are not going to do anything!
		uCmpSeqLen = 0;
	}

	//Calculate length of query sequence to compare
	qCmpSeqLen = q.length() - posQ;
	//Calculate dynamic programming matrix borders
	dim = min(uCmpSeqLen, qCmpSeqLen) + 1;

	//Testing
	//cout << "dim: " << dim << endl;
	//cout << "1" << endl;

	//Make sure that the matrix is not becoming too big
	if(dim > MAX_MATRIX_SIZE){
		//Testing
		//cout << "Do we get here?" << endl;

		dim = MAX_MATRIX_SIZE;
		matShriked = true;
	}

	//Testing
	//cout << "Martixdimension:" << dim << endl;

	//Initialize dynamic programming matrix
	int32_t mat[dim][dim];

	//Initialize upper left cell
	mat[0][0] = maxScore;

	//If one of the two compared sequences is the empty sequence we need to set eMax correctly because we are not entering the main for loop at all
	if(dim == 1) eMax = mat[0][0];

	//Initialize first column
	for(i = 1; i < dim; ++i){
		if(i <= maxGaps){
			//GAP_SCORE is supposed to be negative
			mat[i][0] = maxScore + i * GAP_SCORE;//TODO: This has to be changed as soon as we are not using unit score anymore!
		} else{
			mat[i][0] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
		}
	}

	//Initialize first row
	for(j = 1; j < dim; ++j){
		if(j <= maxGaps){
			mat[0][j] = maxScore + j * GAP_SCORE;//GAP_SCORE is supposed to be negative//TODO: This has to be changed as soon as we are not using unit score anymore!
		} else{
			mat[0][j] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
		}
	}

	//Set borders
	lCalcBorder = 1;
	rCalcBorder = min(dim, maxGaps);

	//Testing
	//cout << "STarte Matrixberechnung" << endl;

	//Fill the matrix
	for(i = 1; i < dim; ++i){
		//Reset abort var
		termCalcs = true;

		//Check if we can adjust left alignment calculation border
		if(maxGaps < i){
			lCalcBorder = i - maxGaps;
		}

		//Adjust right border
		rCalcBorder = min(i + maxGaps, dim);

		for(j = 1; j < dim; ++j){
			//Make sure we calculate a bended alignment if matrix is too big
			if(j < lCalcBorder || j > rCalcBorder){
				mat[i][j] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
			} else{
				//Calculate current cells value and check we don't drop
				if((mat[i][j] = max(mat[i - 1][j - 1] + compUScore(q[posQ + i - 1], uSeq[posU + j - 1]), max(mat[i][j - 1] - 1, mat[i - 1][j] - 1))) > -X){//TODO: This has to be changed as soon as we are not using unit score anymore!
					//Make sure we do not terminate our calculations too soon
					termCalcs = false;

					//Check whether we have a new max score (NOTE: We could check only for better scores here as well. Do it this way will in doubt lead to larger alignments.)
					if(maxScore <= mat[i][j]){
						//Update max score
						maxScore = mat[i][j];
						//Update maximum position
						maxPosQ = posQ + i - 1;
						maxPosU = posU + j - 1;

						//Testing
						//cout << "maxPosQ inside calcSemiGlobAlignment: " << maxPosQ << endl;
					}

					//Testing
					//cout << mat[i][j] << endl;

					//Check whether we are at the right or lower matrix edge and whether we have found a new maximum
					if((i == dim - 1 || j == dim - 1) && eMax <= mat[i][j]){
						//Update edge maximum
						eMax = mat[i][j];
						eMaxPosQ = posQ + i - 1;
						eMaxPosU = posU + j - 1;
					}
				}
			}
		}

		//Check whether calculations may be terminated
		if(termCalcs) break;
	}

	//Testing
	// cout << "Alignment matrix:" << endl;
	// for(i = 0; i < dim; ++i){
	// 	for(j = 0; j < dim; ++j){
	// 		cout << mat[i][j] << " ";
	// 	}
	// 	cout << endl;
	// }

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
bool calcLeftGlobAlignment(const UnitigColorMap<seedlist> &uni, const string &q, uint32_t& posU, uint32_t& posQ, const uint16_t &X, const uint32_t &maxGaps, uint32_t& maxPosQ, uint32_t& maxPosU, int32_t& maxScore, int32_t& eMax, const uint32_t &quorum){//TODO: So far we do not adjust the maximal number of gaps. This has to be changed as soon as we calculate real alignments!//TODO: This function does more than it should!
	bool termCalcs, matShriked = false;
	uint32_t dim, i, j, lCalcBorder, rCalcBorder, eMaxPosU = 0, eMaxPosQ = 0;
	string uSeq;

	//Get unitig's sequence
	uSeq = uni.mappedSequenceToString();

	//Testing
	//cout << "Unitig sequence: " << uSeq << endl;

	//Calculate dynamic programming matrix borders
	dim = min(posU, posQ) + 2;//+2 since we need 1 column and row for the empty string and offset indexes start with 0

	//TODO This is just preliminar and has to be changed a.s.a.p.!
	//Check whether the unitig sequence which has to be compared fulfills the quorum
	if(!quorumFulfilled(uni, 0, posU + 1, quorum)){
		//If the quorum is not fulfilled we are not going to do anything!
		dim = 1;
	}

	//Make sure that the matrix is not becoming too big
	if(dim > MAX_MATRIX_SIZE){
		dim = MAX_MATRIX_SIZE;
		matShriked = true;
	}

	//Initialize dynamic programming matrix
	int32_t mat[dim][dim];
	//Initialize upper left cell
	mat[0][0] = maxScore;
	//Initialize edge maximum
	eMax = -X;

	//Initialize first column
	for(i = 1; i < dim; ++i){
		if(i <= maxGaps){
			mat[i][0] = maxScore + i * GAP_SCORE;//GAP_SCORE is supposed to be negative//TODO: This has to be changed as soon as we are not using unit score anymore!
		} else{
			mat[i][0] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
		}
	}

	//Initialize first row
	for(j = 1; j < dim; ++j){
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
	for(i = 1; i < dim; ++i){
		//Reset abort var
		termCalcs = true;

		//Check if we can adjust left alignment calculation border
		if(maxGaps < i){
			lCalcBorder = i - maxGaps;
		}

		//Adjust right border
		rCalcBorder = min(i + maxGaps, dim);

		for(j = 1; j < dim; ++j){
			//Make sure we calculate a banded alignment if matrix is too big
			if(j < lCalcBorder || j > rCalcBorder){
				mat[i][j] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
			} else{
				//Calculate current cells value and check we don't drop
				if((mat[i][j] = max(mat[i - 1][j - 1] + compUScore(q[posQ - i + 1], uSeq[posU - j + 1]), max(mat[i][j - 1] - 1, mat[i - 1][j] - 1))) > -X){//TODO: This has to be changed as soon as we are not using unit score anymore!
					//Make sure we do not terminate our calculations too soon
					termCalcs = false;

					//Testing
					//cout << "Inside calcLeftGlobAlignment: Characters that have been compared were " << q[posQ - i + 1] << " in the query and " << uSeq[posU - j + 1] << " in the unitig" << endl;

					//Check whether we have a new max score (NOTE: We could check only for better scores here as well. Do it this way will in doubt lead to larger alignments.)//TODO: Maybe we don't what this because this might lead to less significant alignments!
					if(maxScore <= mat[i][j]){
						//Update max score
						maxScore = mat[i][j];
						//Update maximum position
						maxPosQ = posQ - i + 1;
						maxPosU = posU - j + 1;
					}

					//Check whether we are at the right or lower matrix edge and whether we have found a new maximum
					if((i == dim - 1 || j == dim - 1) && eMax <= mat[i][j]){
						//Update edge maximum
						eMax = mat[i][j];
						eMaxPosQ = posQ - i + 1;
						eMaxPosU = posU - j + 1;
					}
				}
			}
		}

		//Check whether calculations may be terminated
		if(termCalcs) break;
	}

	//Testing
	/*cout << "Filled matrix:" << endl;
	for(i = 0; i < dim; ++i){
		for(j = 0; j < dim; ++j){
			cout << mat[i][j] << " ";
		}
		cout << endl;
	}*/

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
bool calcLeftGlobAlignment(const UnitigColorMap<seedlist> &uni, const string &q, uint32_t& posU, uint32_t& posQ, const uint16_t &X, const uint32_t &maxGaps, uint32_t& maxPosQ, uint32_t& maxPosU, int32_t& maxScore, int32_t& eMax, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){//TODO: So far we do not adjust the maximal number of gaps. This has to be changed as soon as we calculate real alignments!//TODO: This function does more than it should!
	bool termCalcs, matShriked = false;
	uint32_t dim, i, j, lCalcBorder, rCalcBorder, eMaxPosU = 0, eMaxPosQ = 0;
	string uSeq;

	//Get unitig's sequence
	uSeq = uni.mappedSequenceToString();

	//Testing
	//cout << "Unitig sequence: " << uSeq << endl;

	//Calculate dynamic programming matrix borders
	dim = min(posU, posQ) + 2;//+2 since we need 1 column and row for the empty string and offset indexes start with 0

	//Testing
	// cout << "1 Option";
	
	//TODO This is just preliminar and has to be changed a.s.a.p.! Especially because it does not work for hits on the reverse complementary strand!
	//Check whether the unitig sequence which has to be compared fulfills the quorum
	if(!quorumFulfilled(uni, 0, posU + 1, quorum, searchSet)){
		//Testing
		// cout << "not ";

		//If the quorum is not fulfilled we are not going to do anything!
		dim = 1;
	}

	//Testing
	// cout << "1" << endl;

	//Make sure that the matrix is not becoming too big
	if(dim > MAX_MATRIX_SIZE){
		dim = MAX_MATRIX_SIZE;
		matShriked = true;
	}

	//Initialize dynamic programming matrix
	int32_t mat[dim][dim];
	//Initialize upper left cell
	mat[0][0] = maxScore;
	//Initialize edge maximum
	eMax = -X;

	//Initialize first column
	for(i = 1; i < dim; ++i){
		if(i <= maxGaps){
			mat[i][0] = maxScore + i * GAP_SCORE;//GAP_SCORE is supposed to be negative//TODO: This has to be changed as soon as we are not using unit score anymore!
		} else{
			mat[i][0] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
		}
	}

	//Initialize first row
	for(j = 1; j < dim; ++j){
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
	for(i = 1; i < dim; ++i){
		//Reset abort var
		termCalcs = true;

		//Check if we can adjust left alignment calculation border
		if(maxGaps < i){
			lCalcBorder = i - maxGaps;
		}

		//Adjust right border
		rCalcBorder = min(i + maxGaps, dim);

		for(j = 1; j < dim; ++j){
			//Make sure we calculate a banded alignment if matrix is too big
			if(j < lCalcBorder || j > rCalcBorder){
				mat[i][j] = INT32_MIN + 1;//TODO: This has to be changed as soon as we are not using unit score anymore!
			} else{
				//Calculate current cells value and check we don't drop
				if((mat[i][j] = max(mat[i - 1][j - 1] + compUScore(q[posQ - i + 1], uSeq[posU - j + 1]), max(mat[i][j - 1] - 1, mat[i - 1][j] - 1))) > -X){//TODO: This has to be changed as soon as we are not using unit score anymore!
					//Make sure we do not terminate our calculations too soon
					termCalcs = false;

					//Testing
					//cout << "Inside calcLeftGlobAlignment: Characters that have been compared were " << q[posQ - i + 1] << " in the query and " << uSeq[posU - j + 1] << " in the unitig" << endl;

					//Check whether we have a new max score (NOTE: We could check only for better scores here as well. Do it this way will in doubt lead to larger alignments.)//TODO: Maybe we don't what this because this might lead to less significant alignments!
					if(maxScore <= mat[i][j]){
						//Update max score
						maxScore = mat[i][j];
						//Update maximum position
						maxPosQ = posQ - i + 1;
						maxPosU = posU - j + 1;
					}

					//Check whether we are at the right or lower matrix edge and whether we have found a new maximum
					if((i == dim - 1 || j == dim - 1) && eMax <= mat[i][j]){
						//Update edge maximum
						eMax = mat[i][j];
						eMaxPosQ = posQ - i + 1;
						eMaxPosU = posU - j + 1;
					}
				}
			}
		}

		//Check whether calculations may be terminated
		if(termCalcs) break;
	}

	//Testing
	/*cout << "Filled matrix:" << endl;
	for(i = 0; i < dim; ++i){
		for(j = 0; j < dim; ++j){
			cout << mat[i][j] << " ";
		}
		cout << endl;
	}*/

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

	if(eMax > -X && matShriked) return false;
	
	return true;
}

//This function calculates the continuation of a gapped alignment on a successive unitig or the next peace of the query considering a quorum
void contRightGappedAlignment(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &posQ, uint32_t &posU, const uint16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum){
	bool explSuc = false, unitigDone;
	uint32_t maxPosQ = posQ, maxPosU = posU;
	int32_t maxBorderScore = -X;
	UnitigColorMap<seedlist> maxUni = uni;

	//
	//++posQ;
	//Calculate gapped alignment
	unitigDone = calcSemiGlobAlignment(uni, q, posU, posQ, X, maxGaps, maxPosQ, maxPosU, score, maxBorderScore, quorum);

	//Testing
	//cout << "posU: " << posU << "\nposQ: " << posQ << endl;
	//cout << "Inside contRightGappedAlignment after calculation of gapped alignment on the current unitig:\nunitigDone: " << unitigDone << "\nmaxBorderScore: " << maxBorderScore << "\nposU: " << posU << "\nposQ: " << posQ << "\nmaxPosQ: " << maxPosQ << "\nmaxPosU: " << maxPosU << "\nscore: " << score << "\nQuery length: " << q.length() << endl;

	//Check outcome of alignment calculation
	if(!unitigDone){
		explSuc = contGappedOnSameUni(uni, q, posQ, posU, X, maxGaps, maxBorderScore, explCount, quorum);
	} else if(maxBorderScore > -X && posQ < q.length() - 1){
		//Continue calculations on the next unitig
		explSuc = contGappedOnSuccUni(uni, q, posQ, posU, X, maxGaps, maxBorderScore, explCount, quorum);
	}

	//Check whether exploration of next unitig was successful (this is only the case if a new maximum score has been found as well)
	if(explSuc && score < maxBorderScore){
		//Score has to be updated if we have found a better one on a successive unitig
		score = maxBorderScore;
	} else{
		//Save new offsets in hit
		posU = maxPosU;
		posQ = maxPosQ;
		//Reset the unitig on which the maximum has been found
		uni = maxUni;
	}
}

//This function calculates the continuation of a gapped alignment on a successive unitig or the next peace of the query considering a quorum and a search color set
void contRightGappedAlignment(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &posQ, uint32_t &posU, const uint16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	bool explSuc = false, unitigDone;
	uint32_t maxPosQ = posQ, maxPosU = posU;
	int32_t maxBorderScore = -X;
	UnitigColorMap<seedlist> maxUni = uni;

	//
	//++posQ;
	//Calculate gapped alignment
	unitigDone = calcSemiGlobAlignment(uni, q, posU, posQ, X, maxGaps, maxPosQ, maxPosU, score, maxBorderScore, quorum, searchSet);

	//Testing
	//cout << "posU: " << posU << "\nposQ: " << posQ << endl;
	//cout << "Inside contRightGappedAlignment after calculation of gapped alignment on the current unitig:\nunitigDone: " << unitigDone << "\nmaxBorderScore: " << maxBorderScore << "\nposU: " << posU << "\nposQ: " << posQ << "\nmaxPosQ: " << maxPosQ << "\nmaxPosU: " << maxPosU << "\nscore: " << score << "\nQuery length: " << q.length() << endl;

	//Check outcome of alignment calculation
	if(!unitigDone){
		explSuc = contGappedOnSameUni(uni, q, posQ, posU, X, maxGaps, maxBorderScore, explCount, quorum, searchSet);
	} else if(maxBorderScore > -X && posQ < q.length() - 1){
		//Continue calculations on the next unitig
		explSuc = contGappedOnSuccUni(uni, q, posQ, posU, X, maxGaps, maxBorderScore, explCount, quorum, searchSet);
	}

	//Check whether exploration of next unitig was successful (this is only the case if a new maximum score has been found as well)
	if(explSuc && score < maxBorderScore){
		//Score has to be updated if we have found a better one on a successive unitig
		score = maxBorderScore;
	} else{
		//Save new offsets in hit
		posU = maxPosU;
		posQ = maxPosQ;
		//Reset the unitig on which the maximum has been found
		uni = maxUni;
	}
}

//This function calculates the continuation of a gapped alignment on a predecessive unitig or the next peace of the query considering a quorum
void contLeftGappedAlignment(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &posQ, uint32_t &posU, const uint16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum){
	bool unitigDone, explSuc = false;
	uint32_t maxPosQ = posQ, maxPosU = posU;
	int32_t maxBorderScore;
	UnitigColorMap<seedlist> maxUni = uni;

	//Calculate gapped alignment
	unitigDone = calcLeftGlobAlignment(uni, q, posU, posQ, X, maxGaps, maxPosQ, maxPosU, score, maxBorderScore, quorum);

	//Testing
	//cout << "Exploration counter: " << explCount << " posQ: " << posQ << endl;

	//Check outcome of alignment calculation
	if(!unitigDone){
		explSuc = contLeftOnSameUni(uni, q, posQ, posU, X, maxGaps, maxBorderScore, explCount, quorum);
	} else if(maxBorderScore > -X && posQ > 0){
		//Testing
		//cout << "Call of contGappedOnPredUni" << endl;

		//Continue calculations on the next unitig
		explSuc = contGappedOnPredUni(uni, q, posQ, posU, X, maxGaps, maxBorderScore, explCount, quorum);
	}

	//Check whether exploration of next unitig was successful (this is only the case if a new maximum score has been found as well)
	if(explSuc && score < maxBorderScore){
		//Score has to be updated if we have found a better one on a successive unitig
		score = maxBorderScore;
	} else{
		//Save new offsets in hit
		posU = maxPosU;
		posQ = maxPosQ;
		//Reset the unitig on which the maximum has been found
		uni = maxUni;
	}
}

//This function calculates the continuation of a gapped alignment on a predecessive unitig or the next peace of the query considering a quorum and a search color set
void contLeftGappedAlignment(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &posQ, uint32_t &posU, const uint16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	bool unitigDone, explSuc = false;
	uint32_t maxPosQ = posQ, maxPosU = posU;
	int32_t maxBorderScore;
	UnitigColorMap<seedlist> maxUni = uni;

	//Calculate gapped alignment
	unitigDone = calcLeftGlobAlignment(uni, q, posU, posQ, X, maxGaps, maxPosQ, maxPosU, score, maxBorderScore, quorum, searchSet);

	//Testing
	//cout << "Exploration counter: " << explCount << " posQ: " << posQ << endl;

	//Check outcome of alignment calculation
	if(!unitigDone){
		explSuc = contLeftOnSameUni(uni, q, posQ, posU, X, maxGaps, maxBorderScore, explCount, quorum, searchSet);
	} else if(maxBorderScore > -X && posQ > 0){
		//Testing
		//cout << "Call of contGappedOnPredUni" << endl;

		//Continue calculations on the next unitig
		explSuc = contGappedOnPredUni(uni, q, posQ, posU, X, maxGaps, maxBorderScore, explCount, quorum, searchSet);
	}

	//Check whether exploration of next unitig was successful (this is only the case if a new maximum score has been found as well)
	if(explSuc && score < maxBorderScore){
		//Score has to be updated if we have found a better one on a successive unitig
		score = maxBorderScore;
	} else{
		//Save new offsets in hit
		posU = maxPosU;
		posQ = maxPosQ;
		//Reset the unitig on which the maximum has been found
		uni = maxUni;
	}
}

//Continues a gapped alignment on the same unitig as before considering a quorum. Returns true if a new maximal score has been found.
bool contGappedOnSameUni(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum){
	uint32_t tmpQOff = qOff + 1, tmpUOff = uOff + 1;
	int32_t tmpScore = score;

	//Testing
	//cout << "Exploration on same unitig" << endl;

	//++qOff;
	//++uOff;
	contRightGappedAlignment(uni, q, tmpQOff, tmpUOff, X, maxGaps, tmpScore, explCount, quorum);

	//Check whether we have found a new maximum
	if(score < tmpScore){
		//Update maximum
		score = tmpScore;
		qOff = tmpQOff;
		uOff = tmpUOff;
		return true;
	}

	return false;
}

//Continues a gapped alignment on the same unitig as before considering a quorum and a search color set. Returns true if a new maximal score has been found.
bool contGappedOnSameUni(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	uint32_t tmpQOff = qOff + 1, tmpUOff = uOff + 1;
	int32_t tmpScore = score;

	//Testing
	//cout << "Exploration on same unitig" << endl;

	//++qOff;
	//++uOff;
	contRightGappedAlignment(uni, q, tmpQOff, tmpUOff, X, maxGaps, tmpScore, explCount, quorum, searchSet);

	//Check whether we have found a new maximum
	if(score < tmpScore){
		//Update maximum
		score = tmpScore;
		qOff = tmpQOff;
		uOff = tmpUOff;
		return true;
	}

	return false;
}

//This function continues a gapped alignment calculation to the left on the same unitig considering a quorum. Returns true if a new maximal score has been found.
bool contLeftOnSameUni(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum){
	uint32_t tmpQOff = qOff - 1, tmpUOff = uOff - 1;
	int32_t tmpScore = score;

	contLeftGappedAlignment(uni, q, tmpQOff, tmpUOff, X, maxGaps, tmpScore, explCount, quorum);

	//Check whether we have found a new maximum
	if(score < tmpScore){
		//Update maximum
		score = tmpScore;
		qOff = tmpQOff;
		uOff = tmpUOff;
		return true;
	}

	return false;
}

//This function continues a gapped alignment calculation to the left on the same unitig considering a quorum and a search color set. Returns true if a new maximal score has been found.
bool contLeftOnSameUni(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	uint32_t tmpQOff = qOff - 1, tmpUOff = uOff - 1;
	int32_t tmpScore = score;

	contLeftGappedAlignment(uni, q, tmpQOff, tmpUOff, X, maxGaps, tmpScore, explCount, quorum, searchSet);

	//Check whether we have found a new maximum
	if(score < tmpScore){
		//Update maximum
		score = tmpScore;
		qOff = tmpQOff;
		uOff = tmpUOff;
		return true;
	}

	return false;
}

//Initiates a gapped alignment calculation on all successive unitigs considering a quorum and returns the alignment's end position in the unitig
bool contGappedOnSuccUni(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum){
	bool suc = false;
	uint32_t tmpPosU, tmpPosQ, maxPosQ, maxPosU;
	int32_t maxScore, tmpScore;
	UnitigColorMap<seedlist> sucUni;
	ForwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> sucIter = uni.getSuccessors();

	//Check whether we have reached the maximum recursion depth of an extension
	if(++explCount > MAXRECURSIONDEPTH){
		//Report this incident//TODO: If this happens the condition that an extended seed cannot be reached anymore is violated. Implement a procedure after the extension that merges such seeds with other seeds that would have reached them!
		//cerr << "Maximum recursion depth reached during extension. Position in q: " << qOff << endl;
		//Terminate this extension
		return suc;
	}

	//Check whether there are successors to continue
	if(sucIter.hasSuccessors()){
		maxScore = score;

		//Explore each neighbor
		for(neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> nI = sucIter.begin(); nI != sucIter.end(); ++nI){
			//Get next unitig
			sucUni = *nI;
			tmpPosQ = qOff + 1;
			tmpPosU = 0;
			tmpScore = score;
			//Calculate gapped alignment on the next unitig
			contRightGappedAlignment(sucUni, q, tmpPosQ, tmpPosU, X, maxGaps, tmpScore, explCount, quorum);

			//Check for new maximum
			if(maxScore < tmpScore){
				//Update maximum
				maxScore = tmpScore;
				maxPosU = tmpPosU;
				maxPosQ = tmpPosQ;
				uni = sucUni;
				suc = true;
			}
		}

		//Only update score and qOff if exploration was successful
		if(suc){
			//Give back best result found
			score = maxScore;
			qOff = maxPosQ;
			uOff = maxPosU;
		}
	}

	return suc;
}

//Initiates a gapped alignment calculation on all successive unitigs considering a quorum and a search color set and returns the alignment's end position in the unitig
bool contGappedOnSuccUni(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	bool suc = false;
	uint32_t tmpPosU, tmpPosQ, maxPosQ, maxPosU;
	int32_t maxScore, tmpScore;
	UnitigColorMap<seedlist> sucUni;
	ForwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> sucIter = uni.getSuccessors();

	//Check whether we have reached the maximum recursion depth of an extension
	if(++explCount > MAXRECURSIONDEPTH){
		//Report this incident//TODO: If this happens the condition that an extended seed cannot be reached anymore is violated. Implement a procedure after the extension that merges such seeds with other seed that would have reached them!
		//cerr << "Maximum recursion depth reached during extension. Position in q: " << qOff << endl;
		//Terminate this extension
		return suc;
	}

	//Check whether there are successors to continue
	if(sucIter.hasSuccessors()){
		maxScore = score;

		//Explore each neighbor
		for(neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> nI = sucIter.begin(); nI != sucIter.end(); ++nI){
			//Get next unitig
			sucUni = *nI;
			tmpPosQ = qOff + 1;
			tmpPosU = 0;
			tmpScore = score;
			//Calculate gapped alignment on the next unitig
			contRightGappedAlignment(sucUni, q, tmpPosQ, tmpPosU, X, maxGaps, tmpScore, explCount, quorum, searchSet);

			//Check for new maximum
			if(maxScore < tmpScore){
				//Update maximum
				maxScore = tmpScore;
				maxPosU = tmpPosU;
				maxPosQ = tmpPosQ;
				uni = sucUni;
				suc = true;
			}
		}

		//Only update score and qOff if exploration was successful
		if(suc){
			//Give back best result found
			score = maxScore;
			qOff = maxPosQ;
			uOff = maxPosU;
		}
	}

	return suc;
}

//This function initiates a gapped alignment calculation on all predecessive unitigs considering a quorum. Returns true if calculations have been successful
bool contGappedOnPredUni(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum){
	bool suc = false;
	uint32_t tmpPosU, tmpPosQ, maxPosQ, maxPosU;
	int32_t maxScore, tmpScore;
	UnitigColorMap<seedlist> predUni;
	BackwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> predIter = uni.getPredecessors();

	//Testing
	//cout << "Last unitig was " << uni.toString() << endl;

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

		//Explore each neighbor
		for(neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> nI = predIter.begin(); nI != predIter.end(); ++nI){
			//Get next unitig
			predUni = *nI;
			tmpPosQ = qOff - 1;
			tmpPosU = predUni.size - predUni.getGraph()->getK();
			tmpScore = score;

			//Testing
			//cout << "\nNext unitig is " << predUni.toString() << endl;

			//Calculate gapped alignment on the predecessive unitig
			contLeftGappedAlignment(predUni, q, tmpPosQ, tmpPosU, X, maxGaps, tmpScore, explCount, quorum);

			//Check for new maximum
			if(maxScore < tmpScore){
				//Update maximum
				maxScore = tmpScore;
				maxPosU = tmpPosU;
				maxPosQ = tmpPosQ;
				uni = predUni;
				suc = true;
			}
		}

		//Only update score and qOff if exploration was successful
		if(suc){
			//Give back best result found
			score = maxScore;
			qOff = maxPosQ;
			uOff = maxPosU;
		}
	}

	return suc;
}

//This function initiates a gapped alignment calculation on all predecessive unitigs considering a quorum and a search color set. Returns true if calculations have been successful
bool contGappedOnPredUni(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	bool suc = false;
	uint32_t tmpPosU, tmpPosQ, maxPosQ, maxPosU;
	int32_t maxScore, tmpScore;
	UnitigColorMap<seedlist> predUni;
	BackwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> predIter = uni.getPredecessors();

	//Testing
	//cout << "Last unitig was " << uni.toString() << endl;

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

		//Explore each neighbor
		for(neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> nI = predIter.begin(); nI != predIter.end(); ++nI){
			//Get next unitig
			predUni = *nI;
			tmpPosQ = qOff - 1;
			tmpPosU = predUni.size - predUni.getGraph()->getK();
			tmpScore = score;

			//Testing
			//cout << "\nNext unitig is " << predUni.toString() << endl;

			//Calculate gapped alignment on the predecessive unitig
			contLeftGappedAlignment(predUni, q, tmpPosQ, tmpPosU, X, maxGaps, tmpScore, explCount, quorum, searchSet);

			//Check for new maximum
			if(maxScore < tmpScore){
				//Update maximum
				maxScore = tmpScore;
				maxPosU = tmpPosU;
				maxPosQ = tmpPosQ;
				uni = predUni;
				suc = true;
			}
		}

		//Only update score and qOff if exploration was successful
		if(suc){
			//Give back best result found
			score = maxScore;
			qOff = maxPosQ;
			uOff = maxPosU;
		}
	}

	return suc;
}

//This function calculates a gapped alignment to the right side of the starting position considering a quorum. ATTENTION: Hit's length attribute will be deprecated after function call!
void startRightGappedAlignment(struct hit *h, const string &q, const uint16_t &X, const uint32_t maxGaps, const uint32_t &quorum){//TODO: Why don't we need an extension pointer here?
	bool explSuc = false, unitigDone = true;
	uint32_t posU = h->rSeedUoff, posQ = h->rSeedQoff, maxPosU = h->rSeedUoff, maxPosQ = h->rSeedQoff, explCount = 0, endBuf = 0;
	int32_t maxScore = 0, maxBorderScore;
	UnitigColorMap<seedlist> currUni = h->rUnitig;

	//Check whether current unitig has a successor
	if(currUni.getSuccessors().hasSuccessors()){
		endBuf = currUni.getGraph()->getK();
	}

	//Check whether we can start the calculation right on this unitig or whether we are already at the end of unitig or query sequence
	if(posQ < q.length() - 1 && posU < currUni.referenceUnitigToString().length() - endBuf){
		//Set alignment start positions (+1, because we want to start right behind the hit)
		++posU;
		++posQ;
		//Calculate gapped alignment
		unitigDone = calcSemiGlobAlignment(currUni, q, posU, posQ, X, maxGaps, maxPosQ, maxPosU, maxScore, maxBorderScore, quorum);
	} else{
		//Set variable for exploration of successive unitigs
		maxBorderScore = 0;
	}

	//Check whether we have reached the end of the current unitig
	if(!unitigDone){
		explSuc = contGappedOnSameUni(h->rUnitig, q, posQ, posU, X, maxGaps, maxBorderScore, explCount, quorum);
	} else if(maxBorderScore > -X && posQ < q.length() - 1){
		//Testing
		//cout << "Gehe zum nächsten unitig" << endl;

		//Continue calculations on the next unitig
		explSuc = contGappedOnSuccUni(h->rUnitig, q, posQ, posU, X, maxGaps, maxBorderScore, explCount, quorum);

		//Testing
		//cout << "Weiter zum Rest" << endl;
	}

	//Check whether exploration of next unitig was successful
	if(explSuc && maxScore < maxBorderScore){
		//Save new offsets in hit
		h->rSeedUoff = posU;
		h->rSeedQoff = posQ;
		//Save the achieved score
		h->score += maxBorderScore;
	} else{
		//Save new offsets in hit
		h->rSeedUoff = maxPosU;
		h->rSeedQoff = maxPosQ;
		//Save the achieved score
		h->score += maxScore;
		//Reset right border unitig
		h->rUnitig = currUni;
	}
}

//This function calculates a gapped alignment to the right side of the starting position considering a quorum and a search color set. ATTENTION: Hit's length attribute will be deprecated after function call!
void startRightGappedAlignment(struct hit *h, const string &q, const uint16_t &X, const uint32_t maxGaps, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	//Testing
	// cout << "We are getting here" << endl;

	bool explSuc = false, unitigDone = true;
	uint32_t posU = h->rSeedUoff, posQ = h->rSeedQoff, maxPosU = h->rSeedUoff, maxPosQ = h->rSeedQoff, explCount = 0, endBuf = 0;
	int32_t maxScore = 0, maxBorderScore;
	UnitigColorMap<seedlist> currUni = h->rUnitig;

	//Testing
	// cout << "2 Option ";

	//Check whether current unitig has a successor
	if(currUni.getSuccessors().hasSuccessors()){
		//Testing
		// cout << "not ";

		endBuf = currUni.getGraph()->getK();
	}

	//Testing
	// cout << "2" << endl;

	//Check whether we can start the calculation right on this unitig or whether we are already at the end of unitig or query sequence
	if(posQ < q.length() - 1 && posU < currUni.referenceUnitigToString().length() - endBuf){
		//Set alignment start positions (+1, because we want to start right behind the hit)
		++posU;
		++posQ;

		//Testing
		// cout << "1 Option 1" << endl;

		//Calculate gapped alignment
		unitigDone = calcSemiGlobAlignment(currUni, q, posU, posQ, X, maxGaps, maxPosQ, maxPosU, maxScore, maxBorderScore, quorum, searchSet);
	} else{
		//Testing
		// cout << "1 Option 2" << endl;

		//Set variable for exploration of successive unitigs
		maxBorderScore = 0;
	}

	//Testing
	//cout << "Inside startRightGappedAlignment after calculation of gapped alignment on the first unitig:\nunitigDone: " << unitigDone << "\nmaxBorderScore: " << maxBorderScore << "\nposU: " << posU << "\nposQ: " << posQ << "\nmaxPosQ: " << maxPosQ << "\nmaxPosU: " << maxPosU << "\nmaxScore: " << maxScore << endl;

	//Check whether we have reached the end of the current unitig
	if(!unitigDone){
		explSuc = contGappedOnSameUni(h->rUnitig, q, posQ, posU, X, maxGaps, maxBorderScore, explCount, quorum, searchSet);
	} else if(maxBorderScore > -X && posQ < q.length() - 1){
		//Testing
		//cout << "Gehe zum nächsten unitig" << endl;

		//Continue calculations on the next unitig
		explSuc = contGappedOnSuccUni(h->rUnitig, q, posQ, posU, X, maxGaps, maxBorderScore, explCount, quorum, searchSet);

		//Testing
		//cout << "Weiter zum Rest" << endl;
	}

	//Check whether exploration of next unitig was successful
	if(explSuc && maxScore < maxBorderScore){
		//Save new offsets in hit
		h->rSeedUoff = posU;
		h->rSeedQoff = posQ;
		//Save the achieved score
		h->score += maxBorderScore;
	} else{
		//Save new offsets in hit
		h->rSeedUoff = maxPosU;
		h->rSeedQoff = maxPosQ;
		//Save the achieved score
		h->score += maxScore;
		//Reset right border unitig
		h->rUnitig = currUni;
	}
}

//This function calculates a gapped alignment to the left side of the starting position considering a quorum. ATTENTION: Hit's length attribute will be deprecated after function call!
void startLeftGappedAlignment(struct hit *h, const string &q, const uint16_t &X, const uint32_t maxGaps, const uint32_t &quorum){
	bool explSuc = false, unitigDone = true;
	uint32_t posU = h->lSeedUoff, posQ = h->lSeedQoff, maxPosU = h->lSeedUoff, maxPosQ = h->lSeedQoff, explCount = 0;
	int32_t maxScore = 0, maxBorderScore;
	UnitigColorMap<seedlist> currUni = h->lUnitig;

	//Check whether we can start the calculation right on this unitig or whether we are already at the begin of unitig or query sequence
	if(posQ > 0 && posU > 0){
		//Testing
		//cout << "Do we get here?" << endl;

		//Set alignment start positions (-1, because we want to start right in front of the hit)
		--posU;
		--posQ;
		//Calculate gapped alignment
		unitigDone = calcLeftGlobAlignment(h->lUnitig, q, posU, posQ, X, maxGaps, maxPosQ, maxPosU, maxScore, maxBorderScore, quorum);
	} else{
		//Set variable for exploration of predecessive unitigs
		maxBorderScore = 0;
	}

	//Check whether we have reached the beginning of the current unitig
	if(!unitigDone){
		explSuc = contLeftOnSameUni(h->lUnitig, q, posQ, posU, X, maxGaps, maxBorderScore, explCount, quorum);
	} else if(maxBorderScore > -X && posQ > 0){
		//Continue calculations on the next unitig
		explSuc = contGappedOnPredUni(h->lUnitig, q, posQ, posU, X, maxGaps, maxBorderScore, explCount, quorum);
	}

	//Check whether exploration of predecessive unitig was successful
	if(explSuc && maxScore < maxBorderScore){
		//Save new offsets in hit
		h->lSeedUoff = posU;
		h->lSeedQoff = posQ;
		//Save the achieved score
		h->score += maxBorderScore;
	} else{
		//Save new offsets in hit
		h->lSeedUoff = maxPosU;
		h->lSeedQoff = maxPosQ;
		//Save the achieved score
		h->score += maxScore;
		//Reset left border unitig
		h->lUnitig = currUni;
	}
}

//This function calculates a gapped alignment to the left side of the starting position considering a quorum and a search color set. ATTENTION: Hit's length attribute will be deprecated after function call!
void startLeftGappedAlignment(struct hit *h, const string &q, const uint16_t &X, const uint32_t maxGaps, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	bool explSuc = false, unitigDone = true;
	uint32_t posU = h->lSeedUoff, posQ = h->lSeedQoff, maxPosU = h->lSeedUoff, maxPosQ = h->lSeedQoff, explCount = 0;
	int32_t maxScore = 0, maxBorderScore;
	UnitigColorMap<seedlist> currUni = h->lUnitig;

	//Check whether we can start the calculation right on this unitig or whether we are already at the begin of unitig or query sequence
	if(posQ > 0 && posU > 0){
		//Testing
		//cout << "Do we get here?" << endl;

		//Set alignment start positions (-1, because we want to start right in front of the hit)
		--posU;
		--posQ;
		//Calculate gapped alignment
		unitigDone = calcLeftGlobAlignment(h->lUnitig, q, posU, posQ, X, maxGaps, maxPosQ, maxPosU, maxScore, maxBorderScore, quorum, searchSet);
	} else{
		//Set variable for exploration of predecessive unitigs
		maxBorderScore = 0;
	}

	//Check whether we have reached the beginning of the current unitig
	if(!unitigDone){
		explSuc = contLeftOnSameUni(h->lUnitig, q, posQ, posU, X, maxGaps, maxBorderScore, explCount, quorum, searchSet);
	} else if(maxBorderScore > -X && posQ > 0){
		//Continue calculations on the next unitig
		explSuc = contGappedOnPredUni(h->lUnitig, q, posQ, posU, X, maxGaps, maxBorderScore, explCount, quorum, searchSet);
	}

	//Check whether exploration of predecessive unitig was successful
	if(explSuc && maxScore < maxBorderScore){
		//Save new offsets in hit
		h->lSeedUoff = posU;
		h->lSeedQoff = posQ;
		//Save the achieved score
		h->score += maxBorderScore;
	} else{
		//Save new offsets in hit
		h->lSeedUoff = maxPosU;
		h->lSeedQoff = maxPosQ;
		//Save the achieved score
		h->score += maxScore;
		//Reset left border unitig
		h->lUnitig = currUni;
	}
}