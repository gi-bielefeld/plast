#ifndef GAPPEDALIGN_HPP
#define GAPPEDALIGN_HPP

#define MAX_MATRIX_SIZE 1000
#define GAP_SCORE -2
#define GAP_RATIO 20
// #define GAP_RATIO 10
#define GAP_SYMB "-"

//This function calculates a semi-global alignment of a unitig and the query sequence considering a quorum and a search color set. Returns true if calculations on unitig could be finished (i.e. the unitig sequence was not to long to be stored inside the edit matrix).
bool calcSemiGlobAlignment(const UnitigColorMap<seedlist> &uni, const string &q, uint32_t& posU, uint32_t& posQ, const uint16_t &X, const uint32_t &maxGaps, uint32_t& maxPosQ, uint32_t& maxPosU, struct Algn &maxAlgn, struct Algn &brdAlgn, int32_t& maxScore, int32_t& eMax, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool srchCritCheck);

//This function calculates a banded alignment of a unitig and the query sequence to the left considering a quorum and a search color set. Returns true if calculations on unitig could be finished (i.e. the unitig sequence was not to long to be stored inside the edit matrix).
bool calcLeftGlobAlignment(const UnitigColorMap<seedlist> &uni, const string &q, uint32_t& posU, uint32_t& posQ, const uint16_t &X, const uint32_t &maxGaps, uint32_t& maxPosQ, uint32_t& maxPosU, struct Algn &maxAlgn, struct Algn &brdAlgn, int32_t& maxScore, int32_t& eMax, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool srchCritCheck);

//This function calculates the continuation of a gapped alignment on a successive unitig or the next peace of the query considering a quorum and a search color set
void contRightGappedAlignment(UnitigColorMap<seedlist> &uni, list<uint16_t> &extPth, const string &q, uint32_t &posQ, uint32_t &posU, const uint16_t &X, const uint32_t &maxGaps, struct Align &align, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function calculates the continuation of a gapped alignment on a predecessive unitig or the next peace of the query considering a quorum and a search color set
void contLeftGappedAlignment(UnitigColorMap<seedlist> &uni, list<uint16_t> &extPth, const string &q, uint32_t &posQ, uint32_t &posU, const uint16_t &X, const uint32_t &maxGaps, struct Algn &algn, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//Continues a gapped alignment on the same unitig as before considering a quorum and a search color set. Returns true if a new maximal score has been found.
bool contGappedOnSameUni(UnitigColorMap<seedlist> &uni, list<uint16_t> &extPth, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, struct Algn &algn, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function continues a gapped alignment calculation to the left on the same unitig considering a quorum and a search color set. Returns true if a new maximal score has been found.
bool contLeftOnSameUni(UnitigColorMap<seedlist> &uni, list<uint16_t> &extPth, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, struct Algn &algn, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//Initiates a gapped alignment calculation on all successive unitigs considering a quorum and a search color set and returns the alignment's end position in the unitig
bool contGappedOnSuccUni(UnitigColorMap<seedlist> &uni, list<uint16_t> &extPth, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, struct Algn &algn, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function initiates a gapped alignment calculation on all predecessive unitigs considering a quorum and a search color set. Returns true if calculations have been successful
bool contGappedOnPredUni(UnitigColorMap<seedlist> &uni, list<uint16_t> &extPth, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, struct Algn &algn, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function calculates a gapped alignment to the right side of the starting position considering a quorum and a search color set. ATTENTION: Hit's length attribute will be deprecated after function call!
void startRightGappedAlignment(hit *h, const string &q, const uint16_t &X, const uint32_t maxGaps, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function calculates a gapped alignment to the left side of the starting position considering a quorum and a search color set. ATTENTION: Hit's length attribute will be deprecated after function call!
void startLeftGappedAlignment(hit *h, const string &q, const uint16_t &X, const uint32_t maxGaps, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function traces back the optimal path of a right gapped extension from a given start position within an edit matrix to get the corresponding alignment
void traceBack(uint32_t &matPosI, uint32_t posQ, const string &query, uint32_t &matPosJ, uint32_t posU, const string &uniSeq, int32_t **mat, struct Algn &algn);

//This function traces back the optimal path of a left gapped extension from a given start position within an edit matrix to get the corresponding alignment
void traceForth(uint32_t &matPosI, uint32_t posQ, const string &query, uint32_t &matPosJ, uint32_t posU, const string &uniSeq, int32_t **mat, struct Algn &algn);

#endif