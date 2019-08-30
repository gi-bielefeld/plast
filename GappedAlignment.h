#ifndef GAPPEDALIGN_HPP
#define GAPPEDALIGN_HPP

#define MAX_MATRIX_SIZE 1000//Testing 10
#define GAP_SCORE -1
#define GAP_RATIO 20

//This function calculates a semi-global alignment of a unitig and the query sequence considering a quorum. Returns true if calculations on unitig could be finished (i.e. the unitig sequence was not to long to be stored inside the edit matrix).
bool calcSemiGlobAlignment(const UnitigMap<seedlist> &uni, const string &q, uint32_t& posU, uint32_t& posQ, const uint16_t &X, const uint32_t &maxGaps, uint32_t& maxPosQ, uint32_t& maxPosU, int32_t& maxScore, int32_t& eMax, const uint32_t &quorum);//TODO Tests for this funtion should be extended for quorum!

//This function calculates a semi-global alignment of a unitig and the query sequence considering a quorum and a search color set. Returns true if calculations on unitig could be finished (i.e. the unitig sequence was not to long to be stored inside the edit matrix).
bool calcSemiGlobAlignment(const UnitigMap<seedlist> &uni, const string &q, uint32_t& posU, uint32_t& posQ, const uint16_t &X, const uint32_t &maxGaps, uint32_t& maxPosQ, uint32_t& maxPosU, int32_t& maxScore, int32_t& eMax, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function calculates a banded alignment of a unitig and the query sequence to the left considering a quorum. Returns true if calculations on unitig could be finished (i.e. the unitig sequence was not to long to be stored inside the edit matrix).
bool calcLeftGlobAlignment(const UnitigColorMap<seedlist> &uni, const string &q, uint32_t& posU, uint32_t& posQ, const uint16_t &X, const uint32_t &maxGaps, uint32_t& maxPosQ, uint32_t& maxPosU, int32_t& maxScore, int32_t& eMax, const uint32_t &quorum);//TODO: So far we do not adjust the maximal number of gaps. This has to be changed as soon as we calculate real alignments!//TODO: This function does more than it should!//TODO Tests have to be extended by quorum!

//This function calculates a banded alignment of a unitig and the query sequence to the left considering a quorum and a search color set. Returns true if calculations on unitig could be finished (i.e. the unitig sequence was not to long to be stored inside the edit matrix).
bool calcLeftGlobAlignment(const UnitigColorMap<seedlist> &uni, const string &q, uint32_t& posU, uint32_t& posQ, const uint16_t &X, const uint32_t &maxGaps, uint32_t& maxPosQ, uint32_t& maxPosU, int32_t& maxScore, int32_t& eMax, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);//TODO: So far we do not adjust the maximal number of gaps. This has to be changed as soon as we calculate real alignments!//TODO: This function does more than it should!

//This function calculates the continuation of a gapped alignment on a successive unitig or the next peace of the query considering a quorum
void contRightGappedAlignment(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &posQ, uint32_t &posU, const uint16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum);//TODO Tests for this funtion should be extended for quorum!

//This function calculates the continuation of a gapped alignment on a successive unitig or the next peace of the query considering a quorum and a search color set
void contRightGappedAlignment(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &posQ, uint32_t &posU, const uint16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function calculates the continuation of a gapped alignment on a predecessive unitig or the next peace of the query considering a quorum
void contLeftGappedAlignment(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &posQ, uint32_t &posU, const uint16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum);//TODO Tests for this function have to be extended by quorum!

//This function calculates the continuation of a gapped alignment on a predecessive unitig or the next peace of the query considering a quorum and a search color set
void contLeftGappedAlignment(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &posQ, uint32_t &posU, const uint16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//Continues a gapped alignment on the same unitig as before considering a quorum. Returns true if a new maximal score has been found.
bool contGappedOnSameUni(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum);//TODO Quorum still has to be included into the tests!

//Continues a gapped alignment on the same unitig as before considering a quorum and a search color set. Returns true if a new maximal score has been found.
bool contGappedOnSameUni(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function continues a gapped alignment calculation to the left on the same unitig considering a quorum. Returns true if a new maximal score has been found.
bool contLeftOnSameUni(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum);//TODO Tests for this function have to be extended by quorum!

//This function continues a gapped alignment calculation to the left on the same unitig considering a quorum and a search color set. Returns true if a new maximal score has been found.
bool contLeftOnSameUni(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//Initiates a gapped alignment calculation on all successive unitigs considering a quorum and returns the alignment's end position in the unitig
bool contGappedOnSuccUni(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum);//TODO Tests for this function have to be extended by quorum!

//Initiates a gapped alignment calculation on all successive unitigs considering a quorum and a search color set and returns the alignment's end position in the unitig
bool contGappedOnSuccUni(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function initiates a gapped alignment calculation on all predecessive unitigs considering a quorum. Returns true if calculations have been successful
bool contGappedOnPredUni(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum);//TODO Tests for this function have to be extended by quorum!

//This function initiates a gapped alignment calculation on all predecessive unitigs considering a quorum and a search color set. Returns true if calculations have been successful
bool contGappedOnPredUni(UnitigColorMap<seedlist> &uni, const string &q, uint32_t &qOff, uint32_t &uOff, const int16_t &X, const uint32_t &maxGaps, int32_t &score, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function calculates a gapped alignment to the right side of the starting position considering a quorum. ATTENTION: Hit's length attribute will be deprecated after function call!
void startRightGappedAlignment(struct hit *h, const string &q, const uint16_t &X, const uint32_t maxGaps, const uint32_t &quorum);//TODO Tests for this function have to be extended by quorum!

//This function calculates a gapped alignment to the right side of the starting position considering a quorum and a search color set. ATTENTION: Hit's length attribute will be deprecated after function call!
void startRightGappedAlignment(struct hit *h, const string &q, const uint16_t &X, const uint32_t maxGaps, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function calculates a gapped alignment to the left side of the starting position considering a quorum. ATTENTION: Hit's length attribute will be deprecated after function call!
void startLeftGappedAlignment(struct hit *h, const string &q, const uint16_t &X, const uint32_t maxGaps, const uint32_t &quorum);//TODO Tests for this functions have to be extended by quorum!

//This function calculates a gapped alignment to the left side of the starting position considering a quorum and a search color set. ATTENTION: Hit's length attribute will be deprecated after function call!
void startLeftGappedAlignment(struct hit *h, const string &q, const uint16_t &X, const uint32_t maxGaps, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function calculates a global, banded alignment of the reverse of a unitig and the query sequence. Returns true if calculations on this unitig could be finished (i.e. the unitig sequence was not to long to be stored inside the edit matrix).
//bool calcLeftGlobAlignment(const UnitigMap<seedlist> &uni, const string &q, uint32_t& startPosU, uint32_t& startPosQ, const uint16_t &X, const uint32_t &maxGaps, uint32_t& maxPosQ, uint32_t& maxPosU, int32_t& maxScore, int32_t& eMax);

//This function calculates a global, banded alignment of the reverse of two sequences. Calculations end as soon as the end of one of the sequences has been reached or the score drops below a specified value X. Returns true if a calculation was successful.
//bool calcLeftGlobAlignment(const string hseq, const string vseq, const uint16_t& X, const uint32_t& bandRadius, uint32_t& maxPos_hseq, uint32_t& maxPos_vseq, int32_t& maxScore, uint32_t& ePos_hseq, uint32_t& ePos_vseq, int32_t& eMaxScore);

#endif