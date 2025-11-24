#ifndef EXTENSION_HPP
#define EXTENSION_HPP

#include "Hit.h"
#include <queue>
#include <tuple>

//Maximum recursion depth we want to allow for an extension
#define MAXRECURSIONDEPTH 1000 //Raising this constant might slightly change the quality of results in very large and complex graphs, but might affect the runtime as well
//Scores for matches and mismatches in unit score system
#define USCORE_MATCH 1
#define USCORE_MISMATCH -1
//Default value for X in X-drop algorithm
#define DEFAULT_X 3

//This function simply calculates the unity score for two bases ((mis)match = (-)1)
inline int32_t compUScore(const char &q, const char &u, const uint16_t &mscore, const int16_t &mmscore){ return (q == u ? mscore : mmscore); }

//The good old X-drop algorithm (extension to the right) for seeds matching the query's reference strand considering quorum and search color set
void startRightX_Drop(Hit* hit, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx, const int16_t extend_modus);

//The good old X-drop algorithm (extension to the right) for seeds matching the query's reverse complement considering quorum and search color set
void startRightX_Drop_OnRevComp(Hit* hit, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

//This function continues an extension to the right on a successive unitig of a seed lying on the query's reference strand considering a quorum and a search color set. Returns the maximum score reached.
int32_t contRightX_Drop(const neighborIterator<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> &sucUnitig, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t &extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastSeedTmpScore, uint32_t uniSeqPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

int32_t contRightX_Drop_BFS(const neighborIterator<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> &sucUnitig, const uint32_t &iniQoff, uint32_t &hitLen, uint32_t &extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, int32_t &lastSeedTmpScore, uint32_t &uniSeqPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx, bool &check);

int32_t contRightX_Drop_BFS_2(const neighborIterator<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> &sucUnitig, const uint32_t &iniQoff, uint32_t &hitLen, uint32_t &extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, int32_t &lastSeedTmpScore, uint32_t &uniSeqPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx, bool &check, int &numOfBases);

//This function continues an extension to the right on a successive unitig of a seed lying on the query's reverse complement considering a quorum and a search color set. Returns the maximum score reached.
int32_t contRightX_Drop_OnRevComp(const neighborIterator<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> &sucUnitig, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t &extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastSeedTmpScore, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

//This function goes through an extension list and incorporates all seeds that are not already part of a seed list into their corresponding seed lists NOTE: This function won't be used at some point anymore!
void processExtPtr(struct Ext_ptr*& extPtr);


using shorterTemp = neighborIterator<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false>;
using pathList = list<uint16_t>;
using shorterContainer = ForwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false>;
using shorterTuple = tuple<shorterTemp, uint32_t, int32_t, pathList, uint32_t, uint32_t, uint32_t>;
using shorterTuple2 = tuple<shorterTemp, uint32_t, int32_t, pathList, uint32_t, uint32_t, uint32_t, int32_t>;
using shorterVector = vector<std::tuple<shorterTemp, uint32_t, int32_t, pathList, uint32_t, uint32_t, uint32_t>>;
using shorterVector2 = vector<std::tuple<shorterTemp, uint32_t, int32_t, pathList, uint32_t, uint32_t, uint32_t, int32_t>>;
using shorterPrioQueue = priority_queue<shorterTuple, vector<shorterTuple>, const bool (*)(const shorterTuple&, const shorterTuple&)>;
using shorterPrioQueue2 = priority_queue<shorterTuple2, vector<shorterTuple2>, const bool (*)(const shorterTuple2&, const shorterTuple2&)>;
//using inputTypes = tuple<queue<shorterTuple2>,uint32_t,string,uint16_t,int16_t,int16_t,list<uint16_t>, uint32_t,uint32_t,list<pair<string, size_t>>,bool,int32_t,int32_t,uint32_t, list<uint16_t>>;
//using outputTypes = tuple<shorterPrioQueue2,int32_t, uint32_t, list<uint16_t>>;
inline const bool prioLongest(const shorterTuple& left, const shorterTuple& right){ return get<1>(left) < get<1>(right); }
inline const bool prioLongest2(const shorterTuple2& left, const shorterTuple2& right){ return get<1>(left) < get<1>(right); }


struct {
    queue<shorterTuple2> extensionQueue;
	uint32_t iniQoff;
	string q;
	uint16_t mscore;
	int16_t mmscore;
	int16_t X;
	uint32_t explCount;
	uint32_t quorum;
	list<pair<string, size_t>> searchSet;
	bool advIdx;
	int32_t maxScore;
	int32_t numOfBases;
	uint32_t hitLen;
	list<uint16_t> bestPath;
} outputTypes;

struct {
    shorterPrioQueue2 bestUnitigsPrioQueue;
    int32_t maxScore;
    uint32_t hitLen;
    list<uint16_t> bestPath;
} inputTypes;

outputTypes calcUnitigsMitBasenberechnung(inputTypes);

queue<shorterTuple2> getBestUnitigs(shorterPrioQueue2,uint);



//This function initiates the extension on all successors of a unitig and returns the best one considering a quorum and a search color set
int32_t extendAtNextUnitig(const ForwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &uniPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

int32_t extendAtNextUnitig_BFS(const ForwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &uniPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

int32_t extendAtNextUnitig_BFS_SMART1(const ForwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &uniPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

int32_t extendAtNextUnitig_BFS_SMART2(const ForwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &uniPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

int32_t extendAtNextUnitig_BFS_SMART3(const ForwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &uniPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

//This function  initiates the extension on all successors of a unitig and returns the best one considering a quorum and a search color set. This function is explicitly designed for seeds lying on the query's reverse complement (considering the overlap between unitigs in sequences' beginning)
int32_t extendAtNextUnitig_OnRevComp(const ForwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

//This function performs an extension on all possible predecessors of a unitig considering a quorum and a search color set and returns the maximum scoring one
int32_t extendAtPrevUnitig(const BackwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> bwIter, uint32_t qPos, uint32_t &hitLen, const uint32_t &tmpExtLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

//This function performs an extension on all possible predecessors of a unitig considering a quorum and a search color set and returns the maximum scoring one. This function is explicitly designed for seeds lying on the query's reverse complement (considering the overlap between unitigs in sequences' beginning)
int32_t extendAtPrevUnitigOnRevComp(const BackwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> bwIter, uint32_t qPos, uint32_t &hitLen, const uint32_t &tmpExtLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, list<uint16_t> &extPth, const uint32_t &lead, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

//This function starts the left extension for seeds lying on the query's reference strand considering a quorum and a search color set
void startLeftX_Drop(Hit* hit, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

//This function starts the left extension for seeds lying on the query's reverse complement considering a quorum and a search color set
void startLeftX_Drop_OnRevComp(Hit* hit, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

//This function continues a left extension on a predecessive unitig of a seed lying on the query's reference strand considering a quorum and a search color set and returns the achieved score
int32_t contLeftX_Drop(const neighborIterator<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> &prevUni, uint32_t qPos, uint32_t &hitLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastSeedTmpScore, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

//This function continues a left extension on a predecessive unitig of a seed lying on the query's reverse complement considering a quorum and a search color set and returns the achieved score
int32_t contLeftX_DropOnRevComp(const neighborIterator<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> &prevUni, uint32_t qPos, uint32_t &hitLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastSeedTmpScore, list<uint16_t> &extPth, uint32_t uPos, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

//This function checks if a hit's start offset for the gapped extension lies inside the overlap of a unitig sequence's end. If so it moves the start position to a unitig where it is not inside the overlap anymore.
void mvStartToValUni(Hit* h, list<uint16_t>& lExtPth);

//This function moves an offset from a given unitig to its successor while keeping left and right extension paths updated. If the offset at the successive unitig lies inside the overlap at the unitig sequence's end the function calls itself recursively
void switUni(uint32_t &offset, UnitigColorMap<UnitigInfo> &currUni, list<uint16_t> &lExtPth, list<uint16_t> &rExtPth);

#endif