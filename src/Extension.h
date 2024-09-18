#ifndef EXTENSION_HPP
#define EXTENSION_HPP

#include "Hit.h"

//Maximum number of unitig switches during a single extension. Raising this constant might slightly change the quality of results in
// very large and complex graphs, but might affect the runtime as well
#define EXPLORED_UNITIGS_MAX 1000
//Scores for matches and mismatches in unit score system
#define USCORE_MATCH 1
#define USCORE_MISMATCH -1
//Default value for X in X-drop algorithm
#define DEFAULT_X 3

//A hit's extension in the graph
class Ext {

	public:
		
		Ext(): offsetU(0), offsetQ(0), tmpQoff(0), score(0), tmpScore(0) {}//TODO: This function still needs to be tested!

		Ext(const int32_t& offU, const int32_t& offQ, const int32_t& tQoff, const int32_t& scr, const int32_t& tScr, const 
		UnitigColorMap<UnitigInfo>& lUni, const list<uint16_t>& pth){//TODO: This function still needs to be tested!
			offsetU = offU;
			offsetQ = offQ;
			tmpQoff = tQoff;
			score = scr;
			tmpScore = tScr;
			ldUni = lUni;
			this->pth = pth;
		}

		//Copy constructor
		Ext(const Ext& e){//TODO: This function still needs to be tested!
			offsetU = e.offsetU;
			offsetQ = e.offsetQ;
			tmpQoff = e.tmpQoff;
			score = e.score;
			tmpScore = e.tmpScore;
			ldUni = e.ldUni;
			pth = e.pth;
		}

		//Offset in leading unitig
		int32_t offsetU;
		//Offset in query
		int32_t offsetQ;
		//Offset in query in the last not yet accepted exploration
		int32_t tmpQoff;
		//The extension's score
		int32_t score;
		//The score of the last not yet accepted exploration
		int32_t tmpScore;
		//The leading unitig
		UnitigColorMap<UnitigInfo> ldUni;
		//The extension's path
		list<uint16_t> pth;
};

//This function encapsulates part of the X-drop algorithm used during an extension. It updates the given extension and exploration 
//length depending on whether the temporary score became positive by adding the given score
inline void updateExtension(Ext& ext, int32_t& posQ , const int32_t& expScr, int32_t& expLen, const bool& isRightExt){
	if((ext.tmpScore += expScr) > 0){
		//Update extension information
		ext.score += ext.tmpScore;

		//If we extend to the right we have to add and otherwise substract the extension length
		if(isRightExt){
			ext.offsetU += expLen;
			ext.offsetQ = posQ + expLen;
		} else{
			ext.offsetU -= expLen;
			ext.offsetQ = posQ - expLen;
		}
		
		//Update posQ
		posQ = ext.offsetQ;
		//Reset exploration length, temporary offset in query and temporary score
		expLen = 0;
		ext.tmpQoff = 0;
		ext.tmpScore = 0;
	}
}

//This function simply calculates the unity score for two bases ((mis)match = (-)1)
inline int32_t compUScore(const char &q, const char &u, const uint16_t &mscore, const int16_t &mmscore){ return (q == u ? mscore : mmscore); }

//The good old X-drop algorithm (extension to the right) for seeds matching the query's reference strand considering quorum and 
//search color set using an iterative approach
void perfRightX_Drop(Hit* hit, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const uint32_t 
	&quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

//The good old X-drop algorithm (extension to the right) for seeds matching the query's reference strand considering quorum and 
//search color set using recursive function calls
void startRightX_Drop(Hit* hit, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

//The good old X-drop algorithm (extension to the right) for seeds matching the query's reverse complement considering quorum and search color set
void startRightX_Drop_OnRevComp(Hit* hit, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

//This function continues an extension to the right on a successive unitig of a seed lying on the query's reference strand considering a quorum and a search color set. Returns the maximum score reached.
int32_t contRightX_Drop(const neighborIterator<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> &sucUnitig, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t &extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastSeedTmpScore, uint32_t uniSeqPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

//This function continues an extension to the right on a successive unitig of a seed lying on the query's reverse complement considering a quorum and a search color set. Returns the maximum score reached.
int32_t contRightX_Drop_OnRevComp(const neighborIterator<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> &sucUnitig, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t &extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastSeedTmpScore, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

//This function goes through an extension list and incorporates all seeds that are not already part of a seed list into their corresponding seed lists NOTE: This function won't be used at some point anymore!
void processExtPtr(struct Ext_ptr*& extPtr);

//This function initiates the extension on all successors of a unitig and returns the best one considering a quorum and a search color set
int32_t extendAtNextUnitig(const ForwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &uniPos, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

//This function  initiates the extension on all successors of a unitig and returns the best one considering a quorum and a search color set. This function is explicitly designed for seeds lying on the query's reverse complement (considering the overlap between unitigs in sequences' beginning)
int32_t extendAtNextUnitig_OnRevComp(const ForwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

//This function performs an extension on all possible predecessors of a unitig considering a quorum and a search color set and returns the maximum scoring one
int32_t extendAtPrevUnitig(const BackwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> bwIter, uint32_t qPos, uint32_t &hitLen, const uint32_t &tmpExtLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, list<uint16_t> &extPth, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

//This function performs an extension on all possible predecessors of a unitig considering a quorum and a search color set and returns the maximum scoring one. This function is explicitly designed for seeds lying on the query's reverse complement (considering the overlap between unitigs in sequences' beginning)
int32_t extendAtPrevUnitigOnRevComp(const BackwardCDBG<DataAccessor<UnitigInfo>, DataStorage<UnitigInfo>, false> bwIter, uint32_t qPos, uint32_t &hitLen, const uint32_t &tmpExtLen, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const int32_t &lastExtSeedTmpScore, list<uint16_t> &extPth, const uint32_t &lead, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

//This function starts the left extension for seeds lying on the query's reference strand considering a quorum and a search color 
//set using an iterative approach
void perfLeftX_Drop(Hit* hit, const string &q, const uint16_t &mscore, const int16_t &mmscore, const int16_t &X, const uint32_t 
	&quorum, const list<pair<string, size_t>> &searchSet, const bool& advIdx);

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