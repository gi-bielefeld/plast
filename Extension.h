#ifndef EXTENSION_HPP
#define EXTENSION_HPP

//Maximum recursion depth we want to allow for an extension
#define MAXRECURSIONDEPTH 2000 //Raising this constant might slightly change the quality of results in very large and complex graphs, but might affect the runtime as well
//Scores for matches and mismatches in unit score system
#define USCORE_MATCH 1
#define USCORE_MISMATCH -1
//Default value for X in X-drop algorithm
#define DEFAULT_X 3

//A pointer to trace back an extension over multiple unitigs
struct Ext_ptr{
	//Reference to the unitig the extension is related to
	UnitigColorMap<seedlist> origUnitig;//TODO Couldn't we actually use the unitig id we use in the seed detection here as well?!
	//Seed that holds the extension seed for this unitig
	struct seed *extSeed;

	// //This flag indicates whether a seed is already part of its unitig's seed list
	// bool inSeedList;
};

//This function simply calculates the unity score for two bases ((mis)match = (-)1)
inline int32_t compUScore(const char &q, const char &u){ return (q == u ? USCORE_MATCH : USCORE_MISMATCH); }

//The good old X-drop algorithm (extension to the right) for seeds matching the query's reference strand considering quorum/*NOTE: For now the seed extension is a little bit redundant. Would be cool if this can be changed in the future! - Is the seed extension still redundant?*/
void startRightX_Drop(struct hit* hit, const string &q, const int16_t &X, const uint32_t &quorum);//TODO This function still needs to be tested!

//The good old X-drop algorithm (extension to the right) for seeds matching the query's reverse complement considering quorum/*NOTE: For now the seed extension is a little bit redundant. Would be cool if this can be changed in the future! - Is the seed extension still redundant?*/
void startRightX_Drop_OnRevComp(struct hit* hit, const string &q, const int16_t &X, const uint32_t &quorum);//TODO This function still needs to be tested!

//The good old X-drop algorithm (extension to the right) for seeds matching the query's reference strand considering quorum and search color set/*NOTE: For now the seed extension is a little bit redundant. Would be cool if this can be changed in the future! - Is the seed extension still redundant?*/
void startRightX_Drop(struct hit* hit, const string &q, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//The good old X-drop algorithm (extension to the right) for seeds matching the query's reverse complement considering quorum and search color set/*NOTE: For now the seed extension is a little bit redundant. Would be cool if this can be changed in the future! - Is the seed extension still redundant?*/
void startRightX_Drop_OnRevComp(struct hit* hit, const string &q, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function continues an extension to the right on a successive unitig of a seed lying on the query's reference strand considering a quorum. Returns the maximum score reached.
int32_t contRightX_Drop(const neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> &sucUnitig, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t &extLen, const string &q, const int16_t &X, struct Ext_ptr*& extPtr, const int32_t &lastSeedTmpScore, uint32_t uniSeqPos, uint32_t &explCount, const uint32_t &quorum);//TODO This function still needs to be tested!

//This function continues an extension to the right on a successive unitig of a seed lying on the query's reverse complement considering a quorum. Returns the maximum score reached.
int32_t contRightX_Drop_OnRevComp(const neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> &sucUnitig, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t &extLen, const string &q, const int16_t &X, struct Ext_ptr*& extPtr, const int32_t &lastUniTmpScore, uint32_t uniSeqPos, uint32_t &explCount, const uint32_t &quorum);//TODO This function still needs to be tested!

//This function continues an extension to the right on a successive unitig of a seed lying on the query's reference strand considering a quorum and a search color set. Returns the maximum score reached.
int32_t contRightX_Drop(const neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> &sucUnitig, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t &extLen, const string &q, const int16_t &X, struct Ext_ptr*& extPtr, const int32_t &lastUniTmpScore, uint32_t uniSeqPos, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function continues an extension to the right on a successive unitig of a seed lying on the query's reverse complement considering a quorum and a search color set. Returns the maximum score reached.
int32_t contRightX_Drop_OnRevComp(const neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> &sucUnitig, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t &extLen, const string &q, const int16_t &X, struct Ext_ptr*& extPtr, const int32_t &lastUniTmpScore, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function deletes an extension list
//NOTE: extList should never be NULL! (What extList?)
inline void delExtPtr(struct Ext_ptr*& extPtr){
	// //Delete the last element's seed
	// delExtPtrSeed(extPtr->inSeedList, extPtr->extSeed, extPtr->origUnitig.getData()->getData(extPtr->origUnitig));

	//Delete extension seed
	free(extPtr->extSeed);
	//Delete extension pointer
	free(extPtr);
}

//This function goes through an extension list and incorporates all seeds that are not already part of a seed list into their corresponding seed lists
void processExtPtr(struct Ext_ptr*& extPtr);

//This function initiates the extension on all successors of a unitig and returns the best one considering a quorum
struct Ext_ptr* extendAtNextUnitig(const ForwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &iniUniPos, uint32_t &explCount, const uint32_t &quorum);//TODO This function still needs to be tested!

//This function initiates the extension on all successors of a unitig and returns the best one considering a quorum. This function is explicitly designed for seeds lying on the query's reverse complement (considering the overlap between unitigs in sequences' beginning)
struct Ext_ptr* extendAtNextUnitig_OnRevComp(const ForwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &iniUniPos, uint32_t &explCount, const uint32_t &quorum);//TODO This function still needs to be tested!

//This function initiates the extension on all successors of a unitig and returns the best one considering a quorum and a search color set
struct Ext_ptr* extendAtNextUnitig(const ForwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &uniPos, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function  initiates the extension on all successors of a unitig and returns the best one considering a quorum and a search color set. This function is explicitly designed for seeds lying on the query's reverse complement (considering the overlap between unitigs in sequences' beginning)
struct Ext_ptr* extendAtNextUnitig_OnRevComp(const ForwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> sucIter, const uint32_t &iniQoff, uint32_t &hitLen, const uint32_t extLen, const string &q, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function performs an extension on all possible predecessors of a unitig consicering a quorum and return the maximum scoring one
struct Ext_ptr* extendAtPrevUnitig(const UnitigColorMap<seedlist> &lastUni, uint32_t qPos, uint32_t &hitLen, const uint32_t &tmpExtLen, const string &q, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &explCount, const uint32_t &quorum);//TODO This function still needs to be tested!

//This function performs an extension on all possible predecessors of a unitig considering a quorum and returns the maximum scoring one. This function is explicitly designed for seeds lying on the query's reverse complement (considering the overlap between unitigs in sequences' beginning)
struct Ext_ptr* extendAtPrevUnitigOnRevComp(const UnitigColorMap<seedlist> &lastUni, uint32_t qPos, uint32_t &hitLen, const uint32_t &tmpExtLen, const string &q, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &explCount, const uint32_t &quorum);//TODO This function still needs to be tested!

//This function performs an extension on all possible predecessors of a unitig considering a quorum and a search color set and returns the maximum scoring one
struct Ext_ptr* extendAtPrevUnitig(const BackwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> bwIter, uint32_t qPos, uint32_t &hitLen, const uint32_t &tmpExtLen, const string &q, const int16_t &X, const int32_t &lastExtSeedTmpScore, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function performs an extension on all possible predecessors of a unitig considering a quorum and a search color set and returns the maximum scoring one. This function is explicitly designed for seeds lying on the query's reverse complement (considering the overlap between unitigs in sequences' beginning)
struct Ext_ptr* extendAtPrevUnitigOnRevComp(const BackwardCDBG<DataAccessor<seedlist>, DataStorage<seedlist>, false> bwIter, uint32_t qPos, uint32_t &hitLen, const uint32_t &tmpExtLen, const string &q, const int16_t &X, const int32_t &lastExtSeedTmpScore, const uint32_t &lead, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function starts the left extension for seeds lying on the query's reference strand considering a quorum
void startLeftX_Drop(struct hit* hit, const string &q, const int16_t &X, const uint32_t &quorum);//TODO This function still needs to be tested!

//This function starts the left extension for seeds lying on the query's reverse complement considering a quorum
void startLeftX_Drop_OnRevComp(struct hit* hit, const string &q, const int16_t &X, const uint32_t &quorum);//TODO This function still needs to be tested!

//This function starts the left extension for seeds lying on the query's reference strand considering a quorum and a search color set
void startLeftX_Drop(struct hit* hit, const string &q, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function starts the left extension for seeds lying on the query's reverse complement considering a quorum and a search color set
void startLeftX_Drop_OnRevComp(struct hit* hit, const string &q, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function continues a left extension on a predecessive unitig of a seed lying on the query's reference strand considering a quorum and returns the achieved score
int32_t contLeftX_Drop(const neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> &prevUni, uint32_t qPos, uint32_t &hitLen, const string &q, const int16_t &X, struct Ext_ptr*& extPtr, const int32_t &lastSeedTmpScore, uint32_t &explCount, const uint32_t &quorum);//TODO This function still needs to be tested!

//This function continues a left extension on a predecessive unitig of a seed lying on the query's reverse complement considering a quorum and returns the achieved score
int32_t contLeftX_DropOnRevComp(const neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> &prevUni, uint32_t qPos, uint32_t &hitLen, const string &q, const int16_t &X, struct Ext_ptr*& extPtr, const int32_t &lastSeedTmpScore, uint32_t &explCount, const uint32_t &quorum);//TODO This function still needs to be tested!

//This function continues a left extension on a predecessive unitig of a seed lying on the query's reference strand considering a quorum and a search color set and returns the achieved score
int32_t contLeftX_Drop(const neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> &prevUni, uint32_t qPos, uint32_t &hitLen, const string &q, const int16_t &X, struct Ext_ptr*& extPtr, const int32_t &lastSeedTmpScore, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function continues a left extension on a predecessive unitig of a seed lying on the query's reverse complement considering a quorum and a search color set and returns the achieved score
int32_t contLeftX_DropOnRevComp(const neighborIterator<DataAccessor<seedlist>, DataStorage<seedlist>, false> &prevUni, uint32_t qPos, uint32_t &hitLen, const string &q, const int16_t &X, struct Ext_ptr*& extPtr, const int32_t &lastSeedTmpScore, uint32_t uPos, uint32_t &explCount, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

#endif