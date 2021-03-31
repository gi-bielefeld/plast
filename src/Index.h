#ifndef INDEX_HPP
#define INDEX_HPP

#include "Smer.h"
#include "UnitigInfo.h"

#define INDEX_VERSION "1.1"
#define ADVANCED_INDEX_FLAG_DEFAULT false

//This function writes PLAST's q-gram profile to disk
void saveProfile(ofstream &fStr, const uint32_t &profSize, uint32_t* qProfile);

//This function loads the q-gram profile from a file and returns it
uint32_t* loadQProf(ifstream &iFile, const uint32_t size);

//This function saves PLAST's index data structures as binary
void saveIndexBin(const char *filename, const int32_t& sdLen, const uint32_t& profSize, uint32_t* qProfile, const uint32_t numUni, UnitigColorMap<UnitigInfo>* uniArr, struct S_mer_pos*& pos, const size_t &seedNum, const int32_t &k, const bool& advIdx);

//This function loads the indices from a binary file
void loadIndexesBin(const char* fName, int32_t& sdLen, uint32_t &qProfSize, uint32_t *&qProf, ColoredCDBG<UnitigInfo>& graph, UnitigColorMap<UnitigInfo>*& uniArr, size_t &posSize, struct S_mer_pos*& pos, bool& advIdx);

//This function loads the indexes from a file
void loadIndexes(const string fName, const uint32_t qSize, uint32_t*& prof, CompactedDBG<UnitigInfo> &cdbg, uint32_t &seedNum, struct S_mer_pos*& lArr);

//This function builds all necessary indexes
void buildIndex(ColoredCDBG<UnitigInfo> &cdbg, const int32_t &minSeedLength, uint32_t* qProf, size_t &numSmers, const uint32_t &profSize, struct S_mer_pos*& linkArr, UnitigColorMap<UnitigInfo>*& uniArray);

#endif