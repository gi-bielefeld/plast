#ifndef INDEX_HPP
#define INDEX_HPP

#include "Smer.h"
#include "Seedlist.h"

//This function writes PLAST's q-gram profile to disk
void saveProfile(ofstream &fStr, const uint32_t &profSize, uint32_t* qProfile);

//This function loads the q-gram profile from a file and returns it
uint32_t* loadQProf(ifstream &iFile, const uint32_t size);

//This function saves PLAST's index data structures as binary
void saveIndexBin(const char *filename, const uint32_t &profSize, uint32_t* qProfile, const uint32_t numUni, UnitigColorMap<seedlist>* uniArr, struct s_mer_pos*& pos, size_t &seedNum, const int32_t &k);

//This function loads the indices from a binary file
void loadIndexesBin(const char* fName, const uint32_t &qProfSize, uint32_t *qProf, ColoredCDBG<seedlist>& graph, UnitigColorMap<seedlist>*& uniArr, size_t &posSize, struct s_mer_pos*& pos);

//This function loads the indexes from a file
void loadIndexes(const string fName, const uint32_t qSize, uint32_t*& prof, CompactedDBG<seedlist> &cdbg, uint32_t &seedNum, struct s_mer_pos*& lArr);

//This function builds all necessary indexes
void buildIndex(ColoredCDBG<seedlist> &cdbg, const int32_t &minSeedLength, uint32_t* qProf, size_t &numSmers, const uint32_t &profSize, struct s_mer_pos*& linkArr, UnitigColorMap<seedlist>*& uniArray);

#endif