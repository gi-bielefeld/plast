#ifndef SEARCH_HPP
#define SEARCH_HPP

#include "GappedAlignment.h"

//Default quorum used for searches
#define DEFAULT_QUORUM 1
//Default number of results to be outputted
#define DEFAULT_NB_RES 250

//This function checks if a quorum is fulfilled for a given offset position in a unitig. Returns the number of checked positions if quorum is fulfilled and 0 otherwise.
inline uint32_t checkSearchCrit(UnitigColorMap<seedlist> unitig, const uint32_t &quorum, const size_t &offset){//TODO: This function could also be modified for a recursive use to verfy the quorum requirement for as much sequence of the unitig as possible. Check whether this might save some time!
	size_t ccount = 0, curID = 0;

	//Testing
	// if(report){
	// 	cerr << "checkSearchCrit" << endl;
	// }

	//If quorum is 1 we do not need to do anything
	if(quorum == 1){
		return unitig.size;
	}
	
	//Calculate used k (we assume here that the uniti's current offset is 0)
	const size_t k = unitig.getGraph()->getK();
	//Set unitig attributes
	unitig.dist = min(offset, unitig.size - k);//make sure that offset exists
	unitig.len = 1;

	//Check all colors that this unitig might have
	while(curID <= unitig.getData()->getUnitigColors(unitig)->colorMax(unitig)){
		//Check whether color is contained
		if(unitig.getData()->getUnitigColors(unitig)->contains(unitig, curID)){
			//Increment counter
			++ccount;

			//Check whether quorum is reached
			if(ccount == quorum) return k;
		}

		//Increment color id
		++curID;
	}

	return 0;
}

//This function checks if a quorum is fulfilled for a given color search set and an offset position in a unitig. Returns k if quorum is fulfilled and 0 otherwise.
inline uint32_t checkSearchCrit(UnitigColorMap<seedlist> unitig, const uint32_t &quorum, const size_t &offset, const list<pair<string, size_t>> &srchColSet){//TODO: This function could also be modified for a recursive use to verfy the quorum requirement for as much sequence of the unitig as possible. Check whether this might save some time!
	size_t ccount = 0, curID;
	//Calculate used k (we assume here that the uniti's current offset is 0)
	const size_t k = unitig.getGraph()->getK();

	//Testing
	//cout << "k: " << k << endl;
	// cout << "Inside checkSearchCrit: Unitig " << unitig.mappedSequenceToString() << " has colors:" << endl;
	// for(UnitigColors::const_iterator j = unitig.getData()->getUnitigColors(unitig)->begin(unitig); j != unitig.getData()->getUnitigColors(unitig)->end(); ++j){
	// 		cout << j.getColorID() << " " << j.getKmerPosition() << endl;
	// }

	//Set unitig attributes
	unitig.dist = min(offset, unitig.size - k);//make sure that offset exists
	unitig.len = 1;

	//Testing
	// cout << "Inside checkSearchCrit: Unitig " << unitig.mappedSequenceToString() << " has colors:" << endl;
	// for(UnitigColors::const_iterator j = unitig.getData()->getUnitigColors(unitig)->begin(unitig); j != unitig.getData()->getUnitigColors(unitig)->end(); ++j){
	// 		cout << j.getColorID() << " " << j.getKmerPosition() << endl;
	// }

	//Check whether a search set is given
	if(srchColSet.empty()){
		//If quorum is 1 we do not need to do anything
		if(quorum == 1){
			return unitig.size;
		}

		curID = 0;

		//Check all colors that this unitig might have
		while(curID <= unitig.getData()->getUnitigColors(unitig)->colorMax(unitig)){
			//Check whether color is contained
			if(unitig.getData()->getUnitigColors(unitig)->contains(unitig, curID)){
				//Increment counter
				++ccount;

				//Check whether quorum is reached
				if(ccount == quorum) return k;
			}

			//Increment color id
			++curID;
		}
	} else{
		//Go through search color set
		for(list<pair<string, size_t>>::const_iterator col = srchColSet.begin(); col != srchColSet.end(); ++col){
			//Check whether the current color is present at our current position
			if(unitig.getData()->getUnitigColors(unitig)->contains(unitig, col->second)){
				//Check whether our quorum is already reached
				if(quorum <= ++ccount) return k;
			}
		}
	}

	return 0;
}

inline bool quorumFulfilled(UnitigColorMap<seedlist> uni, const uint32_t posU, size_t uniLen, const uint32_t &quorum){//TODO This function still needs to be tested!
	size_t ccount = 0, curID = 0;

	//If quorum is 1 there is nothing to do
	if(quorum == 1) return true;

	//Check to which side we want to calculate an alignment
	if(posU > 0){//This can also be false if we want to extend to the right but in that case uniLen is correctly set
		//Adjust  unitig start position
		uni.dist = min((size_t) posU, uni.size - uni.getGraph()->getK());
	} else{
		//Adjust unitig length
		uni.len = uniLen;
	}

	//Check all colors that this unitig might have
	while(curID <= uni.getData()->getUnitigColors(uni)->colorMax(uni)){
		//Check whether color is contained
		if(uni.getData()->getUnitigColors(uni)->contains(uni, curID)){
			//Increment counter
			++ccount;

			//Check whether quorum is reached
			if(ccount == quorum) return true;
		}

		//Increment color id
		++curID;
	}

	return false;
}

//ATTENTION: This function only works correctly if a non-empty search color set is given (Why should that be?)
inline bool quorumFulfilled(UnitigColorMap<seedlist> uni, const uint32_t posU, size_t uniLen, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	size_t ccount = 0, curID = 0;

	//Check to which side we want to calculate an alignment
	if(posU > 0){//This can also be false if we want to extend to the right but in that case uniLen is correctly set
		//Testing
		//cout << "1 Option 2" << endl;

		//If we want to extend to the right adjust unitig start position
		uni.dist = min((size_t) posU, uni.size - uni.getGraph()->getK());
	} else{
		//Testing
		//cout << "1 Option 1" << endl;
		// cout << "uni.len:" << uni.len << endl;
		// cout << "uniLen: " << uniLen << endl;
		
		//If we want to extented to the left adjust unitig length
		uni.len = min(uniLen, uni.len);

		//Testing
		// cout << "uni.len:" << uni.len << endl;
	}

	//Check if we have a search set to consider
	if(searchSet.empty()){
		//If quorum is 1 there is nothing to do
		if(quorum == 1){
			//Testing
			//cout << "2 Option 1" << endl;

			return true;
		}

		//Testing
		//cout << "2 Option 2" << endl;
		// cout << "(Inside quorumFulfilled): Quorum is not 1" << endl << "We control unitig " << uni.mappedSequenceToString() << endl;
		// for(UnitigColors::const_iterator i = uni.getData()->getUnitigColors(uni)->begin(uni); i != uni.getData()->getUnitigColors(uni)->end(); ++i){
		// 	cout << "Position:" << (*i).first << " Color ID:" << (*i).second << endl;
		// }

		//Check all colors that this unitig might have
		while(curID <= uni.getData()->getUnitigColors(uni)->colorMax(uni)){
			//Check whether color is contained
			if(uni.getData()->getUnitigColors(uni)->contains(uni, curID)){
				//Increment counter
				++ccount;

				//Check whether quorum is reached
				if(ccount == quorum){
					//Testing
					//cout << "3 Option 1" << endl;
					// cout << "We can reach the demanded count" << endl;

					return true;
				}
			}

			//Increment color id
			++curID;
		}

		//Testing
		//cout << "3 Option 2" << endl;
	} else{
		//Testing
		//cout << "2 Option " << (quorum == 1 ? "3" : "4") << endl;

		//Go through search color set
		for(list<pair<string, size_t>>::const_iterator col = searchSet.begin(); col != searchSet.end(); ++col){
			//Testing
			//cout << "We get here" << endl;

			//Check whether the current color is present
			if(uni.getData()->getUnitigColors(uni)->contains(uni, col->second)){
				//Testing
				//cout << "We get here too" << endl;
				
				//Check whether our quorum is already reached
				if(quorum <= ++ccount){
					//Testing
					//cout << "3 Option 3" << endl;

					return true;
				}
			}
		}

		//Testing
		//cout << "3 Option 4" << endl;
	}

	return false;
}

//This function calculates the correct offset position depending on which strand we are
inline uint32_t compOffset(const uint32_t &offset, const int32_t &seedLen, const uint32_t &seqLen, const bool &onRefStrand){
	//If we are on the reference strand there is nothing to do
	if(onRefStrand) return offset;

	//Calculate coordinates on reference strand and exit
	return seqLen - offset - seedLen;
}

//This function performs the seed detection if no search color set is given
void detectSeeds(const int32_t &k, const int32_t &minSeedLength, const size_t &numSmers, const uint32_t &quorum, const uint32_t &profileSize, const uint32_t *qProfile, const string &q, const UnitigColorMap<seedlist> *uArr, const struct s_mer_pos *posArray, struct hit *hitArr, const bool isRefSeq);//TODO Tests for this function need to be adjusted!

//This function performs the seed detection if a search color set is given
void detectSeeds(const int32_t &k, const int32_t &minSeedLength, const size_t &numSmers, const uint32_t &quorum, const uint32_t &profileSize, const uint32_t *qProfile, const string &q, const UnitigColorMap<seedlist> *uArr, const struct s_mer_pos *posArray, struct hit *hitArr, const list<pair<string, size_t>> &searchSet, const bool isRefSeq);

//This function extends all seeds found on the queries reference strand considering a quorum
void extendRefSeeds(ColoredCDBG<seedlist> &cdbg, const string &q, const int16_t &X, struct hit *hitArr, const uint32_t &quorum);//TODO This function still needs to be tested!

//This function extends all seeds found on the queries reverse complement considering a quorum
void extendRevCompSeeds(ColoredCDBG<seedlist> &cdbg, const string &q, const int16_t &X, struct hit *hitArr, const uint32_t &quorum);//TODO This function still needs to be tested!

//This function extends all seeds found on the queries reference strand considering a quorum and a search color set
void extendRefSeeds(const ColoredCDBG<seedlist> &cdbg, const string &q, const int16_t &X, struct hit *hitArr, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function extends all seeds found on the queries reverse complement considering a quorum and a search color set
void extendRevCompSeeds(const ColoredCDBG<seedlist> &cdbg, const string &q, const int16_t &X, struct hit *hitArr, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function calculates a banded, semi-global, gapped alignment on a list of results considering a quorum and outputs the result if demanded
void calcGappedAlignment(const list<struct hit*> &resList, const string &q, const int16_t &X, const bool &noOutput, const uint32_t &quorum);//TODO This function still needs to be tested!

//This function calculates a banded, semi-global, gapped alignment on a list of results considering a quorum and a search color set and outputs the result if demanded
void calcGappedAlignment(const list<struct hit*> &resList, const string &q, const int16_t &X, const bool &noOutput, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function performs the actual graph search for a query
void searchQuery(ColoredCDBG<seedlist> &cdbg, const int32_t &kMerLength, const int32_t &minSeedLength, const size_t &numSmers, const uint32_t &quorum, const uint32_t &profileSize, const uint32_t *qProfile, const string &q, const SrchStrd &strand, const UnitigColorMap<seedlist> *uArr, const struct s_mer_pos *posArray, const list<pair<string, size_t>> &searchColors, const int16_t &X, const bool &calcRT, uint16_t nRes);

#endif