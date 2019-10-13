#ifndef SEARCH_HPP
#define SEARCH_HPP

#include "GappedAlignment.h"

//Default quorum used for searches
#define DEFAULT_QUORUM 1
//Default number of results to be outputted
#define DEFAULT_NB_RES 250

//All variations of strands we may perform a search
enum SrchStrd {Plus, Minus, Both};

//This function checks if a quorum is fulfilled for a given offset position in a unitig. Returns the number of checked positions if quorum is fulfilled and 0 otherwise.
inline uint32_t checkSearchCrit(UnitigColorMap<seedlist> unitig, const uint32_t &quorum, const size_t &offset){//TODO: This function could also be modified for a recursive use to verfy the quorum requirement for as much sequence of the unitig as possible. Check whether this might save some time!
	size_t ccount = 0, curID = 0;

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

//This function checks if a unitig fulfills a quorum completely (ATTENTION: This function does not work for a quorum of 1)
inline bool isCovered(const UnitigColorMap<seedlist> &uni, const uint32_t &quorum){//TODO The search criteria tests are now designed in a way that seeds won't be found inside a unitig if search criteria are not fulfilled either from the unitig's beginning or from its end up to the seed's position meaing that we might miss some seeds for which the search criteria are actually fulfilled. Mybe it is a good idea to introduce a parameter which allows to find them as well!
	uint16_t counter = 0;
	size_t curID = SIZE_MAX;

	//Iterate over unitig's colors
	for(UnitigColors::const_iterator i = uni.getData()->getUnitigColors(uni)->begin(uni); curID != i.getColorID(); i.nextColor()){
		//Update current id
		curID = i.getColorID();

		//Testing
		// cout << "Iterate over colors" << endl;
		// cout << "colorMax: " <<  << endl;

		//Check if the color is present on the complete unitig
		if(uni.getData()->getUnitigColors(uni)->contains(uni, curID)){
			//Testing
			// cout << "Color is contained" << endl << "curID: " << curID << endl;

			//Increment counter and check if we are done
			if(++counter == quorum){
				//Testing
				// cout << "Quorum reached" << endl;

				return true;
			}
		}
	}

	return false;
}

//This function checks if a unitig fulfills the search criteria completely
inline bool isCovered(const UnitigColorMap<seedlist> &uni, const uint32_t &quorum, const list<pair<string, size_t>> &srchColSet){
	uint16_t counter = 0;

	//Go through search color set
	for(list<pair<string, size_t>>::const_iterator col = srchColSet.begin(); col != srchColSet.end(); ++col){
		//Testing
		// cout << "Color " << col->second << " is " << (uni.getData()->getUnitigColors(uni)->contains(uni, col->second) ? "" : "not ") << "contained" << endl;
		
		//Check whether the current color is present at our current position
		if(uni.getData()->getUnitigColors(uni)->contains(uni, col->second)){
			//Check whether our quorum is already reached
			if(quorum == ++counter){
				return true;
			}
		}
	}

	return false;
}

//This functions checks up to which offsets search criteria are fulfilled for a unitig and saves the result in the unitig's seedlist. If even the first k-mer from either side is not covered by the search criteria positions are set to -1 and unitig length respectively
inline void calcSrchCritBrds(UnitigColorMap<seedlist> uni, const uint32_t &quorum, const list<pair<string, size_t>> &srchColSet){
	bool lFix, rFix;
	int32_t lBrd, rBrd;

	//Check whether a search set is given
	if(srchColSet.empty()){
		//Testing
		// cout << "1 Option 2" << endl;

		//If quorum is 1 we do not need to do anything
		if(quorum == 1 || isCovered(uni, quorum)){
			//Testing
			// cout << "2 Option 1" << endl;

			//Save offsets
			uni.getData()->getData(uni)->setlBrd(uni.size);
			uni.getData()->getData(uni)->setrBrd(0);
			return;
		}

		//Testing
		// cout << "2 Option 2" << endl;

		//Check if an iterative check makes no sense
		if(uni.len == 1){
			//Testing
			// cout << "3 Option 1" << endl;

			//Save offsets
			uni.getData()->getData(uni)->setlBrd(-1);
			uni.getData()->getData(uni)->setrBrd(uni.size);
			return;
		}

		//Testing
		// cout << "3 Option 2" << endl;

		//Initialize borders for iterative check
		lBrd = 0;
		rBrd = uni.size - uni.getGraph()->getK();
		//Make a copy of our unitig so that we can check both sides at the same time
		UnitigColorMap<seedlist> uniCpy = uni;

		// //Save original unitig length
		// uniSize = uni.size;

		//Set unitig offset
		uniCpy.dist = rBrd;
		//Set unitig's lengths to 1 k-mer
		uni.len = 1;
		uniCpy.len = 1;
		lFix = false;
		rFix = false;

		//Check every position iteratively until either both borders cannot be moved anymore or have reached each other
		while((!lFix || !rFix) && lBrd <= rBrd){
			//Check if we can move the left border
			if(!lFix && isCovered(uni, quorum)){
				//Testing
				// cout << "4 Option 1" << endl;

				//Move it
				++lBrd;
			} else{
				//Testing
				// cout << "4 Option 2" << endl;

				//Report border as fixed
				lFix = true;
			}

			//Testing
			// cout << "isCovered(uniCpy, quorum): " << isCovered(uniCpy, quorum) << endl;
			// cout << "uniCpy: " << uniCpy.mappedSequenceToString() << " quorum: " << quorum << endl;
			// cout << "uniCoy colors:" << endl;
			// for(UnitigColors::const_iterator j = uniCpy.getData()->getUnitigColors(uniCpy)->begin(uniCpy); j != uniCpy.getData()->getUnitigColors(uniCpy)->end(); j.nextColor()) cout << "Pos:" << j.getKmerPosition() << " ID:" << j.getColorID() << endl;

			//Check if we can move the right border
			if(!rFix && isCovered(uniCpy, quorum)){
				//Testing
				// cout << "5 Option 1" << endl;

				//Move it
				--rBrd;
			} else{
				//Testing
				// cout << "5 Option 2" << endl;

				//Report border as fixed
				rFix = true;
			}

			//Testing
			// if(lBrd > rBrd){
			// 	cout << "6 Option 1" << endl;
			// } else{
			// 	// cout << "6 Option 2" << endl;
			// }
			// cout << "lBrd: " << lBrd << endl;

			//Adjust unitigs
			uni.dist = lBrd;
			uniCpy.dist = rBrd;
		}
	} else{
		//Testing
		// cout << "1 Option 1" << endl;

		//Check the complete unitig
		if(isCovered(uni, quorum, srchColSet)){
			//Testing
			// cout << "7 Option 1" << endl;

			//Save offsets
			uni.getData()->getData(uni)->setlBrd(uni.size);
			uni.getData()->getData(uni)->setrBrd(0);
			return;
		}

		//Testing
		// cout << "7 Option 2" << endl;

		//Check if an iterative check makes no sense
		if(uni.len == 1){
			//Testing
			// cout << "8 Option 1" << endl;

			//Save offsets
			uni.getData()->getData(uni)->setlBrd(-1);
			uni.getData()->getData(uni)->setrBrd(uni.size);
			return;
		}

		//Testing
		// cout << "8 Option 2" << endl;

		//Initialize borders for iterative check
		lBrd = 0;
		rBrd = uni.size - uni.getGraph()->getK();
		//Make a copy of our unitig so that we can check both sides at the same time
		UnitigColorMap<seedlist> uniCpy = uni;

		// //Save original unitig length
		// uniSize = uni.size;

		//Set unitig offset
		uniCpy.dist = rBrd;
		//Set unitig's lengths to 1 k-mer
		uni.len = 1;
		uniCpy.len = 1;
		lFix = false;
		rFix = false;

		//Check every position iteratively until either both borders cannot be moved anymore or have reached each other
		while((!lFix || !rFix) && lBrd <= rBrd){
			//Check if we can move the left border
			if(!lFix && isCovered(uni, quorum, srchColSet)){
				//Testing
				// cout << "9 Option 1" << endl;

				//Move it
				++lBrd;
			} else{
				//Testing
				// cout << "9 Option 2" << endl;

				//Report border as fixed
				lFix = true;
			}

			if(!rFix && isCovered(uniCpy, quorum, srchColSet)){
				//Testing
				// cout << "10 Option 1" << endl;

				//Move it
				--rBrd;
			} else{
				//Testing
				// cout << "10 Option 2" << endl;

				//Report border as fixed
				rFix = true;
			}

			//Adjust unitigs
			uni.dist = lBrd;
			uniCpy.dist = rBrd;
		}

		//Testing
		// if(lBrd > rBrd){
		// 	cout << "11 Option 1" << endl;
		// } else{
		// 	// cout << "11 Option 2" << endl;
		// }
	}

	//Testing
	// cout << "rBrd: " << rBrd << "(int32_t) uni.size - uni.getGraph()->getK(): " << (int32_t) uni.size - uni.getGraph()->getK() << endl;

	//Check if left border could be moved at all
	if(lBrd != 0){
		//If the first position is fulfilled the next k-1 are too
		uni.getData()->getData(uni)->setlBrd(lBrd + uni.getGraph()->getK() - 1);
	} else{
		//Save offset
		uni.getData()->getData(uni)->setlBrd(-1);
	}

	//Check if right border could be moved at all
	if(rBrd < (int32_t) uni.size - uni.getGraph()->getK()){
		//Save offset
		uni.getData()->getData(uni)->setrBrd(rBrd);
	} else{
		//This is needed to avoid that left and right border might overlap accidently which would cause false positives
		uni.getData()->getData(uni)->setrBrd(uni.size);
	}
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

//This function looks up the number of positions which fulfill the search criteria depending at which offset we want to start and into which direction we extend with regard to the reference strand
inline int32_t getSrchCritCov(UnitigColorMap<seedlist> uni, const uint32_t &quorum, const list<pair<string, size_t>> &srchCols, const uint32_t &offset, const bool &toRightOnRef){
	int32_t brd;

	//Testing
	// cout << "Search criteria already checked? " << uni.getData()->getData(uni)->srchCritChckd() << endl;
	// cout << "Inside getSrchCritCov" << endl;

	//Check if search criteria have already been checked for this unitig and check if necessary
	if(!uni.getData()->getData(uni)->srchCritChckd()){
		//Testing
		// cout << "1 Option 2" << endl;
		// cout << "Calculate search crit borders" << endl;

		calcSrchCritBrds(uni, quorum, srchCols);
	}// else{
		//Testing
		// cout << "1 Option 1" << endl;
	//}

	//Testing
	// cout << "Done" << endl;

	//Check if search criteria are not fulfilled at all on this unitig
	if(uni.getData()->getData(uni)->getlBrd() < 0 && uni.getData()->getData(uni)->getrBrd() == (int32_t) uni.size){
		//Testing
		// cout << "2" << endl;

		return 0;
	}

	//Testing
	// cout << "Unitig should be covered completely" << endl;

	//Check if unitig is completely covered
	if(uni.getData()->getData(uni)->getlBrd() > uni.getData()->getData(uni)->getrBrd()){
		//Testing
		// cout << "3 Option 1" << endl;
		// cout << "uni.getData()->getData(uni)->getlBrd(): " << uni.getData()->getData(uni)->getlBrd() << " uni.getData()->getData(uni)->getrBrd(): " << uni.getData()->getData(uni)->getrBrd() << endl;
		// cout << "It is" << endl;

		return uni.size;
	}

	//Testing
	// cout << "3 Option 2" << endl;

	//Check into which direction we want to go
	if(toRightOnRef){
		//Testing
		// cout << "4 Option 1" << endl;
		// cout << "lBrd: " << uni.getData()->getData(uni)->getlBrd() << endl;

		//Get right border
		brd = uni.getData()->getData(uni)->getrBrd();

		//Check if the border reaches the first position we have to check
		if((int32_t) offset <= brd){
			//Testing
			// cout << "5 Option 2" << endl;

			//Get left border
			brd = uni.getData()->getData(uni)->getlBrd();

			//Check if our current position is still covered 
			if(brd <= (int32_t) offset){
				//Testing
				// cout << "6 Option 2" << endl;

				return 0;
			} else{
				//Testing
				// cout << "6 Option 1" << endl;

				//Get the number of positions which are still covered
				return brd - (int32_t) offset;
			}
		}

		//Testing
		// cout << "5 Option 1" << endl;

		//Calculate number of covered positions
		return uni.size - (brd + 1);
	} else{
		//Testing
		// cout << "4 Option 2" << endl;
		// cout << "rBrd: " << uni.getData()->getData(uni)->getrBrd() << endl;

		//Get left border
		brd = uni.getData()->getData(uni)->getlBrd();

		//Check if the border reaches the first position we have to check
		if(brd <= (int32_t) offset){
			//Testing
			// cout << "7 Option 2" << endl;
			// cout << "lBrd: " << brd << endl;

			//Get right border
			brd = uni.getData()->getData(uni)->getrBrd();

			//Check if our current position is still covered
			if(brd >= (int32_t) offset){
				//Testing
				// cout << "8 Option 2" << endl;

				return 0;
			} else{
				//Testing
				// cout << "8 Option 1" << endl;
				// cout << "brd: " << brd << endl;
				// cout << "offset: " << offset << endl;
				
				return (int32_t) offset - brd;
			}
		}

		//Testing
		// cout << "7 Option 1" << endl;

		return brd;
	}
}

//This function checks if the search criteria are fulfilled for a region on a unitig starting at some offset and having a certain length
bool vfySrchCrit(UnitigColorMap<seedlist> unitig, const uint32_t &offset, const int32_t &length, const uint32_t &quorum, const list<pair<string, size_t>> &srchColSet);//TODO Would it make sense to make this an inline function?

//This function calculates the correct offset position depending on which strand we are
inline uint32_t compOffset(const uint32_t &offset, const int32_t &seedLen, const uint32_t &seqLen, const bool &onRefStrand){
	//If we are on the reference strand there is nothing to do
	if(onRefStrand) return offset;

	//Calculate coordinates on reference strand and exit
	return seqLen - offset - seedLen;
}

//This function performs the seed detection if a search color set is given
void detectSeeds(const int32_t &k, const int32_t &minSeedLength, const size_t &numSmers, const uint32_t &quorum, const uint32_t &profileSize, const uint32_t *qProfile, const string &q, const UnitigColorMap<seedlist> *uArr, const struct s_mer_pos *posArray, hit *hitArr, const list<pair<string, size_t>> &searchSet, const bool isRefSeq);

// //This function extends all seeds found on the queries reference strand considering a quorum
// void extendRefSeeds(ColoredCDBG<seedlist> &cdbg, const string &q, const int16_t &X, hit *hitArr, const uint32_t &quorum);//TODO This function still needs to be tested!

// //This function extends all seeds found on the queries reverse complement considering a quorum
// void extendRevCompSeeds(ColoredCDBG<seedlist> &cdbg, const string &q, const int16_t &X, hit *hitArr, const uint32_t &quorum);//TODO This function still needs to be tested!

//This function extends all seeds found on the queries reference strand considering a quorum and a search color set
void extendRefSeeds(const ColoredCDBG<seedlist> &cdbg, const string &q, const int16_t &X, hit *hitArr, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function extends all seeds found on the queries reverse complement considering a quorum and a search color set
void extendRevCompSeeds(const ColoredCDBG<seedlist> &cdbg, const string &q, const int16_t &X, hit *hitArr, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

// //This function calculates a banded, semi-global, gapped alignment on a list of results considering a quorum and outputs the result if demanded
// void calcGappedAlignment(const list<hit*> &resList, const string &q, const int16_t &X, const bool &noOutput, const uint32_t &quorum);//TODO This function still needs to be tested!

//This function calculates a banded, semi-global, gapped alignment on a list of results considering a quorum and a search color set and outputs the result if demanded
void calcGappedAlignment(ColoredCDBG<seedlist> &cdbg, list<hit*> &resList, const string &q, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const double &lambda, const double &C);//TODO: Tests for this function need to be adjusted!

//This function performs the actual graph search for a query
void searchQuery(ColoredCDBG<seedlist> &cdbg, const int32_t &kMerLength, const int32_t &minSeedLength, const size_t &numSmers, const uint32_t &quorum, const uint32_t &profileSize, const uint32_t *qProfile, const string &q, const SrchStrd &strand, const UnitigColorMap<seedlist> *uArr, const struct s_mer_pos *posArray, const list<pair<string, size_t>> &searchColors, const int16_t &X, const bool &calcRT, uint16_t nRes, const double &lambda, const double &C, const double &eLim);//TODO Test for this function need to be adjusted!

//This function frees all memory additionally allocated for an hit array
void freeHitArray(hit *arr, uint32_t &arrLen);

#endif