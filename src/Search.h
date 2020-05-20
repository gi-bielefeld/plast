#ifndef SEARCH_HPP
#define SEARCH_HPP

#include "GappedAlignment.h"

//Default quorum used for searches
#define DEFAULT_QUORUM 1
//Default number of results to be outputted
#define DEFAULT_NB_RES 250

//All variations of strands we may perform a search
enum SrchStrd {Plus, Minus, Both};

//This function checks if a unitig fulfills a quorum completely (ATTENTION: This function does not work for a quorum of 1)
inline bool isCovered(const UnitigColorMap<seedlist> &uni, const uint32_t &quorum){
	bool fstId = true;
	uint16_t counter = 0;
	int32_t allwdToMs;
	size_t curID = SIZE_MAX, lstID = 0;

	//Calculate how many colors we are allowed to miss before it is clear that we cannot fulfill the quorum anymore
	allwdToMs = uni.getData()->getUnitigColors(uni)->colorMax(uni) + 1 - quorum;

	//Iterate over unitig's colors
	for(UnitigColors::const_iterator i = uni.getData()->getUnitigColors(uni)->begin(uni); curID != i.getColorID(); i.nextColor()){
		//Update current id
		curID = i.getColorID();

		//If we are dealing with the first color there is no last one
		if(fstId){
			//Decrement by the number of colors we have skipped starting from 0
			allwdToMs -= curID;
			fstId = false;
		} else{
			//Decrement by the number of colors we have skipped between last and current color if they are not consecutive
			allwdToMs -= curID - lstID - 1;
		}

		//Check if number of colors allowed to miss is exceeded
		if(allwdToMs < 0) return false;

		//Check if the color is present on the complete unitig
		if(uni.getData()->getUnitigColors(uni)->contains(uni, curID)){
			//Increment counter and check if we are done
			if(++counter == quorum) return true;
		} else{
			//Decrement number of colors we are still allowed to miss and check if we can stop
			if(--allwdToMs < 0) return false;
		}

		//Update last last id
		lstID = curID;
	}

	return false;
}

//This function checks if a unitig fulfills the search criteria completely
inline bool isCovered(const UnitigColorMap<seedlist> &uni, const uint32_t &quorum, const list<pair<string, size_t>> &srchColSet){
	uint16_t counter = 0;
	uint32_t allwdToMs;

	//Calculate how many colors we are allowed to miss before it is clear that we cannot fulfill the quorum anymore
	allwdToMs = srchColSet.size() - quorum;

	//Go through search color set
	for(list<pair<string, size_t>>::const_iterator col = srchColSet.begin(); col != srchColSet.end(); ++col){
		//Check whether the current color is present at our current position
		if(uni.getData()->getUnitigColors(uni)->contains(uni, col->second)){
			//Check whether our quorum is already reached
			if(quorum == ++counter) return true;
		} else{
			//Decrement number of colors we are still allowed to miss and check if we can stop
			if(allwdToMs-- == 0) return false;
		}
	}

	return false;
}

//This functions checks up to which offsets search criteria are fulfilled for a unitig and saves the result in the unitig's seedlist. If even the first k-mer from either side is not covered by the search criteria positions are set to -1 and unitig length respectively
inline void calcSrchCritBrds(UnitigColorMap<seedlist> uni, const uint32_t &quorum, const list<pair<string, size_t>> &srchColSet){
	bool lFix, rFix;
	int32_t lBrd, rBrd;

	//If the maximum color id (+1, since ids start at 0) is already smaller than our quorum it can never be fulfilled
	if(quorum > uni.getData()->getUnitigColors(uni)->colorMax(uni) + 1){
		//Save offsets
		uni.getData()->getData(uni)->setlBrd(-1);
		uni.getData()->getData(uni)->setrBrd(uni.size);
		return;
	}

	//Check whether a search set is given
	if(srchColSet.empty()){
		//If quorum is 1 we do not need to do anything
		if(quorum == 1 || isCovered(uni, quorum)){
			//Save offsets
			uni.getData()->getData(uni)->setlBrd(uni.size);
			uni.getData()->getData(uni)->setrBrd(0);
			return;
		}

		//Check if an iterative check makes no sense
		if(uni.len == 1){
			//Save offsets
			uni.getData()->getData(uni)->setlBrd(-1);
			uni.getData()->getData(uni)->setrBrd(uni.size);
			return;
		}

		//Initialize borders for iterative check
		lBrd = 0;
		rBrd = uni.size - uni.getGraph()->getK();
		//Make a copy of our unitig so that we can check both sides at the same time
		UnitigColorMap<seedlist> uniCpy = uni;
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
				//Move it
				++lBrd;
			} else{
				//Report border as fixed
				lFix = true;
			}

			//Check if we can move the right border
			if(!rFix && isCovered(uniCpy, quorum)){
				//Move it
				--rBrd;
			} else{
				//Report border as fixed
				rFix = true;
			}

			//Adjust unitigs
			uni.dist = lBrd;
			uniCpy.dist = rBrd;
		}
	} else{
		//Check the complete unitig
		if(isCovered(uni, quorum, srchColSet)){
			//Save offsets
			uni.getData()->getData(uni)->setlBrd(uni.size);
			uni.getData()->getData(uni)->setrBrd(0);
			return;
		}

		//Check if an iterative check makes no sense
		if(uni.len == 1){
			//Save offsets
			uni.getData()->getData(uni)->setlBrd(-1);
			uni.getData()->getData(uni)->setrBrd(uni.size);
			return;
		}

		//Initialize borders for iterative check
		lBrd = 0;
		rBrd = uni.size - uni.getGraph()->getK();
		//Make a copy of our unitig so that we can check both sides at the same time
		UnitigColorMap<seedlist> uniCpy = uni;
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
				//Move it
				++lBrd;
			} else{
				//Report border as fixed
				lFix = true;
			}

			if(!rFix && isCovered(uniCpy, quorum, srchColSet)){
				//Move it
				--rBrd;
			} else{
				//Report border as fixed
				rFix = true;
			}

			//Adjust unitigs
			uni.dist = lBrd;
			uniCpy.dist = rBrd;
		}
	}

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

//ATTENTION: This function only works correctly if a non-empty search color set is given (Why should that be?)
inline bool quorumFulfilled(UnitigColorMap<seedlist> uni, const uint32_t posU, size_t uniLen, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet){
	size_t ccount = 0, curID = 0;

	//Check to which side we want to calculate an alignment
	if(posU > 0){//This can also be false if we want to extend to the right but in that case uniLen is correctly set
		//If we want to extend to the right adjust unitig start position
		uni.dist = min((size_t) posU, uni.size - uni.getGraph()->getK());
	} else{
		//If we want to extented to the left adjust unitig length
		uni.len = min(uniLen, uni.len);
	}

	//Check if we have a search set to consider
	if(searchSet.empty()){
		//If quorum is 1 there is nothing to do
		if(quorum == 1) return true;

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
	} else{
		//Go through search color set
		for(list<pair<string, size_t>>::const_iterator col = searchSet.begin(); col != searchSet.end(); ++col){
			//Check whether the current color is present
			if(uni.getData()->getUnitigColors(uni)->contains(uni, col->second)){
				//Check whether our quorum is already reached
				if(quorum <= ++ccount) return true;
			}
		}
	}

	return false;
}

//This function looks up the number of positions which fulfill the search criteria depending at which offset we want to start and into which direction we extend with regard to the reference strand
inline int32_t getSrchCritCov(UnitigColorMap<seedlist> uni, const uint32_t &quorum, const list<pair<string, size_t>> &srchCols, const uint32_t &offset, const bool &toRightOnRef){
	int32_t brd;

	//Check if search criteria have already been checked for this unitig and check if necessary
	if(!uni.getData()->getData(uni)->srchCritChckd()) calcSrchCritBrds(uni, quorum, srchCols);

	//Check if search criteria are not fulfilled at all on this unitig
	if(uni.getData()->getData(uni)->getlBrd() < 0 && uni.getData()->getData(uni)->getrBrd() == (int32_t) uni.size) return 0;

	//Check if unitig is completely covered
	if(uni.getData()->getData(uni)->getlBrd() > uni.getData()->getData(uni)->getrBrd()) return uni.size;

	//Check into which direction we want to go
	if(toRightOnRef){
		//Get right border
		brd = uni.getData()->getData(uni)->getrBrd();

		//Check if the border reaches the first position we have to check
		if((int32_t) offset <= brd){
			//Get left border
			brd = uni.getData()->getData(uni)->getlBrd();

			//Check if our current position is still covered 
			if(brd <= (int32_t) offset){
				return 0;
			} else{
				//Get the number of positions which are still covered
				return brd - (int32_t) offset;
			}
		}

		//Calculate number of covered positions
		return uni.size - (brd + 1);
	} else{
		//Get left border
		brd = uni.getData()->getData(uni)->getlBrd();

		//Check if the border reaches the first position we have to check
		if(brd <= (int32_t) offset){
			//Get right border
			brd = uni.getData()->getData(uni)->getrBrd();

			//Check if our current position is still covered
			if(brd >= (int32_t) offset){
				return 0;
			} else{
				return (int32_t) offset - brd;
			}
		}

		return brd;
	}
}

//This function checks if the search criteria are fulfilled for a region on a unitig starting at some offset and having a certain length
bool vfySrchCrit(UnitigColorMap<seedlist> unitig, const uint32_t &offset, const int32_t &length, const uint32_t &quorum, const list<pair<string, size_t>> &srchColSet);

//This function calculates the correct offset position depending on which strand we are
inline uint32_t compOffset(const uint32_t &offset, const int32_t &seedLen, const uint32_t &seqLen, const bool &onRefStrand){
	//If we are on the reference strand there is nothing to do
	if(onRefStrand) return offset;

	//Calculate coordinates on reference strand and exit
	return seqLen - offset - seedLen;
}

//This function performs the seed detection if a search color set is given
void detectSeeds(const int32_t &k, const int32_t &minSeedLength, const size_t &numSmers, const uint32_t &quorum, const uint32_t &profileSize, const uint32_t *qProfile, const string &q, const UnitigColorMap<seedlist> *uArr, const struct s_mer_pos *posArray, Hit *hitArr, const list<pair<string, size_t>> &searchSet, const bool isRefSeq);

//This function extends all seeds found on the queries reference strand considering a quorum and a search color set
void extendRefSeeds(const ColoredCDBG<seedlist> &cdbg, const string &q, const int32_t &minSdLen, const int16_t &X, Hit *hitArr, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function extends all seeds found on the queries reverse complement considering a quorum and a search color set
void extendRevCompSeeds(const ColoredCDBG<seedlist> &cdbg, const string &q, const int32_t &minSdLen, const int16_t &X, Hit *hitArr, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet);

//This function calculates a banded, semi-global, gapped alignment on a list of results considering a quorum and a search color set and outputs the result if demanded
void calcGappedAlignment(ColoredCDBG<seedlist> &cdbg, list<Hit*> &resList, const string &q, const int16_t &X, const uint32_t &quorum, const list<pair<string, size_t>> &searchSet, const double &lambda, const double &C);

//This function performs the actual graph search for a query
void searchQuery(ColoredCDBG<seedlist> &cdbg, const int32_t &kMerLength, const int32_t &minSeedLength, const size_t &numSmers, const uint32_t &quorum, const uint32_t &profileSize, const uint32_t *qProfile, const string &q, const SrchStrd &strand, const UnitigColorMap<seedlist> *uArr, const struct s_mer_pos *posArray, const list<pair<string, size_t>> &searchColors, const int16_t &X, const bool &calcRT, uint16_t nRes, const double &lambda, const double &lambdaGap, const double &C, const double &Cgap, const double &eLim, const bool &colOut, const bool &isSim);

//This function frees all memory additionally allocated for an hit array
void freeHitArray(Hit *arr, uint32_t &arrLen);

#endif
