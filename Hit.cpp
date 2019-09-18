#include "Hit.h"
#include "Smer.h"

//This function inserts a hit into the hit list where all hits are ordered increasingly by their left border's q index
void insertHit(hit*& hitToIns, hit* hitList){//TODO This function is apparently not used anymore and I had to modify it to avoid some compilation errors. There is no guarantee that is still does what was initially intended!
	//Check whether our hit list is yet empty
	if(hitList == NULL){
		//Make hitToIns the first element of our hit list
		hitList = hitToIns;

		return;
	}

	//Check if we can insert hitToIns in front of the first element
	if(hitList->offQ > hitToIns->offQ){
		//Put hit as the first element
		hitToIns->nextHit = hitList;
		hitList = hitToIns;

		return;
	}

	hit *currHit = hitList;

	//Go through the hit list
	while(currHit->nextHit != NULL){
		//Check whether we have reached the right place to insert the hit
		if(currHit->nextHit->offQ > hitToIns->offQ){
			//Incorporate the hit into the hit list
			hitToIns->nextHit = currHit->nextHit;
			currHit->nextHit = hitToIns;

			return;
		}

		//Move one hit further
		currHit = currHit->nextHit;
	}

	//If we have reached the hit list's end make hitToIns the new last element
	currHit->nextHit = hitToIns;
}

//This function scans through a hit list decreasingly ordered by score and replaces the worst entry by the given hit if it is not the worst one itself
void replWorseRes(list<hit*> &hList, hit* hit){
	//Go through hit list
	for(list<class hit*>::const_iterator i = hList.begin(); i != hList.end(); ++i){
		//Compare current entry with hit to insert
		if((*i)->score < hit->score){
			//Insert hit
			hList.insert(i, hit);
			//Remove worst hit
			hList.pop_back();
			return;
		}
	}
}

//This function scans through a hit list and inserts the given hit right in front of the first entry which score is smaller or equal to the given hit's score
void insRes(list<hit*> &hList, hit* hit){
	//Go through hit list
	for(list<class hit*>::const_iterator i = hList.begin(); i != hList.end(); ++i){
		//Compare current entry with hit to insert
		if((*i)->score < hit->score){
			//Insert hit
			hList.insert(i, hit);

			//Testing
			/*for(list<struct hit*>::const_iterator j = hList.begin(); j != hList.end(); ++j){
				cout << (*j)->score << endl;
			}*/

			return;
		}
	}

	//Insert hit at the end
	hList.push_back(hit);

	return;
}

//This function retreats an ungapped extension of a seed on the reverse complementary strand if its right border ends within the k-1 overlap in a unitig sequence's end
void makeValidBorders(hit *hit, const string &q){//TODO This function definitely needs to be revised completely! (I have already taken some notes) For the moment we just leaf it as it is (I guess it is okay now as it is)%TODO Now after the merge this function needs a revision again, but maybe we actually do not need it anymore!
	uint32_t k = hit->origUni.getGraph()->getK();
	int32_t overhang;

	//Calculate how much the right border has to be adjusted
	if((overhang = hit->offU - (hit->origUni.size - k)) > 0){
		//Testing
		// cout << "overhang: " << overhang << endl;
		// cout << "offQ: " << hit->offQ << endl;

		//Check if hit is too small to get rid of overhang
		if((uint32_t) overhang >= hit->length || hit->offQ < (uint32_t) overhang){
			//Make sure that we do not consider this hit any further//TODO Maybe there is a nicer way to do so! Of course, we could move the whole seed on all successors, but most likely this will not lead to anything
			hit->score = 0;
			return;
		}

		//Testing
		// cout << "Let's go into the while loop" << endl;

		//Iteratively retreat the extension
		while(overhang > 0){
			//Get unitig sequence
			string uSeq = hit->origUni.mappedSequenceToString();
			//Compare sequences and reduce the score accordingly
			hit->score -= compUScore(uSeq[hit->offU], q[hit->offQ]);
			//Decrease hit's length
			--hit->length;
			//Decrease offsets
			--hit->offQ;
			--hit->offU;
			//Decrease overhang	
			--overhang;

			//Testing
			// cout << "Loop" << endl;
		}
	}

	//Testing
	// cout << "makeValidBorders done" << endl;
}

//This function compresses an extension path of a left or right extension
struct ExtPth cmprExtPth(const list<uint16_t>& extPth){//TODO: Check whether it makes sense to declare this as inline!
	uint16_t nbBytes, i = 0, shifts = 0;
	struct ExtPth ePath;

	//Get number of elements in our path
	ePath.nbElem = extPth.size();
	//Calculate how many bytes we need to save the path
	nbBytes = (ePath.nbElem / SUCCESSORS_PER_BYTE) + ((ePath.nbElem % SUCCESSORS_PER_BYTE) == 0 ? 0 : 1);
	//Allocate space for compressed path
	ePath.path = (unsigned char*) malloc(nbBytes);

	//Walk through extension path
	for(list<uint16_t>::const_iterator j = extPth.begin(); j != extPth.end(); ++j){
		//Check if current byte is full
		if(shifts == SUCCESSORS_PER_BYTE){
			//Move to next byte
			++i;
			//Reset shifts
			shifts = 0;
		} else{
			//Shift all information by 2 bits to the left
			ePath.path[i] <<= BITS_PER_SUCCESSOR;	
		}

		//Add next successor to byte (we have to subtract one to shift to the correct scale 1-4 -> 0-3)
		ePath.path[i] |= ((*j) - 1);
		++shifts;
	}

	//Fill up the last byte if necessary
	if(ePath.nbElem > 0) ePath.path[i] <<= BITS_PER_SUCCESSOR * (SUCCESSORS_PER_BYTE - shifts);	

	return ePath;
}

//This function decompresses an extension path object
const list<uint16_t> decmprExtPth(const struct ExtPth &cmpPth){
	unsigned char successor;
	unsigned char *cmpSucs = cmpPth.path;
	uint16_t nbBytes, shifts, sucLeft = cmpPth.nbElem;
	list<uint16_t> path;

	//Calculate how many bytes our compressed path consists of
	nbBytes = cmpPth.nbElem / SUCCESSORS_PER_BYTE + ((cmpPth.nbElem % SUCCESSORS_PER_BYTE) == 0 ? 0 : 1);

	//Go through compressed path
	for(uint16_t i = 0; i < nbBytes; ++i){
		//We start without shifting
		shifts = 0;

		while(shifts < SUCCESSORS_PER_BYTE){
			//Get rid of information in front of current successor
			successor = *cmpSucs << shifts * BITS_PER_SUCCESSOR;
			//Shift successor to end of byte
			successor >>= 6;
			//Add successor to extension path
			path.push_back(((uint16_t) successor) + 1);

			//Check if we have decompressed the complete path
			if(--sucLeft == 0) break;

			++shifts;
		}

		++cmpSucs;
	}

	return path;
}