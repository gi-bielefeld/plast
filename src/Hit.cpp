#include "Hit.h"
#include "Smer.h"

//This function scans through a hit list decreasingly ordered by score and replaces the worst entry by the given hit if it is not the worst one itself
void replWorseRes(list<Hit*> &hList, Hit* hit){
	//Go through hit list
	for(list<class Hit*>::const_iterator i = hList.begin(); i != hList.end(); ++i){
		//Compare current entry with hit to insert
		if((*i)->score < hit->score){
			//Insert hit
			hList.insert(i, hit);
			//Remove extension paths
			frExtPth(hList.back()->rExt);
			frExtPth(hList.back()->lExt);
			//Remove worst hit
			hList.pop_back();
			return;
		}
	}
}

//This function scans through a hit list and inserts the given hit right in front of the first entry which score is smaller or equal to the given hit's score
void insRes(list<Hit*> &hList, Hit* hit){
	//Go through hit list
	for(list<class Hit*>::const_iterator i = hList.begin(); i != hList.end(); ++i){
		//Compare current entry with hit to insert
		if((*i)->score < hit->score){
			//Insert hit
			hList.insert(i, hit);
			return;
		}
	}

	//Insert hit at the end
	hList.push_back(hit);
	return;
}

//This function compresses an extension path of a left or right extension
struct ExtPth cmprExtPth(const list<uint16_t>& extPth){
	uint16_t nbBytes, i = 0, shifts = 0;
	struct ExtPth ePath;

	//Get number of elements in our path
	ePath.nbElem = extPth.size();

	//Check if we have to allocate memory for path compression
	if(ePath.nbElem){
		//Testing
		// cout << "ePath.nbElem: " << ePath.nbElem << endl;

		//Calculate how many bytes we need to save the path
		nbBytes = (ePath.nbElem / SUCCESSORS_PER_BYTE) + ((ePath.nbElem % SUCCESSORS_PER_BYTE) == 0 ? 0 : 1);
		//Allocate space for compressed path
		ePath.path = (unsigned char*) malloc(nbBytes);
		//Make all bits in the first byte zero
		*ePath.path = 0;
	}

	//Walk through extension path
	for(list<uint16_t>::const_iterator j = extPth.begin(); j != extPth.end(); ++j){
		//Check if current byte is full
		if(shifts == SUCCESSORS_PER_BYTE){
			//Move to next byte
			++i;
			//Reset shifts
			shifts = 0;
			//Make all bits zero
			ePath.path[i] = 0;
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

	//Free memory of compressed path if existing
	if(nbBytes) free(cmpPth.path);
	
	return path;
}