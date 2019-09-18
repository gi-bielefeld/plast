#ifndef HIT_HPP
#define HIT_HPP

// #include "Seedlist.h"

//A unitig may have 4 successors. Thus, we can store 4 successors per byte.
#define SUCCESSORS_PER_BYTE 4
//We need 2 bit to encode a successor
#define BITS_PER_SUCCESSOR 2

//Stores the path an extension takes through the graph
struct ExtPth{
	//Number of elements in the path
	uint16_t nbElem;
	//Actual compressed path
	unsigned char *path;
};

//Stores an alignment
struct Algn{//TODO Instead of storing both sequences it would be fine as well just to store the underlying edit sequence!
	//Alignment sequence from the graph
	string aSeqG;
	//Alignment sequence from query
	string aSeqQ;
};

class hit {

	public:
		
		//Default constructor
		hit () { lExt.path = NULL; rExt.path = NULL; }//TODO Do we really have to set this?

		//Copy constructor
		hit (const hit &h) {
			length = h.length;
			score = h.score;
			offU = h.offU;
			offQ = h.offQ;
			origUni = h.origUni;
			lExt = h.lExt;
			rExt = h.rExt;
			gAlgn = h.gAlgn;
			nextHit = h.nextHit;
		}

		//The hit's length
		uint32_t length;
		//The hit's score
		uint32_t score;
		//The initial seed's offsets
		uint32_t offU, offQ;
		//Initial unitig the extension was started
		UnitigColorMap<seedlist> origUni;
		//Left-side extension path of extended seed
		struct ExtPth lExt;
		//Right-side extension path of extended seed
		struct ExtPth rExt;
		//Gapped alignment of hit
		struct Algn gAlgn;
		//Pointer to the next hit
		hit *nextHit;
};

//This function inserts a hit into the hit list where all hits are order increasingly by their left border's q index
void insertHit(hit*& hitToIns, hit* hitList);

//This function scans through a hit list decreasingly ordered by score and replaces the worst entry by the given hit if it is not the worst one itself
void replWorseRes(list<hit*> &hList, hit* hit);

//This function scans through a hit list and inserts the given hit right in front of the first entry which score is smaller or equal to the given hit's score
void insRes(list<hit*> &hList, hit* hit);

//This function compresses the extension paths of left and right extension and returns a merged extension path object
struct ExtPth cmprExtPth(const list<uint16_t>& extPth);

// 	//Walk through extension path
// 	for(list<uint16_t>::const_iterator j = extPth.begin(); j != extPth.end(); ++j){
// 		//Check if current byte is full
// 		if(shifts == SUCCESSORS_PER_BYTE){
// 			//Move to next byte
// 			++pos;
// 			//Reset shifts
// 			shifts = 0;
// 		} else{
// 			//Shift all information by 2 bits to the left
// 			cmprPth[pos] <<= BITS_PER_SUCCESSOR;
// 			++shifts;	
// 		}

// 		//Testing
// 		//cout << "cmprPth[pos]:" << (uint16_t) cmprPth[pos] << "pos:" << pos << " shifts:" << shifts << endl;

// 		//Add next successor to byte (we have to subtract one to shift to the correct scale 1-4 -> 0-3)
// 		cmprPth[pos] |= ((*j) - 1);
// 	}

// 	//Testing
// 	//cout << "I don't care!" << endl;
// }

//This function decompresses an extension path object
const list<uint16_t> decmprExtPth(const struct ExtPth &cmpPth);

//This function retreats an ungapped extension of a seed on the reverse complementary strand if its right border ends within the k-1 overlap in a unitig sequence's end
void makeValidBorders(hit *hit, const string &q);

#endif