#ifndef HIT_HPP
#define HIT_HPP

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
struct Algn{
	//Alignment sequence from the graph
	string aSeqG;
	//Alignment sequence from query
	string aSeqQ;
};

class Hit {

	public:
		
		//Default constructor
		Hit () { lExt.path = NULL; rExt.path = NULL; }

		//Copy constructor
		Hit (const Hit &h) {
			length = h.length;
			score = h.score;
			offU = h.offU;
			offQ = h.offQ;
			eval = h.eval;
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
		//The hit's e-value
		double eval;
		//Initial unitig the extension was started
		UnitigColorMap<seedlist> origUni;
		//Left-side extension path of extended seed
		struct ExtPth lExt;
		//Right-side extension path of extended seed
		struct ExtPth rExt;
		//Gapped alignment of hit
		struct Algn gAlgn;
		//Pointer to the next hit
		Hit *nextHit;
};

//This function inserts a hit into the hit list where all hits are order increasingly by their left border's q index
void insertHit(Hit*& hitToIns, Hit* hitList);

//This function scans through a hit list decreasingly ordered by score and replaces the worst entry by the given hit if it is not the worst one itself
void replWorseRes(list<Hit*> &hList, Hit* hit);

//This function scans through a hit list and inserts the given hit right in front of the first entry which score is smaller or equal to the given hit's score
void insRes(list<Hit*> &hList, Hit* hit);

//This function retreats an ungapped extension of a seed on the reverse complementary strand if its right border ends within the k-1 overlap in a unitig sequence's end
void makeValidBorders(Hit *hit, const string &q);

//This function compresses the extension paths of left and right extension and returns a merged extension path object
struct ExtPth cmprExtPth(const list<uint16_t>& extPth);

//This function decompresses an extension path object
const list<uint16_t> decmprExtPth(const struct ExtPth &cmpPth);

//This function compares two hits based on their e-values. Returns true if fh is smaller than sh.
inline bool compEvals(const Hit *fh, const Hit *sh){ return (fh->eval < sh->eval); };

#endif