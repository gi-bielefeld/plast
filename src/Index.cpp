#include "Index.h"
#include "Sequence.h"

//This function writes PLAST's q-gram profile to disk
void saveProfile(ofstream &fStr, const uint32_t &profSize, uint32_t* qProfile){
	uint32_t count, shift, i = 0;

	//Iterate over the profile
	do{
		//Get the first profile value
		count = qProfile[i];
		//Initialize the first shift
		shift = 1;
		//Move to the next entry
		++i;

		//Collect all entries with same value
		while(i < profSize && count == qProfile[i]){
			++shift;
			//Move to the next entry
			++i;
		}

		//Write byte array to the graph
		fStr.write((char*)&shift, sizeof(shift));
		fStr.write((char*)&count, sizeof(count));
	}while(i < profSize);
}

//This function loads the q-gram profile from a file and returns it
uint32_t* loadQProf(ifstream &iFile, const uint32_t size){//TODO: Can we make the ifstream const?
	//String to store read in line
	string l;
	//Allocate space for the array
	uint32_t *profile = (uint32_t*) malloc(size * sizeof(uint32_t));

	//Go through the lines of the input file which belong to the index
	for(uint32_t i = 0; i < size; ++i){
		//Try to read next line
		if(!getline(iFile, l)){
			cerr << "An error occured while reading index file. Index might be corrupted!" << endl;
			exit(EXIT_FAILURE);
		}

		//Store the read in value in the index
		profile[i] = stoul(l);
	}



	return profile;
}

//This function saves PLAST's index data structures as binary
void saveIndexBin(const char *filename, const int32_t& sdLen, const uint32_t& profSize, uint32_t* qProfile, const uint32_t numUni, UnitigColorMap<UnitigInfo>* uniArr, struct S_mer_pos*& pos, const size_t &seedNum, const int32_t &k, const bool& advIdx){
	uint16_t qrm;
	size_t i;
	int32_t brd;

	//Open a file to write
	ofstream destFile(filename, ios::out | ios::binary);

	//Check whether the file could be opened
	if(!destFile.is_open()){
		cerr << "ERROR: Destination file to save index could not be opened" << endl;
		exit(EXIT_FAILURE);
	}

	//Save general meta data//

	//Write index version
	destFile.write(INDEX_VERSION, sizeof(INDEX_VERSION));
	//Write minimum seed length used for this index
	destFile.write((char*) &sdLen, sizeof(int32_t));
	//Write if this index contains quorum information
	destFile.write((char*) &advIdx, sizeof(bool));
	//Save the q-gram profile
	saveProfile(destFile, profSize, qProfile);

	//Calculate the number of bytes needed to store a k-mer
	const uint16_t nBytes = (k / BASES_PER_BYTE) + (k % BASES_PER_BYTE == 0 ? 0 : 1);
	//Iterate over the graph
	for(i = 0; i < numUni; ++i){
		//Write the compacted representation of current unitig's starting k-mer to file
		destFile.write(cmprSeq(uniArr[i].referenceUnitigToString().substr(0, k), nBytes), nBytes);

		//Check if we have to save precalculated quorum information
		if(advIdx){
			//Write it to file
			qrm = uniArr[i].getData()->getData(uniArr[i])->getGlobQrm();
			destFile.write((char*) &qrm, sizeof(uint16_t));
			qrm = uniArr[i].getData()->getData(uniArr[i])->getLBrdQrm();
			destFile.write((char*) &qrm, sizeof(uint16_t));
			qrm = uniArr[i].getData()->getData(uniArr[i])->getRBrdQrm();
			destFile.write((char*) &qrm, sizeof(uint16_t));
			brd = uniArr[i].getData()->getData(uniArr[i])->getPrecLBrd();
			destFile.write((char*) &brd, sizeof(int32_t));
			brd = uniArr[i].getData()->getData(uniArr[i])->getPrecRBrd();
			destFile.write((char*) &brd, sizeof(int32_t));
		}
	}

	//Free unitig array
	free(uniArr);

	//Save pos array

	//Walk through pos array
	for(i = 0; i < seedNum; ++i){
		//Write unitig id to file
		destFile.write((char*) &pos[i].unitig_id, sizeof(uint32_t));
		//Write offset to file
		destFile.write((char*) &pos[i].offset, sizeof(uint32_t));
	}

	//Free pos array
	free(pos);

	//Close the file
	destFile.close();
}

//This function loads the indices from a binary file
void loadIndexesBin(const char* fName, int32_t& sdLen, uint32_t &qProfSize, uint32_t *&qProf, ColoredCDBG<UnitigInfo>& graph, UnitigColorMap<UnitigInfo>*& uniArr, size_t &posSize, struct S_mer_pos*& pos, bool& advIdx){
	char *memBlock, *qrmBlock = NULL, *brdBlock = NULL;
	const int32_t kMerSize = graph.getK();
	//The number of bytes needed to store a k-mer
	const uint16_t nBytes = (kMerSize / BASES_PER_BYTE) + (kMerSize % BASES_PER_BYTE == 0 ? 0 : 1);
	uint32_t i, shift;
	//Open file
	ifstream idxFile(fName, ios::in|ios::binary|ios::ate);

	//Check if the file is open
	if(!idxFile.is_open()){
		cerr << "ERROR: Index file could not be opened" << endl;
		exit(EXIT_FAILURE);
	}

	//Get file's size
	streampos fSize = idxFile.tellg();
	//Go to file's beginning
	idxFile.seekg(0, ios::beg);

	//Check meta data//

	//Initialize memory block to load index version number
	memBlock = (char*) malloc(sizeof(INDEX_VERSION));
	//Load index version number
	idxFile.read(memBlock, sizeof(INDEX_VERSION));

	//Make sure this index is not deprecated
	if(strcmp(memBlock, INDEX_VERSION)){
		cerr << "ERROR: Index version is deprecated" << endl << "Please create a new one!" << endl;
		exit(EXIT_FAILURE);
	}

	//Free memory block
	free(memBlock);
	//Initialize memory block to load minimum seed length
	memBlock = (char*) malloc(sizeof(int32_t));
	//Load minimum seed length used to create this index
	idxFile.read(memBlock, sizeof(int32_t));
	sdLen = *((int32_t*) memBlock);
	//Free memory block
	free(memBlock);
	//Initialize memory block to load advanced index flag
	memBlock = (char*) malloc(sizeof(bool));
	//Load flag
	idxFile.read(memBlock, sizeof(bool));
	advIdx = *((bool*) memBlock);
	//Free memory block
	free(memBlock);

	//Load q-gram profile//

	//Calculate profile size
	qProfSize = (uint32_t) pow(SIGMAR, sdLen);
	//Initialize profile
	qProf = (uint32_t*) calloc(qProfSize, sizeof(uint32_t));
	//Initialize memory block to be loaded (We are going to load each shift value pair separately)
	memBlock = (char*) malloc(2 * sizeof(uint32_t));
	i = 0;

	//Check when q-gram profile loading is complete
	while(i < qProfSize){
		//Read next shift value pair
		idxFile.read(memBlock, 2 * sizeof(uint32_t));
		//Get shift
		shift = *((uint32_t*) memBlock);

		//Set profile entries
		for(uint32_t j = 0; j < shift; ++i, ++j){
			qProf[i] = ((uint32_t*) memBlock)[1];
		}
	}

	free(memBlock);

	//Initialize unitig array//

	//Initialize memory block to be loaded
	memBlock = (char*) malloc(nBytes);

	//Check if we have to load quorum infos
	if(advIdx){
		//Initialize memory block to load quorum values
		qrmBlock = (char*) malloc(sizeof(uint16_t));
		//Initialize memory block to load quorum borders
		brdBlock = (char*) malloc(sizeof(int32_t));
	}

	//Allocate memory
	uniArr = (UnitigColorMap<UnitigInfo>*) malloc(graph.size() * sizeof(UnitigColorMap<UnitigInfo>));
	i = 0;

	//Check if all k-mers have been read in
	while(i < graph.size()){
		//Read in next k-mer
		idxFile.read(memBlock, nBytes);
		//Search for k-mer's corresponding unitig
		Kmer k(decmprSeq(memBlock, kMerSize));
		uniArr[i] = graph.find(k, true);

		//Check whether we couldn't find a unitig
		if(uniArr[i].isEmpty){
			cerr << "ERROR: Error while loading unitig array. No corresponding unitig for loaded k-mer found" << endl;
			exit(EXIT_FAILURE);
			
		}
		
		//Assure that objects sequence represents the whole unitig sequence
		uniArr[i].dist = 0;
		uniArr[i].len = uniArr[i].size - kMerSize + 1;

		//Check if we have to load quorum infos
		if(advIdx){
			//Load quorum information
			idxFile.read(qrmBlock, sizeof(uint16_t));
			uniArr[i].getData()->getData(uniArr[i])->setGlobQrm(*((uint16_t*) qrmBlock));
			idxFile.read(qrmBlock, sizeof(uint16_t));
			uniArr[i].getData()->getData(uniArr[i])->setLBrdQrm(*((uint16_t*) qrmBlock));
			idxFile.read(qrmBlock, sizeof(uint16_t));
			uniArr[i].getData()->getData(uniArr[i])->setRBrdQrm(*((uint16_t*) qrmBlock));
			idxFile.read(brdBlock, sizeof(int32_t));
			uniArr[i].getData()->getData(uniArr[i])->setPrecLBrd(*((int32_t*) brdBlock));
			idxFile.read(brdBlock, sizeof(int32_t));
			uniArr[i].getData()->getData(uniArr[i])->setPrecRBrd(*((int32_t*) brdBlock));
		}

		++i;
	}
	
	free(memBlock);

	//Check if we have to load quorum infos
	if(advIdx){
		free(qrmBlock);
		free(brdBlock);
	}

	//Load pos array//

	//Get current position in read file
	streampos curPos = idxFile.tellg();
	//Set pos array size
	posSize = ((size_t) (fSize - curPos)) / sizeof(struct S_mer_pos);
	//Initialize memory block to be loaded
	memBlock = (char*) malloc((size_t) (fSize - curPos));
	//Load rest of file
	idxFile.read(memBlock, fSize - curPos);
	pos = (struct S_mer_pos*) memBlock;
}

//This function loads the indexes from a file
void loadIndexes(const string fName, const uint32_t qSize, uint32_t*& prof, CompactedDBG<UnitigInfo> &cdbg, uint32_t &seedNum, struct S_mer_pos*& lArr){
	//Open the file
	ifstream fStr(fName);

	//Check whether the file is open
	if(!fStr.is_open()){
		cerr << "ERROR: Index file could not be opened" << endl;
		exit(EXIT_FAILURE);
	}

	//Load q-gram index
	prof = loadQProf(fStr, qSize);
	//Close the file
	fStr.close();
}

//This function builds all necessary indexes
void buildIndex(ColoredCDBG<UnitigInfo> &cdbg, const int32_t &minSeedLength, uint32_t* qProf, size_t &numSmers, const uint32_t &profSize, struct S_mer_pos*& linkArr, UnitigColorMap<UnitigInfo>*& uniArray){
	uint32_t seqEnd, i = 0, k = cdbg.getK();
	int32_t pRank;
	//String to store the unitig sequence
	string uSeq;

	//Allocate memory for unitig array
	uniArray = (UnitigColorMap<UnitigInfo>*) malloc(cdbg.size() * sizeof(UnitigMap<UnitigInfo>));

	for(ColoredCDBG<UnitigInfo>::iterator it = cdbg.begin(); it != cdbg.end(); ++it){
		//Store the current unitig inside the unitig array
		uniArray[i] = *it;

		//Check whether the forward or the reverse complementary sequence is stored for the unitig
		if(it->strand){
			uSeq = it->referenceUnitigToString();

			/*NOTE: Actually, we do not consider all s-mers occurring in a unitig's sequence. This is, because the first k - s s-mers are identical between a unitig's and its successors's sequences except if they do not have a successor or predecessor. As a design decision, we say that we do not consider the last k - s s-mers of a unitig's sequences if the unitig has a successor. This decision might be revised in the future.*/
			if(it->getSuccessors().hasSuccessors()){
				seqEnd = it->size - k + 1;
			} else{
				seqEnd = it->size - minSeedLength + 1;
			}
		} else{
			//Save the current unitig
			UnitigColorMap<UnitigInfo> unitig = *it;
			//Change its strand
			unitig.strand = !unitig.strand;
			//Extract unitig's sequence
			uSeq = unitig.referenceUnitigToString();

			//If the unitig represents the reverse complementary sequence we need to check for predecessors instead of successors
			if(it->getPredecessors().hasPredecessors()){
				seqEnd = it->size - k + 1;
			} else{
				seqEnd = it->size - minSeedLength + 1;
			}
		}

		//The rank of the previously calculated q-gram
		pRank = -1;

		for(uint32_t j = 0; j < seqEnd; ++j){
			pRank = compRank(uSeq.substr(j, minSeedLength), pRank);
			++qProf[pRank];
			++numSmers;
		}

		++i;
	}

	for(i  = 1; i < profSize; ++i){
		qProf[i] += qProf[i-1];
	}

	linkArr = (struct S_mer_pos*) malloc(numSmers * sizeof(struct S_mer_pos));
	//Reset unitig counter
	i = 0;

	for(ColoredCDBG<UnitigInfo>::iterator it = cdbg.begin(); it != cdbg.end(); ++it){
		//The rank of the previously calculated q-gram
		pRank = -1;

		if(it->strand){
			uSeq = it->referenceUnitigToString();

			if(it->getSuccessors().hasSuccessors()){
				seqEnd = it->size - k + 1;
			} else{
				seqEnd = it->size - minSeedLength + 1;
			}
		} else{
			UnitigColorMap<UnitigInfo> unitig = *it;
			unitig.strand = !unitig.strand;
			uSeq = unitig.referenceUnitigToString();

			if(it->getPredecessors().hasPredecessors()){
				seqEnd = it->size - k + 1;
			} else{
				seqEnd = it->size - minSeedLength + 1;
			}
		}

		for(uint32_t j = 0; j < seqEnd; ++j){
			pRank = compRank(uSeq.substr(j, minSeedLength), pRank);
			//we need to subtract 1 here to even fill posArray's first entry
			linkArr[qProf[pRank] - 1].unitig_id = i;
			linkArr[qProf[pRank] - 1].offset = j;
			--qProf[pRank];
		}

		//Increment unitig counter
		++i;
	}
}