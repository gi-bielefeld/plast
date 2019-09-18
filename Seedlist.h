#ifndef SEEDLIST_HPP
#define SEEDLIST_HPP

#include "Seed.h"

class seedlist : public CCDBG_Data_t<seedlist> {

	public:

		seedlist(): nextSeedRef(NULL), nextSeedRcomp(NULL), frstSeedRcomp(NULL), lBrd(0), rBrd(0) {}

		//NOTE: Actually, methods join, serialize and sub should be overloaded here. Since the seedlist will only come into play after the graph has already been created and is not used anymore before the graph is written into a file it should be not necessary to overload them.

		struct seed* getSeed(const bool &onRefStrand) const {//get last element of seed list of unextended seeds that lie on unitig's reference strand
			if(onRefStrand) return nextSeedRef;
			return nextSeedRcomp;
		} 

		//struct seed* getLastProcSeed() const { return nextProcSeed; }//get last element of seed list of extended seeds

		//Set a seed to be the first seed in a respective list ATTENTION: Using this for the reverse complementary strand makes the use of pushSeed impossible because frstSeedRcomp is not set correctly!
		void setSeed(struct seed* newSeed, const bool &onRefStrand){ onRefStrand ? nextSeedRef = newSeed : nextSeedRcomp = newSeed; }

		void pushSeed(struct seed* newSeed, const bool &onRefStrand){
			if(onRefStrand){
				//Testing
				// cout << "1 Option 1" << endl;

				//Make the previous seed the new seed's successor
				newSeed->nextSeed = nextSeedRef;
				//Make the new seed the first seed in the seed list
				nextSeedRef = newSeed;
			} else{
				//Testing
				// cout << "1 Option 2" << endl;

				//Make the new seed the last seed in the seed list (to ensure seeds are ordered by decreasing query offset)
				if(!nextSeedRcomp){
					//Testing
					// cout << "2 Option 2" << endl;

					nextSeedRcomp = newSeed;
				} else{
					//Testing
					// cout << "2 Option 1" << endl;

					//Link the previously inserted seed to the new seed
					frstSeedRcomp->nextSeed = newSeed;
				}
				//Store a pointer to the last inserted seed
				frstSeedRcomp = newSeed;
			}
		}

		struct seed* getLastInSd() const { return frstSeedRcomp; }//get the last inserted seed from the seed list of the reverse complementary query

		//void setLastProcSeed(struct seed* newSeed){ nextProcSeed = newSeed; }//insert new element at the seed list's end

		int32_t getlBrd() const { return lBrd; }

		int32_t getrBrd() const { return rBrd; }

		void setlBrd(int32_t offset){ lBrd = offset; }

		void setrBrd(int32_t offset){ rBrd = offset; }

		bool srchCritChckd() const { return lBrd != rBrd; }

	private:

		struct seed *nextSeedRef;
		struct seed *nextSeedRcomp;
		struct seed *frstSeedRcomp;
		//struct seed *nextProcSeed;

		//The lowest offset position in the unitig at which search criteria are not fulfilled anymore
		int32_t lBrd;
		//The highest offset position in the unitig at which search criteria are not fulfilled anymore
		int32_t rBrd;
};

//This function removes a seed from a seed list, i.e. excluding it from the list and deleting it. Returns true if the seed was not the head of its seed list
bool removeSeed(seedlist* sl, struct seed*& seedToRm, struct seed*& prevSeed);

//This function extracts a seed from a seed list of processed seeds. Returns true if there might be another seed to reach
bool extractProcSeed(seedlist* sl, struct seed*& seedToExtr, struct seed*& prevSeed);

//This function deletes a seed from a seed list
inline void delSeed(struct seed*& seedToDel, struct seed*& prevSeed){
	//Check that the seed we want to delete is not the first in the seed list
	if(prevSeed != NULL){
		//Unchain the seed
		prevSeed->nextSeed = seedToDel->nextSeed;
		//Release the seeds memory
		free(seedToDel);
	} else{
		//Release the seeds memory
		free(seedToDel);
	}
}

//This function deletes a seed that has been used during an extension over multiple unitigs
void delExtPtrSeed(const bool& inSeedList, struct seed*& seedToDel, seedlist* seedList);

//This function searches for the closest seed to reach during an extension and returns it if it exists (otherwise NULL)
struct seed* searchRightNeighbor(struct seed* sLSeed, const uint32_t &iniExtSeedOffsQ, const uint32_t &extLen, const uint32_t &curUPos, struct seed*& prevSeed);

//This function searches for the closed seed to reach during an extension and returns it if it exists (otherwise NULL)
struct seed* searchLeftNeighbor(struct seed *seedList, const uint32_t &qOff, const uint32_t &uOff, struct seed*& prevSeed);

#endif