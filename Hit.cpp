#include "Hit.h"

//This function inserts a hit into the hit list where all hits are order increasingly by their left border's q index
void insertHit(struct hit*& hitToIns, struct hit* hitList){
	//Check whether our hit list is yet empty
	if(hitList == NULL){
		//Make hitToIns the first element of our hit list
		hitList = hitToIns;

		return;
	}

	//Check if we can insert hitToIns in front of the first element
	if(hitList->lSeedQoff > hitToIns->lSeedQoff){
		//Put hit as the first element
		hitToIns->nextHit = hitList;
		hitList = hitToIns;

		return;
	}

	struct hit *currHit = hitList;

	//Go through the hit list
	while(currHit->nextHit != NULL){
		//Check whether we have reached the right place to insert the hit
		if(currHit->nextHit->lSeedQoff > hitToIns->lSeedQoff){
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
void replWorseRes(list<struct hit*> &hList, struct hit* hit){
	//Go through hit list
	for(list<struct hit*>::const_iterator i = hList.begin(); i != hList.end(); ++i){
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
void insRes(list<struct hit*> &hList, struct hit* hit){
	//Go through hit list
	for(list<struct hit*>::const_iterator i = hList.begin(); i != hList.end(); ++i){
		//Testing
		//cout << "(*i)->score:" << (*i)->score << "\nhit->score:" << hit->score << endl;

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
void makeValidBorders(struct hit *hit, const string &q){//TODO This function definitely needs to be revised completely! (I have already taken some notes) For the moment we just leaf it as it is (I guess it is okay now as it is)
	uint32_t k = hit->lUnitig.getGraph()->getK();
	int32_t overhang;

	//Calculate how much the right border has to be adjusted
	if((overhang = hit->rSeedUoff - (hit->rUnitig.size - k)) > 0){
		//Testing
		// cout << "1 Option 1" << endl;

		//Check if hit is too small to get rid of overhang
		if((uint32_t) overhang >= hit->length){
			//Testing
			// cout << "2 Option 1" << endl;

			//Make sure that we do not consider this hit any further//TODO Maybe there is a nicer way to do so! Of course, we could move the whole seed on all successors, but most likely this will not lead to anything
			hit->score = 0;
			return;
		}

		//Testing
		//cout << "overhang:" << overhang << endl;
		// cout << "2 Option 2" << endl;

		//Iteratively retreat the extension
		while(overhang > 0){
			//Get unitig sequence
			string uSeq = hit->rUnitig.mappedSequenceToString();

			//Testing
			//cout << "Sequence is " << uSeq << endl << "uoff:" << hit->rSeedUoff << endl << "qoff:" << hit->rSeedQoff << endl << "Hit length:" << hit->length << endl;

			//Compare sequences and reduce the score accordingly
			hit->score -= compUScore(uSeq[hit->rSeedUoff], q[hit->rSeedQoff]);
			//Decrease hit's length
			--hit->length;
			//Decrease offsets
			--hit->rSeedQoff;
			--hit->rSeedUoff;
			//Decrease overhang	
			--overhang;
		}
	} else{
		//Testing
		// cout << "1 Option 2" << endl;
	}

	//Calculate how much the left border has to be adjusted//TODO Do we actually have to care for this case?!
	// if((overhang = hit->lSeedUoff - (hit->lUnitig.size - k)) > 0){
	// 	//Switch back to a previous unitig (doesn't matter which, because we are in the overlap)
	// 	hit->lUnitig = *(hit->lUnitig.getSuccessors().begin());
	// 	//Adjust unitig offset
	// 	hit->lSeedUoff = overhang - 1;
	// }	
}