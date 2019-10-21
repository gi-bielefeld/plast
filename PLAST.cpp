#include <cmath>
#include <chrono>
#include <list>
#include <sys/time.h>   
#include <sys/resource.h>

#include <bifrost/CompactedDBG.hpp>
#include <bifrost/ColoredCDBG.hpp>

//Testing
//static bool report = false;

//TODO: All this includes have to be checked for necessarity to be in this file!
#include "Seed.h"
//#include "Seedlist.h"
//#include "Seedlist.cpp"
// #include "Index.h"
#include "Index.cpp"
#include "Smer.h"
//#include "Smer.cpp"
#include "Extension.h"
//#include "Extension.cpp"
//#include "Hit.h"
//#include "GappedAlignment.h"
//#include "GappedAlignment.cpp"
// #include "Hit.cpp"
#include "IO.cpp"
//#include "Search.h"
#include "Search.cpp"
#include "Graph.h"
#include "Graph.cpp"
#include "Sequence.h"
#include "Statistics.h"

//Buffer size used to write binary files
//const uint8_t bufsize = 2 * sizeof(uint32_t);

//This function seems to be not used anymore
int32_t compSeqs(string &unitigSeq, const string &q, uint32_t &startInUni, uint32_t &startInQ, uint32_t &bpToComp){
	int32_t score = 0;

	for(uint32_t i = 0; i < bpToComp; ++i){
		score += compUScore(unitigSeq[startInUni + i], q[startInQ + i]);
	}

	return score;
}

/*
NOTE:
	-We might miss a seed in the beginning of a unitig's sequence if its reference sequence is the reverse complement of our query and is predecessor's sequence is not, but I guess we can ignore this easily without to miss significant findings
*/

//TODO Would it actually make sense to save k as a static variable in terms of run time?
int main(int argc, char **argv){//TODO This doesn't work yet for the reverse complement!//TODO Try out Guillaumes idea of the color preprocessing per unitig!
	//Staff we need to measure run times
	auto startTime = std::chrono::system_clock::now();
	auto endTime = std::chrono::system_clock::now();
	std::chrono::duration<double> tDiff;

	//A flag indicating whether we need to measure run times
	bool calcRT = RUNTIME_FLAG_DEFAULT;
	//A flag indicating whether color coverage should be outputted for each alignment
	bool repCols = REPORT_COLORS_FLAG_DEFAULT;
	//A flag indicating whether we have to build a graph and an index before actually searching for the query
	int16_t prep = PREPROS_FLAG_DEFAULT;
	//Determine the X-drop parameter
	int16_t X = DEFAULT_X;
	//Determine number of results we want to get in the end
	uint16_t nRes = DEFAULT_NB_RES;//TODO: At some point this should be added to the parameters!
	//Setting k and g
	int32_t kMerLength = DEFAULT_GRAPH_K;
	int32_t miniLength = DEFAULT_GRAPH_G;//TODO This should become an input parameter too!
	//Number of threads to be used for graph building
	int32_t nb_threads = DEFAULT_NB_THREADS;
	//Minimum length of a seed or s of an s-mer
	int32_t minSeedLength = DEFAULT_S;	//to test: 3,5,7,11,15
							/*NOTE: According to Faraz 12 might be a good default parameter (68% of all possible 12-mers occure in the human genome*/
	//Search quorum
	uint32_t quorum = DEFAULT_QUORUM;
	//Query counter
	uint32_t qCounter = 0;
	//Lambda needed for alignment statistic
	double lambda = DEFAULT_LAMBDA;
	//C needed for alignment statistic
	double C = DEFAULT_C;
	//E-value threshold for a hit to be considered
	double eBound = DEFAULT_E_LIMIT;
	//Query file name
	string qFile;
	//List of queries we are searching for
	vector<string> qList;
	//The strand we want to consider during search
	SrchStrd strand = SEARCH_STRAND_DEFAULT;
	//The file prefices of the graph and index files
	string graphFilePref;
	//File name of search color set
	string sColFile;
	//Building options for graph building
	CCDBG_Build_opt bOpt;

	//Check whether input is valid
	/*if(argc < 2){
		displayProgDescr();
		return 1;
	} else{

	}*/

	//Testing
	//cout << "DEFAULT_NB_RES:" << DEFAULT_NB_RES << endl;
	//cout << "Next we are going to parse the arguments" << endl;

	//Parse arguments
	if(!parseArgs(argc, argv, prep, graphFilePref, minSeedLength, kMerLength, miniLength, bOpt, nb_threads, qFile, sColFile, quorum, strand, repCols, X, nRes, lambda, C, eBound)){
		//Display help message
		dispHelp();
		return 1;
	}

	//Testing
	// cout << "Arguments parsed" << "Values are: preprocessing flag: " << prep << " graph file prefix: " << graphFilePref << " minimum seed length: " << minSeedLength << " k: " << kMerLength << " g: " << miniLength << " reference input files: " << endl;
	// for(vector<string>::iterator v = bOpt.filename_ref_in.begin(); v != bOpt.filename_ref_in.end(); ++v) cout << *v << endl;
	// cout << "raw input files:" << endl;
	// for(vector<string>::iterator v = bOpt.filename_seq_in.begin(); v != bOpt.filename_seq_in.end(); ++v) cout << *v << endl;
	// cout << "number of threads: " << nb_threads << " query file: " << qFile << " search color file name: " << sColFile << " quorum: " << quorum << " strand: " << strand << " calcRT: " << (calcRT ? "Yes" : "No") << " X: " << X << " number of results: " << nRes << " lambda: " << lambda << " C: " << C << " e-value threshold: " << eBound << " output colors: " << (repCols ? "Yes" : "No") << endl;
	// return 0;

	const uint32_t profileSize = (uint32_t) pow(SIGMAR, minSeedLength);

	uint32_t *qProfile = (uint32_t*) calloc(profileSize, sizeof(uint32_t));//TODO: Check how performance changes using calloc! Is it better to iterate over the profile separately?

	size_t numSmers = 0;
	struct s_mer_pos *posArray;
	UnitigColorMap<seedlist> *uArr;
	list<pair<string, size_t>> searchColors;

	//Variables needed for the seed detection
	uint32_t nbColors;

	//Initialize the graph
	ColoredCDBG<seedlist> cdbg(kMerLength, miniLength);

	//Check whether a graph has to be built
	if(prep > 0){
		//Build the graph
		genGraph(cdbg, bOpt);

		//Measure and output current runtime if demanded
		if(calcRT){
			//Get current time
			endTime = std::chrono::system_clock::now();
			//Calculate the time difference
			tDiff = endTime - startTime;
			//Output the measured time
			cout << "Getting the graph took " << tDiff.count() << " s" << endl;
			//Update start time for next part
			startTime = std::chrono::system_clock::now();
		}

		//Save the graph
		cdbg.write(graphFilePref, 1, true);
	} else{
		//Store the file name
		//string fName(argv[7]);

		//Testing
		//Get current time
		// endTime = std::chrono::system_clock::now();
		// //Calculate the time difference
		// std::chrono::duration<double> tDiff = endTime - startTime;
		// //Output the measured time
		// cout << "Parsing input and initializing some variables took " << tDiff.count() << " s" << endl;
		// //Update start time for next part
		// startTime = std::chrono::system_clock::now();

		//Load graph
		if(!cdbg.read((graphFilePref + GFA_FILE_ENDING).c_str(), (graphFilePref + COLOR_FILE_ENDING).c_str(), true)){
			cerr << "ERROR: Graph could not be loaded" << endl;
			exit(EXIT_FAILURE);
		}

		//Testing
		//Get current time
		// endTime = std::chrono::system_clock::now();
		// //Calculate the time difference
		// tDiff = endTime - startTime;
		// //Output the measured time
		// cout << "Reading the graph took " << tDiff.count() << " s" << endl;
		// //Update start time for next part
		// startTime = std::chrono::system_clock::now();
	}

	//Testing
	//cout << "Next step is to parse the index" << endl;

	//Check whether an index has to be built
	if(prep > PREPROS_FLAG_DEFAULT){
		/*Building the q-gram index*/
		cout << "Building q-gram index" << endl;

		buildIndex(cdbg, minSeedLength, qProfile, numSmers, profileSize, posArray, uArr);

		//Measure and output current runtime if demanded
		if(calcRT){
			//Get current time
			endTime = std::chrono::system_clock::now();
			//Calculate the time difference
			std::chrono::duration<double> tDiff = endTime - startTime;
			//Output the measured time
			cout << "Index building took " << tDiff.count() << " s" << endl;
			//Update start time for next part
			startTime = std::chrono::system_clock::now();
		}

		//Testing
		// cout << "posArray:" << endl;
		// for(size_t i = 0; i < numSmers; ++i){
		// 	cout << "Unitig id:" << posArray[i].unitig_id << " offset:" << posArray[i].offset << endl;
		// }
		// exit(0);

		//Save indexes
		saveIndexBin((graphFilePref + INDEX_FILE_ENDING).c_str(), profileSize, qProfile, cdbg.size(), uArr, posArray, numSmers, kMerLength);

		return 0;
	} else {
		//Load indexes
		loadIndexesBin((graphFilePref + INDEX_FILE_ENDING).c_str(), profileSize, qProfile, cdbg, uArr, numSmers, posArray);

		//Testing
		// cout << "posArray:" << endl;
		// for(size_t i = 0; i < numSmers; ++i){
		// 	cout << "Unitig id:" << posArray[i].unitig_id << " offset:" << posArray[i].offset << endl;
		// }
		// exit(0);
		//Get current time
		// endTime = std::chrono::system_clock::now();
		// //Calculate the time difference
		// std::chrono::duration<double> tDiff = endTime - startTime;
		// //Output the measured time
		// cout << "Loading indices took " << tDiff.count() << " s" << endl;
		// //Update start time for next part
		// startTime = std::chrono::system_clock::now();
		// return 0;

		//Measure and output current runtime if demanded
		if(calcRT){
			//Get current time
			endTime = std::chrono::system_clock::now();
			//Calculate the time difference
			std::chrono::duration<double> tDiff = endTime - startTime;
			//Output the measured time
			cout << "Loading took " << tDiff.count() << " s" << endl;
			//Update start time for next part
			startTime = std::chrono::system_clock::now();
		}
	}

	//Check whether set quorum exceeds the number of colors in graph
	if(quorum > cdbg.getNbColors()){
		cerr << "WARNING: Quorum exceeds number of colors in the graph. Set to that number instead." << endl;
		quorum = cdbg.getNbColors();
	}

	//Check whether a search color set is given
	if(!sColFile.empty()){
		//Read in color names to search for
		searchColors = loadSearchColors(sColFile.c_str(), nbColors);
		//Map color ids to color names
		mapColorIds(searchColors, cdbg);

		//Check whether set quorum exceeds the number of colors in search set
		if(quorum > nbColors){
			cerr << "WARNING: Quorum exceeds number of colors in search set. Set to that number instead." << endl;
			quorum = nbColors;
		}
	}

	//Testing
	// UnitigColorMap<seedlist> sUni = cdbg.find(Kmer("CTCACCTTGGTTGAAAAATGCTAAAAGTGCC"));
	// cout << "Unitig: " << sUni.mappedSequenceToString() << endl;
	// for(UnitigColors::const_iterator i = sUni.getData()->getUnitigColors(sUni)->begin(sUni); i != sUni.getData()->getUnitigColors(sUni)->end(); ++i){
	// 	cout << "Position:" << (*i).first << " Color:";
	// 	for(list<pair<string, size_t>>::const_iterator cname = searchColors.begin(); cname != searchColors.end(); ++cname){
	// 		if(cname->second == (*i).second) cout << cname->first << endl;
	// 	}
	// }
	// return 0;

	//Load query sequences
	loadQueries(qFile, qList);

	//Search each query
	for(vector<string>::const_iterator q = qList.begin(); q != qList.end(); ++q){
		//Output which query we are working on
		cout << "Query " << ++qCounter << ":" << endl;
		//Search for the current query
		searchQuery(cdbg, kMerLength, minSeedLength, numSmers, quorum, profileSize, qProfile, *q, strand, uArr, posArray, searchColors, X, calcRT, nRes, lambda, C, eBound, repCols);

		//Testing
		// UnitigColorMap<seedlist> sUni = cdbg.find(Kmer("CTCACCTTGGTTGAAAAATGCTAAAAGTGCC"));
		// cout << "Unitig: " << sUni.referenceUnitigToString() << endl;
		// cout << "lBrd: " << sUni.getData()->getData(sUni)->getlBrd() << " rBrd: " << sUni.getData()->getData(sUni)->getrBrd() << endl;
		// return 0;

		//Measure and output current runtime if demanded
		if(calcRT){
			//Get current time
			endTime = std::chrono::system_clock::now();
			//Calculate the time difference
			std::chrono::duration<double> tDiff = endTime - startTime;
			//Output the measured time
			cout << "Freeing the array took " << tDiff.count() << " s" << endl;
			//Update start time for next part
			startTime = std::chrono::system_clock::now();
		}
	}

	// 	while(next != NULL){
	// 		cout << "Seedlist is not NULL" << endl;

	// 		next = next->nextSeed;
	// 	}

	// 	while(next != NULL){
	// 		cout << "Seedlist is not NULL" << endl;

	// //Detect seeds on the reverse complementary query
	// detectSeeds(kMerLength, minSeedLength, numSmers, quorum, profileSize, qProfile, q, uArr, posArray, hitArr, true);
	// if(calcRT){
	// 	//Get current time
	// 	endTime = std::chrono::system_clock::now();
	// 	//Calculate the time difference
	// 	std::chrono::duration<double> tDiff = endTime - startTime;
	// 	//Output the measured time
	// 	cout << "Seed detection took " << tDiff.count() << " s" << endl;
	// 	//Update start time for next part
	// 	startTime = std::chrono::system_clock::now();
	// }

	// for(uint32_t i = 0; i < q.length(); ++i){
	// 	hitArr[i].length = 0;
	// }
	// extendSeeds(cdbg, revQ, X, hitArr, quorum);
	// if(calcRT){
	// 	//Get current time
	// 	endTime = std::chrono::system_clock::now();
	// 	//Calculate the time difference
	// 	std::chrono::duration<double> tDiff = endTime - startTime;
	// 	//Output the measured time
	// 	cout << "Seed extension took " << tDiff.count() << " s" << endl;
	// 	//Update start time for next part
	// 	startTime = std::chrono::system_clock::now();
	// }
	// for(uint32_t i = 0; i < q.length() - minSeedLength + 1; ++i){
	// 	struct hit *hitList = NULL;

	// 	if(hitArr[i].length == 0){
	// 		hitList = NULL;
	// 	} else{
	// 		hitList	= &hitArr[i];
	// 	}

	// 	//Go through all hits of a certain array position
	// 	while(hitList != NULL){
	// 		//Testing
	// 		//cout << "Bis hierher kommen wir schon mal" << endl;

	// 		//Check if we can still just add new results to the result list or need to replace worse existing ones
	// 		if(nRes == 0){
	// 			//Try to replace a worse result
	// 			replWorseRes(resList, hitList);
	// 		} else{
	// 			//Testing
	// 			//cout << "Hierein kommen wir" << endl;

	// 			//Insert a new result
	// 			insRes(resList, hitList);
	// 			--nRes;

	// 			//Testing
	// 			//cout << "Length of resList after insert:" << resList.size() << endl;
	// 		}

	// 		hitList = hitList->nextHit;
	// 	}
	// }
	// calcGappedAlignment(resList, revQ, X, calcRT, quorum);
	// if(calcRT){
	// 	//Get current time
	// 	endTime = std::chrono::system_clock::now();
	// 	//Calculate the time difference
	// 	std::chrono::duration<double> tDiff = endTime - startTime;
	// 	//Output the measured time
	// 	cout << "Gapped extension took " << tDiff.count() << " s" << endl;
	// 	//Update start time for next part
	// 	startTime = std::chrono::system_clock::now();
	// }

	//saveIndex(strcat(argv[5], ".idx"), qProfile, profileSize, posArray, numSeeds);

	return 0;
}
