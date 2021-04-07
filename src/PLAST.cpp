#include <cmath>
#include <chrono>
#include <list>
#include <sys/time.h>   
#include <sys/resource.h>

#include <bifrost/CompactedDBG.hpp>
#include <bifrost/ColoredCDBG.hpp>

#include "Seed.h"
#include "Index.cpp"
#include "Smer.h"
#include "Extension.h"
#include "IO.cpp"
#include "Search.cpp"
#include "Graph.h"
#include "Graph.cpp"
#include "Sequence.h"
#include "Statistics.h"

int main(int argc, char **argv){
	//Staff we need to measure run times
	auto startTime = std::chrono::system_clock::now();
	auto endTime = std::chrono::system_clock::now();
	std::chrono::duration<double> tDiff;

	//A flag indicating whether we need to measure run times
	bool calcRT = RUNTIME_FLAG_DEFAULT;
	//A flag indicating whether color coverage should be outputted for each alignment
	bool repCols = REPORT_COLORS_FLAG_DEFAULT;
	//A flag indicating whether this run is only a simulation run
	bool isSim = SIM_RUN_FLAG_DEFAULT;
	//A flag indicating if quorum information should be added to index during construction
	bool advIdx = ADVANCED_INDEX_FLAG_DEFAULT;
	//A flag indicating whether we have to build a graph and an index before actually searching for the query
	int16_t prep = PREPROS_FLAG_DEFAULT;
	//Determine the X-drop parameter
	int16_t X = DEFAULT_X;
	//Mismatch score
	int16_t mmscr = USCORE_MISMATCH;
	//Determine number of results we want to get in the end
	uint16_t nRes = DEFAULT_NB_RES;
	//Match score
	uint16_t mscr = USCORE_MATCH;
	//Gap open score
	int32_t goscr = GAP_OPEN_SCORE;
	//Gap extension score
	int32_t gescr = GAP_EXT_SCORE;
	//Size of s-mer profile
	uint32_t profileSize;
	//s-mer profile
	uint32_t *qProfile;
	//Setting k and g
	int32_t kMerLength = DEFAULT_GRAPH_K;
	int32_t miniLength = DEFAULT_GRAPH_G;
	//Number of threads to be used for graph building
	int32_t nb_threads = DEFAULT_NB_THREADS;
	//Minimum length of a seed or s of an s-mer
	int32_t minSeedLength = DEFAULT_S;
	//Search quorum
	uint32_t quorum = DEFAULT_QUORUM;
	//Query counter
	uint32_t qCounter = 0;
	//Number of colors in potential search set
	uint32_t nbColors;
	//Number of seeds in the graph
	size_t numSmers = 0;
	//Lambda needed for alignment statistic
	double lambda = DEFAULT_LAMBDA;
	//Gapped lambda needed for alignment statistic
	double lambdaGap = DEFAULT_LAMBDA_G;
	//C needed for alignment statistic
	double C = DEFAULT_C;
	//Gapped C needed for alignment statistic
	double Cgap = DEFAULT_C_G;
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
	//Color id and name mapping
	list<pair<string, size_t>> searchColors;
	//Pos array
	struct S_mer_pos *posArray;
	//Building options for graph building
	CCDBG_Build_opt bOpt;
	//Unitig array
	UnitigColorMap<UnitigInfo> *uArr;

	//Parse arguments
	if(!parseArgs(argc, argv, prep, graphFilePref, minSeedLength, kMerLength, miniLength, bOpt, nb_threads, qFile, sColFile, quorum, strand, repCols, mscr, mmscr, X, goscr, gescr, nRes, lambda, lambdaGap, C, Cgap, eBound, isSim, advIdx)){
		//Display help message
		dispHelp();
		return 1;
	}

	//Initialize the graph
	ColoredCDBG<UnitigInfo> cdbg(kMerLength, miniLength);

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
		cdbg.write(graphFilePref, nb_threads, true);
	} else{
		//Load graph
		if(!cdbg.read(graphFilePref + GFA_FILE_ENDING, graphFilePref + COLOR_FILE_ENDING, nb_threads, false)){
			cerr << "ERROR: Graph could not be loaded" << endl;
			exit(EXIT_FAILURE);
		}
	}

	//Check whether an index has to be built
	if(prep > PREPROS_FLAG_DEFAULT){
		/*Building the q-gram index*/
		cout << "Building q-gram index" << endl;

		//Calculate profile size
		profileSize = (uint32_t) pow(SIGMAR, minSeedLength);
		//Initialize profile
		qProfile = (uint32_t*) calloc(profileSize, sizeof(uint32_t));
		//Build index
		buildIndex(cdbg, minSeedLength, qProfile, numSmers, profileSize, posArray, uArr);

		//Calculate quorum information if needed for the index
		if(advIdx) calcQrms(cdbg);

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

		//Save indexes
		saveIndexBin((graphFilePref + INDEX_FILE_ENDING).c_str(), minSeedLength, profileSize, qProfile, cdbg.size(), uArr, posArray, numSmers, kMerLength, advIdx);

		return 0;
	} else {
		//Load indexes
		loadIndexesBin((graphFilePref + INDEX_FILE_ENDING).c_str(), minSeedLength, profileSize, qProfile, cdbg, uArr, numSmers, posArray, advIdx);

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
		searchColors = loadSearchColors(sColFile.c_str(), cdbg.getColorNames());
		//Get number of colors in search set
		nbColors = searchColors.size();

		//Check whether set quorum exceeds the number of colors in search set
		if(quorum > nbColors){
			cerr << "WARNING: Quorum exceeds number of colors in search set. Set to that number instead." << endl;
			quorum = nbColors;
		}
	}

	//Load query sequences
	loadQueries(qFile, qList);

	//Search each query
	for(vector<string>::const_iterator q = qList.begin(); q != qList.end(); ++q){
		//Output which query we are working on
		cout << "Query " << ++qCounter << ":" << endl;
		//Search for the current query
		searchQuery(cdbg, kMerLength, minSeedLength, numSmers, quorum, profileSize, qProfile, *q, strand, uArr, posArray, searchColors, mscr, mmscr, X, goscr, gescr, calcRT, nRes, lambda, lambdaGap, C, Cgap, eBound, repCols, isSim, advIdx);

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

	//Free pos array
	free(posArray);
	//Free profile
	free(qProfile);
	//Free unitig array
	free(uArr);

	return 0;
}