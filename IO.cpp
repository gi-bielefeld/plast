#include <getopt.h>

#include "IO.h"

//This function parses the program parameters. Returns false if given arguments are not valid.
const bool parseArgs(int& nb_args, char** argList, int16_t& prepros, string& filePref, int32_t& s, int32_t& k, int32_t& g,  CCDBG_Build_opt &gOpt, int32_t& t, string& qFile, string& c, uint32_t& m, SrchStrd& strd, bool& r, int16_t& X, uint16_t &nRes, double &lambda, double &C, double &eValLim){
	int option_index = 0, a;

	//Check wheather arguments are given for anything at all
	if(nb_args < MIN_PARAM_NB) return false;

	static struct option long_options[] = {
        {"graph-prefix",   required_argument,  0, 'i'},
        {"seed-length",    optional_argument,  0, 's'},
        {"kmer-length",    optional_argument,  0, 'k'},
        {"min-length",     optional_argument,  0, 'g'},
        {"raw-input",      optional_argument,  0, 'S'},
        {"ref-input",      optional_argument,  0, 'R'},
        {"threads",        optional_argument,  0, 't'},
        {"query",          optional_argument,  0, 'q'},
        {"colors",         optional_argument,  0, 'c'},
        {"quorum",         optional_argument,  0, 'm'},
        {"X-dropoff",      optional_argument,  0, 'X'},
        {"max-results",    optional_argument,  0, 'n'},
        {"strand",         optional_argument,  0, 'd'},
        {"lambda",         optional_argument,  0, 'l'},
        {"stat-C",         optional_argument,  0, 'C'},
        {"e-value",        optional_argument,  0, 'e'},
        {"report-colors",  no_argument,        0, 'r'},
        {0,                0,                  0,  0 }
    };

    //Try to get the command
	if(!strcmp(argList[1], BUILD_COMMAND)){
		//Set preprocessing flag
		prepros = 0;
	} else if(!strcmp(argList[1], SEARCH_COMMAND)){
		//Set preprocessing flag
		prepros = -1;
	} else{
		cerr << "ERROR: Command unknown" << endl;
		return false;
	}

	//Parse all parameters given
	//while((a = getopt(nb_args, argList, OPTIONS)) != -1){
	while ((a = getopt_long(nb_args, argList, OPTIONS, long_options, &option_index)) != -1){
		//Assign parameter values
		switch(a){
			case 'i':
				//Testing
				//cout << "optarg:" << optarg << endl;

				filePref = optarg;
				break;
			case 'S':
				//Testing
				/*cout << "Why do we not go here?" << endl;*/
				//cout << "optind:" << optind << endl;

				for(--optind; (optind < nb_args) && (*argList[optind] != '-'); ++optind){
					// seqs.push_back(argList[optind]);

					gOpt.filename_seq_in.push_back(argList[optind]);
				}

				break;
			case 'R':
				for(--optind; (optind < nb_args) && (*argList[optind] != '-'); ++optind){
					// seqs.push_back(argList[optind]);

					gOpt.filename_ref_in.push_back(argList[optind]);
				}

				break;
			case 's':
				//Testing
				//cout << "Do we get here?" << endl;

				s = atoi(optarg);
				break;
			case 'k':
				k = atoi(optarg);
				break;
			case 'g':
				g = atoi(optarg);
				break;
			case 't':
				t = atoi(optarg);
				break;
			case 'q':
				qFile = optarg;
				break;
			case 'c':
				c = optarg;
				break;
			case 'm':
				//A quorum has to be positive
				if(atoi(optarg) <= 0){
					cerr << "ERROR: Quorum parameter has to be a positive number" << endl;
					return false;
				} 
				
				m = atoi(optarg);
				break;
			case 'X':
				//Testing
				// cerr << atoi(optarg) << endl;

				//Check if X is applicable
				if(atoi(optarg) < 0 || atoi(optarg) > INT16_MAX){
					//Testing
					// if(atoi(optarg) < 0){
					// 	// cout << "6 Option 1" << endl;
					// } else{
					// 	cout << "6 Option 2" << endl;
					// }
					// if(atoi(optarg) > INT16_MAX){
					// 	cout << "7 Option 1" << endl;
					// } else{
					// 	// cout << "7 Option 2" << endl;
					// }

					cerr << "ERROR: X-dropoff value not applicable" << endl;
					return false;
				}

				X = atoi(optarg);
				break;
			case 'n':
				//Check if value is valid
				if(atoi(optarg) <= 0){
					//Testing
					// cout << "8 Option 1" << endl;

					cerr << "ERROR: Maximum number of alignments should be a positive number" << endl;
					return false;
				} else if(atoi(optarg) > UINT16_MAX){
					//Testing
					// cout << "9 Option 1" << endl;

					cerr << "Maximum number of alignments set to maximum" << endl;
					nRes = UINT16_MAX;
				} else{
					//Testing
					// cout << "9 Option 2" << endl;

					nRes = atoi(optarg);
				}

				//Testing
				// cout << "8 Option 2" << endl;

				break;
			case 'r':
				r = true;
				break;
			case 'd':
				//Testing
				// cout << optarg << endl;

				if(*optarg == PLUS_STRAND){
					strd = Plus;

					//Testing
					// cout << "10 Option 1" << endl;
					// exit(0);
				} else if(*optarg == MINUS_STRAND){
					//Testing
					// cout << "10 Option 2" << endl;

					strd = Minus;
				} else{
					//Testing
					// cout << "10 Option 3" << endl;
					
					cerr << "Unrecognized strand option" <<  endl << "Default is used" << endl;
				}

				break;
			case 'l':
				//Try to read lambda value
				lambda = strtod(optarg, NULL);

				//Testing
				// cout << "11" << endl;

				//Check if the given value is out of range
				if(errno == ERANGE || lambda == 0.0){
					//Testing
					// cout << "12 Option 2" << endl;

					cerr << "ERROR: Invalid lambda value given" << endl;
					return false;
				}

				//Testing
				// cout << "12 Option 1" << endl;

				break;
			case 'C':
				//Try to read C value
				C = strtod(optarg, NULL);

				//Testing
				// cout << "13" << endl;

				//Check if the given value is out of range
				if(errno == ERANGE || C == 0.0){
					//Testing
					// cout << "14 Option 2" << endl;

					cerr << "ERROR: Invalid C value given" << endl;
					return false;
				}

				//Testing
				// cout << "14 Option 1" << endl;

				break;
			case 'e':
				//Try to read e-value threshold
				eValLim = strtod(optarg, NULL);

				//Testing
				// cout << "15" << endl;

				//Check if the given value is out of range
				if(errno == ERANGE || eValLim == 0.0){
					//Testing
					cout << "16 Option 2" << endl;

					cerr << "ERROR: Invalid e-value threshold given" << endl;
					return false;
				}

				//Testing
				// cout << "16 Option 1" << endl;

				break;
			default:
				//Testing
				//cout << "Default: optarg:" << optarg << endl;

				break;
		}
	}

	//Testing
	// cout << "a: " << a << endl;

	//If we have received any input sequences we need to build a graph
	if(!gOpt.filename_seq_in.empty() || !gOpt.filename_ref_in.empty()) ++prepros;

	//Without a file prefix none command can be executed
	if(!strcmp(filePref.c_str(), "")) return false;

	//The search command additionally needs a query sequence
	if(prepros < 0 && !strcmp(qFile.c_str(), "")) return false;

	return true;
}

//This function reads in a file in which colors are stored the search will be based on
const list<pair<string, size_t>> loadSearchColors(const char* filename, uint32_t& nbCols){
	string line;
	pair<string, size_t> color("", -1);//Color id is set to -1 here to avoid false results if color from color search set is not present in the graph//TODO: Avoid this by mapping the colors directly within this function!
	list<pair<string, size_t>> colorlist;

	//Open the file
	ifstream fStr(filename);

	//Check if the file is open
	if(!fStr.is_open()){
		cerr << "ERROR: Color set file could not be opened" << endl;
		exit(EXIT_FAILURE);
	}

	//Reset counter
	nbCols = 0;
	//Read in lines of color set file
	while(getline(fStr, line)){
		//Store color name
		color.first = line;
		//Store each color
		colorlist.push_back(color);
		//Increment counter
		++nbCols;
	}

	//Close the file
	fStr.close();

	return colorlist;
}

//This function loads query sequences from a file and stores it in a list
void loadQueries(const string &filename, vector<string> &qList){
	string line;

	//Open the file
	ifstream fStr(filename);

	//Check if the file is open
	if(!fStr.is_open()){
		cerr << "ERROR: Query sequence file could not be opened" << endl;
		exit(EXIT_FAILURE);
	}

	//Read in queries (one per line)
	while(getline(fStr, line)) qList.push_back(line);
	//Close file
	fStr.close();
}

//This function outputs a hit's alignment
void repAlgn(const hit *res){
	bool largeNb = false;
	uint32_t i = 0, algnCols, gaps = 0, posQ, posAlgn = 0;

	//Output some general info
	cout << "Score: " << res->score << "\tLength: " << res->gAlgn.aSeqQ.length() << "\tE-value: " << res->eval << endl;
	//Initial position in query (we start count from 1 here)
	posQ = res->offQ + 1;

	//Testing
	// if(res->score == 4){
	// 	cerr << "Start outputting alignment having a score of 4" << endl;
	// 	cerr << "origUni: " << res->origUni.mappedSequenceToString() << endl;
	// // 	if(res->origUni.mappedSequenceToString() == "TTTGTTTTCAATTGCTGATGAATGGGGTATGAGTAAACTGAG"){
	// // 		cout << "We have found the interesting unitig and we have to move the starting offset" << endl;
	// // 		exit(0);
	// // 	}
	// 	// exit(0);
	// }

	//Go through the alignment
	while(posAlgn < res->gAlgn.aSeqQ.length()){
		//Check if there is enough alignment to fill another full line
		if(posAlgn + BASES_PER_LINE < res->gAlgn.aSeqQ.length()){
			//We output the maximum amount of bases
			algnCols = BASES_PER_LINE;
		} else{
			//We output only what is left
			algnCols = res->gAlgn.aSeqQ.length() - posAlgn;
		}

		//Upper sequence is the query
		cout << "Query:    ";
		//Output position in query
		cout << posQ << "\t";//We want to start counting from 1 here so we add +1

		if(posQ > 9999999) largeNb = true;

		//Output next query section
		for(i = posAlgn; i < posAlgn + algnCols; ++i){
			//Count non-gap characters
			if(res->gAlgn.aSeqQ[i] == GAP) ++gaps;

			//Output character
			cout << res->gAlgn.aSeqQ[i];
		}

		//Output position in q at the end of this line
		cout << "\t" << (posQ += algnCols - gaps) - 1 << endl;
		//Output match line
		cout << "          " << (largeNb ? "\t\t" : "\t");
		//Reset counter
		i = posAlgn;

		do{
			cout << (res->gAlgn.aSeqQ[i] == res->gAlgn.aSeqG[i] ? "|" : " ");
		} while(++i != posAlgn + algnCols);

		//End the match line
		cout << endl;
		//Lower sequence is the graph
		cout << "Graph:    " << 0 << (largeNb ? "\t\t" : "\t");
		//Output next graph section
		cout << res->gAlgn.aSeqG.substr(posAlgn, algnCols) << "\t0" << endl << endl;//"\t" << res->lSeedUoff + i + algnCols - 1 << endl << endl;
		//Increment position in alignment
		posAlgn += algnCols;
		gaps = 0;
	}
}

//This function outputs the color sets of given result
void outpColSets(ColoredCDBG<seedlist> &cdbg, const hit *res){
	bool identical = false;
	//char *kmerSeq = (char*) malloc(res->origUni.getGraph()->getK());
	string kmerSeq = string(res->origUni.getGraph()->getK(), 'A');

	/*Some additional stuff that Roland needs for an experiment*/
	uint32_t currQpos = res->offQ;
	/*Some additional stuff that Roland needs for an experiment*/

	int32_t kmerSeqPos = 0;
	string color;
	vector<string> colors;
	vector<string>::iterator colIt;
	Kmer currK;
	UnitigColors uniCols;
	UnitigColorMap<seedlist> uni;
	ColSet curColSet = ColSet(0, colors);

	//Testing
	// cout << "gAlgn.aSeqG: " << res->gAlgn.aSeqG << endl;
	// cout << "Seed was length: " << res->length << " offset u: " << res->offU << " q: " << res->offQ << " unitig: " << res->origUni.mappedSequenceToString() << endl;

	//Go through graph alignment sequence
	for(uint32_t i = 0; i < res->gAlgn.aSeqG.length(); ++i){
		//Check if current character in alignment is not a gap
		if(res->gAlgn.aSeqG[i] != GAP){
			//Add base to current k-mer sequence
			kmerSeq[kmerSeqPos] = res->gAlgn.aSeqG[i];
			//Update position in kmer sequence
			++kmerSeqPos;
		}

		/*Some additional stuff that Roland needs for an experiment*/
		//Check if current character in query's alignment sequence is not a gap
		if(res->gAlgn.aSeqQ[i] != GAP) ++currQpos;
		/*Some additional stuff that Roland needs for an experiment*/

		//Testing
		// cout << "i: " << i << endl;
		// cout << "kmerSeq: " << kmerSeq << endl;

		//Check if our k-mer sequence is full
		if(kmerSeqPos == res->origUni.getGraph()->getK()){
			//Testing
			// cout << "We have just found a new k-mer to check" << endl;
			// cout << "kmerSeq: " << kmerSeq << endl;

			//Initialize k-mer
			currK = Kmer(kmerSeq.c_str());
			//Look up k-mer
			uni = res->origUni.getGraph()->find(currK);

			//Testing
			// cout << "Unitig is " << (uni.isEmpty ? "" : "not ") << "empty" << endl;
			// UnitigColorMap<seedlist> bla = uni.getGraph()->find(Kmer("TCTTTACGGCGAAGTTCAGCGCCCTCATAGCC"));
			// cout << "This unitig can " << (bla.isEmpty ? "not " : "") << "be found in the graph" << endl;

			//Get k-mer's color set
			uniCols = *uni.getData()->getUnitigColors(uni);
			//Reset color set
			colors = vector<string>();

			//Testing
			// cout << "Extract colors" << endl;

			//Iterate over k-mer's colors
			for(UnitigColors::const_iterator c = uniCols.begin(uni); c != uniCols.end(); ++c){
				//Get current ID's color name
				color = cdbg.getColorName(c.getColorID());

				//If the color set is empty we have nothing to compare against
				if(colors.empty()){
					colors.push_back(color);
				} else if(colors.back() != color){//Check if we have just seen the current color already
					colors.push_back(color);
				}
			}

			//Testing
			// cout << "We have extracted the color set" << endl;

			//Check if we have dealed with a different k-mer before
			if(!curColSet.colNames.empty()){
				colIt = curColSet.colNames.begin();

				//Iterate over current k-mer's color set
				for(vector<string>::iterator j = colors.begin(); j != colors.end(); ++j){
					//Color sets of different sizes cannot be identical
					if(curColSet.colNames.size() != colors.size()){
						identical = false;
						break;
					}

					//Compare the current colors
					if(*colIt == *j){
						identical = true;
					} else{
						identical = false;
						break;
					}

					++colIt;
				}
			}

			// //Get k-mer's color set
			// colIDs = vector<size_t>();
			// for(UnitigColors::iterator colIt = uni.getData()->getUnitigColors(uni)->begin(uni); colIt != uni.getData()->getUnitigColors(uni)->end(); ++colIt){
			// 	if(colIDs.empty()){
			// 		colIDs.push_back();
			// 	} else if(*j != colIDs.back()){
			// 		colIDs.push_back();
			// 	}
			// }

			// //Compare color sets
			// l = 0;

			// for(vector<size_t>::iterator j = colIDs.begin(); j != colIDs.end(); ++j){
			// 	if(l < curColSet.)
			// }

			// if(*j == colIt->getColorID()){
			// 		identical = true;
			// 	} else{
			// 		identical = false;
			// 		break;
			// 	}

			//Testing
			// cout << "Color set comparison done" << endl;

			//Check outcome of comparison
			if(!identical && !curColSet.colNames.empty()){
				//Output color set
				cout << "Color set ending at alignment position " << curColSet.endPos;

				/*Some additional stuff that Roland needs for an experiment*/
				cout << " (position " << currQpos << " in the query sequence)";
				/*Some additional stuff that Roland needs for an experiment*/

				for(vector<string>::iterator name = curColSet.colNames.begin(); name != curColSet.colNames.end(); ++name){
					cout << " " << *name;
				}

				cout << endl;
			}

			//Testing
			// cout << "kmerSeq:" << kmerSeq << endl;
			// cout << "We have outputted a color set and get here" << endl;

			//Update color set (doesn't do anything if color sets are identical anyways)
			curColSet.colNames = colors;
			//Decrease k-mer sequence length
			--kmerSeqPos;
			//Shift kmer sequence
			for(int32_t j = 0; j < kmerSeqPos; ++j){
				kmerSeq[j] = kmerSeq[j + 1];
			}

			// kmerSeq[kmerSeqPos] <<= sizeof(char);

			//Testing
			// cout << "kmerSeq after shift:" << kmerSeq << endl;
		}

		//Update color set end position
		++curColSet.endPos;
	}

	// //Free sequence
	// free(kmerSeq);
	
	//Output color set
	cout << "Color set ending at alignment position " << curColSet.endPos;

	/*Some additional stuff that Roland needs for an experiment*/
	cout << " (position " << currQpos << " in the query sequence)";
	/*Some additional stuff that Roland needs for an experiment*/

	//Check if color set to be outputted is empty which happens if the graph sequence of our alignment is shorter than k
	if(curColSet.colNames.empty()){
		//Get the unitig the alignment lies on
		uni = res->origUni;
		//Resize unitig so that its color set represents the alignments color set
		uni.dist = min(res->offU, (uint32_t) uni.len - 1);
		//Get color set
		uniCols = *uni.getData()->getUnitigColors(uni);//TODO Put this into a function!
		colors = vector<string>();

		for(UnitigColors::const_iterator c = uniCols.begin(uni); c != uniCols.end(); ++c){
			color = cdbg.getColorName(c.getColorID());

			if(colors.empty()){
				colors.push_back(color);
			} else if(colors.back() != color){
				colors.push_back(color);
			}
		}

		curColSet.colNames = colors;
	}

	for(vector<string>::iterator name = curColSet.colNames.begin(); name != curColSet.colNames.end(); ++name){
		cout << " " << *name;
	}

	cout << endl;
}