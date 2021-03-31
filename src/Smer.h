#ifndef SMER_HPP
#define SMER_HPP

//Default minimal seed length; should always be positive
#define DEFAULT_S 11

//Data structure to store a s-mer's position
struct S_mer_pos{
	//Id that links an s-mer's position to its unitig (Since there are no pointers to unitigs and unitigs have to be stored instead, each unitig is only stored once inside a unitig array and positions are linked to each s-mer)
	uint32_t unitig_id;
	//Position within the unitig sequence at which the s-mer is located
	uint32_t offset;
};

#endif