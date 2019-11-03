#ifndef STATISTICS_HPP
#define STATISTICS_HPP

#define DEFAULT_LAMBDA 0.99
#define DEFAULT_C 21.68
#define DEFAULT_E_LIMIT 10
#define SIM_Q_LENGTH 1000

//This function calculates an e-value for a given score
inline double calcEVal(const uint32_t &score, const double &lambda, const double &C, const uint32_t &qLen){
	return ((double) qLen / SIM_Q_LENGTH) * exp(C) * exp(-lambda * score);
}

#endif