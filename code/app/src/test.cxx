//-----------------------------------------------------------------------------
// Author    : Fearghal O'Donncha, feardonn@ie.ibm.com
//             Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#include <iostream>
#include <random>
//#include "../include/amdados_utils.h"
//#include "../include/matrix.h"

int main(int /*argc*/, char ** /*argv*/)
{/*
	using namespace ::amdados;

	const int    DIM = 3;

	Vector v(DIM);

	for(int i = 0; i < DIM; ++i) {
		std::cerr << v(i) << " ";
	}
	std::cerr << std::endl;
	MakeRandom(v, 'n');
	for(int i = 0; i < DIM; ++i) {
		std::cerr << v(i) << " ";
	}
	std::cerr << std::endl;

	exit(42);*/
//  	  std::mt19937_64 gen(2063);

	std::mersenne_twister_engine<uint64_t, 64, 312, 156, 31, 2842040809, 29, 1431655765, 17, 3987079168, 37, 0, 43, 1284865837 > gen{((uint64_t)2063)};
	        std::normal_distribution<double> distrib;
		auto test = distrib(gen);
		std::cerr << test << std::endl;
		exit(42);
/*	        for (int i = 0; i < size; ++i) { 
			v(i) = distrib(gen); 
		}*/


}

