#include "FCorrelatorSAt.h"

#include <cmath>
//#include "parameters.h"

using namespace SAt;

/////////////////////////////////////////
// Correlator class
/////////////////////////////////////////


Correlator::Correlator(const unsigned int numcorrin, const unsigned int pin, const unsigned int min, unsigned int Nchain) {
	setsize(numcorrin, pin, min, Nchain);
}

Correlator::~Correlator() {

	if (numcorrelators == 0) return;

	delete[] shift;
	delete[] correlation;
	delete[] ncorrelation;
	delete[] accumulator;
	delete[] naccumulator;
	delete[] insertindex;

	delete[] t;
	delete[] f;
}

void Correlator::setsize(const unsigned int numcorrin, const unsigned int pin, const unsigned int min, unsigned int Nchain) {
	numcorrelators = numcorrin;
	p = pin;
	m = min;
	dmin = p / m;

	Nc = Nchain;
	double *w = new double[Nc];


	length = numcorrelators * p;

	shift = new double**[numcorrelators];
	correlation = new double*[numcorrelators];
	ncorrelation = new unsigned long int *[numcorrelators];
	accumulator = new double*[numcorrelators];
	naccumulator = new unsigned int[numcorrelators];
	insertindex = new unsigned int[numcorrelators];

	accval = new double[Nc];

	for (unsigned int j = 0; j < numcorrelators; ++j) {
		shift[j] = new double*[p];
		for (unsigned int i = 0; i < p; ++i) {
			shift[j][i] = new double[Nc];
		}

		accumulator[j] = new double[Nc];
		/* It can be optimized: Apart from correlator 0, correlation and ncorrelation arrays only use p/2 values */
		correlation[j] = new double[p];
		ncorrelation[j] = new unsigned long int[p];
	}

	t = new double[length];
	f = new double[length];
}

void Correlator::initialize() {

	for (unsigned int j = 0; j < numcorrelators; ++j) {
		for (unsigned int i = 0; i < p; ++i) {
			for (unsigned int k = 0; k < Nc; ++k) {
				//VSet(shift[j][i][k], -2E10, -2E10, -2E10);
				shift[j][i][k] = -2E10;
			}
			correlation[j][i] = 0.0;
			ncorrelation[j][i] = 0;
		}
		for (unsigned int a = 0; a < Nc; ++a) {
			//VSet(accumulator[j][a], 0.0, 0.0, 0.0);
			accumulator[j][a] = 0.0;
		}

		naccumulator[j] = 0;
		insertindex[j] = 0;
	}

	for (unsigned int i = 0; i < length; ++i) {
		t[i] = 0.0;
		f[i] = 0.0;
	}

	npcorr = 0;
	kmax = 0;
	for (unsigned int i = 0; i < Nc; ++i) {
		//VSet(accval[i], 0.0, 0.0, 0.0);
		accval[i] = 0.0;
	}

}



void Correlator::add(double *w, const unsigned int k) {
	/// If we exceed the correlator side, the value is discarded
	if (k == numcorrelators) return;
	if (k > kmax) kmax = k;

	/// Insert new value in shift array
	for (unsigned int i = 0; i < Nc; ++i) {
		//VCopy(shift[k][insertindex[k]][i], w[i]);	
		shift[k][insertindex[k]][i] = w[i];
	}

	/// Add to aveage value
	if (k == 0) {
		for (unsigned int a = 0; a < Nc; ++a) {
			//VVAdd(accval[a], w[a]);
			accval[a] += w[a];
		}
	}


	/// Add to accumulator and, if needed, add to next correlator
	/****************************************/
	for (unsigned int a = 0; a < Nc; ++a) {
		//VVAdd(accumulator[k][a], w[a]);
		accumulator[k][a] += w[a];
	}
	/****************************************
	if (naccumulator[k] == 0)
	{
		for (unsigned int a = 0; a < Nc; ++a)
		{
			VCopy(accumulator[k][a], w[a]);
		}
	}
	/****************************************/
	++naccumulator[k];
	if (naccumulator[k] == m) {
		/******************************************************/
		for (unsigned int a = 0; a < Nc; ++a) {
			//VSCopy(accumulator[k][a], 1 / m, accumulator[k][a]);
			//accumulator[k][a].x = accumulator[k][a].x / ((double)m);
			//accumulator[k][a].y = accumulator[k][a].y / ((double)m);
			//accumulator[k][a].z = accumulator[k][a].z / ((double)m);
			accumulator[k][a] = accumulator[k][a] / ((double)m);
		}
		/******************************************************/
		add(accumulator[k], k + 1);
		for (unsigned int a = 0; a < Nc; ++a) {
			//VSet(accumulator[k][a], 0, 0, 0);
			accumulator[k][a] = 0;
		}
		naccumulator[k] = 0;
	}

	/// Calculate correlation function
	unsigned int ind1 = insertindex[k];
	if (k == 0) { /// First correlator is different
		int ind2 = ind1;
		for (unsigned int j = 0; j < p; ++j) {
			for (unsigned int a = 0; a < Nc; ++a) {
				//if (shift[k][ind2][a].x > -1E10&&shift[k][ind2][a].y > -1E10&&shift[k][ind2][a].z > -1E10) 
				if (shift[k][ind2][a]> -1E10)
				{
					//VSub(dr, shift[k][ind2][a], shift[k][ind1][a]);
					//correlation[k][j] += VLenSq(dr);
					//correlation[k][j] += shift[k][ind2][a].x*shift[k][ind1][a].x + shift[k][ind2][a].y*shift[k][ind1][a].y + shift[k][ind2][a].z*shift[k][ind1][a].z;
					correlation[k][j] += shift[k][ind2][a]*shift[k][ind1][a];
					++ncorrelation[k][j];
				}
			}
			--ind2;
			if (ind2 < 0) ind2 += p;
		}
	}
	else {
		int ind2 = ind1 - dmin;
		for (unsigned int j = dmin; j < p; ++j) {
			if (ind2 < 0) ind2 += p;
			for (unsigned int a = 0; a < Nc; ++a) {
				//if (shift[k][ind2][a].x > -1E10&&shift[k][ind2][a].y > -1E10&&shift[k][ind2][a].z > -1E10)
				if (shift[k][ind2][a] > -1E10)
				{
					//VSub(dr, shift[k][ind2][a], shift[k][ind1][a]);
					//correlation[k][j] += VLenSq(dr);
					//correlation[k][j] += shift[k][ind2][a].x*shift[k][ind1][a].x + shift[k][ind2][a].y*shift[k][ind1][a].y + shift[k][ind2][a].z*shift[k][ind1][a].z;
					correlation[k][j] += shift[k][ind2][a]*shift[k][ind1][a];
					++ncorrelation[k][j];
				}
			}
			--ind2;
		}
	}

	++insertindex[k];
	if (insertindex[k] == p) insertindex[k] = 0;
}

void Correlator::evaluate(const bool norm) {
	unsigned int im = 0;

	/***************
	double aux = 0;

	if (norm) {
		aux = (accval.x / ncorrelation[0][0])*(accval.x / ncorrelation[0][0]) + \
			   	  (accval.y / ncorrelation[0][0])*(accval.y / ncorrelation[0][0]) + \
				  (accval.z / ncorrelation[0][0])*(accval.z / ncorrelation[0][0]);
	}
	/***************/

	/// First correlator
	for (unsigned int i = 0; i < p; ++i) {
		if (ncorrelation[0][i] > 0) {
			t[im] = i;
			f[im] = correlation[0][i] / ncorrelation[0][i]; // -aux;
			++im;
		}
	}

	// Subsequent correlators
	/***********************************************/
	for (unsigned int k = 1; k < kmax; ++k) {
		for (unsigned int i = dmin; i < p; ++i) {
			if (ncorrelation[k][i] > 0) {
				t[im] = i * pow((double)m, k);
				f[im] = correlation[k][i] / ncorrelation[k][i];
				++im;
			}
		}
	}
	/***********************************************
	for (int k = 1; k < kmax; ++k) {
		for (int i = dmin; i < p; ++i) {
			if (ncorrelation[k][i] > 0) {
				t[im] = i * pow((double)m, k);
				f[im] = correlation[k][i] / ncorrelation[k][i];
				++im;
			}
		}
	}
	/***********************************************/


	npcorr = im;
}

