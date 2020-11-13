/// /// Definition of correlator classes
#ifndef __FCorrelatorSAt_h
#define __FCorrelatorSAt_h

#include <stdio.h>
//#include "parameters.h"


////////////////////////////////////////////////////
/// Standard Scalar Correlator f(tau)=<A(t)A(t+tau)>
namespace SAt {
	class Correlator {

	protected:
		/** Where the coming values are stored */
		double ***shift; //nMID
		/** Array containing the actual calculated correlation function */
		double **correlation;
		/** Number of values accumulated in cor */
		unsigned long int **ncorrelation;

		/** Accumulator in each correlator */
		double **accumulator; //nMID
		/** Index that controls accumulation in each correlator */
		unsigned int *naccumulator;
		/** Index pointing at the position at which the current value is inserted */
		unsigned int *insertindex;

		/** Number of Correlators */
		unsigned int numcorrelators;

		/** Minimum distance between points for correlators k>0; dmin = p/m */
		unsigned int dmin;

		/*  SCHEMATIC VIEW OF EACH CORRELATOR
		p=N
		<----------------------------------------------->
		_________________________________________________
		|0|1|2|3|.|.|.| | | | | | | | | | | | | | | |N-1|
		-------------------------------------------------
		*/

		/** Lenght of result arrays */
		unsigned int length;
		/** Maximum correlator attained during simulation */
		unsigned int kmax;

	public:
		/** Points per correlator */
		unsigned int p;
		/** Number of points over which to average; RECOMMENDED: p mod m = 0 */
		unsigned int m;
		double *t, *f;
		unsigned int npcorr;

		unsigned int Nc;

		/** Accumulated result of incoming values **/
		//double accval;
		double *accval;
		double dr;

		/** Constructor */
		Correlator() { numcorrelators = 0; };
		Correlator(const unsigned int numcorrin, const unsigned int pin, const unsigned int min, unsigned int Nchain);
		~Correlator();

		/** Set size of correlator */
		void setsize(const unsigned int numcorrin = 32, const unsigned int pin = 16, const unsigned int min = 2, unsigned int Nchain = 0);

		/** Add a scalar/vector to the correlator number k */
		//double *w = new double[Nc];
		void add(double *w, const unsigned int k = 0);

		/** Evaluate the current state of the correlator */
		void evaluate(const bool norm = false);

		/** Initialize all values (current and average) to zero */
		void initialize();

	};

}
#endif