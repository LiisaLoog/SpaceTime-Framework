//---------------------------------------------------------------------------

#ifndef __MERSENNE_TWISTER__
#define __MERSENNE_TWISTER__

#ifndef require
#define require(a) assert(a)
#endif

//if (!(a)) { \
//	fprintf(stderr, "\n%% require( %s ) failed at line %d in file '%s'\n\n", #a, __LINE__, __FILE__ ); \
//		exit(1); \
//}
//#endif

#ifndef uint
typedef unsigned int uint;
#endif

//#include "Util.h"
#include <math.h>
#include <vector>
#include <time.h>
#include <sys/time.h>

//---------------------------------------------------------------------------
// This is the ``Mersenne Twister'' random number generator MT19937, which
// generates pseudorandom integers uniformly distributed in 0..(2^32 - 1)
// starting from any odd seed in 0..(2^32 - 1).  This version is a recode
// by Shawn Cokus (Cokus@math.washington.edu) on March 8, 1998 of a version by
// Takuji Nishimura (who had suggestions from Topher Cooper and Marc Rieffel in
// July-August 1997).
//---------------------------------------------------------------------------
class MT {
	enum {N = 624};				// length of state vector
	static uint    state[N+1];
	static uint   *next;
	static int     left;
	
	static uint Reload();
	
	static double normOld;
	static int    nNormOld;
	
public:

	// computes ln(gamma(z)) for z > 0, from Num.Rec.C ch.6.1 pp. 214
	static inline double LnGamma(double z) {
		static const double coeff[6] = {
			76.18009172947146,     -86.50532032941677,
			24.01409824083091,     -1.231739572450155,
			0.1208650973866179e-2, -0.5395239384953e-5
		};
		
		double x = z, ser = 1.000000000190015;
		for (int j=0; j < 6; j++) ser += coeff[j]/(++x);
		
		x = z + 5.5;
		x -= (z + 0.5)*log(x);
		return log(2.5066282746310005 * ser / z) - x;
	}
	
	static const double pi;


	static void Seed(uint seed);
	
	static uint SeedTime(uint seed = 0) {
		timeval tp;
		if (seed == 0) {
			gettimeofday(&tp, NULL);
			seed = (tp.tv_sec << 16) ^ (tp.tv_sec >> 16) ^ tp.tv_usec;
		}
		Seed(seed);
		return seed;
	}

	static void SaveState(FILE *f);
	static void LoadState(FILE *f);
	
	// random integer in 0 .. 2^32-1
	static inline uint Rand() {    
		if(--left < 0) return(Reload());
		
		uint y  = *next++;
		y ^= (y >> 11);
		y ^= (y <<  7) & 0x9D2C5680U;
		y ^= (y << 15) & 0xEFC60000U;
		return(y ^ (y >> 18));
	}
	
	// Uniform in [0, 1)
	static inline double Uni(){ return Rand() * 2.3283064365386962890625000e-10; }
	
	// Exponential(1)
	static inline double Exp(){
		uint r = Rand();
		return (r == 0) ? 0 : -log(r*2.3283064365386962890625000e-10);
	}

	// geometrical distribution, p(k) = r (1 - r)^k, 0 < r <= 1
	// tested 2006-05-17 ok.
	// Expected value of k = 1/r - 1
	//
	// r = 1/(1 + <k>)
	//
	static inline int Geometric(const double &r) {
		int k = 0;
		if (r > 0.1) {
			const double t = 1 - r;
			double w = t, z = MT::Uni();
			while (w > z) { w *= t; k++; }
		} else {
			k = (int)( log(MT::Uni()) / ( (r < 5e-6) ? -(r + 0.5*r*r) : log(1. - r) ) );
		}
		return k;
	}		
	
	// Gamma distributed numbers
	// from Num.Rec.C ch.7.3 pp.292-293
	//
	// Returns a deviate distributed as a gamma distribution of integer order ia, 
	// i.e., a waiting time to the ia:th event in a Poisson process of unit mean, 
	// using ran1(idum) as the source of uniform deviates. 
	//
	static double Gamma(const int ia) 
	{ 
		double x;
		if (ia < 6) { // Use direct method, adding waiting times.
			x = 1.0; 
			for (int j = 0; j < ia; j++) x *= Uni(); x = -log(x); 
		} else { // Use rejection method.
			const double am = ia - 1; 
			const double s  = sqrt((double)(2*ia - 1));
			double y, e;
			
			do {
				do {
					double v1, v2;
					do { 
						// These four lines generate the tangent of a random angle, 
						// i.e., they are equivalent to y = tan(¹* ran1(idum)).
						v1 = Uni(); 
						v2 = 2.0*Uni() - 1.0; 
					} while (v1*v1 + v2*v2 > 1.0); 
					y = v2/v1;
					x = s*y + am; // We decide whether to reject x:
				} while (x <= 0.0); // Reject in region of zero probability.
				e = (1.0 + y*y)*exp(am*log(x/am) - s*y); // Ratio of prob. fn. to comparison fn.
			} while (Uni() > e); // Reject on basis of a second uniform deviate
		}
		return x;
	}
	
	static double GammaRnd(const double ia) 
	{ 
		double x;

		// Use rejection method.
		const double am = ia - 1; 
		const double s  = sqrt((double)(2*ia - 1));
		double y, e;
		
		do {
			do {
				double v1, v2;
				do { 
					// These four lines generate the tangent of a random angle, 
					// i.e., they are equivalent to y = tan(¹* ran1(idum)).
					v1 = Uni(); 
					v2 = 2.0*Uni() - 1.0; 
				} while (v1*v1 + v2*v2 > 1.0); 
				y = v2/v1;
				x = s*y + am; // We decide whether to reject x:
			} while (x <= 0.0); // Reject in region of zero probability.
			e = (1.0 + y*y)*exp(am*log(x/am) - s*y); // Ratio of prob. fn. to comparison fn.
		} while (Uni() > e); // Reject on basis of a second uniform deviate
		
		return x;
	}
		
	// poisson distributed numbers
	// from Num.Rec.C ch.7.3 pp.294-295 
	// static uint Poisson(double lambda);
	static inline int Poisson(double lambda)
	{
		static const double PiFactor = 3.1415926535897932384626433832795028842/4294967296.0;
		double g = exp(-lambda), t = 1;
		int    n = -1;
		
		if (lambda < 12.0) {
			do {
				++n;
				t *= Uni();
			} while (t > g);
		} else {
			double sq = sqrt(2.0 * lambda);
			double alxm = log(lambda);
			g = lambda * alxm - LnGamma(lambda + 1);
			
			// sample until not rejected
			double y,z;
			do {
				do {
					// y is deviate from Lorentzian
					y = tan(PiFactor * Rand());
					// comparison function
					z = sq * y + lambda;
				} while (z < 0.0);
				z = floor(z); // convert to integer
				t = 0.9 * (1 + y*y) * exp(z*alxm - LnGamma(z + 1) - g);
			} while (Uni() > t);
			// set n to accepted value
			n = (uint)z;
		}
		
		return n;
	}
	

	// binomially distributed numbers
	// from Num.Rec.C ch.7.3 pp.294-295 
	// implementation tested ok anders 2005-09-13
	static int Binomial(const int n, const double &pp); 
	
	// N(0,1) using Box-Müller method w. rejection test
	static inline double Norm() {
		if (--nNormOld == 0) return normOld;
		nNormOld = 1;

		double x1, x2, w;

		do {
			x1 = Rand() * 4.6566128730773925781250000e-10 - 1.0;
			x2 = Rand() * 4.6566128730773925781250000e-10 - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w > 1.0 || w == 0.0);

		w = sqrt(((-2.0) * log(w)) / w);
		normOld = x1 * w;
		return x2 * w;
	}

  // Do a Kolmogorov-Smirnov test against a uniform distribution.
	// from Numerical Recipies in C, pp. 623 - 627.
	static double KolmogorovSmirnovTest(int n, int data[], double scale);
	
	
};

//---------------------------------------------------------------------------


template<class T> void RandomPermutation(T v[], const int n)
{
	if (n < 2) return;

	for (int i = 1; i < n; i++) std::swap(v[i], v[MT::Rand()%(i+1)]);
}

#endif
