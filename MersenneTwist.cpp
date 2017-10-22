//---------------------------------------------------------------------------
#include <math.h>
#include <algorithm>
#include <stdio.h>

typedef unsigned int uint;

#include <assert.h>
#include "MersenneTwist.h"

//---------------------------------------------------------------------------
// uint must be an unsigned integer type capable of holding at least 32
// bits; exactly 32 should be fastest, but 64 is better on an Alpha with
// GCC at -O3 optimization so try your options and see what's best for you
//---------------------------------------------------------------------------
#define M              (397)                 // a period parameter
#define K              (0x9908B0DFU)         // a magic constant
#define hiBit(u)       ((u) & 0x80000000U)   // mask all but highest   bit of u
#define loBit(u)       ((u) & 0x00000001U)   // mask all but lowest    bit of u
#define loBits(u)      ((u) & 0x7FFFFFFFU)   // mask     the highest   bit of u
#define mixBits(u, v)  (hiBit(u)|loBits(v))  // move hi bit of u to hi bit of v

uint    MT::state[N+1];     // state vector + 1 extra to not violate ANSI C
uint   *MT::next;           // next random value is computed from here
int     MT::left = -1;      // can *next++ this many times before reloading


const double MT::pi = 3.14159265358979323846264338328;

//---------------------------------------------------------------------------
double MT::normOld  = 0;
int    MT::nNormOld = 0;

//---------------------------------------------------------------------------

void MT::SaveState(FILE *f)
{
#define SAVE(a) fwrite(&a, sizeof(a), 1, f);
	require(f != NULL);
	fwrite(&state[0], sizeof(uint), N+1, f);
	int k = next - state;
	SAVE(k);
	SAVE(left);
	SAVE(normOld);
	SAVE(nNormOld);
#undef SAVE
}

void MT::LoadState(FILE *f)
{
#define LOAD(a) require(fread(&a, sizeof(a), 1, f) == 1);
	require(f != NULL);
	require(fread(&state[0], sizeof(uint), N+1, f) == N+1);
	int k; LOAD(k); next = state + k;
	LOAD(left);
	LOAD(normOld);
	LOAD(nNormOld);
#undef LOAD
}

//---------------------------------------------------------------------------
void MT::Seed(uint seed)
{
	// We initialize state[0..(N-1)] via the generator
	//
	//   x_new = (69069 * x_old) mod 2^32
	//
	// from Line 15 of Table 1, p. 106, Sec. 3.3.4 of Knuth's
	// _The Art of Computer Programming_, Volume 2, 3rd ed.
	//
	// Notes (SJC): I do not know what the initial state requirements
	// of the Mersenne Twister are, but it seems this seeding generator
	// could be better.  It achieves the maximum period for its modulus
	// (2^30) iff x_initial is odd (p. 20-21, Sec. 3.2.1.2, Knuth); if
	// x_initial can be even, you have sequences like 0, 0, 0, ...;
	// 2^31, 2^31, 2^31, ...; 2^30, 2^30, 2^30, ...; 2^29, 2^29 + 2^31,
	// 2^29, 2^29 + 2^31, ..., etc. so I force seed to be odd below.
	//
	// Even if x_initial is odd, if x_initial is 1 mod 4 then
	//
	//   the          lowest bit of x is always 1,
	//   the  next-to-lowest bit of x is always 0,
	//   the 2nd-from-lowest bit of x alternates      ... 0 1 0 1 0 1 0 1 ... ,
	//   the 3rd-from-lowest bit of x 4-cycles        ... 0 1 1 0 0 1 1 0 ... ,
	//   the 4th-from-lowest bit of x has the 8-cycle ... 0 0 0 1 1 1 1 0 ... ,
	//    ...
	//
	// and if x_initial is 3 mod 4 then
	//
	//   the          lowest bit of x is always 1,
	//   the  next-to-lowest bit of x is always 1,
	//   the 2nd-from-lowest bit of x alternates      ... 0 1 0 1 0 1 0 1 ... ,
	//   the 3rd-from-lowest bit of x 4-cycles        ... 0 0 1 1 0 0 1 1 ... ,
	//   the 4th-from-lowest bit of x has the 8-cycle ... 0 0 1 1 1 1 0 0 ... ,
	//    ...
	//
	// The generator's potency (min. s>=0 with (69069-1)^s = 0 mod 2^32) is
	// 16, which seems to be alright by p. 25, Sec. 3.2.1.3 of Knuth.  It
	// also does well in the dimension 2..5 spectral tests, but it could be
	// better in dimension 6 (Line 15, Table 1, p. 106, Sec. 3.3.4, Knuth).
	//
	// Note that the random number user does not see the values generated
	// here directly since reloadMT() will always munge them first, so maybe
	// none of all of this matters.  In fact, the seed values made here could
	// even be extra-special desirable if the Mersenne Twister theory says
	// so-- that's why the only change I made is to restrict to odd seeds.

	register uint x = (seed | 1U) & 0xFFFFFFFFU, *s = state;
	register int    j;

	for(left=0, *s++=x, j=N; --j;	) *s++ = (x*=69069U) & 0xFFFFFFFFU;
}

//---------------------------------------------------------------------------
uint MT::Reload()
{
	register uint *p0=state, *p2=state+2, *pM=state+M, s0, s1;
	register int    j;

	if(left < -1) MT::Seed(4357U);

	left=N-1, next=state+1;

	for(s0=state[0], s1=state[1], j=N-M+1; --j; s0=s1, s1=*p2++)
		*p0++ = *pM++ ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);

	for(pM=state, j=M; --j; s0=s1, s1=*p2++)
		*p0++ = *pM++ ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);

	s1=state[0], *p0 = *pM ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);
	s1 ^= (s1 >> 11);
	s1 ^= (s1 <<  7) & 0x9D2C5680U;
	s1 ^= (s1 << 15) & 0xEFC60000U;
	return(s1 ^ (s1 >> 18));
}

//---------------------------------------------------------------------------
//	uint MT::Poisson(double lambda)
//	{
//		static const double PiFactor = 3.1415926535897932384626433832795028842/4294967296.0;
//		double g = exp(-lambda), t = 1;
//		int    n = -1;
//
//		if (lambda < 12.0) {
//			do {
//				++n;
//				t *= Uni();
//			} while (t > g);
//		} else {
//			double sq = sqrt(2.0 * lambda);
//			double alxm = log(lambda);
//			g = lambda * alxm - LnGamma(lambda + 1);
//
//			// sample until not rejected
//			double y,z;
//			do {
//				do {
//					// y is deviate from Lorentzian
//					y = tan(PiFactor * Rand());
//					// comparison function
//					z = sq * y + lambda;
//				} while (z < 0.0);
//				z = floor(z); // convert to integer
//				t = 0.9 * (1 + y*y) * exp(z*alxm - LnGamma(z + 1) - g);
//			} while (Uni() > t);
//			// set n to accepted value
//			n = (uint)z;
//		}
//
//		return n;
//	}


//---------------------------------------------------------------------------

int MT::Binomial(const int n, const double &pp) 
{
// 	const double p = (pp <= 0.5 ? pp : 1.0 - pp); // The binomial distribution is invariant under changing pp to 1-pp, 
// 												// if we also change the answer to n minus it self; weÕll remember to 
// 												// do this below. 
// 	
// 	int bnl = 0; 
// 	for(int j = 0; j < n; j++) if( Uni() < p) ++bnl;
// 	
// 	if( p != pp ) bnl = n - bnl; // Remember to undo the symmetry transformation.	
// 	return bnl; 

#define PI 3.1415926535897932385 
// Returns as a floating-point number an integer value that is a 
//	random deviate drawn from a binomial distribution of n trials 
// each of probability pp, using Uni() as a source of uniform 
	// random deviates. 
	int j, bnl; 
	static int nold = -1; 
	static double pold = (-1.0), pc, plog, pclog, en, oldg; 
	const double p = (pp <= 0.5 ? pp : 1.0 - pp); // The binomial distribution is invariant under changing pp to 1-pp, 
												// if we also change the answer to n minus it self; weÕll remember to 
												// do this below. 
	const double am = n*p; // This is the mean of the deviate to be produced. 
	if ( n < 25 ) { // Use the direct method while n is not too large. This can require up to 25 calls to ran1.
		bnl = 0; 
		for(int j = 1; j <= n; j++) if( Uni() < p) ++bnl; 
	} else if (am < 1.0) { // If fewer than one event is expected out of 25 or more trials, 
								  // then the distribution is quite accurately Poisson. Use direct Poisson method. 
		const double g = exp(-am); 
		double t = 1.0;
		for(j = 0; j <= n; j++) { 
			t *= Uni(); 
			if(t < g) break; 
		} 
		bnl = (j <= n ? j : n); 
		
	} else { // Use the rejection method.
		
		if(n != nold) { // If n has changed, then compute useful quantities.
			en = n; oldg = LnGamma(en+1.0); nold = n; 
		}
		
		if(p != pold){ // If p has changed, then compute useful quantities.
			pc = 1.0 - p; plog = log(p); pclog = log(pc); pold = p; 
		} 
		
		const double sq = sqrt(2.0*am*pc); // The following code should by now seem familiar: rejection method with a Lorentzian comparison function. 
		double t, y, em;
		do { 
			do { 
				const double angle = PI*Uni(); 
				y = tan(angle); 
				em = sq*y + am; 
			} while(em < 0.0 || em >= (en + 1.0));  // Reject. 
			
			em = floor(em); // Trick for integer-valued distribution. 
			t = 1.2 * sq * (1.0 + y*y) * exp(oldg - LnGamma(em + 1.0) - LnGamma(en - em + 1.0) + em*plog + (en - em)*pclog); 
			
		} while( Uni() > t); // Reject. This happens about 1.5 times per deviate, on average.
		
		bnl = (int)em; 
	} 
	
	if( p != pp ) bnl = n - bnl; // Remember to undo the symmetry transformation.
	
	return bnl; 
}

//---------------------------------------------------------------------------

static double ProbKS(double x)
{
	double a2,fac=2.0, sum=0.0, term, termbf=0.0;

	a2 = -2.0 * x * x;
	for (int j=1; j <=100; j++) {
		term = exp(a2 * (j * j));
		sum += fac * term;
		if (term <= 1.0e-10 * sum) return sum;
		fac    = -fac;
		termbf = term;
	}
	return 1.0;
}


double MT::KolmogorovSmirnovTest(int n, int data[], double scale)
{
	double d = 0, dt, en, ff, fn, fo =0;

	std::sort(data, data + n);
	
	en = 1./(double)n;
	for (int j=0; j < n; j++) {
		fn = j*en; ff = data[j]*scale;
		dt = fabs(fo-ff);
		{ 
			double yy = fabs(fn-ff);
			if (yy > dt) dt = yy;
		}
		if (dt > d) d = dt;
		fo = fn;
	}

	return ProbKS((sqrt(n) + 0.12 + 0.11 * sqrt(en)) * d);
}

//---------------------------------------------------------------------------
