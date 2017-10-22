//---------------------------------------------------------------------------
// (c) Anders Eriksson 2017
//
// To compile, run
// g++ -Wall -O3 -o SimChangeMigrate SimChangeMigrate.cpp MersenneTwist.cpp
//
//---------------------------------------------------------------------------


#include "Util.h"
#include "MersenneTwist.h"

//---------------------------------------------------------------------------
// Main function
//---------------------------------------------------------------------------

int main(int ac, char *av[])
{
	// first print general info
	printf("%% "); { for (int i=0; i < ac; i++) printf("%s ", av[i]); } printf("\n\n");
	printf("%% seed = %d;\n", MT::SeedTime());
	
	
	//---------------------------------------------------------------------------
	// read parameters for population dynamics
	//---------------------------------------------------------------------------

	Parameter(const int, nTraits, 1000); // Number of traits to simulate
	Parameter(const int, N, 200); // Population size
	Parameter(const int, nTrans,10000); // Number of transient steps
	Parameter(const double, pDeath,0.1); // Live on average 1/pDeath time steps
	Parameter(const double, mu, 0.01); // Mutation rate
	Parameter(const int, tSigma1, 2000); // Time of first episode of constant mobility
	Parameter(const int, tSigma2, tSigma1); // Time of second mobility episode
	Parameter(const int, tSigma3, tSigma1); // Time of third mobility episode
	
	Parameter(const double, sigma1, 0.001); // Mobility rate in first mobility episode
	Parameter(const double, sigma2, 0.001); // Mobility rate in second mobility episode
	Parameter(const double, sigma3, 0.001); // Mobility rate in third mobility episode

	const int nMeasureSteps = tSigma1+tSigma2+tSigma3;
	
	Parameter(const int, nMeasureInds,  1000); // Number of individuals to sample
	
	fflush(stdout);
	
	double * const trait = new double[N*nTraits];
	double * const xPos  = new double[N];
	double * const yPos  = new double[N];
	
	// init positions
	for (int i = 0; i < N; i++) {
		xPos[i] = MT::Uni();
		yPos[i] = MT::Uni();
		for (int k = 0; k < nTraits; k++) trait[i*nTraits+k] = 0;
	}
	
	const double sqrt_of_3 = sqrt(3);
	
	// transient
	fprintf(stderr,"Transient");
	for (int t = 0; t < nTrans; t++) {
		
		double sigma = sigma1;
		
		if ((t+1) % 500 == 0) fprintf(stderr," %d", t+1);
		
		for (int nChange = MT::Poisson(pDeath*N); nChange > 0; --nChange) {
			const int i = MT::Rand() % N, j = MT::Rand() % N;
			xPos[j] = xPos[i];
			yPos[j] = yPos[i];
			for (int k = 0; k < nTraits; k++) trait[j*nTraits+k] = trait[i*nTraits+k];
		}
		
		for (int i = 0; i < N; i++) {
			xPos[i] += sigma*(sqrt_of_3*(2*MT::Uni()-1));
			if (xPos[i] < 0) xPos[i] = -xPos[i];
			if (xPos[i] > 1) xPos[i] = 2-xPos[i];
			assert(0 <= xPos[i] && xPos[i] <= 1);
			
			yPos[i] += sigma*(sqrt_of_3*(2*MT::Uni()-1));
			if (yPos[i] < 0) yPos[i] = -yPos[i];
			if (yPos[i] > 1) yPos[i] = 2-yPos[i];
			assert(0 <= yPos[i] && yPos[i] <= 1);

			for (int k = 0; k < nTraits; k++) trait[i*nTraits+k] += mu*(sqrt_of_3*(2*MT::Uni()-1));
		}
	}
	
	fprintf(stderr,"\nMeasure");
	
	int tMeasure[nMeasureInds];
	for (int i = 0; i < nMeasureInds; i++) tMeasure[i] = MT::Rand() % nMeasureSteps;
	std::sort(&tMeasure[0],&tMeasure[nMeasureInds]);
	for (int i = 1; i < nMeasureInds; i++) assert(tMeasure[i-1] <= tMeasure[i]);
	
	
	for (int t = 0, iActiveInd = 0; t < nMeasureSteps; t++) {
		
		double sigma = sigma1;
		if (t >= tSigma1) {
			sigma = sigma2;
			if (t >= tSigma1+tSigma2) sigma = sigma3;
		}
		
		if ((t+1) % 200 == 0) fprintf(stderr," %d", t+1);

		if (iActiveInd < nMeasureInds) {
			
			assert(tMeasure[iActiveInd] >= t);
			while (iActiveInd < nMeasureInds && tMeasure[iActiveInd] == t) {
				const int i = MT::Rand() % N;
				printf("%d %.6lg %.6lg", t, xPos[i],yPos[i]);
				for (int k = 0; k < nTraits; k++) printf(" %.6lg", trait[i*nTraits+k]);
				printf("\n");
				iActiveInd++;
			}
		}
		
		for (int nChange = MT::Poisson(pDeath*N); nChange > 0; --nChange) {
			const int i = MT::Rand() % N, j = MT::Rand() % N;
			xPos[j] = xPos[i];
			yPos[j] = yPos[i];
			for (int k = 0; k < nTraits; k++) trait[j*nTraits+k] = trait[i*nTraits+k];
		}
		
		for (int i = 0; i < N; i++) {
			xPos[i] += sigma*(sqrt_of_3*(2*MT::Uni()-1));
			if (xPos[i] < 0) xPos[i] = -xPos[i];
			if (xPos[i] > 1) xPos[i] = 2-xPos[i];
			assert(0 <= xPos[i] && xPos[i] <= 1);
			
			yPos[i] += sigma*(sqrt_of_3*(2*MT::Uni()-1));
			if (yPos[i] < 0) yPos[i] = -yPos[i];
			if (yPos[i] > 1) yPos[i] = 2-yPos[i];
			assert(0 <= yPos[i] && yPos[i] <= 1);
			
			for (int k = 0; k < nTraits; k++) trait[i*nTraits+k] += mu*(sqrt_of_3*(2*MT::Uni()-1));
		}
	}
	
	fprintf(stderr,"\n");

	return 0;
}
