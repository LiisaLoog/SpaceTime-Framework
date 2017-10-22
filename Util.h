#ifndef __UTIL_H__
#define __UTIL_H__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <assert.h>

using std::vector;


#ifndef require
#define require(a) assert(a)
#endif


//#ifndef require
//#define require(a) \
//if (!(a)) { \
//	fprintf(stderr, "\n%% require( %s ) failed at line %d in file '%s'\n\n", #a, __LINE__, __FILE__ ); \
//		exit(1); \
//}
//#endif

template<class T>
bool ReadParameter(const int ac, char *av[], 
						 const char *str, T &var)
{
	require(str != 0);
	const int nMax = 100;
	int n = strlen(str); require(n > 0); require(n + 1 < nMax);
	require(str[n] == 0);
	char buf[nMax]; strcpy(buf, str); buf[n] = '='; buf[++n] = 0;
	require(buf[n] == 0);
	
	bool set = false;
	for (int i = 0; i < ac; i++) {
		if (strncmp(buf, av[i], n) == 0) {
			require(av[i][n-1] == '=');
			var = (T)atof(av[i] + n);
#ifndef READPARAMETER_IS_QUIET
			printf("%% parameter %s set to %lg\n", str, (double)var);
#endif
			set = true;
		}
	}
	
#ifndef READPARAMETER_IS_QUIET
	if (!set) printf("%% parameter %s is %lg\n", str, (double)var);
#endif
	return set;
}



//---------------------------------------------------------------	
// Parameter loading for all numerical types
//---------------------------------------------------------------	

template<class T>
T ReadParameter2(const int ac, char *av[], 
					  const char *str, const T defaultValue)
{
	require(str != 0);
	const int nMax = 100;
	int n = strlen(str); require(n > 0); require(n + 1 < nMax);
	require(str[n] == 0);
	char buf[nMax]; strcpy(buf, str); buf[n] = '='; buf[++n] = 0;
	require(buf[n] == 0);
	
	T var = defaultValue;
	
	bool set = false;
	for (int i = 0; i < ac; i++) {
		if (strncmp(buf, av[i], n) == 0) {
			require(av[i][n-1] == '=');
			var = (T)atof(av[i] + n);
#ifndef READPARAMETER_IS_QUIET
#ifdef READPARAMETER_M_FILE
			printf("%s = %lg; %% parameter set\n", str, (double)var);
#else
			printf("%% %s = %lg; %% parameter set\n", str, (double)var);
			// printf("%% parameter %s set to %lg\n", str, (double)var);
#endif
#endif
			set = true;
		}
	}
	
#ifndef READPARAMETER_IS_QUIET
#ifdef READPARAMETER_M_FILE
	if (!set) printf("%s = %lg; %% default value\n", str, (double)var);
#else
	if (!set) printf("%% %s = %lg; %% default value\n", str, (double)var);
	// if (!set) printf("%% parameter %s is %lg\n", str, (double)var);
#endif
#endif
	
	return var;
}

//---------------------------------------------------------------	
// Template specialization for strings
//---------------------------------------------------------------	

template<>
const char * ReadParameter2(const int ac, char *av[], 
		  		  	             const char *str, const char *defaultValue)
{
	require(str != 0);
	const int nMax = 100;
	int n = strlen(str); require(n > 0); require(n + 1 < nMax); require(str[n] == 0);
	char buf[nMax]; strcpy(buf, str); buf[n] = '='; buf[++n] = 0; require(buf[n] == 0);
	
	const char *var = defaultValue;
	
	bool set = false;
	for (int i = 0; i < ac; i++) {
		if (strncmp(buf, av[i], n) == 0) {
			require(av[i][n-1] == '=');
			var = av[i] + n;
#ifndef READPARAMETER_IS_QUIET
#ifdef READPARAMETER_M_FILE
			printf("%s = '%s'; %% parameter set\n", str, var);
#else
			printf("%% %s = '%s'; %% parameter set\n", str, var);
			// printf("%% parameter %s set to %s\n", str, var);
#endif
#endif
			set = true;
		}
	}
	
#ifndef READPARAMETER_IS_QUIET
#ifdef READPARAMETER_M_FILE
	if (!set) printf("%s = '%s'; %% default value\n", str, var);
#else
	if (!set) printf("%% %s = '%s'; %% default value\n", str, var);
	// if (!set) printf("%% parameter %s is '%s'\n", str, var);
#endif
#endif
	
	return var;
}

//---------------------------------------------------------------	
// Template specialization for bools
//---------------------------------------------------------------	

template<>
bool ReadParameter2(const int ac, char *av[], 
						  const char *str, bool defaultValue)
{
	require(str != 0);
	const int nMax = 100;
	int n = strlen(str); require(n > 0); require(n + 1 < nMax); require(str[n] == 0);
	char buf[nMax]; strcpy(buf, str); buf[n] = '='; buf[++n] = 0; require(buf[n] == 0);
	// fprintf(stderr,"'%s',%d\n", buf, n);
	
	bool var = defaultValue;
	
	bool set = false;
	for (int i = 0; i < ac; i++) {
		if (strncmp(buf, av[i], n) == 0) {
			require(av[i][n-1] == '=');
			set = true;
			switch(av[i][n]) {
				case 't': case 'T': var = true;  break;
				case 'f': case 'F': var = false; break;
				default: set = false;
			}
			
			if (set) {
#ifndef READPARAMETER_IS_QUIET
#ifdef READPARAMETER_M_FILE
				printf("%s = %s; %% parameter set\n", str, (var ? "true" : "false"));
#else
				printf("%% %s = %s; %% parameter set\n", str, (var ? "true" : "false"));
				// printf("%% parameter %s set to %s\n", str,  (var ? "true" : "false"));
#endif
#endif
			}
		}
	}
	
#ifndef READPARAMETER_IS_QUIET
#ifdef READPARAMETER_M_FILE
	if (!set) printf("%s = %s; %% default value\n", str,  (var ? "true" :"false"));
#else
	if (!set) printf("%% %s = %s; %% default value\n", str,  (var ? "true" :"false"));
	//	if (!set) printf("%% parameter %s is %s\n", str,  (var ? "true" :"false"));
#endif
#endif
	
	return var;
}


//---------------------------------------------------------------	
// Loadnumbers
//
// given a buffer and a vector of a numeric type, read numbers from the buffer into the vector.
// return # of numbers read
//
// given a file pointer and a vector of a numeric type, read  numbers from the next line of the file into the vector.
// return # of numbers read (>= 0). if end of file, return -1.
//
//---------------------------------------------------------------	

template<class T>
int LoadNumbers(const char buf[], vector<T> &f, bool clear_f = true)
{
	if (clear_f) f.clear();
	
	int n_col = 0, offs = 0, n_read = 0;
	double x;
	while (sscanf(buf + offs, "%lg%n", &x, &n_read) == 1) {
		f.push_back((T)x);
		n_col++;
		offs += n_read;
	}
	
	return n_col;
}

template<class T>
int LoadNumbers(FILE *file, vector<T> &f, bool clear_f = true)
{
	const int n_buf = 10240;
	char buf[n_buf];
	
	// load a line
	buf[0] = 0; fgets(buf, n_buf, file); if (feof(file)) return -1;
	
	LoadNumbers(buf, f, clear_f);
}

//---------------------------------------------------------------	
// Print vector and matrix
//---------------------------------------------------------------	

void PrintMatrix(const int A[], const int m, const int n) 
{
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) printf(" %2d", A[n*i + j]);
		printf("\n");
	}
}

void PrintMatrix(const double A[], const int m, const int n) 
{
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) printf(" %.12lg", A[n*i + j]);
		printf("\n");
	}
}


#define PrintVec(a,n)   { printf( #a " = [\n"); PrintMatrix(&(a)[0], 1, n);    printf("];\n"); }
#define PrintMat(a,m,n) { printf( #a " = [\n"); PrintMatrix(&(a)[0][0], m, n); printf("];\n"); }

//---------------------------------------------------------------	
// Short macros
//---------------------------------------------------------------	

#define RP(param) ReadParameter(ac, av, #param, param)

#define RPC(param, defaultValue) ReadParameter2(ac, av, #param, defaultValue)

#define Parameter(type, param, defaultValue) type param = ReadParameter2(ac, av, #param, (type)defaultValue)

typedef unsigned int uint;

#define REP(param) printf(#param "=%.12lg; ", (double)param)
#define ENDL printf("\n")
#define REPL(param) printf(#param " = %.12lg;\n\n", (double)param)


#define REPE(param) fprintf(stderr,#param "=%.12lg; ", (double)param)
#define ENDLE fprintf(stderr,"\n")

#endif

