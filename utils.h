#ifndef UTILS_H
#define UTILS_H
#include <time.h>

typedef struct {
	time_t tv_sec; // seconds 
	long tv_nsec;  // nanoseconds 
}timespect;


typedef enum {
  TRUE = 1, 
  FALSE = 0
} boolean;


//boolean checkerror(const double *resp, const double *ress, int dim);
boolean checkerror(const float *resp, const float *ress, int dim);

//void getmul(const double *val, const double *vec, const int *rIndex,
//	const int *cIndex, int nz, double *res);

void getmul(const float *val, const float *vec, const int *rIndex,
	const int *cIndex, int nz, float *res);


void sortVals(float *val, const int M, const int N, int *rIndex, int *cIndex, 
              const int nz, float *vec, int **rsIdx, int **reIdx, int **rCount);
 

//void quicksort(float* a, double* vindex, int* rindex, int* cindex, int n);
void quicksort(float* a, long* vindex, int* rindex, int* cindex, int n);


void quickSortRowCount(int left, int right, int *rCount, int *rsIndex, int *reIndex);

void printSortedVal(const int M, const int N, const int nz, const int blockSize,
                    const int *rIndex,  const int *cIndex,
                    const int *rsIndex, const int *reIndex, const int *rCount,
                    float *val, float *vec);

void printPrimeArrays(const int nnz, const int nChunk, const int threadsPerBlock, 
                      const int *rIdxPrime, const int *cIdxPrime, const float *valPrime);

void printResult(const int M, float *res);
void printInput(const int M, const int N, const int nz, const int *rIndex, const int *cIndex, float *val, float *vec);
void printOutput(const int M, float *res);

void print_p(const int M, const int N, const int nz, const int *rIndex, float *p);
void print_p2(const int M, const int N, const int nz, const int blocksize, const int *rIndex,  float *p);
void print_p3(const int M, const int N, const int nz, const int blocksize, const int *rIndex,  float *p);
void print_p4(const int M, const int N, const int nz, const int blocksize, const int *rIndex,  float *p);
void print_val(const int M, const int N, const int nz, const int *rIndex,  float *p, float *v);



#endif
