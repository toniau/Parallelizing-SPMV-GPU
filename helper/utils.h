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
              const int nz, float *vec, int **rsIndex, int **reIndex); 

void quicksort(float* a, double* vindex, int* rindex, int* cindex, int n);

void printSortedVal(const int M, const int N, const int nz, 
                    const int *rIndex,  const int *cIndex,
                    const int *rsIndex, const int *reIndex, 
                    float *val, float *vec);

void printResult(const int M, float *res);
void printOutput(const int M, float *res);

void print_p(const int M, const int N, const int nz, const int *rIndex, float *p);
void print_p2(const int M, const int N, const int nz, const int blocksize, const int *rIndex,  float *p);
void print_p3(const int M, const int N, const int nz, const int blocksize, const int *rIndex,  float *p);
void print_p4(const int M, const int N, const int nz, const int blocksize, const int *rIndex,  float *p);
void print_val(const int M, const int N, const int nz, const int *rIndex,  float *p, float *v);



#endif
