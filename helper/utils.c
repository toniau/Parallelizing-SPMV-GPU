#include "utils.h"
#include <stdio.h>

//void getmul(const double* val, const double* vec, const int* rIndex, const int*cIndex, int nz, double* res)
void getmul(const float* val, const float* vec, const int* rIndex, const int*cIndex, int nz, float* res)
{
	int i; 
	for (i = 0; i < nz; i++)
	{
		int rInd = rIndex[i];
		int cInd = cIndex[i];
		res[rInd] += val[i] * vec[cInd];

/*
           if (rInd==0) printf("getmul: val[%d]=%lg, vec[%d]=%lg, res[%d]=%lg\n",
                     i, val[i], cInd, vec[cInd], rInd, res[rInd]);

                double dval = (double)val[i];
                double dvec = (double)vec[cInd];
                double dres = dval * dvec;

		res[rInd] += (float)dres;
*/
	}
}




//boolean checkerror(const double* resp, const double* ress, int dim)
boolean checkerror(const float* resp, const float* ress, int dim)
{
	//int errCount = 10;
	int i;

	for (i = 0; i < dim; i++)
	{
		double diff;
		diff = (double)(resp[i] - ress[i]);
		if ( diff < 0 ) diff = diff * (-1);

		//if (resp[i] != ress[i]) {
		//if ( diff > 0.00000000001 ) {

		if ( diff > 0.0001 ) {
			printf("checkerror: resp[%d]= %f, ress[%d]= %f, diff= %20.18lg\n", 
                                i, resp[i], i, ress[i], diff );

			//if ( --errCount == 0 ) return FALSE;
			return FALSE;
		}
	}

	return TRUE;

}
