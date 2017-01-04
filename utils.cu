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



void sortVals(float *val, const int M, const int N, int *rIndex, int *cIndex, 
              const int nz, float *vec, int **rsIdx, int **reIdx, int **rCnt) 
{
	int i;

	//preprocess the dataset to make the calculation can be parallelized
	//
	long *vIndex = (long*)malloc(nz*sizeof(long));
	memset(vIndex, 0, nz*sizeof(long));
	for (i = 0; i < nz; i++)
	{
		vIndex[i] = (long)rIndex[i] * N + cIndex[i];
		if (vIndex[i] < 0)
		{	
	           printf("Error!\n");
                   printf("i= %d, rIndex[%d]= %d, N= %d, cIndex[%d]= %d, vIndex[%d]= %lg\n", 
                           i, i, rIndex[i], N, i, cIndex[i], i, vIndex[i]);
	           exit(1);
	        }
	}

	quicksort(val, vIndex, rIndex, cIndex, nz);


	//we use rsIndex/reIndex to keep the start/end position of each row. The intial values are 
	//-1 for all entries.  rsIndex[i] indicates the start poistion of the i-th row. Hence 
	//the position index of the i-th row is from rsIndex[i] to reIndex[i]
	//
     
	int *rsIndex = (int*)malloc(M*sizeof(int)); //start/end position of each row
	memset(rsIndex, -1, M*sizeof(int));

	int *reIndex = (int*)malloc(M*sizeof(int));
	memset(reIndex, -1, M*sizeof(int));

	int *rCount = (int*)malloc(M*sizeof(int));
	memset(rCount, 0, M*sizeof(int));


        int curRowIndex = 0;

	for (i = 0; i<nz; i++)
	{
            curRowIndex = rIndex[i];

            rsIndex[curRowIndex] = i;
            while ((i < nz) && (rIndex[i] == rIndex[i + 1])) {
                // next element belongs to the same row
                i++;
            }; 

            // next element belongs to the next row
            reIndex[curRowIndex] = i;
            rCount[curRowIndex]  = reIndex[curRowIndex] - rsIndex[curRowIndex] + 1;
        }



/*
	for (i = 0; i<nz; i++)
	{
		int tmp = (int)(vIndex[i] / N);
if ( (tmp < 0) || (tmp > nz-1)) printf("invalid tmp= %d, vIndex[%d]= %lg, N= %d\n", tmp, i, vIndex[i],N);

		if (rsIndex[tmp] == -1)
		{
			rsIndex[tmp] = i;
			reIndex[tmp] = i;
                        rCount[tmp]  = 1;
		}
		else {
			reIndex[tmp] = i;
                        rCount[tmp]++; 
                }
	}
*/


        *rsIdx = rsIndex;
        *reIdx = reIndex;
        *rCnt  = rCount;

        free(vIndex);
}




//sorting according to the index
//void quicksort(float* a, double* vindex, int* rindex, int* cindex, int n)
void quicksort(float* a, long* vindex, int* rindex, int* cindex, int n)
{
	int i, j, m;

	double p, s;
	float  t;

	if (n < 2)
		return;
	p = vindex[n / 2];

	for (i = 0, j = n - 1;; i++, j--) {
		while (vindex[i]<p)
			i++;
		while (p<vindex[j])
			j--;
		if (i >= j)
			break;
		t = a[i];
		a[i] = a[j];
		a[j] = t;

		s = vindex[i];
		vindex[i] = vindex[j];
		vindex[j] = s;

		m = rindex[i];
		rindex[i] = rindex[j];
		rindex[j] = m;

		m = cindex[i];
		cindex[i] = cindex[j];
		cindex[j] = m;
	}
	quicksort(a, vindex, rindex, cindex, i);
	quicksort(a + i, vindex + i, rindex + i, cindex + i, n - i);
}



void swap(int num1, int num2, int *rCount, int *rsIndex, int *reIndex) 
{
   int temp = rCount[num1];
   rCount[num1] = rCount[num2];
   rCount[num2] = temp;

   temp = rsIndex[num1];
   rsIndex[num1] = rsIndex[num2];
   rsIndex[num2] = temp;

   temp = reIndex[num1];
   reIndex[num1] = reIndex[num2];
   reIndex[num2] = temp;

}



int partition(int left, int right, int pivot, int *rCount, int *rsIndex, int *reIndex) 
{
   int leftPointer = left -1;
   int rightPointer = right;

   while(true) {
      //while(rCount[++leftPointer] < pivot) {   //Sort low to high
      while(rCount[++leftPointer] > pivot) {     //sort high to low
         //do nothing
      }
		
      //while(rightPointer > 0 && rCount[--rightPointer] > pivot) { //Sort low to high
      while(rightPointer > 0 && rCount[--rightPointer] < pivot) {   //sort high to low
         //do nothing
      }

      if(leftPointer >= rightPointer) {
         break;
      } else {
         //printf(" item swapped :%d,%d\n", rCount[leftPointer],rCount[rightPointer]);
         swap(leftPointer, rightPointer, rCount, rsIndex, reIndex);
      }
   }
	
   //printf(" pivot swapped :%d,%d\n", rCount[leftPointer],rCount[right]);
   swap(leftPointer, right, rCount, rsIndex, reIndex);
   //printf("Updated Array: "); 
   return leftPointer;
}


void quickSortRowCount(int left, int right, int *rCount, int *rsIndex, int *reIndex)
{

   if(right-left <= 0) {
      return;   
   } else {
      int pivot = rCount[right];
      int partitionPoint = partition(left, right, pivot, rCount, rsIndex, reIndex);
      quickSortRowCount(left, partitionPoint-1, rCount, rsIndex, reIndex);
      quickSortRowCount(partitionPoint+1, right, rCount, rsIndex, reIndex);
   }        
}



void printSortedVal(const int M, const int N, const int nz, const int blockSize,
                    const int *rIndex,  const int *cIndex,
                    const int *rsIndex, const int *reIndex, const int *rCount,
                    float *val, float *vec)
{
	FILE *f;
	int i;

	f = fopen("sortedval.txt", "w");
        fprintf(f, "M= %d, N= %d, nz= %d\n\n", M, N, nz);
	fprintf(f, "row\trsIndex\treIndex\trCount\trIndex\tcIndex\t\tval\t\tvec\tproduct\n");
	fprintf(f, "------------------------------------------------------------------------------------------\n");

	for (i=0; i < nz; i++) 
	{
		if ( i<M) {
                   if ( (i % blockSize) == 0 ) fprintf(f, "\n");

		   //fprintf(f,"%d\t%d\t%d\t%d\t%d\t%d\t%12.7lg\t%12.7lg\t%12.7lg\n", i, 
		   //	rsIndex[i], reIndex[i], reIndex[i]-rsIndex[i]+1,
		   //	rIndex[i], cIndex[i], val[i], vec[i], val[i]*vec[cIndex[i]]);

		   fprintf(f,"%d\t%d\t%d\t%d\t%d\t%d\t%12.7lg\t%12.7lg\t%12.7lg\n", i, 
			rsIndex[i], reIndex[i], rCount[i],
			rIndex[i], cIndex[i], val[i], vec[i], val[i]*vec[cIndex[i]]);
		}
		else {
		   fprintf(f,"%d\t-- \t-- \t-- \t%d\t%d\t%12.7lg\t--\t--\n", i, 
			rIndex[i], cIndex[i], val[i]);
		}
	}

	fclose(f);
}



void printPrimeArrays(const int nnz, const int nChunk, const int threadsPerBlock, 
                      const int *rIdxPrime, const int *cIdxPrime, const float *valPrime)
{
	FILE *f;
	int i;

	f = fopen("prime.txt", "w");
        fprintf(f, "nnz= %d, nChunk= %d, threadsPerBlock= %d\n\n", nnz, nChunk, threadsPerBlock);
	fprintf(f, "row\trIdxPrime\tcIdxPrime\tvalPrime\n");
	fprintf(f, "--------------------------------------------------------\n");

	for (i=0; i < nnz; i++) 
	{
            if ( (i % threadsPerBlock) == 0 ) fprintf(f, "\n");

	    fprintf(f,"%d\t%d\t\t%d\t\t%lg\n", i, rIdxPrime[i], cIdxPrime[i], valPrime[i]); 
	}

	fclose(f);
}





void printInput(const int M, const int N, const int nz, const int *rIndex, const int *cIndex, float *val, float *vec)
{
	FILE *f;
	int i;

	f = fopen("input.txt", "w");
        fprintf(f, "M= %d, N= %d, nz= %d\n\n", M, N, nz);
	fprintf(f, "row\trIndex\tcIndex\t\tval\t\tvec\n");
	fprintf(f, "------------------------------------------------------------------------------------------\n");

	for (i=0; i < nz; i++) 
	{
		if ( i<M) {
		   fprintf(f,"%d\t%d\t%d\t%12.7lg\t%12.7lg\n", 
                        i, rIndex[i], cIndex[i], val[i], vec[i]);
		}
		else {
		   fprintf(f,"%d\t-- \t-- \t%12.7lg\t--\n", 
                        i, val[i]);
		}
	}

	fclose(f);
}




void print_p(const int M, const int N, const int nz, const int *rIndex,  float *p)
{
	FILE *f;
	int i;

	f = fopen("p.txt", "w");
        fprintf(f, "M= %d, N= %d, nz= %d\n\n", M, N, nz);
	fprintf(f, "row\trIndex\t\tproduct\n");
	fprintf(f, "------------------------------------------------------------------------------\n");

	for (i=0; i < nz; i++) 
	{
	   fprintf(f,"%d\t%d\t%12.7lg\n", i, rIndex[i], p[i]);
	}

	fclose(f);
}



void print_p2(const int M, const int N, const int nz, const int blocksize, const int *rIndex,  float *p)
{
	FILE *f;
	int i;

	f = fopen("p2.txt", "w");
        fprintf(f, "M= %d, N= %d, nz= %d\n\n", M, N, nz);
	fprintf(f, "row\trIndex\t\tsum after scan\n");
	fprintf(f, "------------------------------------------------------------------------------\n");

	for (i=0; i < nz; i++) 
	{
           if ( (i % blocksize) == 0 ) fprintf(f, "\n");
	   fprintf(f,"%d\t%d\t%12.7lg\n", i, rIndex[i], p[i]);
	}

	fclose(f);
}



void print_p3(const int M, const int N, const int nz, const int blocksize, const int *rIndex,  float *p)
{
	FILE *f;
	int i;

	f = fopen("p3.txt", "w");
        fprintf(f, "M= %d, N= %d, nz= %d\n\n", M, N, nz);
	fprintf(f, "row\trIndex\t\tpartial sum before scan\n");
	fprintf(f, "------------------------------------------------------------------------------\n");

	for (i=0; i < nz; i++) 
	{
           if ( (i % blocksize) == 0 ) fprintf(f, "\n");
	   fprintf(f,"%d\t%d\t%12.7lg\n", i, rIndex[i], p[i]);
	}

	fclose(f);
}



void print_p4(const int M, const int N, const int nz, const int blocksize, const int *rIndex,  float *p)
{
	FILE *f;
	int i;

	f = fopen("p4.txt", "w");
        fprintf(f, "M= %d, N= %d, nz= %d\n\n", M, N, nz);
	fprintf(f, "row\trIndex\t\tsum after scan\n");
	fprintf(f, "------------------------------------------------------------------------------\n");

	for (i=0; i < nz; i++) 
	{
           if ( (i % blocksize) == 0 ) fprintf(f, "\n");
	   fprintf(f,"%d\t%d\t%12.7lg\n", i, rIndex[i], p[i]);
	}

	fclose(f);
}



void print_val(const int M, const int N, const int nz, const int *rIndex,  float *p, float *v)
{
	FILE *f;
	int i;

	f = fopen("val.txt", "w");
        fprintf(f, "M= %d, N= %d, nz= %d\n\n", M, N, nz);
	fprintf(f, "row\trIndex\t\tproduct\t\tval\n");
	fprintf(f, "------------------------------------------------------------------------------\n");

	for (i=0; i < nz; i++) 
	{
	   fprintf(f,"%d\t%d\t%12.7lg\t%12.7lg\n", i, rIndex[i], p[i], v[i]);
	}

	fclose(f);
}



void printResult(const int M, float *res)
{
	FILE *f;
	int i;

	f = fopen("answer.txt", "w");
	fprintf(f, "row\tres\n");
	fprintf(f, "------------------------------------\n");

	for (i=0; i < M; i++) 
	{
		fprintf(f,"%d\t%lg\n", i, res[i]);
	}

	fclose(f);
}


void printOutput(const int M, float *res)
{
	FILE *f;
	int i;

	f = fopen("output.txt", "w");
	fprintf(f, "row\tres\n");
	fprintf(f, "------------------------------------\n");

	for (i=0; i < M; i++) 
	{
		fprintf(f,"%d\t%12.7lg\n", i, res[i]);
	}

	fclose(f);
}

