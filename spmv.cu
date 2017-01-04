/*
*********************************************
*  314 Principles of Programming Languages  *
*  Fall 2016                                *
*********************************************
*
*   Read a real (non-complex) sparse matrix from a Matrix Market (v. 2.0) file
*   and a vector from a txt file, perform matrix multiplication and store the
*   result to output.txt. 
*
*
*
*   NOTES:
*
*   1) Matrix Market files are always 1-based, i.e. the index of the first
*      element of a matrix is (1,1), not (0,0) as in C.  ADJUST THESE
*      OFFSETS ACCORDINGLY offsets accordingly when reading and writing
*      to files.
*
*   2) ANSI C requires one to use the "l" format modifier when reading
*      double precision floating point numbers in scanf() and
*      its variants.  For example, use "%lf", "%lg", or "%le"
*      when reading doubles, otherwise errors will occur.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "mmio.h"
#include "utils.h"




///////////////////////////////////////////////////////////////////////////

typedef enum {
    ALG_ATOM = 0,
    ALG_SCAN,
    ALG_OPTSCAN,
    ALG_CACHE,
    ALG_DESIGN
} ALG_TYPE;



///////////////////////////////////////////////////////////////////////////


/*
__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                                          (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed, 
                        __double_as_longlong(val + 
                        __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}
*/


///////////////////////////////////////////////////////////////////////////


__global__ void spmv_atomic_kernel( const int nnz,
                                    const int *coord_row,
                                    const int *coord_col,
                                    const float *A,
                                    const float *x,
                                          float *y )
{
    int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
    int thread_num = blockDim.x * gridDim.x;
    int iter = nnz % thread_num ? nnz/thread_num + 1 : nnz/thread_num ;

    //
    //  y[rIndex[i]] += val[i] * vec[cIndex[i]] 
    //
    //
    // dataid    = offset to the "val" array
    // A[dataid] = value of the "val" array at row "dataid"
    // coord_row = rIndex array
    // coord_col = cIndex array
    // x[col]    = vector value for column "col"
    // y[row]    = result value for row "row"
    //
    for (int i = 0; i < iter; i++)
    {
        int dataid = thread_id + i * thread_num;
        if ( dataid < nnz ) {
            float data = A[ dataid ];
            int row = coord_row[ dataid ];
            int col = coord_col[ dataid ];
            float temp = data * x[ col ];
            atomicAdd( &y[row], temp );
        }
    }

    __syncthreads();
}



///////////////////////////////////////////////////////////////////////////


__device__ void segmented_scan( const int nz, int threadsPerBlock, const int *rows, float *vals)
{

    int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
/*
    int thread_num = blockDim.x * gridDim.x;
    int iter = nnz % thread_num ? nnz/thread_num + 1 : nnz/thread_num ;
*/

    if ( thread_id >= nz ) return;

    if ( threadsPerBlock > 32)  threadsPerBlock = 32;

    const int lane = threadIdx.x % threadsPerBlock;


    //
    // lane = thread offset in the thread warp
    // vals[] = multiplication results of one element 
    // rows[] = row indices
    //
    if (lane >= 1 && rows[thread_id] == rows[thread_id - 1] )
            vals[thread_id] += vals[thread_id - 1];

    if (lane >= 2 && rows[thread_id] == rows[thread_id - 2] )
            vals[thread_id] += vals[thread_id - 2];

    if (lane >= 4 && rows[thread_id] == rows[thread_id - 4] )
            vals[thread_id] += vals[thread_id - 4];

    if (lane >= 8 && rows[thread_id] == rows[thread_id - 8] )
            vals[thread_id] += vals[thread_id - 8];

    if (lane >= 16 && rows[thread_id] == rows[thread_id - 16] )
            vals[thread_id] += vals[thread_id - 16];

}



///////////////////////////////////////////////////////////////////////////


__global__ void add_partial_sum_kernel( const int nnz, const int threadsPerBlock,
                                 const int *rIndex, float *val, float *res)
{

    int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
/*
    int thread_num = blockDim.x * gridDim.x;
    int iter = nnz % thread_num ? nnz/thread_num + 1 : nnz/thread_num ;
*/


    if ( thread_id >= nnz ) return;

    //
    // Step 1: Perform segmented scan
    //
    segmented_scan(nnz, threadsPerBlock, rIndex, val);


    //
    // Step 2: Save all non-partial results
    //
    int row = rIndex[thread_id];

    if ( thread_id == nnz -1 ) {
        // This is the last element
        res[row] = val[thread_id];
    }
    else {
        if ( rIndex[thread_id] != rIndex[thread_id + 1] ) {
            // This is not a partial sum.  Save it
            res[row] = val[thread_id];
        }
    }
}



///////////////////////////////////////////////////////////////////////////


__global__ void spmv_scan_kernel( const int nnz, const int threadsPerBlock,
                                 const int *rIndex, const int *cIndex, 
                                 float *val, float *vec, float *res)
{

    int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
/*
    int thread_num = blockDim.x * gridDim.x;
    int iter = nnz % thread_num ? nnz/thread_num + 1 : nnz/thread_num ;
*/


    if ( thread_id >= nnz ) return;

    //
    // Step 1: Each thread calculate the product for one element 
    //
    int row = rIndex[thread_id];
    int col = cIndex[thread_id];

    val[thread_id] = val[thread_id] * vec[col];


    //
    // Step 2: Perform segmented scan
    //
    segmented_scan(nnz, threadsPerBlock, rIndex, val);


    //
    // Step 3: Save all non-partial results
    //
    if ( thread_id == nnz -1 ) {
        // This is the last element
        //
        //res[row] += val[thread_id];
        atomicAdd( &res[row], val[thread_id] );
    }
    else {
        if ( rIndex[thread_id] != rIndex[thread_id + 1] ) {
            // This is not a partial sum.  Save it
            //
            //res[row] += val[thread_id];
            atomicAdd( &res[row], val[thread_id] );
        }
    }
    //__syncthreads();
}



///////////////////////////////////////////////////////////////////////////


void getParams(int argc, char **argv, char **mfile, char **vfile, ALG_TYPE *alg, int *blocksize)
{
    int i;

    // Get all command line parameters
    //
    for ( i= 1; i < argc; i++) 
    {
        if ( strcmp(argv[i], "-mat") == 0 ) {
            // The next parameter is the matrix filename
            //
            *mfile = argv[++i];
            //printf("spmv: matrixfile= %s\n", *mfile);
            continue;
        } // End if "-mat"


        if ( strcmp(argv[i], "-ivec") == 0 ) {
            // The next parameter is the vector filename
            //
            *vfile = argv[++i];
            //printf("spmv: vectorfile= %s\n", *vfile);
            continue;
        } // End if "-ivec"


        if ( strcmp(argv[i], "-alg") == 0 ) {
            // The next parameter is the selected algorithm
            //
            i++;
            if ( strcmp(argv[i], "atom") == 0 ) {
                *alg = ALG_ATOM;
                continue;
            }
            else if ( strcmp(argv[i], "scan") == 0 ) {
                *alg = ALG_SCAN;
                continue;
            }
            else if ( strcmp(argv[i], "optscan") == 0 ) {
                *alg = ALG_OPTSCAN;
                continue;
            }
            else if ( strcmp(argv[i], "cache") == 0 ) {
                *alg = ALG_CACHE;
                printf("spmv: algorithm \"cache\" not implemented\n");
		exit(1);
                continue;
            }
            else {
                //
                // ALG_DESIGN is not implemented.  Default to atomic add.
                //

                *alg = ALG_DESIGN;
                //printf("spmv: algorithm \"design\" not implemented\n");
		//exit(1);

                *alg = ALG_ATOM;
                continue;
            }
        } // End if "-alg"


        if ( strcmp(argv[i], "-blocksize") == 0 ) {
            // The next parameter is the thread block size
            //
            *blocksize = atoi(argv[++i]);
            if ( *blocksize > 1024) {
                printf("spmv: blocksize cannot be larger than 1024\n");
                *blocksize = 1024;
            } 

            continue;
        } // End if "-blocksize"

    }

}


///////////////////////////////////////////////////////////////////////////


int Rollup_scan(const int nz, int threadsPerBlock, 
                int **prIndex, float **pval) 
{
    if ( nz <= 0 ) return 0;

    if ( threadsPerBlock > 32)  threadsPerBlock = 32;


    int new_nz = 0; // Partial sum count

    int *rIndex = *prIndex;
    float *val = *pval;

    int *ps_rIndex = (int *)malloc( nz * sizeof(int));
    float *ps_val  = (float *)malloc( nz * sizeof(float));

    //
    // There is a potential partial sum at every 
    // boundary of "threadsPerBlock".  For example,
    // if "threadsPerBlock" is 32, there is a potential
    // partial sum between row "31" and row "32+".
    //
    for ( int iBoundary= threadsPerBlock -1; iBoundary < nz; iBoundary += threadsPerBlock )
    {
        if ( rIndex[iBoundary] == rIndex[iBoundary + 1] ) {
            // This is the first partial sum value
            ps_rIndex[new_nz] = rIndex[iBoundary];
            ps_val[new_nz] = val[iBoundary];
            new_nz++;          


            // Look for the last partial sum value
            do
            {
                iBoundary++;
                if ( (iBoundary % threadsPerBlock) == threadsPerBlock - 1 ) 
                {
                    if (rIndex[iBoundary] == rIndex[iBoundary + 1]) {
                        // Another intermediate partial sum value
                        ps_rIndex[new_nz] = rIndex[iBoundary];
                        ps_val[new_nz] = val[iBoundary];
                        new_nz++;          
                    }
                }
            } while (rIndex[iBoundary] == rIndex[iBoundary + 1]);

            // This is the last partial sum value
            ps_rIndex[new_nz] = rIndex[iBoundary];
            ps_val[new_nz] = val[iBoundary];
            new_nz++;          

            // Re-align the iBoundary counter to the next segment
            iBoundary = iBoundary - (iBoundary % threadsPerBlock) - 1;
        }
    }

    free( rIndex );
    rIndex = (int *)realloc( ps_rIndex, new_nz * sizeof(int));
    *prIndex = rIndex;

    free( val );
    val = (float *)realloc( ps_val, new_nz * sizeof(float));
    *pval = val;

    return new_nz;
}



///////////////////////////////////////////////////////////////////////////


void createPrimeArrays(const int M, const int N, const int nz, 
                       const int threadsPerBlock, const int T,
                       int *rsIndex, int *reIndex, int *rCount,
                       int *rIndex, int *cIndex, 
                       float *val, float *vec,
                       int *nnz, int **p_rIdxPrime, int **p_cIdxPrime, float **p_valPrime,
                       int *nnz2, int **p_rIdxPrime2, int **p_cIdxPrime2, float **p_valPrime2)
{
    int i = 0;
    int tmp_nnz = 0;
    int nChunk = 32;

    int *rIdxPrime = NULL;
    int *cIdxPrime = NULL;
    float *valPrime = NULL;

    rIdxPrime = *p_rIdxPrime;
    cIdxPrime = *p_cIdxPrime;
    valPrime  = *p_valPrime;

    //
    // Step 1: prepare row, column and val prime arrays to be 
    //         used for segment scan.
    //
    while (nChunk >= T) {

        quickSortRowCount(0, M-1, rCount, rsIndex, reIndex);
        //printSortedVal(M, N, nz, threadsPerBlock, rIndex, cIndex, rsIndex, reIndex, rCount, val, vec);
    
        for (i=0; i<M; i++) {
            if ((rCount[i] > 0) && (rCount[i] >= nChunk)) {
                if ( rIdxPrime == NULL ) rIdxPrime = (int *)malloc(nz * sizeof(int));
                if ( cIdxPrime == NULL ) cIdxPrime = (int *)malloc(nz * sizeof(int));
                if ( valPrime == NULL ) valPrime = (float *)malloc(nz * sizeof(float));
    
                // move this many rows to the prime arrays
                int moveCount = (rCount[i] / nChunk) * nChunk;
    
                //printf("createPrimeArrays: i= %d, rCount[%d]= %d, nChunk= %d, moveCount= %d, nnz= %d\n",
                //           i, i, rCount[i], nChunk, moveCount, nnz); 
    
                for (int j= 1; j <= moveCount; j++) {
                    int srcStartIdx = rsIndex[i];
                    rIdxPrime[tmp_nnz] = rIndex[srcStartIdx];
                    cIdxPrime[tmp_nnz] = cIndex[srcStartIdx];
                    valPrime[tmp_nnz]  = val[srcStartIdx];
    
                    rsIndex[i] += 1;
                    rCount[i]  -= 1; 
                    tmp_nnz++;
                }
            }
        }

        nChunk = nChunk / 2;

    } // End of while


    if (rIdxPrime != NULL) rIdxPrime = (int *)realloc( rIdxPrime, tmp_nnz * sizeof(int));
    if (cIdxPrime != NULL) cIdxPrime = (int *)realloc( cIdxPrime, tmp_nnz * sizeof(int));
    if (valPrime  != NULL) valPrime = (float *)realloc( valPrime, tmp_nnz * sizeof(float));

    *p_rIdxPrime = rIdxPrime;
    *p_cIdxPrime = cIdxPrime;
    *p_valPrime  = valPrime;
    *nnz         = tmp_nnz;

    //printPrimeArrays(*nnz, nChunk, threadsPerBlock, rIdxPrime, cIdxPrime, valPrime);
    //printf("createPrimeArrays: nnz= %d\n", *nnz);


    //
    // Step 2: prepare row, column and val in another set of 
    //         prime arrays to be used for atomic add.
    //

    tmp_nnz = 0;
    rIdxPrime = *p_rIdxPrime2;
    cIdxPrime = *p_cIdxPrime2;
    valPrime  = *p_valPrime2;


    quickSortRowCount(0, M-1, rCount, rsIndex, reIndex);
    //printSortedVal(M, N, nz, threadsPerBlock, rIndex, cIndex, rsIndex, reIndex, rCount, val, vec);

    for (i=0; i<M; i++) {
        if (rCount[i] > 0) {
            if ( rIdxPrime == NULL ) rIdxPrime = (int *)malloc(nz * sizeof(int));
            if ( cIdxPrime == NULL ) cIdxPrime = (int *)malloc(nz * sizeof(int));
            if ( valPrime == NULL ) valPrime = (float *)malloc(nz * sizeof(float));

            // move this many rows to the prime arrays
            int moveCount = rCount[i];

            //printf("createPrimeArrays: i= %d, rCount[%d]= %d, nChunk= %d, moveCount= %d, nnz= %d\n",
            //           i, i, rCount[i], nChunk, moveCount, nnz); 

            for (int j= 1; j <= moveCount; j++) {
                int srcStartIdx = rsIndex[i];
                rIdxPrime[tmp_nnz] = rIndex[srcStartIdx];
                cIdxPrime[tmp_nnz] = cIndex[srcStartIdx];
                valPrime[tmp_nnz]  = val[srcStartIdx];

                rsIndex[i] += 1;
                rCount[i]  -= 1; 
                tmp_nnz++;
            }
        }
    }

    if (rIdxPrime != NULL) rIdxPrime = (int *)realloc( rIdxPrime, tmp_nnz * sizeof(int));
    if (cIdxPrime != NULL) cIdxPrime = (int *)realloc( cIdxPrime, tmp_nnz * sizeof(int));
    if (valPrime  != NULL) valPrime = (float *)realloc( valPrime, tmp_nnz * sizeof(float));

    *p_rIdxPrime2 = rIdxPrime;
    *p_cIdxPrime2 = cIdxPrime;
    *p_valPrime2  = valPrime;
    *nnz2         = tmp_nnz;

    //printf("createPrimeArrays: nnz2= %d\n", *nnz2);
    //printPrimeArrays(*nnz2, 32, threadsPerBlock, rIdxPrime, cIdxPrime, valPrime);
}

///////////////////////////////////////////////////////////////////////////


int main(int argc, char *argv[])
{
    char *matrixfile = NULL;
    char *vectorfile = NULL;
    ALG_TYPE  alg = ALG_DESIGN;
    int threadBlockSize = 128;

    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;   //M is row number, N is column number and nz is the number of non-zero entries
    int i, vecdim, *rIndex, *cIndex;
    float *val, *res, *vec;



    if ( (argc < 5) || ((argc % 2) != 1) ) 
    {
    	fprintf(stderr, "Usage: %s -mat [matrixfile] -ivec [vectorfile] -alg [atom|scan|optscan|cache|design] -blocksize [threadBlockSize]\n", argv[0]);
    	exit(1);
    }

    //
    // Get all the command line parameters
    //
    getParams(argc, argv, &matrixfile, &vectorfile, &alg, &threadBlockSize);


    //
    //Open and load the input matrix file  
    //
    //printf("\nOpening input matrix file: %s\n", matrixfile);
    if ((f = fopen(matrixfile, "r")) == NULL)
    {
    	printf("Fail to open the input matrix file!\n");
    	exit(1);
    }
    if (mm_read_banner(f, &matcode) != 0)
    {
    	printf("Could not process Matrix Market banner.\n");
    	exit(1);
    }

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */
    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
    	mm_is_sparse(matcode))
    {
    	printf("Sorry, this application does not support ");
    	printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
    	exit(1);
    }

    /* find out size of sparse matrix .... */
    ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz);
    if (ret_code != 0) exit(1);


    /* reserve memory for matrices */
    rIndex = (int *)malloc(nz * sizeof(int));
    cIndex = (int *)malloc(nz * sizeof(int));
    val = (float *)malloc(nz * sizeof(float));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */
    for (i = 0; i<nz; i++)
    {
    	fscanf(f, "%d %d %f\n", &rIndex[i], &cIndex[i], &val[i]);
    	rIndex[i]--;  /* adjust from 1-based to 0-based */
    	cIndex[i]--;
    }

    if (f != stdin) fclose(f);


    //
    //Open and load the input vector file 
    //
    //printf("Opening input vector file: %s\n", vectorfile);

    if ((f = fopen(vectorfile, "r")) == NULL)
    {
    	printf("Fail to open the input vector file!\n");
    	exit(1);
    }
    fscanf(f, "%d\n", &vecdim);
    if (vecdim != N)
    {
    	printf("dimension mismatch!\n");
    	exit(1);
    }
    vec = (float*)malloc(vecdim * sizeof(float));
    for (i = 0; i<vecdim; i++)
    {
    	fscanf(f, "%f\n", &vec[i]);
    }
    if (f != stdin) fclose(f);


    //
    // Sort the input maxtrix by row and then by column
    //
    int *rsIndex = NULL; 
    int *reIndex = NULL;
    int *rCount  = NULL;
    sortVals(val, M, N, rIndex, cIndex, nz, vec, &rsIndex, &reIndex, &rCount); 

    //printSortedVal(M, N, nz, threadBlockSize, rIndex, cIndex, rsIndex, reIndex, rCount, val, vec);



    res = (float*)malloc(M*sizeof(float));
    memset(res, 0, M*sizeof(float));
    

    //
    // Load the set of arrays into GPU device memory
    //
    int    *d_rIndex, *d_cIndex;
    float  *d_val, *d_vec, *d_res;

    cudaMalloc(&d_rIndex, nz * sizeof(int));
    cudaMalloc(&d_cIndex, nz * sizeof(int));
    cudaMalloc(&d_val, nz * sizeof(float));
    cudaMalloc(&d_vec, vecdim * sizeof(float));
    cudaMalloc(&d_res, M * sizeof(float));

    cudaMemcpy( d_rIndex, rIndex, nz*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy( d_cIndex, cIndex, nz*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy( d_val, val, nz*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy( d_vec, vec, vecdim*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy( d_res, res, M*sizeof(float), cudaMemcpyHostToDevice);


    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    
    cudaDeviceProp deviceProp;
    int device;
    for (device = 0; device < deviceCount; ++device) {
        cudaGetDeviceProperties(&deviceProp, device);

        //printf("Device %d, \"%s\", has compute capability %d.%d.\n",
        //       device, deviceProp.name, deviceProp.major, deviceProp.minor);
        //printf("Device %d, multiProcessorCount=  %d\n", device, deviceProp.multiProcessorCount);
        //printf("Device %d, isMultiGpuBoard=  %d\n", device, deviceProp.isMultiGpuBoard);
        //printf("Device %d, warpSize=  %d\n", device, deviceProp.warpSize);
        //printf("Device %d, maxThreadsPerBlock=  %d\n", device, deviceProp.maxThreadsPerBlock);
        //printf("Device %d, maxThreadsDim[0]=  %d\n", device, deviceProp.maxThreadsDim[0]);
        //printf("Device %d, maxGridSize[0]=  %d\n", device, deviceProp.maxGridSize[0]);
    }

    if ( cudaSetDevice(0) != cudaSuccess ) {
        printf("spmv: failed to set current GPU to device 0\n");
        exit(1);
    }



    //
    // Set the GPU launch configuration
    //
    int threadsPerBlock = threadBlockSize;
    int blocksPerGrid = (nz + threadsPerBlock - 1) / threadsPerBlock; 
    

    // Configuration settings used by opt-scan
    int T = 8;

    int nnz = 0;
    int *rIdxPrime = NULL;
    int *cIdxPrime = NULL;
    float *valPrime  = NULL;

    int nnz2 = 0;
    int *rIdxPrime2 = NULL;
    int *cIdxPrime2= NULL;
    float *valPrime2  = NULL;

    if ( alg == ALG_OPTSCAN ) {
        //
        // Prepare two separate set of prime arrays to be
        // used for calculation.  The first set of arrays
        // is used for "segment scan".  The second set of
        // arrays is used for "atomic add".
        //

        // Configure threshold to stop segment scan after groups of "T" 
        T = 16;

        createPrimeArrays(M, N, nz, threadsPerBlock, T,
                          rsIndex, reIndex, rCount,
                           rIndex, cIndex, val, vec,
                           &nnz, &rIdxPrime, &cIdxPrime, &valPrime,
                           &nnz2, &rIdxPrime2, &cIdxPrime2, &valPrime2);

        //printf("main: nz= %d, nnz= %d, nnz2= %d\n", nz, nnz, nnz2);
        //printPrimeArrays(nnz, 32, threadsPerBlock, rIdxPrime, cIdxPrime, valPrime);
    }

    // Prepare cuda events for time measurements
    cudaEvent_t startTime, stopTime;
    cudaEventCreate(&startTime);
    cudaEventCreate(&stopTime);
    float elapsedTime = 0;

    //
    // Start recording time measurement
    //
    cudaEventRecord(startTime, 0);



    /**********************************************************/
    /* Start the spmarse matrix vector multiplication         */
    /**********************************************************/

    switch (alg) {
        case ALG_ATOM:

            //printf(" spmv_atomic_kernel: M= %d, N= %d, nz= %d, threadsPerBlock= %d, blocksPerGrid= %d\n",
            //          M, N, nz, threadsPerBlock, blocksPerGrid);

            spmv_atomic_kernel<<< blocksPerGrid, threadsPerBlock >>> (nz, d_rIndex, d_cIndex, d_val, d_vec, d_res);

            cudaMemcpy( res, d_res, M*sizeof(float), cudaMemcpyDeviceToHost);
            break;

        case ALG_SCAN:

            //printf(" spmv_scan_kernel: M= %d, N= %d, nz= %d, threadsPerBlock= %d, blocksPerGrid= %d\n",
            //          M, N, nz, threadsPerBlock, blocksPerGrid);

            spmv_scan_kernel<<<blocksPerGrid, threadsPerBlock>>>
                     (nz, threadsPerBlock, d_rIndex, d_cIndex, d_val, d_vec, d_res);

            cudaMemcpy( val, d_val, nz*sizeof(float), cudaMemcpyDeviceToHost);

            nz = Rollup_scan(nz, threadsPerBlock, &rIndex, &val);

            while (nz > 0 ) 
            {
                cudaFree(&d_rIndex);
                cudaMalloc(&d_rIndex, nz * sizeof(int));
                cudaMemcpy( d_rIndex, rIndex, nz*sizeof(int), cudaMemcpyHostToDevice);

                cudaFree(&d_val);
                cudaMalloc(&d_val, nz * sizeof(float));
                cudaMemcpy( d_val, val, nz*sizeof(float), cudaMemcpyHostToDevice);

                add_partial_sum_kernel<<<blocksPerGrid, threadsPerBlock>>>
                            ( nz, threadsPerBlock, d_rIndex, d_val, d_res);

                cudaMemcpy( val, d_val, nz*sizeof(float), cudaMemcpyDeviceToHost);

                nz = Rollup_scan(nz, threadsPerBlock, &rIndex, &val);
            }
 
            cudaMemcpy( res, d_res, M*sizeof(float), cudaMemcpyDeviceToHost);
            break;

        case ALG_OPTSCAN:

            if ( nnz > 0 ) {
                //printf("main: nnz= %d, do segment scan\n",nnz);

                // Do the segment scan with the reduced prime arrays
                cudaFree(&d_rIndex);
                cudaMalloc(&d_rIndex, nnz * sizeof(int));
                cudaMemcpy( d_rIndex, rIdxPrime, nnz*sizeof(int), cudaMemcpyHostToDevice);

                cudaFree(&d_cIndex);
                cudaMalloc(&d_cIndex, nnz * sizeof(int));
                cudaMemcpy( d_cIndex, cIdxPrime, nnz*sizeof(int), cudaMemcpyHostToDevice);

                cudaFree(&d_val);
                cudaMalloc(&d_val, nnz * sizeof(float));
                cudaMemcpy( d_val, valPrime, nnz*sizeof(float), cudaMemcpyHostToDevice);

                spmv_scan_kernel<<<blocksPerGrid, threadsPerBlock>>>
                         (nnz, threadsPerBlock, d_rIndex, d_cIndex, d_val, d_vec, d_res);

                if (rIdxPrime != NULL) { free(rIdxPrime); rIdxPrime= NULL; }
                if (cIdxPrime != NULL) { free(cIdxPrime); cIdxPrime= NULL; }
                if (valPrime  != NULL) { free(valPrime);  valPrime = NULL; }

                //cudaMemcpy( res, d_res, M*sizeof(float), cudaMemcpyDeviceToHost);
                //printResult(M, res);
            }


            if ( nnz2 > 0 ) {
                //printf("main: nnz2= %d, do atomic add\n",nnz2);

                // Do the atomic add scan with the reduced prime arrays
                cudaFree(&d_rIndex);
                cudaMalloc(&d_rIndex, nnz2 * sizeof(int));
                cudaMemcpy( d_rIndex, rIdxPrime2, nnz2 * sizeof(int), cudaMemcpyHostToDevice);

                cudaFree(&d_cIndex);
                cudaMalloc(&d_cIndex, nnz2 * sizeof(int));
                cudaMemcpy( d_cIndex, cIdxPrime2, nnz2 * sizeof(int), cudaMemcpyHostToDevice);

                cudaFree(&d_val);
                cudaMalloc(&d_val, nnz2 * sizeof(float));
                cudaMemcpy( d_val, valPrime2, nnz2 * sizeof(float), cudaMemcpyHostToDevice);

                spmv_atomic_kernel<<< blocksPerGrid, threadsPerBlock >>> 
                           (nnz2, d_rIndex, d_cIndex, d_val, d_vec, d_res);

                if (rIdxPrime2 != NULL) { free(rIdxPrime2); rIdxPrime2= NULL; }
                if (cIdxPrime2 != NULL) { free(cIdxPrime2); cIdxPrime2= NULL; }
                if (valPrime2  != NULL) { free(valPrime2);  valPrime2 = NULL; }

                //cudaMemcpy( res, d_res, M*sizeof(float), cudaMemcpyDeviceToHost);
                //printResult(M, res);
            }

            cudaMemcpy( res, d_res, M*sizeof(float), cudaMemcpyDeviceToHost);
            //printResult(M, res);

            break;

        case ALG_CACHE:
            break;

        case ALG_DESIGN:
        default:
            break;
    }


    // Stop event
    cudaEventRecord(stopTime, 0);
    cudaEventSynchronize(stopTime);

    cudaEventElapsedTime(&elapsedTime, startTime, stopTime); // that's our time!
    printf(" The total kernel running time on GPU [%s] is %8.3f milli-seconds\n", deviceProp.name, elapsedTime);


    //
    // save the result in an "output.txt" file
    //
    if ((f = fopen("output.txt", "w")) == NULL)
    {
    	printf("Fail to open the output file!\n");
    	exit(1);
    }
    for (i = 0; i<M; i++)
    {
    	fprintf(f, "%f\n", res[i]);
    }
    fclose(f);

    free(res);
    free(vec);
    free(rIndex);
    free(cIndex);
    free(val);
    free(rsIndex);
    free(reIndex);
    free(rCount);

    cudaFree(&d_rIndex);
    cudaFree(&d_cIndex);
    cudaFree(&d_val);
    cudaFree(&d_vec);
    cudaFree(&d_res);

    // Clean up:
    cudaEventDestroy(startTime);
    cudaEventDestroy(stopTime);

    return 0;
}


