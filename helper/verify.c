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
#include <math.h>
#include <float.h>
#include <sys/time.h>
#include "mmio.h"
#include "utils.h"




///////////////////////////////////////////////////////////////////////////

int getRowsInMaxtrix(const char *matrixfile)
{
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;   //M is row number, N is column number and nz is the number of non-zero entries


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
        fclose(f);
    	exit(1);
    }

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */
    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
    	mm_is_sparse(matcode))
    {
    	printf("Sorry, this application does not support ");
    	printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        fclose(f);
    	exit(1);
    }

    /* find out size of sparse matrix .... */
    mm_read_mtx_crd_size(f, &M, &N, &nz);

    //printf("getRowsInMaxtrix: M= %d, N= %d, nz= %d\n", M, N, nz);
    fclose(f);
    return M;
}

///////////////////////////////////////////////////////////////////////////


void createAnswerFile( const char *matrixfile, const char *answerfile ) 
{
    char vectorfile[128] = "";

    MM_typecode matcode;
    FILE *f;
    int M, N, nz;   //M is row number, N is column number and nz is the number of non-zero entries
    int i, vecdim, *rIndex, *cIndex;
    float *val, *vec;


    i = strlen(matrixfile) - 4;  // Remove the extension, ".mtx"

    strcpy( vectorfile, matrixfile);
    vectorfile[i] = '\0';  // Remove the ".mtx" extension
    strcat( vectorfile, "_vector.txt");

    //
    //Open and load the input matrix file  
    //
    //printf("createAnswerFile: Opening input matrix file: %s\n", matrixfile);
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
    if (  mm_read_mtx_crd_size(f, &M, &N, &nz) != 0 )  exit(1); 

    //printf("createAnswerFile: M= %d, N= %d, nz= %d\n", M, N, nz);



    /* reserve memory for matrices */
    rIndex = (int *)malloc(nz * sizeof(int));
    cIndex = (int *)malloc(nz * sizeof(int));
    val = (float *)malloc(nz * sizeof(float));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */
    for (i = 0; i<nz; i++)
    {
    	//fscanf(f, "%d %d %lg\n", &rIndex[i], &cIndex[i], &val[i]);
    	fscanf(f, "%d %d %f\n", &rIndex[i], &cIndex[i], &val[i]);
    	rIndex[i]--;  /* adjust from 1-based to 0-based */
    	cIndex[i]--;
    }

    if (f != stdin) fclose(f);  //Close the matrix file
    //printf("createAnswerFile: Closed input matrix file: %s\n", matrixfile);


    //
    //Open and load the input vector file 
    //
    //printf("createAnswerFile: Opening input vector file: %s\n", vectorfile);

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
    //printf("createAnswerFile: Closed input vector file: %s\n", vectorfile);




    //
    // The original calculation result.  This is
    // the correct answer to the output vector.
    //
    float * res_seq = (float*)malloc(M*sizeof(float));
    memset(res_seq, 0, M*sizeof(float));

    getmul(val, vec, rIndex, cIndex, nz, res_seq);


    //
    // Save the answer to the answer file.
    //
    //
    if ((f = fopen(answerfile, "w")) == NULL)
    {
    	printf("Fail to open the answer file \"%s\" for write\n", answerfile);
    	exit(1);
    }

    for (i = 0; i<M; i++)
    {
    	fprintf(f, "%f\n", res_seq[i]);
        //if (i==0) printf("createAnswerFile: wrote answer[0]= %lg\n", res_seq[0]);
    }

    if (f != stdin) fclose(f);

    free(res_seq);
}


///////////////////////////////////////////////////////////////////////////

static boolean nearlyEqual(float a, float b, float epsilon) {

    float absA = fabs(a);
    float absB = fabs(b);
    float diff = fabs(a - b);


    if (a == b) { // shortcut, handles infinities
        return TRUE;

    } else if (a == 0 || b == 0 || diff < FLT_MIN) {
        // a or b is zero or both are extremely close to it
        // relative error is less meaningful here
        //
        if ( diff < (epsilon * FLT_MIN )) {
            return TRUE;

        } else { 
            //printf("diff >= (epsilon * FLT_MIN)\n");
            //printf("epsilon= %lg, FLT_MIN= %lg, absA= %f, absB= %f, diff= %f\n",
            //         FLT_EPSILON, FLT_MIN, absA, absB, diff);

            return FALSE; 
        }

    } else { 
         // use relative error
        if ( diff / (absA + absB) < epsilon ) {
            return TRUE;
        } else { 
            //printf("diff / (absA + absB) >= epsilon\n");
            //printf("epsilon= %lg, FLT_MIN= %lg, absA= %lg, absB= %lg, diff= %lg, absA+absB= %lg, diff/(absA+absB)= %lg\n",
            //FLT_EPSILON, FLT_MIN, absA, absB, diff, absA+absB, diff/(absA+absB));

            return FALSE; 
        }
    }

}


///////////////////////////////////////////////////////////////////////////


int main(int argc, char *argv[])
{
    char matrixfile[128] = "";
    char answerfile[128] = "";
    char outputfile[128] = ""; 

    int i;
    FILE *f_output, *f_answer;


    if (argc < 3) 
    {
    	fprintf(stderr, "Usage: %s [matrixfile] [outputfile]\n", argv[0]);
    	exit(1);
    }

    //
    // Get all the command line parameters
    //
    strcpy( matrixfile, argv[1] );
    strcpy( outputfile, argv[2] );

    i = strlen(matrixfile) - 4;  
    strcpy( answerfile, matrixfile);
    answerfile[i] = '\0';  // Remove the ".mtx" extension
    strcat( answerfile, "_answer.txt");


    // Check the results.  Compare the "output.txt"
    // file with the correct answer. 
    //
    if ((f_output = fopen(outputfile, "r")) == NULL)
    {
    	printf("Fail to open the output file \"%s\"\n", outputfile);
    	exit(1);
    }
    //printf("verify: outputfile= %s\n",outputfile);


    //printf("verify: trying to open answerfile, %s\n",answerfile);
    if ((f_answer = fopen(answerfile, "r")) == NULL)
    {
        createAnswerFile( matrixfile, answerfile );
    	//printf("Created an answer file \"%s\"\n", answerfile);
    	
        if ((f_answer = fopen(answerfile, "r")) == NULL)
        {
    	    printf("Fail to open the answer file \"%s\"\n", answerfile);
    	    exit(1);
        }
    }
    //printf("verify: answerfile= %s\n",answerfile);

    int M = getRowsInMaxtrix(matrixfile);


    float *val = (float *)malloc( M * sizeof(float));
    float *ans = (float *)malloc( M * sizeof(float));

    for (i=0; i<M; i++) {
        fscanf(f_output, "%f\n", &val[i]);
        fscanf(f_answer, "%f\n", &ans[i]);
        //if (i==0) printf("main: read from answer file - ans[0]= %lg\n", ans[0]);
    }
        
    fclose(f_output);
    fclose(f_answer);
 


    float epsilon = 0.005;

    for (i=0; i<M; i++) {
        //if (i==0) printf("verify: val[0]= %lg, ans[0]= %lg\n", val[0], ans[0]);

        if ( nearlyEqual(ans[i], val[i], epsilon) == FALSE ) {
            float diff;
            diff = ans[i] - val[i];
            if ( diff < 0 ) diff = diff * (-1);

            printf("Calculation Error: calculated value[%d]= %lg, correct answer[%d]= %lg, diff= %lg\n",
                       i, val[i], i, ans[i], diff); 
            break;
        }
    }

    free(val);
    free(ans);
}


