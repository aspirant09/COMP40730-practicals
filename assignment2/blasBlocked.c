/**
cache blocked blas implmentation of matrix multiplication
*/
#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <sys/time.h>

struct timeval start_time, end_time;
double elapsed_time;

void initMat( int M, int N, double mat[], double val )
{
        int i, j;
        for (i= 0; i< M; i++)
                for (j= 0; j< N; j++)
                        mat[i*N+j] = val;
}

void cacheBlockedMatrixMultiplication(int n, double* A, double* B, double* C, int block_size) {
    int i, j, k;
    for (i = 0; i < n/block_size; i ++) {
        for (j = 0; j < n/block_size; j ++) {
            for (k = 0; k < n/block_size; k ++) {
              cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, block_size, block_size, block_size, 1.0, &A[(i * block_size*  n) + (k*block_size)], n, &B[(k * n * block_size)+ (j * block_size)], n, 1.0, &C[(i * n * block_size) + (j * block_size)], n);

            }
        }
    }
}

/*
void cacheBlockedMatrixMultiplication(int n, double *A, double *B, double *C,int b) {
  int i, j, l;

  // Divide the matrices into blocks.
  for (i = 0; i < n; i += b) {
    for (j = 0; j < n; j += b) {
      for (l = 0; l < n; l += b) {
        // Call blas_dgemm to multiply each pair of blocks.
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, b, b,
                   b, 1.0, &A[i * n + l], n, &B[l * n + j], n,
                   1.0, &C[i * n + j], n);
      }
    }
  }
}
*/
int main(void)
{
        double ALPHA=1.0, BETA=1.0;
        int n = 2048;
        double* A;
        double* B;
        double* C;
        int block_size=2;
        int i=0;
        int j=0;

        printf ("Please enter matrix dimension n : ");
        scanf("%d", &n);

        printf("Please enter cache block size : ");
        scanf("%d", &block_size);

        // allocate memory for the matrices
        // allocate memory for the matrices
        A = malloc (n*n*sizeof(double));
        B = malloc (n*n*sizeof(double));
        C = malloc (n*n*sizeof(double));
        if (!A || !B || !C)
        {
                printf( "Out of memory, reduce N value.\n");
                exit(-1);
        }

        // initialise matrices A and B
        initMat(n, n, A, 125);
        initMat(n, n, B, 64);

        gettimeofday(&start_time, NULL);
        //multuply matrices
        cacheBlockedMatrixMultiplication(n, A, B, C, block_size);
        gettimeofday(&end_time, NULL);
        elapsed_time = (end_time.tv_sec - start_time.tv_sec) +
               (end_time.tv_usec - start_time.tv_usec) / 1.0e6;
        printf("Done... \n");
        printf("Elapsed time: %lf seconds\n", elapsed_time);

        //deallocate memory
        free(A);
        free(B);
        free(C);

        return 0;
}
