#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <sys/time.h>

struct timeval start_time, end_time;
double elapsed_time;

void initMat( int M, int N, double mat[], double val )
{
        int     i, j;
        for (i= 0; i< M; i++)
                for (j= 0; j< N; j++)
                        mat[i*N+j] = val;
}

int main(void)
{
        double ALPHA=1.0, BETA=1.0;
        int n = 2048;
        double* A;
        double* B;
        double* C;
        int numreps = 5;

        int i=0;
        int j=0;

        printf ("Please enter matrix dimension n : ");
        scanf("%d", &n);

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
        cblas_dgemm( CblasRowMajor,  CblasNoTrans, CblasNoTrans, n, n, n,
                                        ALPHA, A, n, B, n, BETA, C, n );
        gettimeofday(&end_time, NULL);
        elapsed_time = (end_time.tv_sec - start_time.tv_sec) +
               (end_time.tv_usec - start_time.tv_usec) / 1.0e6;
        printf("Done ...\n");


        for(i=0;i<n;i++){
                printf("\n");
                for(j=0;j<n;j++){
                        printf("%lf ",C[i*n+j]);
                }
        }

        printf("Elapsed time: %lf seconds\n", elapsed_time);

        //deallocate memory
        free(A);
        free(B);
        free(C);

        return 0;
}
