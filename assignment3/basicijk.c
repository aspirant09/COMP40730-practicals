/**
basic ijk implementation of matrix multiplication without any cache blocking.
*/
#include <stdio.h>
#include <stdlib.h>
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

void Multiply(int n, double* a, double* b, double* c)
{
        int i=0;
        int j=0;
        int r=0;
        int x=0,k=0,sum=0;



        for(i = 0; i < n; i++) {
                for( j = 0; j < n; j++){
                        sum=0;
                        for( r = 0; r < n; r++)
                                sum += a[i * n + r] * b[r * n + j];
                        c[i * n + j]= sum;
                }
        }

        /*
        int x=0;
        int k=0;
        for (k=0;k<n;k++){
                for (j=0;j<n;j++){
                        x= b[k * n + j];
                        for(i=0;i<n;i++)
                                c[i * n + j] = a[i * n + k] * x;
                }
        }



        for(k=0;k<n;k++){
                for(i=0;i<n;i++){
                        x=a[i * n + k];
                        for(j=0;j<n;j++)
                                c[i * n + j]= x + b[j * n + k];
                }
        }
        */



}

int main(void)
{
        int n=1000;
        double* A;
        double* B;
        double* C;
        double *R;

        int i=0;
        int j=0;

        printf ("Please enter matrix dimension n : ");
        scanf("%d", &n);


        // allocate memory for the matrices
        A = malloc (n*n*sizeof(double));
        B = malloc (n*n*sizeof(double));
        C = malloc (n*n*sizeof(double));
        R = malloc (n*sizeof(double));
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
        Multiply(n,A,B,C);
        gettimeofday(&end_time, NULL);
        for (int i = 0; i < n; i++)
            R[i] = 0.0;
        for(j=0;j<n;j++){
                for(i=0;i<n;i++){
                        R[i]= R[i]>C[i*n+j]?R[i]:C[i*n+j];
                }
        }
        printf("Result.. \n");
        for(int i=0;i<n;i++){
            printf("%lf ", R[i]);
        }
        elapsed_time = (end_time.tv_sec - start_time.tv_sec) +
               (end_time.tv_usec - start_time.tv_usec) / 1.0e6;
        printf("Elapsed time: %lf seconds\n", elapsed_time);
        


        //deallocate memory
        free(A);
        free(B);
        free(C);
        free(R);
        return 0;
}
