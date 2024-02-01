/**
 * OPENMP implementation of matrix column sum norm
 * 
*/
#include<stdio.h>
#include<stdlib.h>
#include<omp.h>
#include <sys/time.h>
#include <math.h>


#define MAXTHRDS 124

struct timeval start_time, end_time;
double elapsed_time;

double *A;
double *B;
double *R;

void max_col_sum_norm_of_product(int num_threads, int n) {
    #pragma omp parallel num_threads(num_threads) reduction(max:R[:n])
    {
        int tid = omp_get_thread_num();
        int num_rows = n / num_threads;
        int start_row = tid * num_rows;
        int end_row = (tid == num_threads - 1) ? n : start_row + num_rows;

        for (int i = start_row; i < end_row; i++) {
            for (int j = 0; j < n; j++) {
                double sum = 0.0;
                for (int k = 0; k < n; k++) {
                    sum += A[i * n + k] * B[k * n + j];
                }
                R[j] = fmax(R[j], fabs(sum));
                // #pragma omp critical
                //     R[j]=R[j]>sum?R[j]:sum;
            }
            

        }

    }
}

int main() {
    int num_threads;
    int n;

    printf("Enter number of threads: ");
    if (scanf("%d", &num_threads) < 1 || num_threads > MAXTHRDS) {
        printf("Check input for the number of processors. Bye.\n");
        return -1;
    }

    printf("Enter matrix dimension: ");
    scanf("%d", &n);

    A = malloc(n * n * sizeof(double));
    B = malloc(n * n * sizeof(double));
    R = malloc(n * sizeof(double));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i * n + j] = rand();
            B[i * n + j] = rand();
        }
    }

    for (int i = 0; i < n; i++) {
        R[i] = 0.0;
    }

    
    // printf("Printing A.. \n");
    // for(int i=0;i<n;i++){
    //     for(int j=0;j<n;j++){
    //         printf("%lf ",A[i*n+j]);
    //     }
    //     printf("\n");
    // }
    // printf("Printing B.. \n");
    // for(int i=0;i<n;i++){
    //     for(int j=0;j<n;j++){
    //         printf("%lf ",B[i*n+j]);
    //     }
    //     printf("\n");
    // }

    gettimeofday(&start_time, NULL);

    max_col_sum_norm_of_product(num_threads, n);

    gettimeofday(&end_time, NULL);

    printf("Result:\n");
    for (int i = 0; i < n; i++) {
        printf("%lf ", R[i]);
    }
    
    elapsed_time = (end_time.tv_sec - start_time.tv_sec) +
                   (end_time.tv_usec - start_time.tv_usec) / 1.0e6;
    printf("\nElapsed time: %lf seconds\n", elapsed_time);

    free(A);
    free(B);
    free(R);

    return 0;
}
