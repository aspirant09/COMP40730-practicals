#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#define MAXTHRDS 124

struct timeval start_time, end_time;
double elapsed_time;

double *A;
double *B;
double *R;

void max_col_sum_norm_of_product_serial(int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += A[i * n + k] * B[k * n + j];
            }
            R[j] = fmax(R[j], fabs(sum));
        }
    }
}

int main() {
    int n;

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

    gettimeofday(&start_time, NULL);

    max_col_sum_norm_of_product_serial(n);

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
