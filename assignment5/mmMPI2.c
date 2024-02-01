#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <mpi.h>
#include<math.h>


int n, p;

int main(int argc, char **argv) {
    int myn, myrank;
    double *local_a, *local_b, *local_c, *allB, start, sum, *allC, sumdiag;
    int i, j, k;

    n = atoi(argv[1]);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    myn = n / sqrt(p);
    local_a = malloc(myn * n * sizeof(double));
    local_b = malloc(myn * n * sizeof(double));
    local_c = malloc(myn * n * sizeof(double));
    allB = malloc(n * n * sizeof(double));
    allA = malloc(n * n * sizeof(double));

    for (i = 0; i < myn * n; i++) {
        for(int j=0;j<myn;j++){
            local_a[i*n+j] = 1.0;
            local_b[i*n+j] = 2.0;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (myrank == 0)
        start = MPI_Wtime();

    // for (i = 0; i < p; i++){
    //     MPI_Gather(local_b, myn * n, MPI_DOUBLE, allB, myn * n, MPI_DOUBLE, i, MPI_COMM_WORLD);
    //     MPI_Gather(a, myn * n, MPI_DOUBLE, allB, myn * n, MPI_DOUBLE, i, MPI_COMM_WORLD);
    // }
        

    for (i = 0; i < myn; i++)
    {
        // MPI_Gather(local_b, myn * n, MPI_DOUBLE, allB, myn * n, MPI_DOUBLE, i, MPI_COMM_WORLD);
        // MPI_Gather(a, myn * n, MPI_DOUBLE, allA, myn * n, MPI_DOUBLE, i, MPI_COMM_WORLD);
        MPI_Sendrecv(localData, localRows * localCols, MPI_DOUBLE, dest, 0, receivedData, localRows * localCols, MPI_DOUBLE, source, 0, comm, MPI_STATUS_IGNORE);
        for (j = 0; j < n; j++) {
            sum = 0.0;
            for (k = 0; k < n; k++)
                sum += local_a[i * n + k] * allB[k * n + j];
            local_c[i * n + j] = sum;
        }

    }
        

    free(allB);

    MPI_Barrier(MPI_COMM_WORLD);

    if (myrank == 0)
        printf("It took %f seconds to multiply 2 %dx%d matrices.\n", MPI_Wtime() - start, n, n);

    if (myrank == 0)
        allC = malloc(n * n * sizeof(double));

    MPI_Gather(local_c, myn * n, MPI_DOUBLE, allC, myn * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (myrank == 0) {
        for (i = 0; i < n; i++){
            for(int j=0;j < n;j++)
                printf(" %lf",allc[i*n+j])
        }
        printf("The trace of the resulting matrix is %f\n", sumdiag);
    }

    if (myrank == 0)
        free(allC);

    MPI_Finalize();
    free(local_a);
    free(local_b);
    free(local_c);

    return 0;
}


// void exchangeData(int *localA, int *localB, int *receivedA, int *receivedB, int localRows, int rank, int size, MPI_Comm comm) {
//     int source, dest;

//     // Exchange rows with the processor below
//     dest = (rank + 1) % size;
//     source = (rank - 1 + size) % size;

//     MPI_Sendrecv(localA, localRows * n, MPI_INT, dest, 0, receivedA, localRows * n, MPI_INT, source, 0, comm, MPI_STATUS_IGNORE);

//     // Exchange columns with the processor to the right
//     dest = (rank + 1) % size;
//     source = (rank - 1 + size) % size;

//     MPI_Sendrecv(localB, localRows * n, MPI_INT, dest, 0, receivedB, localRows * n, MPI_INT, source, 0, comm, MPI_STATUS_IGNORE);
// }


void exchangeData(double *localData, double *receivedData, int localRows, int localCols, int rank, int size, MPI_Comm comm, int isColumn) {
    int dest;
    int source;

    if (isColumn) {
        // Column-wise communication
        dest = (rank + 1) % sqrt(size);        // Send to the processor to the right
        source = (rank - 1 + sqrt(size)) % sqrt(size); // Receive from the processor to the left

        MPI_Sendrecv(localData, localRows * localCols, MPI_DOUBLE, dest, 0, receivedData, localRows * localCols, MPI_DOUBLE, source, 0, comm, MPI_STATUS_IGNORE);
    } else {
        // Row-wise communication
        dest = (rank + sqrt(size)) % size;     // Send to the processor below
        source = (rank - sqrt(size) + size) % size; // Receive from the processor above

        MPI_Sendrecv(localData, localRows * localCols, MPI_DOUBLE, dest, 0, receivedData, localRows * localCols, MPI_DOUBLE, source, 0, comm, MPI_STATUS_IGNORE);
    }
}

