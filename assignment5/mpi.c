#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N 64

#define MASTER 0
#define SLAVE 1

MPI_Status status;

double a[N][N];
double b[N][N];
double c[N][N];
double b_transpose[N][N];
double start_t;

void output() {
    
}

int main(int argc, char **argv) {
    int prank, np, destination, source, row_offset1, offset2, row_offset2;
    int rows, columns; // columns and rows per partition
    int mtype;
    double start_t;
    int i, j, k;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &prank);

    if (prank == 0) {
        int i, j;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                a[i][j] = 3.0;
                b[i][j] = 2.0;
                c[i][j] = 0.0;
                b_transpose[j][i] = b[i][j];
            }
        }
        start_t = MPI_Wtime();
        mtype = MASTER;
        rows = (N / np) * 2;
        row_offset1 = rows;
        columns = N / 2;
        offset2 = N / 2;
        row_offset2 = 0;

/**
 * Sending for firs half of matrices. They come below the root node in the grid
 * offset 0 belongs to root
*/
        for (destination = 1; destination < np / 2; destination++) {
            MPI_Send(&row_offset1, 1, MPI_INT, destination, mtype, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_INT, destination, mtype, MPI_COMM_WORLD);
            MPI_Send(&columns, 1, MPI_INT, destination, mtype, MPI_COMM_WORLD);
            MPI_Send(&a[row_offset1][0], rows * N, MPI_DOUBLE, destination, mtype, MPI_COMM_WORLD);
            MPI_Send(&b_transpose[0][0], columns * N, MPI_DOUBLE, destination, mtype, MPI_COMM_WORLD);
            row_offset1 = row_offset1 + rows;
        }

/**
 * Sending offsets for the second half of cluster. They come on the right side of the root node.
*/
        for (destination = np / 2; destination < np; destination++) {
            MPI_Send(&row_offset2, 1, MPI_INT, destination, mtype, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_INT, destination, mtype, MPI_COMM_WORLD);
            MPI_Send(&offset2, 1, MPI_INT, destination, mtype, MPI_COMM_WORLD);
            MPI_Send(&columns, 1, MPI_INT, destination, mtype, MPI_COMM_WORLD);
            MPI_Send(&a[row_offset2][0], rows * N, MPI_DOUBLE, destination, mtype, MPI_COMM_WORLD);
            MPI_Send(&b_transpose[offset2][0], columns * N, MPI_DOUBLE, destination, mtype, MPI_COMM_WORLD);
            row_offset2 = row_offset2 + rows;
        }
//master multiplication
        for (i = 0; i < rows; i++) {
            for (j = 0; j < columns; j++) {
                for (k = 0; k < N; k++) {
                    c[i][j] += a[i][k] * b_transpose[j][k];
                }
            }
        }

//receiving /gathering response from the worker nodes
/**
 * Receives from first set of nodes. data is identified based on the offset. The same wil be replicatd below in reverse
 * when worker nodes send the data to master
*/
        mtype = SLAVE;
        for (source = 1; source < np / 2; source++) {
            MPI_Recv(&row_offset1, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&columns, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
            for (i = row_offset1; i < row_offset1 + rows; i++) {
                MPI_Recv(&c[i][0], columns, MPI_DOUBLE, source, mtype, MPI_COMM_WORLD, &status);
            }
        }

/**
 * Receicing from submatrices on the RHS
*/
        mtype = SLAVE;
        for (source = np / 2; source < np; source++) {
            MPI_Recv(&row_offset2, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&offset2, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&columns, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
            for (i = row_offset2; i < row_offset2 + rows; i++) {
                MPI_Recv(&c[i][offset2], columns, MPI_DOUBLE, source, mtype, MPI_COMM_WORLD, &status);
            }
        }
    } else if (prank > 0 && prank < np / 2) {   //Worker nodes logic(worker in the first half)
    // SAME CODE FROM MASTER . WHERE MASTER SENT, WORKER RECEIVES
        mtype = MASTER;
        MPI_Recv(&row_offset1, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&columns, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&a[row_offset1][0], rows * N, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&b_transpose[0][0], columns * N, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &status);

        
        for (i = row_offset1; i < row_offset1 + rows; i++)
            for (j = 0; j < columns; j++)
                for (k = 0; k < N; k++)
                    c[i][j] += a[i][k] * b_transpose[j][k];

        //SENDING OFFSETS
        mtype = SLAVE;
        MPI_Send(&row_offset1, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
        MPI_Send(&columns, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
        //Send all rows of the submatrix of C the worker processed.
        for (i = row_offset1; i < row_offset1 + rows; i++) {
            MPI_Send(&c[i][0], columns, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD);
        }
    } else if (prank >= np / 2) { //WORKER ON THE SECOND HALF OF THE VERTICAL SPLIT
        mtype = MASTER;
        MPI_Recv(&row_offset2, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&offset2, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&columns, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&a[row_offset2][0], rows * N, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&b_transpose[offset2][0], columns * N, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &status);

        for (i = row_offset2; i < row_offset2 + rows; i++)
            for (j = offset2; j < N; j++)
                for (k = 0; k < N; k++)
                    c[i][j] += a[i][k] * b_transpose[j][k];

        mtype = SLAVE;
        MPI_Send(&row_offset2, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
        MPI_Send(&offset2, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
        MPI_Send(&columns, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
        for (i = row_offset2; i < row_offset2 + rows; i++) {
            MPI_Send(&c[i][offset2], columns, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();
    if (prank == 0) {
        int i, j;
    printf("\n Final Output is:\n");
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("%7.2f", c[i][j]);
        }
        printf("\n");
    }
        printf("\n Execution time on %d nodes is: %f ", np, MPI_Wtime() - start_t);
    }
    return 0;
}
