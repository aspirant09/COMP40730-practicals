/****
 * 2D-Matrix-Matrix Multiplication using MPI
 * Task-1
 ****/

#include <stdio.h>
#include<stdlib.h>
#include<mpi.h>

#define SIZE 4                 /* SIZE should be a multiple of number of nodes*/
#define MASTER 1
#define SLAVE 2

MPI_Status status;

double a[N][N];
double b[N][N];
double c[N][N];
double b_transpose[N][N];

void initialization()
{
        int i,j;
        int count =1;

        for(i = 0;i<N;i++)
        {
                for(j = 0;j<N;j++)
                {

            a[i][j] = count;
            b[i][j] = count;
            count++;

                        c[i][j] = 0.0; /*Initially*/
                        b_transpose[j][i] = b[i][j]; /*Since it is easy to send rows than columns */
                }
        }

        printf("A is.... ");
        for(i = 0;i<N;i++)
        {
                for(j = 0;j<N;j++)
                {

            printf(" %lf ",a[i][j]);
                }
                printf("\n");
        }
        printf("B is.... ");
        for(i = 0;i<N;i++)
        {
                for(j = 0;j<N;j++)
                {

            printf(" %lf ",b[i][j]);
                }
                printf("\n");
        }
}

static void output()
{
        int i, j;
        printf("\n Final Output is:\n");
        for(i = 0;i<N;i++)
        {
                for(j = 0;j<N;j++)
                {
                        printf("%7.2f", c[i][j]);
                }
                printf("\n");
        }
}

int main(int argc, char **argv)
{
        int prank, p;
        int cols;                                                                                                                                                               /* number of columns per worker */
        int rows;                                                                                                                                                               /* number of rows per worker */
        int mtype;                                                                                                                                                              /* message type */
        int dest, src, offset1, offset2, new_offset1;
        double start_t, end_t;
        int i, j ,k;

        MPI_Init(&argc,&argv);
        MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &prank);

        if(p>=4)  /*Number of nodes greater than or equal to 4*/
        {
                if (prank == 0)                                                                                 /*Master Job */
                {
                        initialization();
                        printf("\nSIZE = %d, number of nodes = %d\n", N, p);
                        start_t = MPI_Wtime();
                        mtype = MASTER;
                        rows = (N/p)*2;                                                                                                                                     /* Number of rows in matrix A for each node*/
                        offset1 = rows;
                        cols = N/2;                                                                                                                                          /* Number of cols in matrix B which is same as number of rows in B-Transpose for each node*/
                        offset2 = cols;
                        new_offset1 = 0;

                        for(dest = 1; dest <p/2; dest++) /*Sending data for first set of slaves*/
                        {
                                MPI_Send(&offset1, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
                                MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
                                MPI_Send(&cols, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
                                MPI_Send(&a[offset1][0], rows*N, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);     /* STart element in A matrix */
                                MPI_Send(&b_transpose[0][0], cols*N, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD); /* number of cols in matrix b is equal to number of rows in b_transpose */
                                offset1 = offset1 + rows;
                        }

                        for(dest = p/2;dest<p;dest++) /*Sending data for second set of slaves*/
                        {
                                MPI_Send(&new_offset1, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
                                MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
                                MPI_Send(&offset2, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
                                MPI_Send(&cols, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
                                MPI_Send(&a[new_offset1][0], rows*N, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD); /* Starting point in A matrix */
                                MPI_Send(&b_transpose[offset2][0], cols*N, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD); /* number of cols in matrix b is equal to number of rows in b_transpose */
                                new_offset1 = new_offset1 + rows;
                        }

                        /* Master Calculation Part */
                        printf("$$$$");
                        printf("Current node rank is %d \n",prank);
                        for(i = 0; i<rows;i++)
                        {
                                
                                for(j=0;j< cols;j++)
                                {
                                                printf("%lf",a[i][j]);
                                }
                                printf("\n");
                        }

                        for(i = 0; i<rows; i++)
                        {
                                for(j = 0; j< cols; j++)
                                {
                                        for(k = 0; k <N; k++)
                                        {
                                                c[i][j] += a[i][k] * b_transpose[j][k];
                                        }
                                }
                        }


                        /* collect part-1 results from slaves */

                        mtype = SLAVE;
                        for(src = 1; src < p/2; src++)
                        {
                                MPI_Recv(&offset1, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
                                MPI_Recv(&rows, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
                                MPI_Recv(&cols, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
                                for(i = offset1;i<offset1+rows;i++)
                                {
                                        MPI_Recv(&c[i][0], cols, MPI_DOUBLE, src, mtype, MPI_COMM_WORLD, &status);
                                }
                        }
                        
                        mtype = SLAVE;
                        for(src = p/2;src<p;src++)
                        {
                                MPI_Recv(&new_offset1, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
                                MPI_Recv(&rows, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
                                MPI_Recv(&offset2, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
                                MPI_Recv(&cols, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
                                for(i = new_offset1;i<new_offset1+rows;i++)
                                {
                                        MPI_Recv(&c[i][offset2], cols, MPI_DOUBLE, src, mtype, MPI_COMM_WORLD, &status);
                                }
                        }


                        end_t = MPI_Wtime();
                }

                else if(prank>0 && prank<p/2) /* First set of slaves */
                {                                                                                                                                                                                               /* Receive data from Master */

                        mtype = MASTER;
                        MPI_Recv(&offset1, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
                        MPI_Recv(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
                        MPI_Recv(&cols, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
                        MPI_Recv(&a[offset1][0], rows*N, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &status);
                        MPI_Recv(&b_transpose[0][0], cols*N, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &status);

                        printf("*****");
                        printf("Current node rank is %d \n",prank);
                        for(i = offset1; i<offset1+rows;i++)
                        {
                                
                                for(j=0;j< cols;j++)
                                {
                                        
                                                printf("%lf",a[i][j]);
                                }
                                printf("\n");
                        }
                        /* slaves calculation part-1 */

                        for(i = offset1; i<offset1+rows;i++)
                        {
                                for(j=0;j< cols;j++)
                                {
                                        for(k=0;k< N;k++)
                                        {
                                                c[i][j] += a[i][k] * b_transpose[j][k];
                                        }
                                }
                        }

                        
                        /*Now sending  part-1 of the results back to the master */

                        mtype = SLAVE;
                        MPI_Send(&offset1, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
                        MPI_Send(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
                        MPI_Send(&cols, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
                        for(i = offset1;i<offset1+rows;i++)
                        {
                                MPI_Send(&c[i][0], cols, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD);
                        }
                        }

                else if(prank>=p/2) /*Second set of slaves */
                {
                        /*Receiving data from master*/

                        mtype = MASTER;
                        MPI_Recv(&new_offset1, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
                        MPI_Recv(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
                        MPI_Recv(&offset2, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
                        MPI_Recv(&cols, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
                        MPI_Recv(&a[new_offset1][0], rows*N, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &status);
                        MPI_Recv(&b_transpose[offset2][0], cols*N, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &status);


                        printf("###");
                        printf("Current node rank is %d \n",prank);
                        for(i = new_offset1; i<new_offset1+rows;i++)
                        {
                                
                                for(j=offset2;j< N;j++)
                                {
                                        
                                                printf("%lf",a[i][j]);
                                }
                                printf("\n");
                        }
                        /* Slaves calculations part-2 */

                        for(i=new_offset1;i<new_offset1 + rows;i++)
                        {
                                for(j=offset2;j<N;j++)
                                {
                                        for(k = 0; k<N;k++)
                                        {
                                                c[i][j] +=a[i][k] * b_transpose[j][k];
                                        }
                                }
                        }

                        /* Sending Part-2 of the results back to the master */

                        mtype = SLAVE;
                        MPI_Send(&new_offset1, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
                        MPI_Send(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
                        MPI_Send(&offset2, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
                        MPI_Send(&cols, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
                        for(i = new_offset1;i<new_offset1+rows;i++)
                        {
                                MPI_Send(&c[i][offset2], cols, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD);
                        }
                }
        }

        MPI_Finalize();
        if(prank==0){
            
            output();

            printf("\n Execution time on %d nodes is: %f ", p, end_t-start_t);
        }
        return 0;

}