/**
 * The following code implements a matrix multiplication 
 * of 2 matrices and then finding the absolute max column sum of the resultant
 * matrix.
 * 
 * This is implemented using pthreads.
 * Each thread receives a row of the first matrix as input, calculates 
 * the corresponding row of the resultant product matrix.
 * 
 * Each thread compares the resultant row with a global result which has 
 * absolute max column sum norm.
 * 
 * In this approach each thread receives a set of rows of the left matrix
 * at a time. The threads synchronize themselves by acquiring a lock over the globally
 * shared result set R
 * 
*/

#include<stdio.h>
#include<stdlib.h>
#include<pthread.h>
#include <sys/time.h>
#define MAXTHRDS 124
struct timeval start_time, end_time;
double elapsed_time;

typedef struct{
    double *sub_A;
    double *result;
    pthread_mutex_t *mutex;
    int num_rows;
    int length;
} subResult;

double *A;
double *B;
double *R;
pthread_mutex_t *mutex;
void *max_col_sum_norm_of_product( void *args){
    
    subResult *s;
    s=args;
    
    int n=s->length;
    for(int i=0;i<s->num_rows;i++){
        for(int j=0;j<n;j++){
                double sum=0.0;
                for(int k=0;k<n;k++){
                    sum+= s->sub_A[i*n+k]*B[k*n+j];
                }
                pthread_mutex_lock(s->mutex);
                R[j]= R[j]>sum?R[j]:sum;
                pthread_mutex_unlock(s->mutex);
                
        }
 
    }      
    pthread_exit(NULL);
}

int main(){
    int num_threads;
    int n;
    int num_rows;
    pthread_t *worker_thread;
    pthread_mutex_t *main_mutex;
    void *status;
    subResult *sub;

    printf("Enter number of threads: ");
    if(scanf("%d", &num_threads) < 1 || num_threads > MAXTHRDS) {
        printf("Check input for number of processors. Bye.\n");
        return -1;
    }
    printf("Enter matrix dimension: ");
    scanf("%d",&n);

    num_rows=n/num_threads;
    A=malloc(n*n*sizeof(double));
    B=malloc(n*n*sizeof(double));
    R=malloc(n*sizeof(double));
    worker_thread=malloc(num_threads*sizeof(pthread_t));
    main_mutex=malloc(sizeof(pthread_mutex_t));
    sub=malloc(num_threads*sizeof(subResult));
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            A[i*n +j]=rand();
            B[i *n + j]=rand();
        }
    }
    // A[2*n+3]=60000.0;
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
    for (int i = 0; i < n; i++) {
    R[i] = 0.0;
    }
    pthread_mutex_init(main_mutex,NULL);
    int k=0;
    gettimeofday(&start_time, NULL);

    for(int i=0;i<num_threads;i++){
            
            sub[i].sub_A= A+i*num_rows*n;
            sub[i].mutex=main_mutex;
            sub[i].num_rows=(i==num_threads-1)? 
                n%(num_threads)+num_rows: num_rows;
            sub[i].length=n;
            pthread_create(&worker_thread[i],NULL, max_col_sum_norm_of_product,
            (void*)&sub[i]);
 
    }
        

    for(int i=0;i<num_threads;i++)
            pthread_join(worker_thread[i],&status);

    gettimeofday(&end_time, NULL);
    
    printf("Result.. \n");
    for(int i=0;i<n;i++){
         printf("%lf ", R[i]);
    }
    elapsed_time = (end_time.tv_sec - start_time.tv_sec) +
               (end_time.tv_usec - start_time.tv_usec) / 1.0e6;
    printf("Elapsed time: %lf seconds\n", elapsed_time);
    free(A);
    free(B);
    free(worker_thread);
    free(sub);
    pthread_mutex_destroy(main_mutex);
}