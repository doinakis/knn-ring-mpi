/*
*Implementation of kNN sequential
*Doinakis Michail
*e-mail: doinakis@ece.auth.gr
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <sys/times.h>
#include <math.h>
#include <cblas.h>

#include "../inc/knnring.h"
#define N 1423// points
#define M 1423 // query data points
#define D 37 // dimensions
#define K 13 // neighbors

// helper variables to count the time that the code needs to execute
struct timeval startwtime,endwtime;
double p_time;






int main(){

    double *X=(double *)malloc(N*D*sizeof(double));
    double *Y=(double *)malloc(M*D*sizeof(double));

    //printf("--Corpus set X--\n");
    for(int i=0; i<N; i++){
      for(int j=0; j<D; j++){
        *(X+i*D+j)=(double) (rand())  / (double) RAND_MAX; //creating random double points;
      }
    }
    // printf("--Corpus set Y--\n");
     for(int i=0; i<M; i++){
       for(int j=0; j<D; j++){
         *(Y+i*D+j)=(rand() % 10001) / 10000.0; //creating random double points
       }
     }
    gettimeofday(&startwtime,NULL);
    knnresult a = kNN(X,Y,N,M,D,K);
    gettimeofday(&endwtime,NULL);
    p_time = (double)((endwtime.tv_usec-startwtime.tv_usec)/1.0e6+endwtime.tv_sec-startwtime.tv_sec);
    printf("Time %f\n",p_time);

    free(a.ndist);
    free(X);
    free(Y);
    return 0;
}
