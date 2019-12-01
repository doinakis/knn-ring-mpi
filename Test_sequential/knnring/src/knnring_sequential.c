/*
*Sequential of kNN
*Doinakis Michail
*e-mail: doinakis@ece.auth.gr
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#include "../inc/knnring.h"

/*
  Version 0 Sequential Implementation
  Finds for each point in a query set Y the k  nearest neighbors in the corpus set X
*/
knnresult kNN(double * X, double * Y, int n , int m , int d , int k){
  knnresult result;                                 //the variable where we store the result

  result.nidx = (int *)malloc(m*k*sizeof(int));     // allocating space for the ids and the diastances of the neighbors of its point
  result.ndist = (double *)malloc(m*k*sizeof(double));
  result.m = m;                                     // setting the number of query points and the number of neighbors
  result.k = k;
  int *ids = (int *)malloc(n*sizeof(int));          //allocating space for the ids

  double *D = (double *)malloc(n*m*sizeof(double)); // m-by-n matrix for the distances


  double temp=0.0;                                  // helper variable

  // variables we are gonna need to call the blas routine
  double alpha , beta;
  int ld = d;
  alpha = -2.0;
  beta = 1.0;


  for(int i=0; i<n; i++){                           // initialising the ids
    ids[i] = i;
  }
  /*we can actually calculate and put the values of sum(Y.^2,2) and sum(X.^2,2).'
  in the same D array with the use of the temp variable we basically end up with Y*X' at the end of those calculations*/
  for(int i=0; i<m; i++){
    for(int j=0; j<d; j++){                        // calculating the power of every dimension of the corresponding point of Y
      temp = temp + pow(*(Y+i*d+j),2);             //sum(Y.^2,2) this part of matlab code
    }
    for(int k=0; k < n; k++){
      *(D+i*n+k) = temp;                           //copy the first column n times so we can do the summation of the arrays [m-by-n array]
    }
  temp = 0.0;
  }

  for(int i=0; i<n; i++){                          // calculating the power of every dimension of the corresponding point of X
    for(int j=0; j<d; j++){
      temp = temp + pow(*(X+i*d+j),2);             //sum(X.^2,2).' this part of the matlab code
    }
    for(int k=0; k < m; k++){
      *(D+k*n+i) = *(D+k*n+i) +temp;               //copy the first column m times so we can do the summation of the arrays [m-by-n array]
    }
    temp = 0.0;
  }

  // calling the cblas routine to calculate D<-alpha*Y*B + beta*D
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,m,n,d,alpha,Y,ld,X,ld,beta,D,n);//m-by-n!!
  //calculating the squere root of every element(distance) in the D array
  for(int i=0; i<m; i++){
      for(int j=0; j<n; j++){
        if(*(D+i*n+j) < 0.00000001)               // if the distacne is smaller tha 10^(-8) we set it to be 0
          *(D+i*n+j) = 0;
        *(D+i*n+j) = sqrt(*(D+i*n+j));
      }
  }
  /*for every raw of the arrays result.ndist we call the quickselect function to sort the elements of each raw
    and set the correct ids to each element */
  for(int i=0; i<m; i++){
    for (int j=0; j<k; j++){
      *(result.ndist + i*k + j) = quickselect((D+i*n),ids,n,j);
      *(result.nidx + i*k + j) = ids[j];
    }
    for(int j=0; j<n; j++){                       //reinitialize the ids for the next iteration
      ids[j]=j;
    }
  }
  //free the arrays we used for the calculation and return the result
  free(ids);
  free(D);
  return result;
}

/*The quickselect function where arr is the array with the distances of the points
ids is an array that helps us keep track of the ids of the points, cause it switches
the Indices at the same time we switch the distances in the quick select
n is the number of points and k is wich neighbor we want to calculate
 */
double quickselect(double *arr,int *ids, int n, int k) // where arr is the array n is the number of points and k is the k-th smallest element
{
  unsigned long i,ir,j,l,mid;
  double a,temp,b;

  l=0;
  ir=n-1;
  for(;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
	SWAP(arr[l],arr[ir]);
  SWAP(ids[l],ids[ir]);
      }
      return arr[k];
    }
    else {
      mid=(l+ir) >> 1;  //binary shift by 1; devides by 2 keeps the lowest int ,the other parts are like quicksort
      SWAP(arr[mid],arr[l+1]);
      SWAP(ids[mid],ids[l+1]);
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir]);
  SWAP(ids[l],ids[ir]);
      }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1],arr[ir]);
  SWAP(ids[l+1],ids[ir]);
      }
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l],arr[l+1]);
  SWAP(ids[l],ids[l+1]);
      }
      i=l+1;
      j=ir;
      a=arr[l+1];
      b=ids[l+1];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	SWAP(arr[i],arr[j]);
  SWAP(ids[i],ids[j]);
      }
      arr[l+1]=arr[j];
      ids[l+1]=ids[j];
      arr[j]=a;
      ids[j]=b;
      if (j >= k) ir=j-1;
      if (j <= k) l=i;
    }
  }
}
