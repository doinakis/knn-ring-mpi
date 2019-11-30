/*
*Implementation of kNN
*Doinakis Michail
*e-mail: doinakis@ece.auth.gr
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include "mpi.h"
#include "../inc/knnring.h"


/*
  V2. Asynchronous Implementation
  Computes distributed all-kNN of points in X
*/
knnresult distrAllkNN(double *X,int n,int d,int k){
int p,id;                                                                       // p is the number of processes and id is the id of the process
knnresult result;                                                               //Variable where we are gonna store the result

MPI_Comm_rank(MPI_COMM_WORLD,&id);                                              //Initialising the id of the process
MPI_Comm_size(MPI_COMM_WORLD,&p);                                               //Initialising the number of processes
MPI_Status Stat_rec,Stat_send;
int tag = 1;
MPI_Request mpireq_rec,mpireq_send;
/*if the receive tag is less than zero then the process has id=0 and need to send to the process with p-1 id
if the receive tag is more than the number of processes then this process sends to the process with 0 id*/

int rcv = id - 1;
int dst = id + 1;
if(dst > p-1) dst = 0;
if(rcv < 0) rcv = p - 1;
double *corpus = (double *)malloc(n*d*sizeof(double));                          //Allocating memory for the corpus we are going to receive
double *corpustosend = (double *)malloc(n*d*sizeof(double));                    //Allocating memory for the corpus we are goind to send

/*Using the MPI_Isend and MPI_Irecv in order to do hide the cost
of communications by transering the data while we do the computations*/
MPI_Isend(X,n*d,MPI_DOUBLE,dst,tag,MPI_COMM_WORLD,&mpireq_send);
MPI_Irecv(corpus,n*d,MPI_DOUBLE,rcv,tag,MPI_COMM_WORLD,&mpireq_rec);

result = kNN(X,X,n,n,d,k);                                                      //Store the result of the calculation

/* Every process corresponds to a part of a big array
we map the ids to point at the whole array
the elements of process 0 goes last */
for(int i=0; i<n; i++){
  for(int j=0; j<k; j++){
      *(result.nidx+i*k+j) = *(result.nidx+i*k+j) + ((id-1+p)%p)*n;
  }
}

//wait for the send and the receive to end before we procceed
MPI_Wait(&mpireq_send,&Stat_send);
MPI_Wait(&mpireq_rec,&Stat_rec);

/*The processes with even number of id first send and then receive
while the processes with odd number of id first receive and then send
we are gonna need to make p-1 iterations to have all the corpuses moved to everyone */

for(int iter=0; iter<p-1; iter++){
  for(int i=0; i<n; i++){
    for(int j=0; j<d; j++){
      *(corpustosend+i*d+j) = *(corpus+i*d+j);
    }
  }
  /*since we have already completed the first send and receive we only need to send
  p-2 times more*/

  if(iter != p-2){
    MPI_Isend(corpustosend,n*d,MPI_DOUBLE,dst,tag,MPI_COMM_WORLD,&mpireq_send);
    MPI_Irecv(corpus,n*d,MPI_DOUBLE,rcv,tag,MPI_COMM_WORLD,&mpireq_rec);
  }
  //starting the computations without waiting

  knnresult comps = kNN(corpustosend,X,n,n,d,k);

  //setting the correct ids depending on which proccess it was and how many times it received(iter)

  for(int i=0; i<n; i++){
    for(int j=0; j<k; j++){
        *(comps.nidx+i*k+j) = *(comps.nidx+i*k+j) + ((id-2-iter+p)%p)*n;
    }
  }
  /*Checking for every element in the i row of the comps.dist if its bigger than the last element of the result.ndist
  if it is then there is no other element to check since the comps.dist are sorted so we break the loop
  if its smaller then we set the last element to be the smallest and then sort again the array using quickselect*/

  for(int i=0; i<n; i++){
    for(int j=0; j<k; j++){
        if(*(comps.ndist+i*k+j) >= *(result.ndist+i*k+k-1)){
            break;
        }else{
            *(result.ndist+i*k+k-1) = *(comps.ndist+i*k+j);
            *(result.nidx+i*k+k-1) = *(comps.nidx+i*k+j);
            for(int h=0;h<k;h++){
            quickselect((result.ndist+i*k),(result.nidx+i*k),k,h);
            }
          }
        }
      }
  //we need to wait before we proceed to the next send and receive

  MPI_Wait(&mpireq_send,&Stat_send);
  MPI_Wait(&mpireq_rec,&Stat_rec);
}
free(corpustosend);
free(corpus);
//Global reductions and all-to-all
//
// double min_distance,max_distance;                                               //variables that every process have
// min_distance = RAND_MAX;                                                        //setting the minimum value to the biggest random number
// max_distance = 0;                                                               //setting the maximum to 0 cause we have distances that are less or equal to 0
// /*each process compares its results and find the minimum
// and the maximum value*/
// for(int i=0; i<n; i++){
//   for(int j=0; j<k; j++){
//     if(min_distance > *(result.ndist + i*k + j)){
//       min_distance = *(result.ndist + i*k + j);
//     }
//     if(max_distance < *(result.ndist + i*k + j)){
//       max_distance = *(result.ndist + i*k + j);
//     }
//   }
// }
// //variables where we want to store the global minimum and maximum
// double max,min;
// /*call the mpi reduce for the min_distance and max_distance and we store
// its rerurn values in the min max variables*/
// MPI_Reduce(&min_distance,&min,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
// MPI_Reduce(&max_distance,&max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
// // the MPI_Reduce returns the values to the process with id 0
// if( id==0 ){
//   printf("MAX : %f \n MIN : %f \n",max,min ); // and then we can print them
// }
return result;
}



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

  for(int i=0; i<m; i++){
    for(int j=0; j<d; j++){                        // calculating the power of every dimension of the corresponding point of Y
      temp = temp + pow(*(Y+i*d+j),2);
    }
    for(int k=0; k < n; k++){
      *(D+i*n+k) = temp;
    }
  temp = 0.0;
  }

  for(int i=0; i<n; i++){                          // calculating the power of every dimension of the corresponding point of X
    for(int j=0; j<d; j++){
      temp = temp + pow(*(X+i*d+j),2);
    }
    for(int k=0; k < m; k++){
      *(D+k*n+i) = *(D+k*n+i) +temp;
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

/*The quickselect function where arr is the array with the distances of the pints
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
