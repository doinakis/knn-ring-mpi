/*
*knn.h header with the accessors
*Doinakis Michail
*e-mail: doinakis@ece.auth.gr
*/
#ifndef KNN_H
#define KNN_H
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp; // swaps the values of a and b

//! Computes k nearest neighbors of each point in X [n-by-d]

/*
\param ndix  Indices (0-based) of nearest neighbors [m-by-k]
\param  ndist   Distance of nearest neighbors[m-by-k]
\param  m       Number of query points[scalar]
\param  k       Number of nearest neighbors[scalar]
*/

typedef struct knnresult{
  int *nidx;
  double *ndist;
  int m;
  int k;
}knnresult;

//! Computes k nearest neighbors of each point in X [n-by-d]
/*
\param X  Corpus data points[n-by-d]
\param  Y   Query data points[m-by-d]
\param  n   Number of corpus points[scalar]
\param  m   Number of query points[scalar]
\param  d   Number of dimensions[scalar]
\param  k   Number of neighbors[scalar]
\return The kNN result*/

knnresult kNN(double * X, double * Y, int n , int m , int d , int k);

//! Computes distributed all-kNN of points in X
/*!\param X  Data points[n-by-d]
\param  n Number of data points[scalar]
\param  d Number of dimensions[scalar]
\param  k Number of neighbors[scalar]
\return The kNN result*/

knnresult distrAllkNN(double *X,int n,int d,int k);

/*Modified quickselect from the previous assignment
\param  arr Array with the points we want to find the kth smallest element
\param  ids Array that contains the ids of the elements
\param  n   Number of elements in the array
\param  k the kth smallest element we want to find
\return the kth smallest element*/

double quickselect(double *arr,int * ids, int n, int k);
#endif
