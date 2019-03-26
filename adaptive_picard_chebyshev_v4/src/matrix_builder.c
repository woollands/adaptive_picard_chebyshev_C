/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com) (based on a similar code by Brent Macomber)
*  DATE WRITTEN:     May 2017
*  LAST MODIFIED:    May 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      One time build & store constant matrices required for the Adaptive Picard-Chebyshev numerical integration method
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "clenshaw_curtis_ivpII.h"
#include "c_functions.h"
#include "matrix_loader.h"

int N, M;

// Initialize Storage Arrays
double MAT_T2[(Nmax-Nmin+1)][(Nmax+1)*(Nmax+1)];
double MAT_P2[(Nmax-Nmin+1)][(Nmax+1)*(Nmax+1)];
double MAT_T1[(Nmax-Nmin+1)][(Nmax+1)*(Nmax+1)];
double MAT_P1[(Nmax-Nmin+1)][(Nmax+1)*(Nmax+1)];
double MAT_Ta[(Nmax-Nmin+1)][(Nmax+1)*(Nmax+1)];
double MAT_A[(Nmax-Nmin+1)][(Nmax+1)*(Nmax+1)];

int main(){

  printf("Building constant matrices...\n");

  for (int i=Nmin; i<=Nmax; i++){
    N = i;
    M = N;
    printf("N = %i\n",N);

    // Compute Clenshaw-Curtis Quadrature Constant Matrices
    double T2[(M+1)*(N+1)];
    memset( T2, 0.0, ((M+1)*(N+1)*sizeof(double)));
    double P2[(N+1)*N];
    memset( P2, 0.0, ((N+1)*N*sizeof(double)));
    double T1[(M+1)*N];
    memset( T1, 0.0, ((M+1)*N*sizeof(double)));
    double P1[N*(N-1)];
    memset( P1, 0.0, (N*(N-1)*sizeof(double)));
    double Ta[(M+1)*(N-1)];
    memset( Ta, 0.0, ((M+1)*(N-1)*sizeof(double)));
    double A[(N-1)*(M+1)];
    memset( A, 0.0, ((N-1)*(M+1)*sizeof(double)));
    clenshaw_curtis_ivpII(N,M,T2,P2,T1,P1,Ta,A);

    // Build & Store Arrays
    for (int j=1; j<=M+1; j++){
      for (int k=1; k<=N+1; k++){
        MAT_T2[i-Nmin][ID2(j,k,Nmax+1)] = T2[ID2(j,k,M+1)];
      }
    }
    for (int j=1; j<=N+1; j++){
      for (int k=1; k<=N; k++){
        MAT_P2[i-Nmin][ID2(j,k,Nmax+1)] = P2[ID2(j,k,N+1)];
      }
    }
    for (int j=1; j<=M+1; j++){
      for (int k=1; k<=N; k++){
        MAT_T1[i-Nmin][ID2(j,k,Nmax+1)] = T1[ID2(j,k,M+1)];
      }
    }
    for (int j=1; j<=N; j++){
      for (int k=1; k<=N-1; k++){
        MAT_P1[i-Nmin][ID2(j,k,Nmax+1)] = P1[ID2(j,k,N)];
      }
    }
    for (int j=1; j<=M+1; j++){
      for (int k=1; k<=N-1; k++){
        MAT_Ta[i-Nmin][ID2(j,k,Nmax+1)] = Ta[ID2(j,k,M+1)];
      }
    }
    for (int j=1; j<=N-1; j++){
      for (int k=1; k<=M+1; k++){
        MAT_A[i-Nmin][ID2(j,k,Nmax+1)] = A[ID2(j,k,N-1)];
      }
    }

  }

  // Open Files
  FILE* fT2 = fopen("../matrices/T2_matrices.bin","wb");
  FILE* fP2 = fopen("../matrices/P2_matrices.bin","wb");
  FILE* fT1 = fopen("../matrices/T1_matrices.bin","wb");
  FILE* fP1 = fopen("../matrices/P1_matrices.bin","wb");
  FILE* fTa = fopen("../matrices/Ta_matrices.bin","wb");
  FILE* fA  = fopen("../matrices/A_matrices.bin","wb");

  // Confirm Opening
  if ( !fT2 ){
    printf("Failure to open fT2 for binary write: CHECK PATH\n");
  }
  if ( !fP2 ){
    printf("Failure to open fP2 for binary write: CHECK PATH\n");
  }
  if ( !fT1 ){
    printf("Failure to open fT1 for binary write: CHECK PATH\n");
  }
  if ( !fT1 ){
    printf("Failure to open fT1 for binary write: CHECK PATH\n");
  }
  if ( !fTa ){
    printf("Failure to open fTa for binary write: CHECK PATH\n");
  }
  if ( !fA ){
    printf("Failure to open fA for binary write: CHECK PATH\n");
  }

  printf("Saving constant matrices... ");

  // Write Binary Data
  fwrite( MAT_T2, sizeof(double), (Nmax-Nmin+1)*(Nmax+1)*(Nmax+1), fT2 );
  fwrite( MAT_P2, sizeof(double), (Nmax-Nmin+1)*(Nmax+1)*(Nmax+1), fP2 );
  fwrite( MAT_T1, sizeof(double), (Nmax-Nmin+1)*(Nmax+1)*(Nmax+1), fT1 );
  fwrite( MAT_P1, sizeof(double), (Nmax-Nmin+1)*(Nmax+1)*(Nmax+1), fP1 );
  fwrite( MAT_Ta, sizeof(double), (Nmax-Nmin+1)*(Nmax+1)*(Nmax+1), fTa );
  fwrite( MAT_A, sizeof(double), (Nmax-Nmin+1)*(Nmax+1)*(Nmax+1), fA );

  // Close Files
  fclose( fT2 );
  fclose( fP2 );
  fclose( fT1 );
  fclose( fP1 );
  fclose( fTa );
  fclose( fA );

  printf("Complete!\n");

}
