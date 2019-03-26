/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     May 2017
*  LAST MODIFIED:    May 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Propagates the initial conditions
*
* INPUT:
*    r0            -- Initial position vector (km)
*    v0            -- Initial velocity vector (km/s)
*    t0            -- Initial time (s)
*    t_final       -- Final time (s)
*    deg           -- Gravity Degree (max 100)
*    tol           -- Tolerance
*    Period        -- Period (s)
*    tvec          -- Segment start and end times (s)
*    t_orig        -- Segment start and end times for first segment (s)
*    seg           -- Segments per orbit
*    N             -- Polynomial degree
*    M             -- Sample points
*    prep_HS       -- Hot start switch function
*    coeff_size    -- Length of coefficient array
*    soln_size     -- Size of solution array
*    P1            -- First integration operator
*    P2            -- Second integration operator
*    T1            -- Chebyshev position matrix
*    T2            -- Chebyshev velocity matrix
*    A             -- Least squares operator
*    Ta            -- Chebyshev acceleration matrix
*    W1            -- Time scale factor 1
*    W2            -- Time scale factor 2
*    Feval         -- Function evaluation counter
*
* OUTPUTS:
*    total_seg     -- Total segments
*    ALPHA         -- Position coefficients
*    BETA          -- Velocity coefficients
*    segment_times -- Array of segment start and end times
*/

#include "picard_chebyshev_propagator.h"
#include "prepare_propagator.h"
#include "picard_iteration.h"
#include "FandG.h"
#include "reosc_perigee.h"
#include "c_functions.h"

void picard_chebyshev_propagator(double* r0, double* v0, double t0, double t_final,double deg, double tol, double Period,
   double* tvec, double* t_orig, int seg, int N, int M, int* prep_HS, int coeff_size, int soln_size, int* total_seg,
   double* P1, double* P2, double* T1, double* T2, double* A, double* Ta, double* W1, double* W2, double* Feval,
   double* ALPHA, double* BETA, double* segment_times){

  int loop    = 0;      // Break loop condition
  int k       = 0;      // Counter: segments per orbit
  int hot     = 0;      // Hot start switch
  int seg_cnt = 0;      // Counter: total segments
  int sz      = (int) ceil(1.1*t_final/Period)*seg;
  double w1, w2, tf;

  double HotX[seg*(M+1)*3];
  memset( HotX, 0.0, (seg*(M+1)*3*sizeof(double)));
  double HotV[seg*(M+1)*3];
  memset( HotV, 0.0, (seg*(M+1)*3*sizeof(double)));

  // PROPAGATION
  while (loop == 0){

    // Compute cosine time vector for a given segment
    t0 = tvec[k];
    tf = tvec[k+1];
    while (tf == 0.0){
      k = k+1;
      t0 = tvec[k];
      tf = tvec[k+1];
    }
    if (tf > t_final){
      tf = t_final;
    }
    w1 = (tf + t0)/2.0;
    w2 = (tf - t0)/2.0;
    W1[seg_cnt] = w1;
    W2[seg_cnt] = w2;

    double z0[6] = {0.0};
    double z[6] = {0.0};
    z0[0] = r0[0]; z0[1] = r0[1]; z0[2] = r0[2];
    z0[3] = v0[0]; z0[4] = v0[1]; z0[5] = v0[2];

    double tau[M+1];
    memset( tau, 0.0, ((M+1)*sizeof(double)));
    double times[M+1];
    memset( times, 0.0, ((M+1)*sizeof(double)));
    double X[(M+1)*3];
    memset( X, 0.0, ((M+1)*3*sizeof(double)));
    double V[(M+1)*3];
    memset( V, 0.0, ((M+1)*3*sizeof(double)));
    double Beta[N*3];
    memset( Beta, 0.0, (N*3*sizeof(double)));
    double Alpha[(N+1)*3];
    memset( Alpha, 0.0, ((N+1)*3*sizeof(double)));

    // KEPLERIAN WARM START
    for (int cnt=0; cnt<=M; cnt++){
      tau[cnt]   = -cos(cnt*C_PI/M);
      times[cnt] = tau[cnt]*w2 + w1;
      FandG(z0,z,times[cnt]-t0);
      X[ID2(cnt+1,1,M+1)] = z[0];
      X[ID2(cnt+1,2,M+1)] = z[1];
      X[ID2(cnt+1,3,M+1)] = z[2];
      V[ID2(cnt+1,1,M+1)] = z[3];
      V[ID2(cnt+1,2,M+1)] = z[4];
      V[ID2(cnt+1,3,M+1)] = z[5];
    }
    // Warm Start
    double WSX[(M+1)*3];
    memset( WSX, 0.0, ((M+1)*3*sizeof(double)));
    double WSV[(M+1)*3];
    memset( WSV, 0.0, ((M+1)*3*sizeof(double)));
    memcpy(WSX,X,(M+1)*3*sizeof(double));
    memcpy(WSV,V,(M+1)*3*sizeof(double));

    // HOT START (after 1+ orbits)
    // if (hot == 1){
    //   for (int i = 1; i<=M+1; i++){
    //     X[ID2(i,1,M+1)] = X[ID2(i,1,M+1)] + HotX[ID2(i+(k*(M+1)),1,seg*(M+1))];
    //     X[ID2(i,2,M+1)] = X[ID2(i,2,M+1)] + HotX[ID2(i+(k*(M+1)),2,seg*(M+1))];
    //     X[ID2(i,3,M+1)] = X[ID2(i,3,M+1)] + HotX[ID2(i+(k*(M+1)),3,seg*(M+1))];
    //     V[ID2(i,1,M+1)] = V[ID2(i,1,M+1)] + HotV[ID2(i+(k*(M+1)),1,seg*(M+1))];
    //     V[ID2(i,2,M+1)] = V[ID2(i,2,M+1)] + HotV[ID2(i+(k*(M+1)),2,seg*(M+1))];
    //     V[ID2(i,3,M+1)] = V[ID2(i,3,M+1)] + HotV[ID2(i+(k*(M+1)),3,seg*(M+1))];
    //   }
    // }

    // PICARD ITERATION
    picard_iteration(r0,v0,X,V,times,N,M,deg,hot,tol,P1,P2,T1,T2,A,Feval,Alpha,Beta);

    // Loop exit condition
    if (fabs(tf - t_final)/tf < 1e-12){
      loop = 1;
    }

    // Prepare Hot Start
    // if (*prep_HS == -1){
    //   for (int i=1; i<=M+1; i++){
    //     HotX[ID2(i+(k*(M+1)),1,seg*(M+1))] = X[ID2(i,1,M+1)] - WSX[ID2(i,1,M+1)];
    //     HotX[ID2(i+(k*(M+1)),2,seg*(M+1))] = X[ID2(i,2,M+1)] - WSX[ID2(i,2,M+1)];
    //     HotX[ID2(i+(k*(M+1)),3,seg*(M+1))] = X[ID2(i,3,M+1)] - WSX[ID2(i,3,M+1)];
    //     HotV[ID2(i+(k*(M+1)),1,seg*(M+1))] = V[ID2(i,1,M+1)] - WSV[ID2(i,1,M+1)];
    //     HotV[ID2(i+(k*(M+1)),2,seg*(M+1))] = V[ID2(i,2,M+1)] - WSV[ID2(i,2,M+1)];
    //     HotV[ID2(i+(k*(M+1)),3,seg*(M+1))] = V[ID2(i,3,M+1)] - WSV[ID2(i,3,M+1)];
    //   }
    //
    //   if (k == seg-1){
    //     hot = 1;
    //   }
    // }

    // Assign new initial Conditions (next segment)
    for (int j=1; j<=3; j++){
      r0[j-1] = X[ID2(M+1,j,M+1)];
      v0[j-1] = V[ID2(M+1,j,M+1)];
    }

    /* REOSCULATE PERIGEE
    Reosculate Keplerian perigee after each orbit propagation. If this is not
    done then the precomputed segment times no longer align with true segment
    breaks (due to perturbations) and the precomputed segment scheme fails to
    produce a solution that satisfies the required tolerance. This effect
    increases with increasing eccentricity. */
    double orb_end = 0.0;
    reosc_perigee(X,V,times,Alpha,Beta,tf,t_final,t_orig,N,M,&k,seg,prep_HS,tol,&orb_end,tvec,r0,v0);

    // Segments per orbit counter
    k = k+1;

    // STORE TRAJECTORY COEFFICIENTS
    for (int i=1; i<=N; i++){
      BETA[ID2(i+(seg_cnt*N),1,coeff_size)] = Beta[ID2(i,1,N)];
      BETA[ID2(i+(seg_cnt*N),2,coeff_size)] = Beta[ID2(i,2,N)];
      BETA[ID2(i+(seg_cnt*N),3,coeff_size)] = Beta[ID2(i,3,N)];
    }
    for (int i=1; i<=N+1; i++){
      ALPHA[ID2(i+seg_cnt*(N+1),1,coeff_size)] = Alpha[ID2(i,1,N+1)];
      ALPHA[ID2(i+seg_cnt*(N+1),2,coeff_size)] = Alpha[ID2(i,2,N+1)];
      ALPHA[ID2(i+seg_cnt*(N+1),3,coeff_size)] = Alpha[ID2(i,3,N+1)];
    }

    segment_times[seg_cnt+1] = tf;
    if (orb_end != 0.0){
      segment_times[seg_cnt+1] = orb_end;
    }

    // Total segments counter
    seg_cnt = seg_cnt + 1;

  }
  *total_seg = seg_cnt;

}
