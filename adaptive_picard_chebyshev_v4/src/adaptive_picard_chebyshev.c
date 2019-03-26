/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     May 2017
*  LAST MODIFIED:    May 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Adaptive Picard-Chebyshev Numerical Integration
*
* INPUT:
*    r0        -- Initial position vector (km)
*    v0        -- Initial velocity vector (km/s)
*    t0        -- Initial time (s)
*    tf        -- Final time (sec)
*    dt        -- Solution Output Time Interval (s)
*    deg       -- Gravity Degree (max 100)
*    tol       -- Tolerance
*    soln_size -- Length of solution array
*    Feval     -- Function evaluation counter
*
* OUTPUTS:
*    Soln      -- [Position (km), Velocity (km/s)] solution at user specified times
*
* REFERENCES:
* 1. Junkins, J.L., and Woollands, R., "Nonlinear Differential Equations Solvers via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics",
*    AAS/AIAA Astrodynamics Specialist Conference, Stevenson, WA, 2017.
* 2. Junkins, J.L., and Woollands, R., "Nonlinear Differential Equation Solvers via Adaptve Picard-Chebyshev Iteration: Applications in Astrodynamics",
*    JGCD, submitted 2017.
*
* COMMENTS:
*
*/

#include "adaptive_picard_chebyshev.h"
#include "polydegree_segments.h"
#include "prepare_propagator.h"
#include "picard_chebyshev_propagator.h"
#include "interpolate.h"
#include "c_functions.h"

void adaptive_picard_chebyshev(double* r0,double* v0, double t0, double tf, double dt, double deg, double tol, int soln_size, double* Feval, double* Soln){

  /* 1. DETERMINE DEGREE/SEGMENTATION SCHEME
  Compute the polynomial degree and number of segments per orbit that will
  result in a solution that satisfies the user specified tolerance. */
  int seg, N;
  double tp, Period;
  polydegree_segments(r0,v0,deg,tol,Feval,&seg,&N,&tp,&Period);
  // printf("N %i\t",N);
  // printf("seg %i\n",seg);

  // Array size for coefficients and solution
  int coeff_size;
  coeff_size = (int) (tf/Period + 1.0)*(seg+2.0)*(N+1);

  /* 2. PREPARE PROPAGATOR
  Compute and store the begin and end times for each segment (based on true
  anomaly segmentation) and load the constant matrices corresponding to N. */
  // Initialize Arrays
  int M = N;                // # sample points = polynomial degree
  int prep_HS = -1;         // Hot start switch condition
  double T2[(M+1)*(N+1)];   // [(M+1)x(N+1)]
  memset( T2, 0.0, ((M+1)*(N+1)*sizeof(double)));
  double P2[(N+1)*N];       // [(N+1)xN]
  memset( P2, 0.0, ((N+1)*N*sizeof(double)));
  double T1[(M+1)*N];       // [(M+1)xN]
  memset( T1, 0.0, ((M+1)*N*sizeof(double)));
  double P1[N*(N-1)];       // [Nx(N-1)]
  memset( P1, 0.0, (N*(N-1)*sizeof(double)));
  double Ta[(M+1)*(N-1)];   // [(M+1)x(N-1)]
  memset( Ta, 0.0, ((M+1)*(N-1)*sizeof(double)));
  double A[(N-1)*(M+1)];    // [(N-1)x(M+1)]
  memset( A, 0.0, ((N-1)*(M+1)*sizeof(double)));
  double t_orig[seg+1];
  memset( t_orig, 0.0, ((seg+1)*sizeof(double)));
  double tvec[seg+1];
  memset( tvec, 0.0, ((seg+1)*sizeof(double)));
  prepare_propagator(r0,v0,t0,tf,dt,tp,tol,N,M,seg,&prep_HS,t_orig,tvec,P1,P2,T1,T2,A,Ta);

  /* 3. PICARD-CHEBYSHEV PROPAGATOR
  Propagate from t0 to tf, iterating on each segment (Picard Iteration), until
  completion. */
  double *ALPHA;
  ALPHA = calloc((coeff_size*3),sizeof(double));
  double *BETA;
  BETA = calloc((coeff_size*3),sizeof(double));
  int total_seg = 0;
  int sz = (int) ceil(1.1*tf/Period)*seg;
  double segment_times[sz];
  memset( segment_times, 0.0, (sz*sizeof(double)));
  double W1[sz];
  memset( W1, 0.0, (sz*sizeof(double)));
  double W2[sz];
  memset( W2, 0.0, (sz*sizeof(double)));
  picard_chebyshev_propagator(r0,v0,t0,tf,deg,tol,Period,tvec,t_orig,seg,N,M,&prep_HS,coeff_size,soln_size,&total_seg,
    P1,P2,T1,T2,A,Ta,W1,W2,Feval,ALPHA,BETA,segment_times);

  // /* 4. INTERPOLATE SOLUTION
  // The Chebyshev coefficients from each of the orbit segments are used to compute
  // the solution (position & velocity) at the user specified times. */
  interpolate(ALPHA,BETA,soln_size,coeff_size,N,segment_times,W1,W2,t0,tf,dt,total_seg,Soln);

  free(ALPHA);
  free(BETA);
}
