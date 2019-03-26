/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     May 2017
*  LAST MODIFIED:    May 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Reosculate to Keplerian perigee after each orbit
*
* INPUT:
*    X       -- Position array for current segment (km)
*    V       -- Velocity array for current segment (km/s)
*    times   -- Time array for current segment (s)
*    Alpha   -- Position coefficients
*    Beta    -- Velocity coefficients
*    tf      -- End time of current segment (s)
*    t_final -- Final time (s)
*    t_orig  -- Segment time intervals for the first orbit (s)
*    N       -- Degree of the Chebyshev polynomial
*    M       -- Number of sample points
*    k       -- Counter for number of segments per orbit
*    seg     -- Number of segments per orbit
*    prep_HS -- Hot start switch function
*    tol     -- Tolerance
*    orb_end -- Time at end of current orbit (s)
*
* OUTPUTS:
*    tvec -- Array of segment start and end times for the next orbit
*    r0   -- Initial position for the start of the next orbit segment (km)
*    v0   -- Initial velocity for the start of the next orbit segment (km/s)
*/

#include "reosc_perigee.h"
#include "c_functions.h"
#include "rv2elm.h"

void reosc_perigee(double* X, double* V, double* times, double* Alpha, double* Beta,
  double tf, double t_final, double* t_orig, int N, int M, int* k, int seg, int* prep_HS,
  double tol, double* orb_end, double* tvec, double* r0, double* v0){

  // Initialization
  int peri_check = 0;
  double w1, w2, t1, t2, e, TAU_old, TAU, f_old, f_new, TAU_new, df_dtau, err;
  int itrf;
  double TA[N+1];
  memset( TA, 0.0, ((N+1)*sizeof(double)));
  double TB[N];
  memset( TB, 0.0, ((N)*sizeof(double)));

  // Compute True Anomaly (along trajectory)
  double r[3]      = {0.0};
  double v[3]      = {0.0};
  double elm[10]   = {0.0};
  double fvec[M+1];
  memset( fvec, 0.0, ((M+1)*sizeof(double)));
  for (int i=1; i<=M+1; i++){
    for (int j=1; j<=3; j++){
      r[j-1] = X[ID2(i,j,M+1)];
      v[j-1] = V[ID2(i,j,M+1)];
    }
    rv2elm(r,v,tol,elm);
    fvec[i-1] = elm[6];
    e = elm[2];

  }

  if (*k == seg-1 || (*k == *prep_HS && fabs(tf - t_final)/tf > tol)){
    if (fabs(e) > 1e-6){    // Skip for zero eccentricity (no need to re-osculate as perigee is undefined)
      for (int i=1; i<=M+1; i++){
        // If passing through perigee
        if ((i > 2) && (fvec[i-1] < fvec[i-2])){
          // Prepare for secant method
          t1        = times[i-2];
          t2        = times[i-1];
          w1        = (times[M] + times[0])/2.0;
          w2        = (times[M] - times[0])/2.0;
          TAU_old   = (t1 - w1)/w2;
          TAU       = (t2 - w1)/w2;
          f_old     = fvec[i-2];
          f_new     = fvec[i-1];

          // Initialize
          itrf = 0;
          err  = 10.0;

          // Secant method to find tau at f0
          while (err > tol){
            for (int kk=0; kk<=N; kk++){
              TA[kk] = cos(kk*acos(TAU));
            }
            for (int kk=0; kk<=N-1; kk++){
              TB[kk] = cos(kk*acos(TAU));
            }
            r0[0] = 0.0; r0[1] = 0.0; r0[2] = 0.0;
            v0[0] = 0.0; v0[1] = 0.0; v0[2] = 0.0;
            for(int j=0; j<=N; j++){
              r0[0] = r0[0] + TA[j]*Alpha[ID2(j+1,1,N+1)];
              r0[1] = r0[1] + TA[j]*Alpha[ID2(j+1,2,N+1)];
              r0[2] = r0[2] + TA[j]*Alpha[ID2(j+1,3,N+1)];
            }
            for (int j=0; j<=N-1; j++){
              v0[0] = v0[0] + TB[j]*Beta[ID2(j+1,1,N)];
              v0[1] = v0[1] + TB[j]*Beta[ID2(j+1,2,N)];
              v0[2] = v0[2] + TB[j]*Beta[ID2(j+1,3,N)];
            }
            tf = TAU*w2 + w1;

            rv2elm(r0,v0,tol,elm);
            f_new = elm[6];

            if (f_new > C_PI){
              f_new = f_new - 2.0*C_PI;
            }

            // Compute new TAU
            df_dtau = (f_new - f_old)/(TAU - TAU_old);
            if (df_dtau == 0){
              break;
            }
            TAU_new = TAU + (0 - f_new)/df_dtau;
            if (TAU_new > 1.0){
              TAU_new = 1.0;
              break;
            }

            // Check error
            err     = fabs(f_new);

            // Update
            f_old   = f_new;
            TAU_old = TAU;
            TAU     = TAU_new;
            itrf = itrf + 1;
            if (itrf > 10){
              break;
            }
          }

          // Compute time vector for next orbit
          for (int j=0; j<=seg; j++){
            tvec[j] = t_orig[j] + tf;
          }
          *k    = -1;
          *orb_end    = tvec[0];
          peri_check = 1;
          break;
        }
      }
    }

    // Compute time vector for next orbit
    if (peri_check == 0){
      for (int j=0; j<=seg; j++){
        tvec[j] = t_orig[j] + tf;
      }
    }
    *orb_end = tvec[0];

    // Reset counters
    *prep_HS = -1;
    *k       = -1;
  }

}
