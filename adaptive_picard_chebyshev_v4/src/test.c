/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     May 2017
*  LAST MODIFIED:    May 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Set up an Adaptive-Picard-Chebyshev integration test case
*  REFERENCE:        Woollands, R., and Junkins, J., "Nonlinear Differential Equation Solvers
*                    via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics", JGCD, 2016.
*/

#include "adaptive_picard_chebyshev.h"
#include "c_functions.h"
#include "EGM2008.h"
#include <time.h>

FILE *fID;

int main(){

  // Initialize Input Variables
  // LEO
  double r0[3] = {7000.0, 0.0, 0.0};                // Initial Position (km)
  double v0[3] = {0.0, 8.003798178945150, 0.0};     // Initial Velocity (km/s)
  double t0    = 0.0;                               // Initial Times (s)
  double tf    = 3.0*7.121081577578024e+03;         // Final Time (s)
  // MEO
  // double r0[3] = {9000.0, 0.0, 0.0};                                // Initial Position (km)
  // double v0[3] = {0.0, 6.7419845635570, 1.806509319188210};         // Initial Velocity (km/s)
  // double t0    = 0.0;                                               // Initial Times (s)
  // double tf    = 3.0*9.952014050491189e+03;                         // Final Time (s)
  // GEO
  // double r0[3] = {42000, 0.0, 0.0};                              // Initial Position (km)
  // double v0[3] = {0.0, 3.080663355435613, 0.0};                  // Initial Velocity (km/s)
  // double t0    = 0.0;                                            // Initial Times (s)
  // double tf    = 3.0*8.566135031791795e+04;                      // Final Time (s)
  // GTO
  // double r0[3] = {8064, 0.0, 0.0};                               // Initial Position (km)
  // double v0[3] = {0.0, 9.112725097814229, 0.0};                  // Initial Velocity (km/s)
  // double t0    = 0.0;                                            // Initial Times (s)
  // double tf    = 3.0*3.981179798339227e+04;                      // Final Time (s)
  // Molniya
  // double r0[3] = {7435.12, 0.0, 0.0};                            // Initial Position (km)
  // double v0[3] = {0.0, 4.299654205302486, 8.586211043023614};    // Initial Velocity (km/s)
  // double t0    = 0.0;                                            // Initial Times (s)
  // double tf    = 5.0*4.306316113361824e+04;                      // Final Time (s)

  double dt    = 30.0;                             // Soution Output Time Interval (s)
  double deg   = 70.0;                             // Gravity Degree (max 100)
  double tol   = 1.0e-15;                          // Tolerance

  // Initialize Output Variables
  int soln_size = (int) (1.1*(tf/dt));
  if (soln_size == 1){
    soln_size = 2;
  }
  double *Soln;
  Soln = calloc(soln_size*6,sizeof(double));       // Position (km) & Velocity (km/s)

  double Feval[2] = {0.0};

  // Call Adaptive Picard Chebyshev Integrator
  clock_t startTime = clock();
  for (int tt=0; tt<=1; tt++){
    adaptive_picard_chebyshev(r0,v0,t0,tf,dt,deg,tol,soln_size,Feval,Soln);
  }
  clock_t endTime = clock();
  float elapsedTime = ((float) (endTime - startTime))/CLOCKS_PER_SEC/(1.0);
  printf("Elapsed time: %f s\t",elapsedTime);

  // Number of function evaluations
  int total;
  total = (int) ceil(Feval[0] + Feval[1]*pow(6.0,2)/pow(deg,2));
  printf("Func Evals: %i\t",total);

  // Save as text file [time r v H]
  double state[6] = {0.0};
  double H    = 0.0;
  double H0   = 0.0;
  double Hmax = 0.0;
  fID = fopen((const char*)"./output.txt",(const char*)"w");
  double t_curr = t0;
  for (int i=1; i<=soln_size; i++){
    fprintf(fID,"%1.16E\t", t_curr);
    for (int j=1; j<=6; j++){
      fprintf(fID,"%1.16E\t",Soln[ID2(i,j,soln_size)]);
      state[j-1] = Soln[ID2(i,j,soln_size)];
    }
    jacobiIntegral(t_curr,state,&H,deg);
    if (i == 1){
      H0 = H;
    }
    if (fabs((H-H0)/H0) > Hmax){
      Hmax = fabs((H-H0)/H0);
    }
    fprintf(fID,"%1.16E\t",fabs((H-H0)/H0));
    t_curr = t_curr + dt;
    if (t_curr > tf){
      break;
    }
    fprintf(fID,"\n");
  }
  printf("Hmax %1.16E\n",Hmax);
  fclose(fID);
  free(Soln);

}
