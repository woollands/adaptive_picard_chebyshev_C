/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Feb 2016
*  LAST MODIFIED:    Feb 2016
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Generates linear correction term for acceleration
*
* INPUT:
*    X     -- Position (km)
*    del_X -- Position error (km)
*
* OUTPUTS:
*    del_A -- Acceleration correction (km/s^2)
*
* REFERENCES:
* 1. Junkins, J.L., and Woollands, R., "Nonlinear Differential Equations Solvers via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics",
*    AAS/AIAA Astrodynamics Specialist Conference, Stevenson, WA, 2017.
* 2. Junkins, J.L., and Woollands, R., "Adaptive-Picard-Chebyshev for Propagating Perturbed Two-Body Orbits",
*    JGCD, submitted 2017.
*/

#include "picard_error_feedback.h"

void picard_error_feedback(double* X, double* del_X, double* del_a){

  double R3, R5;
  R3 = pow(sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]),3);
  R5 = pow(sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]),5);

  double J[3][3] = {0.0};

  J[0][0] = 3.0*C_MU*pow(X[0],2)/R5 - C_MU/R3;
  J[0][1] = 3.0*C_MU*X[0]*X[1]/R5;
  J[0][2] = 3.0*C_MU*X[0]*X[2]/R5;
  J[1][0] = 3.0*C_MU*X[1]*X[0]/R5;
  J[1][1] = 3.0*C_MU*pow(X[1],2)/R5 - C_MU/R3;
  J[1][2] = 3.0*C_MU*X[1]*X[2]/R5;
  J[2][0] = 3.0*C_MU*X[2]*X[0]/R5;
  J[2][1] = 3.0*C_MU*X[2]*X[1]/R5;
  J[2][2] = 3.0*C_MU*pow(X[2],2)/R5 - C_MU/R3;

  // Acceleration Error Feedback
  del_a[0] = J[0][0]*del_X[0] + J[1][0]*del_X[1] + J[2][0]*del_X[2];
  del_a[1] = J[0][1]*del_X[0] + J[1][1]*del_X[1] + J[2][1]*del_X[2];
  del_a[2] = J[0][2]*del_X[0] + J[1][2]*del_X[1] + J[2][2]*del_X[2];

}
