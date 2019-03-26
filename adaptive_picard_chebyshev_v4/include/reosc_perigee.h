/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     May 2017
*  LAST MODIFIED:    May 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __PER__
#define __PER__

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include "const.h"

void reosc_perigee(double* X, double* V, double* times, double* Alpha, double* Beta,
  double tf, double t_final, double* t_orig, int N, int M, int* k, int seg, int* prep_HS,
  double tol, double* orb_end, double* tvec, double* r0, double* v0);

#endif
