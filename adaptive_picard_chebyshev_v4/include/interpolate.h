/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Feb 2017
*  LAST MODIFIED:    Feb 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __INT__
#define __INT__

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include "const.h"

void interpolate(double* ALPHA, double* BETA, int soln_size, int coeff_size, int N, double* segment_times,
  double* W1, double* W2, double t0, double tf, double dt, int total_seg, double* Soln);

#endif
