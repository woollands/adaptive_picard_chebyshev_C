/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Feb 2017
*  LAST MODIFIED:    Feb 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __PROP__
#define __PROP__

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include "const.h"

void  picard_chebyshev_propagator(double* r0, double* v0, double t0, double t_final,double deg, double tol, double Period,
   double* tvec, double* t_orig, int seg, int N, int M, int* prep_HS, int coeff_size, int soln_size, int* total_seg,
   double* P1, double* P2, double* T1, double* T2, double* A, double* Ta, double* W1, double* W2, double* Feval,
   double* ALPHA, double* BETA, double* TimeVEC);

#endif
