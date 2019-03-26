
void adaptive_picard(double* r0, double* v0, double t0, double t_final, double dt,
  int var_fid, int rad_grav, int quasi_lin, int hot_start, int atm_drag,
  double Deg, double tol, double* Tout, double* Xout, double* Vout, double* seg, int* N, double* Feval){




    // Tolerance Computation
    double coeff_tol = tol/100.0;
    double fit_tol   = tol/10.0;

    // Function Evaluations Counters
    double Jeval = 0.0;
    double Feval = 0.0;

    // PART 2: Determine Segmentation Scheme

    // Compute Keplerian Orbit Period
    double a, e, Period;
    rv2elm(r0,v0,tol,elm);
    a      = elm[1];
    e      = elm[2];
    Period = 2.0*C_PI*sqrt(pow(a,3)/C_MU);

    // Compute Time Vector for 1 Orbit Period
    double w1, w2;
    w1 = Period/2.0;
    w2 = Period/2.0;

    // More stuff here

    int fit_check = 0;
    double coeff;

    seg   = 3.0;
    coeff = 3.0;

    while (fit_check <= coeff){
      //
      d = 360.0/seg;
      f[0] = 0.0;
      for (int i=1; i<=seg+1; i++){
        f[i] = f[i-1] + d;    // True anomaly
      }
      E = 2.0*atan2(tan(0.5*f(2))*sqrt(1-e),sqrt(1+e)); // Eccentric anomaly
      if (E < 0.0){
        E = 2.0*C_PI + E;
      }
      Mf = E - e*sin(E);      // Mean anomaly
      n  = 2.0*C_PI/Period;   // Mean motion
      t0 = 0.0;               // Initial Time
      tf = Mf/n;              // Final Time
      // Scaling parameter
      w1 = (tf+t0)/2.0;
      w2 = (tf-t0)/2.0;

      Nint = 10;
      jmax = 3;

      // Loop through different values for N
      for (int j=0; j<=jmax; j++){
        // INitialize indices
        pcnt = 1;
        gcnt = 1;
        // clear X V G tau time
        // Generate discrete cosine nodes
        if (j == 0){
          N = Nint*(j+1);
          tau = ;
          time = ;
        }
        if (j > 0){
          N = Nint*(j+1);
          tau = ;
          time = ;
        }
        // Compute F&G analytical solution

        // Compute Gravity Radial Adaptation

        // Compute Acceleration

        // Biild Acceleration Vector (include new points as computed)
        if (j == 0){
          for (int i=0; i<=N; i++){
              Gprev[] = G[];
          }
        }



      }












    }











































}
