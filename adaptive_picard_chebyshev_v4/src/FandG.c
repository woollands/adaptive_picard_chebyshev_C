#include "FandG.h"

/*! \brief Calculates the analytic solution to the unperturbed orbit problem using F&G solution
 *   Used as a "warm-start" for MCPI.
 *  \param[in] z0  Initial state vector [r0,v0] in units of [e.r] and [e.r/ctu]
 *  \param[out] z  Output state vector [rf,vf] in units of [e.r] and [e.r/ctu]
 *  \param[in] dt  Desired time of output in units of [ctu]
*/

//Define F and G magic numbers.
const int newtonMaxIt = 300;
const double newtonTol = 1E-13;

void FandG( const double* z0, double* zf, const double dt)
{
    int kk;

    double rMag = sqrt( z0[0]*z0[0] + z0[1]*z0[1] + z0[2]*z0[2] );
    double vSq  = z0[3]*z0[3] + z0[4]*z0[4] + z0[5]*z0[5];
    double a    = 1.0 / ( 2.0/rMag - vSq/C_MU );
    double sig0 = (z0[0]*z0[3] + z0[1]*z0[4] + z0[2]*z0[5])/sqrt(C_MU);
    double EHat = newtonFandG( a, dt, rMag, sig0, newtonTol );
    double r    = a + (rMag-a)*cos(EHat) + sqrt(a)*sig0*sin(EHat);
    double F    = 1.0 - a/rMag*( 1 - cos(EHat) );
    double G    = dt + sqrt( a*a*a/C_MU )*( sin(EHat) - EHat )  ;
    double Fdot = -sqrt(a*C_MU)*sin(EHat)/(r*rMag);
    double Gdot = 1.0 - a/r*(1.0-cos(EHat));

    // Current state is a linear combination of the initial condition vectors
    int ii;
    for ( ii=0; ii<3; ii++ ) {
        zf[ii]   = F*z0[ii]    + G*z0[ii+3];
        zf[ii+3] = Fdot*z0[ii] + Gdot*z0[ii+3];
    }
}

/*! \brief Newton solver to solve Kepler's equation in the F&G solution.
 *  \param[in] a     Semimajor axis [e.r]
 *  \param[in] dt    Desired time of solution [ctu]
 *  \param[in] rMag  Magnitude of the radius vector [e.r.]
 *  \param[in] sig0  Measure of orthogonality between instantaneous positiona and velocity vector
 *  \param[in] tol   Convergence threshold
 *  \param[out] Ehat Change in eccentric anomaly
*/

double newtonFandG( const double a, const double dt, const double rMag,
                    const double sig0, const double tol )
{

    // Initial guess
    //double EHat = C_PI;
    double EHat = C_PI;
    int ctr = 0;
    double fx, dfx;
    double dE = 1.0;

    while ( ( fabs(dE/EHat) > tol ) && ( ctr <= newtonMaxIt ) ) {

        // Calculate dE correction
        fx = sqrt(C_MU)/(a*sqrt(a))*dt  - EHat + (1.0-rMag/a)*sin(EHat) +  \
            sig0/sqrt(a)*(cos(EHat)-1.0);
        dfx = -1.0 + ((1.0-rMag/a)*cos(EHat) - sig0/sqrt(a)*sin(EHat));
        dE  = -1.0*fx/dfx;
        if ( ctr == newtonMaxIt ) {
        }

        // Add correction and increment ctr
        EHat += dE;
        ctr++;

    }

    return EHat;
}
