              ADAPTIVE PICARD-CHEBYSHEV ITERATION
                  (Perturbed Two-body Problem)
                      TEXAS A&M UNIVERSITY
                      (Version 3, Mar 2018)

            Robyn Woollands (robyn.woollands@gmail.com)


1. Edit "test.c" to specify initial conditions, propagation time, output time-step,
   spherical harmonic gravity degree and solution tolerance.

2. Compile the matrix from within the folder /src
   >> make matrix_builder

3. Perform the one time build of the constant Picard-Chebyshev matrices from within the folder /src
   >> ./matrix_builder

4. Compile the propagator from within the folder /src
   >> make

5. Propagate the test case orbit from within the folder /src
   >> ./test

6. The output solution is saved in output.txt

NOTES:
1. Jan 2018: Added the Hamiltonian check and fixed a few bugs that were present in Version 1 (v1: Dec 2017).
2. Mar 2018: Fixed a bug that I found when interpolating position for the Molniya orbit in Version 2 (v2: Jan 2018).
