# Eigenvalue Decomposition of the Orr-Sommerfeld operator
This repo contains the MATLAB files for the eigenvalue analysis of the Orr-Sommerfeld matrix and is adapted from the MATLAB code found in "Appendix A.6: MATLAB codes for hydrodynamic stability" in "Stability and Transition in Shear Flows" by Peter J. Schmid and Dan S. Henningson (2001).

The Orr-Sommerfeld matrix describes the normal velocity and vorticity fluctutions in a plane Poiseuille flow and is derived form the linearized Navier-Stokes equation. In this problem, we investigate the Orr-Sommerfeld matrix with streamwise wavenumber (kx) = spanwise wavenumber (kz) = 1 and Re = 5,000. The test.m file will derive the Orr-Sommerfeld matrix and plot the eigenspectrum and eigenfunction of the problem.

While this analysis produces the correct eigenspectrum, the eigenfunctions produced from this script does not accurately protray what is presented in the textbook in Fig. 3.2 (refer to eigenfunction.png for image). 

