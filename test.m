%test.m

%osmat.m
%
% Program to compute the Orr-Sommerfeld matrix for three
% dimensional Poiseuille or Couette flows and to compute
% energy matrix
%
%
% INPUT
%
% nosmod = number of Orr-Sommerfeld modes
% R = Reynolds number
% alp = alpha (streamwise wave number)
% beta = beta (spanwise wave number)
% iflow = type of flow (Poiseuille=1, Couette=2)
% nosmod = total number of modes for normal velocity
% iflag = flag
%        iflag = 1: compute the maximum growth and initial condition in
%        time and interval [0,T]
%        iflag = 2: compute the initial disturbance yielding maximum growth
%        and time T
%
% OUTPUT
% d = 3D Orr-Sommerfeld matrix
% M = energy matrix
%
zi=sqrt(-1);
% input data
iflow = 1; %Poiseuille flow
nosmod= 200; %number of OS modes
R= 10000; %Reynolds Number
alp= 1; %input alpha value
beta= 0; %input for beta value

% generate Chebyshev differentiation matrices
[D0,D1,D2,D4]=Dmat(nosmod);

% set up Orr-Sommerfeld matrices A and B
[A,B]=pois(nosmod,alp,beta,R,D0,D1,D2,D4);

% generate energy weight matrix
ak2=alp^2+beta^2;
M=energy(nosmod+1,nosmod+1,ak2);

% compute the Orr-Sommerfeld matrix (by inverting B)
d=B\A;

% Phase 1: Compute eigenvalues and eigenfunctions of
% Orr-Sommerfeld matrix and sort in order of descending
% imaginary part. The function nlize normalizes the
% eigenfunctions with respect to the weight matrix M.
[xs,es]=iord2(d);
% 
% %plot the spectrum of the Orr-Sommerfeld equations
% figure
% plot(es,'o')
% ylim([-1 0.1])
% xlabel('Re(\lambda)')
% ylabel('Im(\lambda)')
% title('Orr-Sommerfeld Spectrum')

% Phase 2: compute the eigenvalues of the perturbed matrix and plot the
% pseudospectrum of the Orr-Sommerfeld equations

all = [];
s = nosmod*2+2;
for j=1:10
   w = rand(s,s)-1/2;
   v = rand(s,s)-1/2;
   E = w+1i*v;
   E = 1e-7/norm(E)*E;
   P = d + E;
   [~,ep]=iord2(P);
   all = [all ep];
end
% %first make the perturbation matrix E
% s = nosmod*2+2;
% E = rand(s,s)+1i*rand(s,s);
% E = 1e-6/norm(E)*E; 
% %now, compute the eigenvalues of the perturbed matrix P
% P = d + E;
% [~,ep]=iord2(P);

%plot the pseudospectrum of the Orr-Sommerfeld equations
figure
plot(all,'.')
ylim([-1 0.1])
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
title('Orr-Sommerfeld Pseudospectrum (\epsilon = 10^{-7})')



