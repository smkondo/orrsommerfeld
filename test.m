% test.m
%
% Program to compute the Orr-Sommerfeld matrix for three
% dimensional Poiseuille flows and to compute
% energy matrix
%
%
% INPUT
%
% n = number of chebyshev modes
% Re = Reynolds number
% kx = streamwise wave number
% kz = spanwise wave number
%
% OUTPUT
% d = 3D Orr-Sommerfeld matrix
% M = weighted energy matrix (for normalization)
%
% input data
n= 100; %number of OS modes
ny = 1000; %number of Chebyshev discretization points
Re = 5000; %Reynolds Number
kx = 1; %input alpha value
kz = 1; %input for beta value

% generate Chebyshev differentiation matrices
[D0,D1,D2,D4]=Dmat(n);

% set up Orr-Sommerfeld matrices A and B
[A,B]=pois(n,kx,kz,Re,D0,D1,D2,D4);

% generate energy weight matrix
ak2=kx^2+kz^2;
M=energy(n+1,n+1,ak2);

% compute the Orr-Sommerfeld matrix (by inverting B)
d=B\A;

% solve the eigenvalue problem 
[x,e] = eigs(d,98,'smallestabs');

% plot the eigenspectrum
Cr = diag(real(e)); %real component of the eigenvalues
Ci = diag(imag(e)); %imaginary component of the eigenvalues
figure(1)
plot(Cr,Ci,'.k','MarkerSize',16)
ylim([-1 0.1]);
xlim([0 1]);
xlabel('Cr', 'FontSize', 18)
ylabel('Ci','FontSize', 18)
title('Eigenspectrum of Orr-Sommerfeld Problem', 'FontSize', 18)

% plot the eigenfunction 
xn = nlize(x,M);

% A-branch
xtest = xn(1:101,25);
vec=(0:ny)';
yj = cos(pi*vec/ny);
xeigfun = abs(cheb_expansion_soln(yj,xtest));

data5 = importdata('Pv.csv');
Pvy = data5(:,1); Pvv = data5(:,2);

figure(2)
plot(yj,xeigfun/0.0021399,'-g')
hold on
plot(Pvy, Pvv,'or')
hold off
title('Eigenfunction at Branch A (v)')
xlabel('y')
ylabel('v')


