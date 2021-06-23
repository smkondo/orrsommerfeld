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
[x,e] = eig(d);

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
xtest = abs(x(:,4));
figure(2)
plot(xtest,'o')

function val = cheb_expansion_soln(y,a)
    %function that evaluates the Chebyshev expansion for given coeffs
    % y = points to evaluate Chebyshev expansion at
    % a = vector of coeffs for the Chebyshev expansion
    val = zeros(length(y),1);
    for i=1:length(y)
        for j=1:length(a)
            val(i) = val(i) + a(j)*cheb_basis(y(i),j-1);
        end
    end
end

function val = cheb_basis(y,n)
    %Returns value of Chebyshev polynomial of order n evaluated at y
    val = cos(n*acos(y));
end


