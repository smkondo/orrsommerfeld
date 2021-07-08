%test4.m
%conducting resolvent analysis on the Orr-Sommerfeld operator
n= 100; %number of OS modes
ny = 1000; %number of Chebyshev discretization points
Re = 2000; %Reynolds Number
kx = 1; %input alpha value
kz = 1; %input for beta value
% om = 0.2; 

% generate Chebyshev differentiation matrices
[D0,D1,D2,D4]=Dmat(n);

% set up Orr-Sommerfeld matrices A and B
[A,B]=pois(n,kx,kz,Re,D0,D1,D2,D4);

% generate energy weight matrix
ak2=kx^2+kz^2;
M=energy(n+1,n+1,ak2);

% sig1 = zeros(1,2/0.1+1);
% sig2 = zeros(1,2/0.1+1);
samp = 50;
omegas = linspace(-2,0,samp);

sig1 = [];
% sig2 = [];
for j=1:samp
    om = omegas(j);
%     count = 1;
    H = inv(1i*om*eye(2*n+2)+B\A);
%     Hnorm = norm(H);
    [~,s,~] = svds(H,1,'largest');
    sig1 = [sig1 s];
%     sig1(count) = s(1,1);
%     sig2(count) = s(2,2);
%     count = count+1;
end

% x = 0:0.01:2;
semilogy(omegas,real(sig1),'ko')
% hold on
% plot(x,sig2,'bo')
% hold off

% % compute the resolvent operator of the Orr-Sommerfeld matrix
% H = inv(-1i*om*eye(2*n+2)+B\A);
% 
% [u,s,v] = svds(H,2);
% 
% sing = real(diag(s));
% 
% plot(sing,'o')

