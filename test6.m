%test6.m
%writing the script using Chebyshev discretization method instead of using
%the Chebfun operators

%   Example 3: computing singular values for frequency responses
%   linearized equations for 3D channel flow. This code reproduces Figure 4.10 in
%   Stability and transition in shear flows by Schmid and Henningson.
%
%   I've used parfor to use all cores in the computer. If you don't like
%   parallel, change parfor to for.
%   Written by Gokul Hariharan, harih020@umn.edu


%% Calculate
% clear;
n= 100; %number of OS modes
Re = 2000;
kx = 1;
kz = 1;
k2 = kx*kx + kz*kz;
k4 = k2*k2;
om = 0.2;

% generate Chebyshev differentiation matrices
[D0,D1,D2,D4]=Dmat(n);

% set up Orr-Sommerfeld matrices A and B
[A,E,F]=pois2(n,kx,kz,Re,D0,D1,D2,D4);

samp = 50;
omegas = linspace(-2,0,samp);
sig1 = [];
for j=1:samp
    om = omegas(j);
    RA = inv(1i*om*eye(2*n+2)-A);
    H = F*RA*E;
    Hnorm = norm(H);
%     [~,s,~] = svds(H,1,'largest');
    sig1 = [sig1 Hnorm];
%     sig1(count) = s(1,1);
%     sig2(count) = s(2,2);
%     count = count+1;
end

semilogy(omegas,real(sig1),'ko')



%{
samp = 50;
omegas = linspace(-2,0,samp);

svals = zeros(samp,1);

for i = 1:samp

omega = omegas(i);
A = chebop([-1,1]);
B = chebop([-1,1]);
C = chebop([-1,1]);
y = chebfun('y',[-1,1]);

U = 1 - y*y;
Uy = diff(U);
Uyy = diff(Uy);
A.op = @(x,v,n)([ii*omega*(diff(v,2) - k2*v) - (diff(v,4)-2*k2*diff(v,2) ...
        + k4*v)/Re - ii*kx*Uyy*v  + ii*kx*U*(diff(v,2) - k2*v);...
        ii*omega*n + ii*kz*Uy*v + ii*kx*U*n - (diff(n,2) - k2*n)/Re]);

A.lbc = @(v,n)[diff(v);v;n];
A.rbc = @(v,n)[diff(v);v;n];

B.op = @(x,dx,dy,dz) ([-ii*kx*diff(dx) - k2*dy - ii*kz*diff(dz);...
                        ii*kz*dx - ii*kx*dz]);
C.op = @(x,v,n)([ii*kx*diff(v)/k2 - ii*kz*n/k2;...
                v ; ...
                ii*kz*diff(v)/k2 + ii*kx*n/k2]);



%cheboppref.setDefaults('minDimension',511);
cheboppref.setDefaults('maxDimension',1000);

% Choose one of Chebfun's discretization, ultraS works much better
cheboppref.setDefaults('discretization',@ultraS)
%cheboppref.setDefaults('discretization',@chebcolloc2)
lam = svdfr(A,B,C,1,'LR');
svals(i) = lam;
disp(i)
end


%% Plot
set(0,'DefaultTextInterpreter', 'latex');
semilogy(omegas,real(svals),'-k');
ylabel('$\sigma_i$');
xlabel('$\omega$');
% hh = legend('$Re = 4000$','$Re = 2000$', '$Re = 1000$','$Re = 500$','location','northwest');
% hh.Interpreter = 'latex';
ax = gca;
ax.YTick = [1 10 100];
% print('-painters','-dsvg','docs/pics/schmidSols');
%} 


