%test4.m
%conducting resolvent analysis on the Orr-Sommerfeld operator
n= 100; %number of OS modes
ny = 1000; %number of Chebyshev discretization points
Res = [4000,2000,1000,500];
% Re = 2000; %Reynolds Number
kx = 1; %input alpha value
kz = 1; %input for beta value
% om = 0.2; 

% generate Chebyshev differentiation matrices
[D0,D1,D2,D4]=Dmat(n);
sig1 = zeros(samp,1);
ak2=kx^2+kz^2;
M=energy(n+1,n+1,ak2);

samp = 100;
omegas = linspace(0,1.2,samp);
for i = 1:4
    Re = Res(i);

    % set up Orr-Sommerfeld matrices A and B
    [A,B]=pois(n,kx,kz,Re,D0,D1,D2,D4);

    for j=1:samp
        om = omegas(j);
        H = inv(B\A-om*eye(2*n+2));
        [~,s,~] = svds(H,1,'largest');
        sig1(j,i)=s;
    end
end

% data = importdata('Rvalue.csv');
% som = data(:,1); sR = data(:,2)/0.1326196199;
semilogy(omegas,real(sig1(:,1)),'-k')
hold on
semilogy(omegas,real(sig1(:,2)),'--k')
semilogy(omegas,real(sig1(:,3)),'-.k')
semilogy(omegas,real(sig1(:,4)),':k')
hold off
ylabel('\sigma_i');
xlabel('\omega');
legend('Re=4000','Re=2000','Re=1000','Re=500','location','northwest')
title('Resolvent norm vs. Frequency')


