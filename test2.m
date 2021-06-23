
% test2.m

%Chebyshev
%Parameters
kx = 1; kz = 1; Re = 5000; n = 98;

ak2=kx^2+kz^2;

dom = [-1 1]; %Define the Domain
% uX = chebfun('1-x.^2',n);
uX = chebfun('1-x.^2');
uXp = diff(uX, 1); %Define the first derivative of U
uXpp = diff(uX,2); %Define second derivative of U


A = chebop(dom); %Define the a chebop structure's domain
%Define the LHS of the OS equation
A.op = @(x,u,eta)([kx*uX(x).*(diff(u,2)-ak2*u)-kx*uXpp(x).*u-1/Re/1i...
    *(diff(u,4)-2*ak2*diff(u,2)+ak2^2*u); kz*uXp(x)*u + kx*uX(x)*eta-1/Re/1i*(diff(eta,2)-ak2*eta)]); 
       
B = chebop(dom); %Define the a chebop structure's domain
% B.op = @(x,u) diff(u,2)-ak2*u;  %Define the LHS of the OS equation
B.op = @(x,u,eta) ([diff(u,2)-ak2*u ; eta]);


% A.lbc = @(u) [diff(u),u]; %Implement the boundary conditions A
% A.rbc = @(u) [diff(u),u];
% A.lbc = 0;
% A.rbc = 0;
A.lbc = @(u,eta) [diff(u),u,eta]; %Implement the boundary conditions A
A.rbc = @(u,eta) [diff(u),u,eta];

[v, e] = eigs(A,B,98,'SM');

% Cr = diag(real(e));
% Ci = diag(imag(e));
% figure(1)
% plot(Cr,Ci,'o')
% % ylim([-1 0.1]);
% % xlim([0 1]);

normalv = chebfun(v(1,:));
normaleta = chebfun(v(2,:));

Av = normalv(:,3); An = normaleta(:,3);
Sv = normalv(:,69); Sn = normaleta(:,69);
Pv = normalv(:,84); Pn = normaleta(:,84);

% test = abs(normalv(1,:));
% figure(1)
% plot(test,'o')

% data = importdata('spectrum.csv');  
% sCr = data(:,1); sCi = data(:,2);

data1 = importdata('Av.csv');
Avy = data1(:,1); Avv = data1(:,2);

data2 = importdata('An.csv');
Any = data2(:,1); Ann = data2(:,2);

data3 = importdata('Sv.csv');
Svy = data3(:,1); Svv = data3(:,2);

data4 = importdata('Sn.csv');
Sny = data4(:,1); Snn = data4(:,2);

data5 = importdata('Pv.csv');
Pvy = data5(:,1); Pvv = data5(:,2);

data6 = importdata('Pn.csv');
Pny = data6(:,1); Pnn = data6(:,2);


figure(2)
plot(abs(Av)/15.401,'-k')
hold on
plot(Avy, Avv,'or')
hold off
title('Eigenfunction at Branch A (v)')
xlabel('y')
ylabel('v')


figure(3)
plot(abs(An)/16.016,'-k')
hold on
plot(Any, Ann,'or')
hold off
title('Eigenfunction at Branch A (\eta)')
xlabel('y')
ylabel('\eta')

figure(4)
plot(abs(Sv)/0.087456,'-k')
% plot(abs(Sv),'-k')
hold on
plot(Svy, Svv ,'or')
hold off
title('Eigenfunction at Branch S (v)')
% xlim([-1.1,1])
xlabel('y')
ylabel('v [10^4]')

figure(5)
plot(abs(Sn)/0.87135,'-k')
hold on
plot(Sny, Snn ,'or')
hold off
title('Eigenfunction at Branch S (\eta)')
xlabel('y')
ylabel('\eta [10^3]')

figure(6)
plot(abs(Pv)/0.006279,'-k')
hold on
plot(Pvy, Pvv ,'or')
hold off
title('Eigenfunction at Branch P (v)')
xlabel('y')
ylabel('v [10^3]')

figure(7)
plot(abs(Pn)/0.18158,'-k')
hold on
plot(Pny, Pnn ,'or')
hold off
title('Eigenfunction at Branch P (\eta)')
xlabel('y')
ylabel('\eta [10^2]')

