%test3.m
%from test.m
n= 100; ny = 1000; Re = 5000; kx = 1; kz = 1; 
[D0,D1,D2,D4]=Dmat(n);
[A,B]=pois(n,kx,kz,Re,D0,D1,D2,D4);
ak2=kx^2+kz^2;
M=energy(n+1,n+1,ak2);
d=B\A;
[x,~] = eigs(d,98,'smallestabs');
xn = nlize(x,M);

%from test2.m
dom = [-1 1]; %Define the Domain
uX = chebfun('1-x.^2');
uXp = diff(uX, 1); %Define the first derivative of U
uXpp = diff(uX,2); %Define second derivative of U

A = chebop(dom);
A.op = @(x,u,eta)([kx*uX(x).*(diff(u,2)-ak2*u)-kx*uXpp(x).*u-1/Re/1i...
    *(diff(u,4)-2*ak2*diff(u,2)+ak2^2*u); kz*uXp(x)*u + kx*uX(x)*eta-1/Re/1i*(diff(eta,2)-ak2*eta)]);        
B = chebop(dom);
B.op = @(x,u,eta) ([diff(u,2)-ak2*u ; eta]);

A.lbc = @(u,eta) [diff(u),u,eta]; %Implement the boundary conditions A
A.rbc = @(u,eta) [diff(u),u,eta];
[v, ~] = eigs(A,B,98,'SM');
normalv = chebfun(v(1,:));
normaleta = chebfun(v(2,:));

vec=(0:ny)';
yj = cos(pi*vec/ny);



% Branch A
data1 = importdata('Av.csv');
Avy = data1(:,1); Avv = data1(:,2);

data2 = importdata('An.csv');
Any = data2(:,1); Ann = data2(:,2);

Av = normalv(:,3); An = normaleta(:,3);

xav = xn(1:101,3);
Av2 = abs(cheb_expansion_soln(yj,xav));
xan = xn(102:202,3);
An2 = abs(cheb_expansion_soln(yj,xan));

figure(1)
plot(yj,xeigfun/14.8239,'-g','LineWidth',4)
hold on
plot(abs(Av)/15.401,'-k','LineWidth',2)
plot(Avy, Avv,'or','LineWidth', 2, 'MarkerSize', 4)
hold off
title('Eigenfunction at Branch A (v)')
xlabel('y')
ylabel('v')

figure(2)
plot(yj,xeigfun/14.924,'-g','LineWidth',4)
hold on
plot(abs(An)/16.016,'-k','LineWidth',2)
plot(Any, Ann,'or','LineWidth', 2, 'MarkerSize', 4)
hold off
title('Eigenfunction at Branch A (\eta)')
xlabel('y')
ylabel('\eta')

% S branch
data3 = importdata('Sv.csv');
Svy = data3(:,1); Svv = data3(:,2);

data4 = importdata('Sn.csv');
Sny = data4(:,1); Snn = data4(:,2);

Sv = normalv(:,69); Sn = normaleta(:,69);
xsv = xn(1:101,69);
Sv2 = abs(cheb_expansion_soln(yj,xsv));
xsn = xn(102:202,69);
Sn2 = abs(cheb_expansion_soln(yj,xsn));

figure(3)
plot(yj,Sv2/0.063222,'-g','LineWidth',4)
hold on
plot(abs(Sv)/0.087456,'-k','LineWidth',2)
plot(Svy, Svv ,'or','LineWidth', 2, 'MarkerSize', 4)
hold off
title('Eigenfunction at Branch S (v)')
xlabel('y')
ylabel('v [10^4]')

figure(4)
plot(yj,Sn2/0.633453,'-g','LineWidth',4)
hold on
plot(abs(Sn)/0.87135,'-k','LineWidth',2)
plot(Sny, Snn ,'or','LineWidth', 2, 'MarkerSize', 4)
hold off
title('Eigenfunction at Branch S (\eta)')
xlabel('y')
ylabel('\eta [10^3]')


%P branch
data5 = importdata('Pv.csv');
Pvy = data5(:,1); Pvv = data5(:,2);

data6 = importdata('Pn.csv');
Pny = data6(:,1); Pnn = data6(:,2);

Pv = normalv(:,27); Pn = normaleta(:,38);
xpv = xn(1:101,25);
Pv2 = abs(cheb_expansion_soln(yj,xpv));
xpn = xn(102:202,38);
Pn2 = abs(cheb_expansion_soln(yj,xpn));


figure(5)
plot(yj,Pv2/0.0021399,'-g','LineWidth',4)
hold on
plot(abs(Pv)/0.003318,'-k','LineWidth',2)
plot(Pvy, Pvv ,'or','LineWidth', 2, 'MarkerSize', 4)
hold off
title('Eigenfunction at Branch P (v)')
xlabel('y')
ylabel('v [10^3]')

figure(6)
plot(yj,Pn2/0.1237,'-g','LineWidth',4)
hold on
plot(abs(Pn)/0.211676,'-k','LineWidth',2)
plot(Pny, Pnn ,'or','LineWidth', 2, 'MarkerSize', 4)
hold off
title('Eigenfunction at Branch P (\eta)')
xlabel('y')
ylabel('\eta [10^2]')

