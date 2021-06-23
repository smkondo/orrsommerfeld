
% test4.m

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
A.op = @(x,u,eta)([1i*kx*uX(x).*(diff(u,2)-ak2*u)-1i*kx*uXpp(x).*u-1/Re...
    *(diff(u,4)-2*ak2*diff(u,2)+ak2^2*u); 1i*kz*uXp(x)*u + 1i*kx*uX(x)*eta-1/Re*(diff(eta,2)-ak2*eta)]); 
       
B = chebop(dom); %Define the a chebop structure's domain
% B.op = @(x,u) diff(u,2)-kx^2*u;  %Define the LHS of the OS equation
B.op = @(x,u,eta) ([diff(u,2)-kx^2*u ; eta]);


% A.lbc = @(u) [diff(u),u]; %Implement the boundary conditions A
% A.rbc = @(u) [diff(u),u];
% A.lbc = [0,0,0];
% A.rbc = [0,0,0];
A.lbc = @(u,eta) [u,diff(u),eta]; %Implement the boundary conditions A
A.rbc = @(u,eta) [u,diff(u),eta];

[v, e] = eigs(A,B);

% Cr = diag(real(e));
% Ci = diag(imag(e));
% figure(1)
% plot(Ci,Cr,'o')
% % ylim([-1 0.1]);
% % xlim([0 1]);


normalv = chebfun(v(1,:));
normaleta = chebfun(v(2,:));

Av = normalv(:,2); An = normaleta(:,2);
% Sv = normalv(:,72); Sn = normaleta(:,72);
% Pv = normalv(:,94); Pn = normaleta(:,94);

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
plot(abs(Av)/9.681,'-k')
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

% figure(4)
% plot(abs(Sv)/0.02748,'-k')
% hold on
% plot(Svy, Svv ,'or')
% hold off
% title('Eigenfunction at Branch S (v)')
% xlabel('y')
% ylabel('v [10^4]')
% 
% figure(5)
% plot(abs(Sn),'-k')
% hold on
% plot(Sny, Snn ,'or')
% hold off
% title('Eigenfunction at Branch S (\eta)')
% xlabel('y')
% ylabel('\eta [10^3]')
% 
% figure(6)
% plot(abs(Pv)/0.14706,'-k')
% hold on
% plot(Pvy, Pvv ,'or')
% hold off
% title('Eigenfunction at Branch P (v)')
% xlabel('y')
% ylabel('v [10^3]')
% 
% figure(7)
% plot(abs(Pn)/0.281,'-k')
% hold on
% plot(Pny, Pnn ,'or')
% hold off
% title('Eigenfunction at Branch P (\eta)')
% xlabel('y')
% ylabel('\eta [10^2]')

a=3;
p=72;
sbest =94;
s=5170769478;



% xtest1 = abs(normalv(:,s));
% xtest2 = real(v(:,s));
% xtest3 = imag(v(:,s));
% figure(2)
% plot(xtest1,'o')
% hold on
% % plot(xtest2, '-k')
% % plot(xtest3, '-g')
% hold off

%{
imagy = [-1.00000, -0.97595, -0.95591, -0.91984, -0.87575, -0.83567,...
    -0.79559, -0.76353, -0.72745, -0.69539, -0.64329, -0.59118, -0.53507,...
    -0.48297, -0.41483, -0.33467, -0.26253, -0.17435, -0.07415, 0.00601, ...
    0.10621, 0.20641, 0.32265, 0.42285, 0.52705, 0.63126, 0.73948, 0.83166,...
    0.94389, 0.90782, 0.86774, 0.77555, 0.96393, 0.99599];

imagv = [0.00002, 0.00012, 0.00023, 0.00033, 0.00031, 0.00029, 0.00033,...
    0.00044, 0.00056, 0.00067, 0.00081, 0.00092, 0.00102, 0.00108, 0.00115,...
    0.00123, 0.00127, 0.00131, 0.00134, 0.00134, 0.00134, 0.00131, 0.00125,...
    0.00115, 0.00104, 0.00087, 0.00054, 0.00031, 0.00033, 0.00037, 0.00033,...
    0.00043, 0.00020, 0.00004];

magy = [-1.00000, -0.97595, -0.96393, -0.94790, -0.93186, -0.91984, -0.90381,...
    -0.88778, -0.86373, -0.83968, -0.81162, -0.77555, -0.73146, -0.68337,...
    -0.61122, -0.54709, -0.45892, -0.32665, -0.21443, -0.04609, 0.08617,...
    0.21042, 0.38277, 0.49900, 0.57114, 0.65130, 0.71944,0.77154, 0.82766,...
    0.88778, 0.94790, 0.92385, 0.86774, 0.79960, 0.75150, 0.68337, 0.90381,...
    0.99599, 0.97595, 0.96393];
magv = [0.00002, 0.00021, 0.00043, 0.00066, 0.00098, 0.00129, 0.00159,...
    0.00203, 0.00255, 0.00313, 0.00370, 0.00445, 0.00514, 0.00567, 0.00632,...
    0.00676, 0.00724, 0.00774, 0.00803, 0.00816, 0.00814, 0.00797, 0.00755,...
    0.00701, 0.00655, 0.00592, 0.00516, 0.00439, 0.00326, 0.00200, 0.00064,...
    0.00108, 0.00244, 0.00391, 0.00471, 0.00556,  0.00156, 0.00004, 0.00025,...
    0.00044]; 

realy = [-1.00000, -0.97595, -0.96393, -0.94790, -0.93186, -0.91984,...
    -0.90381, -0.88778, -0.86373, -0.83968, -0.81162, -0.77555, -0.73146,...
    -0.68337, -0.61122, 0.57114, 0.65130, 0.71944, 0.77154, 0.82766, 0.88778,...
    0.94790, 0.92385, 0.86774, 0.79960, 0.75150, 0.68337, 0.90381, 0.99599,...
    0.97595, 0.96393, -0.51904, -0.39880, -0.30661, -0.19439, -0.09018, 0.03808,...
    0.20641, 0.37074, 0.48697];
realv = [0.00002, 0.00021, 0.00043, 0.00066, 0.00098, 0.00129, 0.00159,...
    0.00203, 0.00255, 0.00313, 0.00370, 0.00445, 0.00514, 0.00567, 0.00632,...
    0.00655, 0.00592, 0.00516, 0.00439, 0.00326, 0.00200, 0.00064, 0.00108,...
    0.00244, 0.00391, 0.00471, 0.00556, 0.00156, 0.00004, 0.00025, 0.00044,...
    0.00686, 0.00740, 0.00768, 0.00791, 0.00803, 0.00805, 0.00789, 0.00747, 0.00697];


plot(vmag, 'k', 'LineWidth', 5)
hold on 
plot(vr, 'r','LineWidth', 2)
plot(vi, 'b','LineWidth', 2)
plot(magy,magv,'.k','MarkerSize',13)
plot(realy,realv,'.r','MarkerSize',13)
plot(imagy,imagv, '.b','MarkerSize',13)
hold off
xlabel('y','FontSize', 18)
ylabel('v','FontSize', 18)
ylim([0 0.012])
legend('magnitude','real','imaginary', 'magnitude(SH)', 'real(SH)', 'imaginary(SH)','FontSize', 14)
%}

%{
n=200; kx = 1; kz = 1; Re = 10000;
[D0, D1, D2, D4] = Dmat(n); 
[A,B] = pois(n,kx,kz,Re,D0,D1,D2,D4);


[v,e] = eigs(A,B,1,'largestimag');

yd = chebpts(201);

vhat = v(1:201);

vhati = imag(vhat);
vhatr = real(vhat);
vhatmag = abs(vhati)+abs(vhatr);

plot(yd, vhati)
hold on
plot(yd, vhatr)
plot(yd, vhatmag)
hold off



% plot(yd,vhat)
% Cr = diag(real(e));
% Ci = diag(imag(e));
% 
% plot(Cr, Ci, 'o')
% ylim([-1 0.1]) 
%}
