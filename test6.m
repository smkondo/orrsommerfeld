
%test6.m
n=100;
Re = 2000;

[D0,D1,D2,D4]=Dmat(n);

samp = 100;
samp2 = 50;
alpha = linspace(10^-4,1,samp);
beta = linspace(10^-2,10,samp);

sigmav = zeros(samp);

for k=1:samp
    for j=1:samp
        kx = alpha(k);
        kz = beta(j);
        
        [A,B]=pois(n,kx,kz,Re,D0,D1,D2,D4);
        
        omegas = linspace(-2,0,samp2);
        svals = zeros(samp2,1);
        for i = 1:samp2
           
            om = omegas(i);
            H = inv(B\A-om*eye(2*n+2));
            [~,s,~] = svds(H,1,'largest');
            svals(i)=real(s);
        end
        gamma = max(svals);
        sigmav(k,j) = gamma;
    end 
    disp(k)
end




[X,Z] = meshgrid(beta,alpha);
contourf(X,Z,log10(sigmav))
colorbar()
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('k_x')
ylabel('k_z')
title('log_{10}(||H||_2)') 
colormap(turbo)

% %% Plot
% set(0,'DefaultTextInterpreter', 'latex');
% semilogy(omegas,real(svals),'-k');
% ylabel('$\sigma_i$');
% xlabel('$\omega$');
% % hh = legend('$Re = 4000$','$Re = 2000$', '$Re = 1000$','$Re = 500$','location','northwest');
% % hh.Interpreter = 'latex';
% ax = gca;
% ax.YTick = [1 10 100];
% % print('-painters','-dsvg','docs/pics/schmidSols');
