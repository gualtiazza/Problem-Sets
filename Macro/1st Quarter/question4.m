%Question4
%set parameters
clear
addpath(genpath('mfe_toolbox'));

rho_list=[0 0.95 0.99 0.999];
N_list=[2 5 10 20];
T=6500;
sigmaz=1;
k=1;

for i=1:size(rho_list,2)
    rho=rho_list(i);
    p=(1+rho)/2;
    q=p;
    
for j=1:size(N_list,2)
    
    N=N_list(j);
    psi=(N-1)^(1/2)*sigmaz;
    s0=floor((N+1)/2);
    
    [Y P]=arprocess(p,q,psi,N);
    
    [X_approx,s]=markov(P,T+1,s0,Y);X_approx=X_approx';
    [X_filter,e]=armaxfilter_simulate(T,0,1,rho); X_filter = X_filter*sqrt(1-rho^2)
    %% Plot
plots=figure;    
subplot(2,2,1)
h1 = histogram(X_approx(1500:end),N);
h1.FaceColor='r'
p1 = findobj(gca,'Type','patch');
set(p1,'EdgeColor','w','facealpha',0.5);
hold on;
h2 = histogram(X_filter(1500:end),h1.NumBins);
h2.FaceColor='w'
p2 = findobj(gca,'Type','patch');
set(p2,'facealpha',0.5);
legend({'Approx.','Exact'},'FontSize',8); 
legend boxoff
title(['\rho = ', num2str(rho),'; N = ',num2str(N)],'FontSize',16);

subplot(2,2,2)
scatter(0:10,autocov(X_approx(1500:end),10),'rx');
hold on
[acorr,sig2] = acf(rho,0,10);
scatter(0:10,acorr,'bo');
legend({'Approx.','Exact'},'FontSize',8); 
legend boxoff
title(['\rho = ', num2str(rho),'; N = ',num2str(N)],'FontSize',16);

subplot(2,2,3)
plot(1500:T,X_approx(1500:end),'r')
hold on
plot(1500:T,X_filter(1500:end),'b-.')
legend({'Approx.','Exact'},'FontSize',8); 
legend boxoff
title(['\rho = ', num2str(rho),'; N = ',num2str(N)],'FontSize',16);

end
end
for i=1:16
saveas(figure(i),['figure_',num2str(i),'.png'])
end

