%% Problem Set 3 code - Gualtiero Azzalini %%

%% Exercise 1

addpath('/Users/Gualtiero/Dropbox/A NYU/Metrics/crump/psets/3')
addpath('/Users/Gualtiero/Documents/MATLAB/Add-Ons/mitbook/utilities_tom')

clear; clc; close all;

% Set Parameters 

beta0   = 1;    beta1   = 1;    beta2   = 1;
gamma0  = 0.1;  gamma1  = 0.1;  gamma2  = 0.1;
rho     = 0.5;
T       = 200;
N       = 1e4;
k       = 3;
 
% Set model
e       = randn(T,N);
x1      = sqrt(12)*rand(T,N);
x2      = sqrt(1-rho^2)*sqrt(12)*rand(T,N) + rho*x1;
sigma   = exp(gamma0 + gamma1*x1 + gamma2*x2);
y       = beta0 + beta1*x1 + beta2*x2 + e;
i       = ones(T,1);
a       = [0; 1; 0];

% Containers
beta        = NaN(k,N);
s2          = NaN(1,N);
Qxx         = NaN(k,k,N);
Omg_loo     = NaN(k,k,N); 
Omg_bar     = NaN(k,k,N);
Vbeta_0     = NaN(k,k,N);  
Vbeta_loo   = NaN(k,k,N); 
Vbeta_bar   = NaN(k,k,N);
tstat_0     = NaN(1,N);     
tstat_loo   = NaN(1,N);   
tstat_bar   = NaN(1,N);

for n = 1:N

    X = [i, x1(:,n), x2(:,n)];
    Y = y(:,n);

    beta(:,n)       = X \ Y; 
    e_hat           = Y - X*beta(:,n);
    hii             = diag(X*inv(X'*X)*X');
    e_loo           = e_hat./(1-hii);
    e_bar           = e_hat./sqrt(1-hii);
    s2(:,n)         = (T/(T-k))*mean(e_hat.^2);
    Qxx(:,:,n)      = (1/T)*X'*X;
    Omg_loo(:,:,n)  = (1/T)*(X.*e_loo)'*(X.*e_loo);
    Omg_bar(:,:,n)  = (1/T)*(X.*e_bar)'*(X.*e_bar);

    Vbeta_0(:,:,n)  = inv(Qxx(:,:,n))*inv(Qxx(:,:,n));
    Vbeta_loo(:,:,n)= inv(Qxx(:,:,n))*Omg_loo(:,:,n)*inv(Qxx(:,:,n));
    Vbeta_bar(:,:,n)= inv(Qxx(:,:,n))*Omg_bar(:,:,n)*inv(Qxx(:,:,n));

    

    tstat_0(:,n)    = (beta(2,n)-beta1)/sqrt((1/T)*s2(n)*a'*Vbeta_0(:,:,n)*a);
    tstat_loo(:,n)  = (beta(2,n)-beta1)/sqrt((1/T)*s2(n)*a'*Vbeta_loo(:,:,n)*a);
    tstat_bar(:,n)  = (beta(2,n)-beta1)/sqrt((1/T)*s2(n)*a'*Vbeta_bar(:,:,n)*a);
end

 

% Plot t-statistics and compare to quantiles of the Normal distribution

figure(1)
subplot(1,2,1)
hist(tstat_0,50)
title('Histogram: $$ t_n^0 $$','interpreter','latex')
subplot(1,2,2)
h = qqplot(tstat_0);
set(h(1),'marker','^','markersize',8,'markeredgecolor',[0.6 0 0]);
set(h(2),'linewidth',2,'color',[0 0 0]);
set(h(3),'linewidth',2,'color',[0 0 0]);
title('QQ Plot: $$ t_n^0 $$','interpreter','latex')
h=gcf;
set(h,'PaperPositionMode','auto');         
set(h,'PaperOrientation','landscape');
set(h,'Position',[50 50 1200 800]);
saveas(figure(1),'t0.png')
 
figure(2)
subplot(1,2,1)
hist(tstat_loo,50)
title('Histogram: $$ \tilde{t}_n $$','interpreter','latex')
subplot(1,2,2)
h=qqplot(tstat_loo);
set(h(1),'marker','^','markersize',8,'markeredgecolor',[0.6 0 0]);
set(h(2),'linewidth',2,'color',[0 0 0]);
set(h(3),'linewidth',2,'color',[0 0 0]);
title('QQ Plot: $$ \tilde{t}_n $$','interpreter','latex')
h=gcf;
set(h,'PaperPositionMode','auto');         
set(h,'PaperOrientation','landscape');
set(h,'Position',[50 50 1200 800]);
saveas(figure(2),'t_tilde.png')

figure(3)
subplot(1,2,1)
hist(tstat_bar,50)
title('Histogram: $$ \bar{t}_n $$','interpreter','latex')
subplot(1,2,2)
h=qqplot(tstat_bar);
set(h(1),'marker','^','markersize',8,'markeredgecolor',[0.6 0 0]);
set(h(2),'linewidth',2,'color',[0 0 0]);
set(h(3),'linewidth',2,'color',[0 0 0]);
title('QQ Plot: $$ \bar{t}_n $$','interpreter','latex')
h=gcf;
set(h,'PaperPositionMode','auto');         
set(h,'PaperOrientation','landscape');
set(h,'Position',[50 50 1200 800]);
saveas(figure(3),'t_bar.png')
 

% Quantile table

qtiles = [.01 .025 .05 .1 .2 .3 .4 .5 .6 .7 .8 .9 .95 .975 .99]';
[qtiles, quantile(tstat_0,qtiles), quantile(tstat_loo,qtiles), quantile(tstat_bar,qtiles)]

 

%% Exercise #2

% Set parameters

alpha_y = 0.98;
alpha_x = 0.97;
y0      = 0;
x0      = 0;
T       = 200;
N       = 1e4;
k       = 2;
beta1   = 0;

% Simulate T draws N times

u       = randn(T,N);
v       = randn(T,N);
y       = zeros(T+1,N);
x       = zeros(T+1,N);
i       = ones(T,1);
a       = [0; 1;];

for t   = 2:(T+1)
    y(t,:) = alpha_y*y(t-1,:) + u(t-1,:);
    x(t,:) = alpha_x*x(t-1,:) + v(t-1,:);    
end

% Set containers
beta        = NaN(k,N);
s2          = NaN(1,N);
Qxx         = NaN(k,k,N);
Vbeta_0     = NaN(k,k,N); 
tstat_0     = NaN(1,N);   

for n = 1:N

    X = [i, x(2:end,n)];
    Y = y(2:end,n);
    beta(:,n)       = X \ Y; 
    e_hat           = Y - X*beta(:,n);
    s2(:,n)         = (T/(T-k))*mean(e_hat.^2);
    Qxx(:,:,n)      = (1/T)*X'*X;

    Vbeta_0(:,:,n)  = inv(Qxx(:,:,n))*inv(Qxx(:,:,n));    
    tstat_0(:,n)    = (beta(2,n)-beta1)/sqrt((1/T)*s2(n)*a'*Vbeta_0(:,:,n)*a);
end

figure(4)
subplot(1,2,1)
hist(tstat_0,50)
title('Histogram: t_n^0')
subplot(1,2,2)
h = qqplot(tstat_0);
set(h(1),'marker','^','markersize',8,'markeredgecolor',[0.6 0 0]);
set(h(2),'linewidth',2,'color',[0 0 0]);
set(h(3),'linewidth',2,'color',[0 0 0]);
title('QQ Plot: t_n^0')
h=gcf;
set(h,'PaperPositionMode','auto');         
set(h,'PaperOrientation','landscape');
set(h,'Position',[50 50 1200 800]);
saveas(figure(4),'t0_exercise2.png')

% Quantile table

qtiles = [.01 .025 .05 .1 .2 .3 .4 .5 .6 .7 .8 .9 .95 .975 .99]';
[qtiles, quantile(tstat_0,qtiles)]

