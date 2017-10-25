%% Econometrics - Part II %%
%  Problem Set 5
%  Gualtiero Azzalini

addpath('/Users/Gualtiero/Dropbox/A NYU/Metrics/christensen/psets/5')
addpath('/Users/Gualtiero/Documents/MATLAB/Add-Ons/mfe_toolbox/utility')
clear all;

%% Exercise 2

n = 10;
k = 1000;
reps = 3000;
alpha=0.05;
theta0 = 1;
tic
for j=1:k
    x = randn(n,1).*sqrt(6);
    theta = exp(mean(x));
    bs=randi([1 10],n,reps);
    x_star(:,:)=x(bs(:,:));
    for i=1:reps
        Tn(i,1) = (sqrt(n)*(exp(mean(x_star(:,i)))-theta)*(std(x_star(:,i))*theta)^(-1))';
    end
    Tn_abs = abs(Tn);
    z(j,1) = quantile(Tn_abs,1-alpha);
    CI_bs(j,1) = theta-z(j,1)*(std(x)*theta)*sqrt(n)^(-1);
    CI_asy(j,1) = theta-1.96*(std(x)*theta)*sqrt(n)^(-1);
    CI_bs(j,2) = theta+z(j,1)*(std(x)*theta)*sqrt(n)^(-1);
    CI_asy(j,2) = theta+1.96*(std(x)*theta)*sqrt(n)^(-1);
    if rem(j/200,1) == 0
        j
    end
    if theta0<CI_bs(j,1) | theta0>CI_bs(j,2)
        dummy_bs(j,1)=0;
    else
        dummy_bs(j,1)=1;
    end
        
    if theta0<CI_asy(j,1) | theta0>CI_asy(j,2)
        dummy_asy(j,1)=0;
    else
        dummy_asy(j,1)=1;
    end
end
toc

%Percentage theta0 is in the interval

perc_bs = 100*sum(dummy_bs)/k; 
perc_asy = 100*sum(dummy_asy)/k; 

summary = table(perc_bs,perc_asy)


