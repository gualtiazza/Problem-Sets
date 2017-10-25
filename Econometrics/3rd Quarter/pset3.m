%% Econometrics - Part II %%
%  Problem Set 3
%  Gualtiero Azzalini

addpath('/Users/Gualtiero/Dropbox/A NYU/Metrics/christensen/psets/3')
addpath('/Users/Gualtiero/Documents/MATLAB/Add-Ons/mfe_toolbox/utility')
clear all;

%% Exercise 1 - First specification for z
clear;
% Import the .csv file
M = csvread('hsdata.csv',1,0);
t = csvread('hsdata.csv',1,0,[1,0,250,0])';
c = csvread('hsdata.csv',1,1,[1,1,250,1])';
d = csvread('hsdata.csv',1,2,[1,2,250,2])';
r = csvread('hsdata.csv',1,3,[1,3,250,3])';

% Balance the observations

for i=1:(length(c)-1)
c_1(i)=c(i+1);
r_1(i)=r(i+1);
d_1(i)=d(i+1);
end
r(250)=[];c(250)=[];d(250)=[]; 
T = length(c);
z=[ones(1,T); exp(r); (1+c)];

% Compute GMM - first step 
W_1=eye(3);
fun1 = @(x)objective(x, W_1, z, T, c_1, r_1, r, c);
x0=[0.9 0.9];
[theta1_gmm, Qval1] = fminsearch(fun1, x0);

% Compute GMM - second step 

% Euler conditions at parameters estimated in step 1
eul = ee(theta1_gmm, z, T, c_1, r_1, r, c);
% 2nd-step weighting matrix - the long-run variance, with lags (b/c dependent data)
W_2=inv(longrunW(eul, T, 12));
fun2 = @(x)objective(x, W_2, z, T, c_1, r_1, r, c);
[theta2_gmm, Qval2]= fminsearch(fun2, theta1_gmm);

% Compute Asymptotic variance (G'S^(-1)G)^(-1) using consistent estimators
% G_hat and W_hat - look at notes for computation of G

G = [sum(exp(-theta2_gmm(2).*c_1+r_1))/T, sum(-c_1.*theta2_gmm(1).*exp(-theta2_gmm(2).*c_1+r_1))/T;
    sum(exp(-theta2_gmm(2).*c_1+r_1+r))/T, sum(-c_1.*theta2_gmm(1).*exp(-theta2_gmm(2).*c_1+r_1+r))/T;
    sum(exp(-theta2_gmm(2).*c_1+r_1).*(1+c))/T, sum(-c_1.*theta2_gmm(1).*exp(-theta2_gmm(2).*c_1+r_1).*(1+c))/T];

eul2 = ee(theta2_gmm, z, T, c_1, r_1, r, c);
W_2=inv(longrunW(eul2, T, 0));
AVar = inv(G'*W_2*G);

delta_CI = [-1.96*sqrt(AVar(1,1)/(T))+theta2_gmm(1) 1.96*sqrt(AVar(1,1)/(T))+theta2_gmm(1)];
gamma_CI = [-1.96*sqrt(AVar(2,2)/(T))+theta2_gmm(2) 1.96*sqrt(AVar(2,2)/(T))+theta2_gmm(2)];

% Hansen-Sargan Test for over-identification
J = 2*T*Qval2;
pvalue_J = 1-chi2cdf(J,1); 

z1_3inst = table(theta1_gmm,theta2_gmm,delta_CI,gamma_CI,J,pvalue_J)

%% Exercise 1 - Second specification for z

clear;

% Import the .csv file
M = csvread('hsdata.csv',1,0);
t = csvread('hsdata.csv',1,0,[1,0,250,0])';
c = csvread('hsdata.csv',1,1,[1,1,250,1])';
d = csvread('hsdata.csv',1,2,[1,2,250,2])';
r = csvread('hsdata.csv',1,3,[1,3,250,3])';

% Balance the observations

for i=1:(length(c)-2)
c_1(i)=c(i+1);
r_1(i)=r(i+1);
d_1(i)=d(i+1);
c_2(i)=c(i+2);
r_2(i)=r(i+2);
d_2(i)=d(i+2);
end
r(250)=[];c(250)=[];d(250)=[]; 
r(249)=[];c(249)=[];d(249)=[]; 
T = length(c);
z=[ones(1,T); exp(r_1); (1+d_1); (1+c_1);
    exp(r); (1+d); (1+c)];

% Compute GMM - first step 
W_1=eye(7);
fun1 = @(x)objective(x, W_1, z, T, c_2, r_2, r_1, c_1);
x0=[0.9 0.9];
[theta1_gmm, Qval1] = fminsearch(fun1, x0);

% Compute GMM - second step 

% Euler conditions at parameters estimated in step 1
eul = ee(theta1_gmm, z, T, c_2, r_2, r_1, c_1);
% 2nd-step weighting matrix - the long-run variance, with lags (b/c dependent data)
W_2=inv(longrunW(eul, T, 12));
fun2 = @(x)objective(x, W_2, z, T, c_2, r_2, r_1, c_1);
[theta2_gmm, Qval2]= fminsearch(fun2, theta1_gmm);

% Compute Asymptotic variance (G'S^(-1)G)^(-1) using consistent estimators
% G_hat and W_hat - look at notes for computation of G

G = [sum(exp(-theta2_gmm(2).*c_2+r_2))/T, sum(-c_2.*theta2_gmm(1).*exp(-theta2_gmm(2).*c_2+r_2))/T;
    sum(exp(-theta2_gmm(2).*c_2+r_2+r_1))/T, sum(-c_2.*theta2_gmm(1).*exp(-theta2_gmm(2).*c_2+r_2+r_1))/T;
    sum(exp(-theta2_gmm(2).*c_2+r_2).*(1+d_1))/T, sum(-c_2.*theta2_gmm(1).*exp(-theta2_gmm(2).*c_2+r_2).*(1+d_1))/T;
    sum(exp(-theta2_gmm(2).*c_2+r_2).*(1+c_1))/T, sum(-c_2.*theta2_gmm(1).*exp(-theta2_gmm(2).*c_2+r_2).*(1+c_1))/T;
    sum(exp(-theta2_gmm(2).*c_2+r_2+r))/T, sum(-c_2.*theta2_gmm(1).*exp(-theta2_gmm(2).*c_2+r_2+r))/T;
    sum(exp(-theta2_gmm(2).*c_2+r_2).*(1+d))/T, sum(-c_2.*theta2_gmm(1).*exp(-theta2_gmm(2).*c_2+r_2).*(1+d))/T;
    sum(exp(-theta2_gmm(2).*c_2+r_2).*(1+c))/T, sum(-c_2.*theta2_gmm(1).*exp(-theta2_gmm(2).*c_2+r_2).*(1+c))/T];

eul2 = ee(theta2_gmm, z, T, c_1, r_1, r, c);
W_2=inv(longrunW(eul2, T, 0));
AVar = inv(G'*W_2*G);

delta_CI = [-1.96*sqrt(AVar(1,1)/(T))+theta2_gmm(1) 1.96*sqrt(AVar(1,1)/(T))+theta2_gmm(1)];
gamma_CI = [-1.96*sqrt(AVar(2,2)/(T))+theta2_gmm(2) 1.96*sqrt(AVar(2,2)/(T))+theta2_gmm(2)];

% Hansen-Sargan Test for over-identification
J = 2*T*Qval2;
pvalue_J = 1-chi2cdf(J,5); 

z2_7inst = table(theta1_gmm,theta2_gmm,delta_CI,gamma_CI,J,pvalue_J)

%% Exercise 2 

clear;
% Import the .csv file
M = csvread('hsdata.csv',1,0);
t = csvread('hsdata.csv',1,0,[1,0,250,0])';
c = csvread('hsdata.csv',1,1,[1,1,250,1])';
d = csvread('hsdata.csv',1,2,[1,2,250,2])';
r = csvread('hsdata.csv',1,3,[1,3,250,3])';

% Balance the observations

for i=1:(length(c)-1)
c_1(i)=c(i+1);
r_1(i)=r(i+1);
end
r(250)=[];c(250)=[];d(250)=[]; 
T = length(c);


% Min the criterion function

%theta=(delta,gamma,alpha,beta,mu_c,sigma_c,sigma_r,rho)
data=[c' c_1' r' r_1'];
b=regress(c_1',[c' r']); %find alpha and beta starting values
x0=[0.9 3 b(1) b(2) mean(c) std(c) std(r) corr(c',r')];
options = optimset; 
[theta_ml] = fminunc(@like_bivnorm,x0,options,data);

delta=theta_ml(1);gamma=theta_ml(2);alpha=theta_ml(3);beta=theta_ml(4);
mu_c=theta_ml(5);sigma_c=theta_ml(6);sigma_r=theta_ml(7);rho=theta_ml(8);

% Compute Hessian and Gradient
H = HessMp(@like_bivnorm,theta_ml,data) % this returns the negative of the Hessian
% Asymptotic variance
AVar = inv(H); 
% Confidence intervals
delta_CI = [theta_ml(1,1)-1.96*sqrt(AVar(1,1)/(T)) theta_ml(1,1)+1.96*sqrt(AVar(1,1)/(T))];
gamma_CI = [theta_ml(1,2)-1.96*sqrt(AVar(2,2)/(T)) theta_ml(1,2)+1.96*sqrt(AVar(2,2)/(T))];

mle_summary = table(delta,gamma,alpha,beta,mu_c,sigma_c,sigma_r,rho,delta_CI,gamma_CI)

