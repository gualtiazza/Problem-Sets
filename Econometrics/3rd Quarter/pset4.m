%% Econometrics - Part II %%
%  Problem Set 4
%  Gualtiero Azzalini

addpath('/Users/Gualtiero/Dropbox/A NYU/Metrics/christensen/psets/4')
addpath('/Users/Gualtiero/Documents/MATLAB/Add-Ons/mfe_toolbox/utility')
clear all;

%% Exercise 2 
clear;
% Import the .csv file
M = csvread('hsdata.csv',1,0);
t = csvread('hsdata.csv',1,0,[1,0,250,0])';
c_or = csvread('hsdata.csv',1,1,[1,1,250,1])';
d_or = csvread('hsdata.csv',1,2,[1,2,250,2])';
r_or = csvread('hsdata.csv',1,3,[1,3,250,3])';
data = [c_or;d_or;r_or];

% Balance the observations

for i=1:(length(c_or)-1)    
c_1(i)=c_or(i+1);
r_1(i)=r_or(i+1);
d_1(i)=d_or(i+1);
end
c=c_or;d=d_or;r=r_or;
r(250)=[];c(250)=[];d(250)=[]; 
T = length(c);
z=[ones(1,T); exp(r); (1+c)];

% Compute GMM - first step 
W_1=eye(3);
fun1 = @(x)objective(x, W_1, z, T, c_1, r_1);
x0=[0.9 0.9];
[theta1_gmm, Qval1] = fminsearch(fun1, x0);

% Compute GMM - second step 

% Euler conditions at parameters estimated in step 1
eul = ee(theta1_gmm, z, T, c_1, r_1);
% 2nd-step weighting matrix - the long-run variance, without lags (mds)
W_2=inv(longrunW(eul, T, 0));
fun2 = @(x)objective(x, W_2, z, T, c_1, r_1);
[theta2_gmm, Qval2]= fminsearch(fun2, theta1_gmm);

% Compute Asymptotic variance (G'S^(-1)G)^(-1) using consistent estimators
% G_hat and W_hat - look at notes for computation of G

G = [sum(exp(-theta2_gmm(2).*c_1+r_1))/T, sum(-c_1.*theta2_gmm(1).*exp(-theta2_gmm(2).*c_1+r_1))/T;
    sum(exp(-theta2_gmm(2).*c_1+r_1+r))/T, sum(-c_1.*theta2_gmm(1).*exp(-theta2_gmm(2).*c_1+r_1+r))/T;
    sum(exp(-theta2_gmm(2).*c_1+r_1).*(1+c))/T, sum(-c_1.*theta2_gmm(1).*exp(-theta2_gmm(2).*c_1+r_1).*(1+c))/T];

eul2 = ee(theta2_gmm, z, T, c_1, r_1);
W_2_2=inv(longrunW(eul2, T, 0));
AVar = inv(G'*W_2_2*G); 

delta_CI = [-1.96*sqrt(AVar(1,1)/(T))+theta2_gmm(1) 1.96*sqrt(AVar(1,1)/(T))+theta2_gmm(1)];
gamma_CI = [-1.96*sqrt(AVar(2,2)/(T))+theta2_gmm(2) 1.96*sqrt(AVar(2,2)/(T))+theta2_gmm(2)];

% Test1 H0: delta=0.99, gamma=4 %
null_1 = [0.99 4];
A_theta_1 = [1 0;0 1];

% WALD
wald_1 = T*(theta2_gmm-null_1)*inv(A_theta_1*AVar*A_theta_1')*(theta2_gmm-null_1)';
pvalue_wald_1 = 1-chi2cdf(wald_1,2);

% LM %

theta2_cons_1 = null_1;
% Euler conditions at parameters estimated in step 1
eul_c_1 = ee(theta2_cons_1, z, T, c_1, r_1);
% 2nd-step weighting matrix - the long-run variance, without lags (mds)
W_2_cons_1=inv(longrunW(eul_c_1, T, 0));
Qval_cons2_1 = objective(theta2_cons_1, W_2_cons_1, z, T, c_1, r_1);

% Compute Asymptotic variance (G'S^(-1)G)^(-1) for constrained 

G_cons_1 = [sum(exp(-theta2_cons_1(2).*c_1+r_1))/T, sum(-c_1.*theta2_cons_1(1).*exp(-theta2_cons_1(2).*c_1+r_1))/T;
    sum(exp(-theta2_cons_1(2).*c_1+r_1+r))/T, sum(-c_1.*theta2_cons_1(1).*exp(-theta2_cons_1(2).*c_1+r_1+r))/T;
    sum(exp(-theta2_cons_1(2).*c_1+r_1).*(1+c))/T, sum(-c_1.*theta2_cons_1(1).*exp(-theta2_cons_1(2).*c_1+r_1).*(1+c))/T];

AVar_cons_1 = inv(G_cons_1'*W_2_cons_1*G_cons_1);
eul_c_1 = ee(theta2_cons_1, z, T, c_1, r_1);

lm_1 = T.*(-(T^(-1)).*G_cons_1'*W_2_cons_1*transpose(sum(eul_c_1')))'*AVar_cons_1*...
    (-(T^(-1)).*G_cons_1'*W_2_cons_1*transpose(sum(eul_c_1')));
pvalue_lm_1 = 1-chi2cdf(lm_1,2);
 
% QLR %

qlr_1 = -2*T*(Qval2-Qval_cons2_1);
pvalue_qlr_1 = 1-chi2cdf(qlr_1,2);

tests_3inst_1 = table(wald_1,pvalue_wald_1,lm_1,pvalue_lm_1,qlr_1,pvalue_qlr_1)

% Test2 H0: delta=0.99 %

null_2 = 0.99;
A_theta_2 = [1 0]';

% WALD %
wald_2 = T*(theta2_gmm*A_theta_2-null_2)*inv(A_theta_2'*AVar*A_theta_2)*(theta2_gmm*A_theta_2-null_2)';
pvalue_wald_2 = 1-chi2cdf(wald_2,1);

% LM %
% Compute constrained GMM - first step 
W_1=eye(3);
fun3 = @(x)objective([0.99 x], W_1, z, T, c_1, r_1);
x0=[0.9];
gamma1_cons = fminsearch(fun3, x0);
theta1_cons_1=[0.99 gamma1_cons];

% Compute GMM constrained - second step 

% Euler conditions at parameters estimated in step 1
eul_c_2 = ee(theta1_cons_1, z, T, c_1, r_1);
% 2nd-step weighting matrix - the long-run variance, without lags (mds)
W_2_cons_2=inv(longrunW(eul_c_2, T, 0));
fun4 = @(x)objective([0.99 x], W_2_cons_2, z, T, c_1, r_1);
[gamma2_cons Qval_cons2_2]= fminsearch(fun4, gamma1_cons);
theta2_cons_2=[0.99 gamma2_cons];

% Compute Asymptotic variance (G'S^(-1)G)^(-1) for constrained 

G_cons_2 = [sum(exp(-theta2_cons_2(2).*c_1+r_1))/T, sum(-c_1.*theta2_cons_2(1).*exp(-theta2_cons_2(2).*c_1+r_1))/T;
    sum(exp(-theta2_cons_2(2).*c_1+r_1+r))/T, sum(-c_1.*theta2_cons_2(1).*exp(-theta2_cons_2(2).*c_1+r_1+r))/T;
    sum(exp(-theta2_cons_2(2).*c_1+r_1).*(1+c))/T, sum(-c_1.*theta2_cons_2(1).*exp(-theta2_cons_2(2).*c_1+r_1).*(1+c))/T];

eul_c_2 = ee(theta2_cons_2, z, T, c_1, r_1);
W_2_cons_2=inv(longrunW(eul_c_2, T, 0));
AVar_cons_2 = inv(G_cons_2'*W_2_cons_2*G_cons_2);


lm_2 = T.*(-(T^(-1)).*G_cons_2'*W_2_cons_2*transpose(sum(eul_c_2')))'*AVar_cons_2*...
    (-(T^(-1)).*G_cons_2'*W_2_cons_2*transpose(sum(eul_c_2')));
pvalue_lm_2 = 1-chi2cdf(lm_2,1);

% QLR %

qlr_2 = -2*T*(Qval2-Qval_cons2_2);
pvalue_qlr_2 = 1-chi2cdf(qlr_2,1);

tests_3inst_2 = table(wald_2,pvalue_wald_2,lm_2,pvalue_lm_2,qlr_2,pvalue_qlr_2)

%% Exercise 3

% Part a: reporting values from previous pset

part_a = table(theta2_gmm,delta_CI,gamma_CI)

% Part b: Block Bootstrap 

% Routine elements

reps = 5000; %number of bootstrap repetitions
g_center=(sum(eul2')/T)';
reps_done=0;
tic
for k=[5 2 10];
    for j=1:reps
        
        data_b=block(data,k);
        c_b = data_b(1,:);
        d_b = data_b(2,:);
        r_b = data_b(3,:);

        % Balance the observations
            c = c_b;
            r = r_b;
            d = d_b;
            c_1=c(2:end);
            r_1=r(2:end);
            d_1=d(2:end);
        r(250)=[];c(250)=[];d(250)=[]; 
        T = length(c);
        z=[ones(1,T); exp(r); (1+c)];

        % Compute GMM - first step 
        W_1=eye(3);
        fun1_b = @(x)objectivebs(x,g_center, W_1, z, T, c_1, r_1);
        x0=[0.9 0.9];
        options = optimset('MaxFunEvals',1000);
        [theta1_gmm_b, Qval1_b] = fminsearch(fun1_b, x0,options);

        % Compute GMM - second step 

        % Euler conditions at parameters estimated in step 1
        eul_b = eebs(theta1_gmm_b,g_center, z, T, c_1, r_1);
        % 2nd-step weighting matrix - the long-run variance, without lags (mds)
        W_2_b=inv(longrunW(eul_b, T, 0));
        fun2_b = @(x)objectivebs(x,g_center, W_2_b, z, T, c_1, r_1);
        [theta2_gmm_b, Qval2_b]= fminsearch(fun2_b, theta1_gmm_b,options);

        % Compute Asymptotic variance (G'S^(-1)G)^(-1) using consistent estimators
        % G_hat and W_hat - look at notes for computation of G

        G_b = [sum(exp(-theta2_gmm_b(2).*c_1+r_1))/T, sum(-c_1.*theta2_gmm_b(1).*exp(-theta2_gmm_b(2).*c_1+r_1))/T;
            sum(exp(-theta2_gmm_b(2).*c_1+r_1+r))/T, sum(-c_1.*theta2_gmm_b(1).*exp(-theta2_gmm_b(2).*c_1+r_1+r))/T;
            sum(exp(-theta2_gmm_b(2).*c_1+r_1).*(1+c))/T, sum(-c_1.*theta2_gmm_b(1).*exp(-theta2_gmm_b(2).*c_1+r_1).*(1+c))/T];

        eul2_b = eebs(theta2_gmm_b,g_center, z, T, c_1, r_1);
        W_2_2 = inv(longrunW(eul2_b, T, 0));
        AVar_b = inv(G_b'*W_2_2*G_b); 
        

        theta_bs(j,:)=theta2_gmm_b; 
        Tn_delta_bs(j,:) = (theta2_gmm_b(1)-theta2_gmm(1))/sqrt(AVar_b(1,1)/T);
        Tn_gamma_bs(j,:) = (theta2_gmm_b(2)-theta2_gmm(2))/sqrt(AVar_b(2,2)/T);
        if rem(j/200,1) == 0
            j
        end
    end
    
    eval(['theta_bs_',num2str(k),'=theta_bs;']);
    eval(['Tn_delta_bs_',num2str(k),'=Tn_delta_bs;']);
    eval(['Tn_gamma_bs_',num2str(k),'=Tn_gamma_bs;']);
       
end
toc

figure(1)
subplot(3,2,1);hist(Tn_delta_bs_5,50);title('Tndelta5')
subplot(3,2,2);hist(Tn_gamma_bs_5,50);title('Tngamma5')
subplot(3,2,3);hist(Tn_delta_bs_2,50);title('Tndelta2')
subplot(3,2,4);hist(Tn_gamma_bs_2,50);title('Tngamma2')
subplot(3,2,5);hist(Tn_delta_bs_10,50);title('Tndelta10')
subplot(3,2,6);hist(Tn_gamma_bs_10,50);title('Tngamma10')
saveas(figure(1),'hist.png');


%Bootstrap Confidence intervals
alpha=0.10;
delta_CI_bs_5=[theta2_gmm(1)-quantile(abs(Tn_delta_bs_5),1-alpha/2)*sqrt(AVar(1,1)/T) ...
                theta2_gmm(1)+quantile(abs(Tn_delta_bs_5),1-alpha/2)*sqrt(AVar(1,1)/T)];
delta_CI_bs_2=[theta2_gmm(1)-quantile(abs(Tn_delta_bs_2),1-alpha/2)*sqrt(AVar(1,1)/T) ...
                theta2_gmm(1)+quantile(abs(Tn_delta_bs_2),1-alpha/2)*sqrt(AVar(1,1)/T)];
delta_CI_bs_10=[theta2_gmm(1)-quantile(abs(Tn_delta_bs_10),1-alpha/2)*sqrt(AVar(1,1)/T) ...
                theta2_gmm(1)+quantile(abs(Tn_delta_bs_10),1-alpha/2)*sqrt(AVar(1,1)/T)];

gamma_CI_bs_5=[theta2_gmm(2)-quantile(abs(Tn_gamma_bs_5),1-alpha/2)*sqrt(AVar(2,2)/T) ...
                theta2_gmm(2)+quantile(abs(Tn_gamma_bs_5),1-alpha/2)*sqrt(AVar(2,2)/T)];
gamma_CI_bs_2=[theta2_gmm(2)-quantile(abs(Tn_gamma_bs_2),1-alpha/2)*sqrt(AVar(2,2)/T) ...
                theta2_gmm(2)+quantile(abs(Tn_gamma_bs_2),1-alpha/2)*sqrt(AVar(2,2)/T)];
gamma_CI_bs_10=[theta2_gmm(2)-quantile(abs(Tn_gamma_bs_10),1-alpha/2)*sqrt(AVar(2,2)/T) ...
                theta2_gmm(2)+quantile(abs(Tn_gamma_bs_10),1-alpha/2)*sqrt(AVar(2,2)/T)];    
            
delta_CI_bs=table(delta_CI_bs_5,delta_CI_bs_2,delta_CI_bs_10)
gamma_CI_bs=table(gamma_CI_bs_5,gamma_CI_bs_2,gamma_CI_bs_10)

quantiles_delta=table(quantile(abs(Tn_delta_bs_5),1-alpha/2),quantile(abs(Tn_delta_bs_2),1-alpha/2),quantile(abs(Tn_delta_bs_10),1-alpha/2))
quantiles_gamma=table(quantile(abs(Tn_gamma_bs_5),1-alpha/2),quantile(abs(Tn_gamma_bs_2),1-alpha/2),quantile(abs(Tn_gamma_bs_10),1-alpha/2))




