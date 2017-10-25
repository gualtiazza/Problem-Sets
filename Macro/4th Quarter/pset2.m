%% Macroeconomics II - Part II %%
%  Problem Set 2
%  Gualtiero Azzalini

addpath('/Users/Gualtiero/Dropbox/A NYU/Macro/gilchrist/psets/2/mycode')
addpath('/Users/Gualtiero/Documents/MATLAB/Add-Ons/mfe_toolbox/utility')
addpath('/Users/Gualtiero/Documents/MATLAB/Add-Ons/mfe_toolbox/timeseries')
addpath('/Users/Gualtiero/Documents/MATLAB/Add-Ons/habit_homework_cheb')

clear all;

% Import Dataset

data    = xlsread('data_homework_nyu.xls');
years   = data(:,1);	
C_ND  	= data(:,2);
P_ND	= data(:,3);
exr     = data(:,4)/100;
smb     = data(:,5)/100;
hml     = data(:,6)/100;
rf      = data(:,7)/100;
s1      = data(:,8)/100;
s2      = data(:,9)/100;
s3      = data(:,10)/100;
b1      = data(:,11)/100;
b2      = data(:,12)/100;
b3      = data(:,13)/100;

T       = length(years);

%% Question 1

r = exr+rf; % Market Return

% Compute C growth and Pi 

c_g = log(C_ND(2:end,1)./C_ND(1:end-1,1));
pi  = log(P_ND(2:end,1)./P_ND(1:end-1,1));


% Real returns

rf_r = rf(1:end-1) - pi;
r_r  = r(1:end-1) - pi;
s1_r = s1(1:end-1) - pi;
s2_r = s2(1:end-1) - pi;
s3_r = s3(1:end-1) - pi;
b1_r = b1(1:end-1) - pi;
b2_r = b2(1:end-1) - pi;
b3_r = b3(1:end-1) - pi;

% Mean, sd, autocorrelation
means = [mean(c_g),mean(rf_r),mean(r_r),mean(s1_r),mean(s2_r),...
                    mean(s3_r),mean(b1_r),mean(b2_r),mean(b3_r)];
stds  = [std(c_g),std(rf_r),std(r_r),std(s1_r),std(s2_r),...
                    std(s3_r),std(b1_r),std(b2_r),std(b3_r)];
autocorrs = [sacf(c_g,1),sacf(rf_r,1),sacf(r_r,1),...
                    sacf(s1_r,1),sacf(s2_r,1),sacf(s3_r,1),...
                    sacf(b1_r,1),sacf(b2_r,1),sacf(b3_r,1)];  

Q1Part2_mean = table(means);
Q1Part2_std  = table(stds);
Q1Part2_autocorr = table(autocorrs);
close all
                
% Mean, sd of Log-returns

lmeans = [mean(log(1+rf_r)),mean(log(1+r_r)),mean(log(1+s1_r)),...
        mean(log(1+s2_r)),mean(log(1+s3_r)),mean(log(1+b1_r)),...
        mean(log(1+b2_r)),mean(log(1+b3_r))];
lstds  = [std(log(1+rf_r)),std(log(1+r_r)),std(log(1+s1_r)),...
         std(log(1+s2_r)),std(log(1+s3_r)),std(log(1+b1_r)),...
         std(log(1+b2_r)),std(log(1+b3_r))];

Q1Part3_mean = table(lmeans);
Q1Part3_std  = table(lstds);

% Computing the Betas
portfolio = [s1_r s2_r s3_r b1_r b2_r b3_r];
port_mean = [mean(s1_r-rf_r),mean(s2_r-rf_r),mean(s3_r-rf_r),...
            mean(b1_r-rf_r),mean(b2_r-rf_r),mean(b3_r-rf_r)]';
X1 = [ones(length(r_r),1) r_r];
X2 = [ones(length(c_g),1) c_g];
for j=1:6
    beta_market(j,:) = regress(portfolio(:,j),X1)';
    beta_cons(j,:)   = regress(portfolio(:,j),X2)';
end

beta_cons(7,:) = regress(r_r,X2)';
port_mean_c = [port_mean; mean(r_r)];

figure(1)
subplot(1,2,1)
scatter(beta_market(:,2),port_mean,'filled');title('CAPM');
ylabel('Mean Excess Return');xlabel('\beta_{mkt}');lsline;
axis([min(beta_market(:,2)) max(beta_market(:,2)) 0 0.3]);
subplot(1,2,2)
scatter(beta_cons(:,2)/beta_cons(7,2),port_mean_c,'filled');title('CCAPM');
ylabel('Mean Excess Return');xlabel('\beta_{cons}');lsline;
saveas(figure(1),'figure1.png');

%% Question 2

% H-J Bounds

R1 = [(rf_r+1) (r_r+1)];
ER1 = mean(R1)';
VR1 = cov(R1);

R2 = [(s1_r+1) (b3_r+1) (rf_r+1) (r_r+1) (s2_r+1) (b2_r+1) (s3_r+1) (b1_r+1)];
ER2 = mean(R2)';
VR2 = cov(R2);

p=1;
for j=0:.01:3
VM1(p,1) = sqrt((ones(length(ER1),1)-j*ER1)'*inv(VR1)*(ones(length(ER1),1)-j*ER1));
VM2(p,1) = sqrt((ones(length(ER2),1)-j*ER2)'*inv(VR2)*(ones(length(ER2),1)-j*ER2));
p=p+1;
end

figure(2)
plot(0:.01:3,VM1);
hold on
plot(0:.01:3,VM2,'r--');
ylabel('\sigma(M)');xlabel('E(M)');
legend('Rf&R','All Assets');
saveas(figure(2),'figure2.png');
hold off

% CRRA
beta = 0.998;
gamma = [0.5 1 2 4 6 8 10 15 20 30 34 35 36 40 44 45 46 47 50 60 70 90 100 150 200]';

for jj = 1:length(gamma)
    for tt=1:T-1
    M(tt,1) = beta*(C_ND(tt+1,1)/C_ND(tt,1))^(-gamma(jj,1));
    end
    E_SDF(jj,1) = mean(M);
    DEV_SDF(jj,1) = std(M);
end

figure(3)
plot(0:.01:3,VM1);
hold on
plot(0:.01:3,VM2,'r--');
ylabel('\sigma(M)');xlabel('E(M)');
hold on 
scatter(E_SDF,DEV_SDF);
labelpoints(E_SDF, DEV_SDF, gamma);
plot(E_SDF,DEV_SDF,'-.');
axis([0 2 0 40]);
legend('Rf&R','All Assets','gamma','E(M)/\sigma(M)');
saveas(figure(3),'figure3.png');
hold off

% Habit model

theta = [0.9 0.95]';
delta = [0.1 0.01]';
number = 1;
H(1,1) = C_ND(1,1);
gamma1 = [0.5 1 2 3 4 5 6 7 8]';
for kk = 1:length(delta)
    for k = 1:length(theta)
        for ttt=2:T
            H(ttt,1) = (1-delta(kk,1))*H(ttt-1,1)+delta(kk,1)*C_ND(ttt-1,1);
        end
        for jj = 1:length(gamma1)
            for tt=1:T-1
                M_h(tt,1) = beta*theta(k,1)*((C_ND(tt+1,1)-theta(k,1)*H(tt+1,1))/...
                    (C_ND(tt,1)-theta(k,1)*H(tt,1)))^(-gamma1(jj,1));
            end
            E_SDF_h(jj,1) = mean(M_h);
            DEV_SDF_h(jj,1) = std(M_h);
        end
        figure(4)
        subplot(2,2,number)
        plot(0:.01:3,VM1);
        hold on
        plot(0:.01:3,VM2,'r--');
        ylabel('\sigma(M)');xlabel('E(M)');
        hold on 
        scatter(E_SDF_h,DEV_SDF_h);
        labelpoints(E_SDF_h, DEV_SDF_h, gamma1);
        title(['\delta=',num2str(delta(kk,1)),',','\theta=',num2str(theta(k,1))]);
        plot(E_SDF_h,DEV_SDF_h,'-.');
        axis([0 2 0 40]);
        legend('Rf&R','All Assets','gamma','E(M)/\sigma(M)');
        number = number+1;
    end
end
saveas(figure(4),'figure4.png');
hold off

%% Question 3

beta  = 0.998;
ngamma = 50;
gamma  = linspace(0,ngamma,ngamma)';

ret_1 = exr(1:end-1,1);
ret_2 = r_r;
ret_3 = rf_r;

% Only one return - parts 1,2,3

for d = 1:3
    if d ==1
        cons = 0;
    else
        cons = 1;
    end
    eval(['ret','=ret_',num2str(d),';']);
    for rr = 1:length(gamma)
        [Q, g] = obj(1,beta,C_ND,rr,ret,cons);
        eval(['Q_',num2str(d),'(rr,1)','=Q',';']);
    end    
    gamma_hat = fminsearch(@(x)obj(1,beta,C_ND,x,ret,cons),0);
    [value, err] = obj(1,beta,C_ND,gamma_hat,ret,cons);
    eval(['gamma_hat_',num2str(d),'=gamma_hat',';']);
    eval(['value_',num2str(d),'=value',';']);
    eval(['err_',num2str(d),'=err',';']);
    
    % Asy Variance (G'S^(-1)G)^(-1)
    if d==1
    for j=1:size(ret,2)
    G(j,1)  = mean(-beta.*(C_ND(2:end,1)./C_ND(1:end-1,1)).^(-gamma_hat)...
            .*(ret(:,j)).*log(C_ND(2:end,1)./C_ND(1:end-1,1))); 
    end
    else
    for j=1:size(ret,2)
    G(j,1)  = mean(-beta.*(C_ND(2:end,1)./C_ND(1:end-1,1)).^(-gamma_hat)...
            .*(1+ret(:,j)).*log(C_ND(2:end,1)./C_ND(1:end-1,1))); 
    end
    end

    AVar = inv(G'*1*G);
    gamma_CI = [-1.96*sqrt(AVar(1,1)/length(ret))+gamma_hat ...
                1.96*sqrt(AVar(1,1)/length(ret))+gamma_hat];
    eval(['gamma_CI_',num2str(d),'=gamma_CI',';']);
    
    figure(5)
    subplot(2,3,d)
    plot(gamma,eval(['Q_',num2str(d)]));
    title(['Part',num2str(d)]);
    xlim([0 ngamma]);
    ylabel('Q');xlabel('\gamma');
    hold on
end
gammas_1  = table(gamma_hat_1,gamma_hat_2,gamma_hat_3)
errs_1    = table(err_1,err_2,err_3)

% More than 1 return
% Part 4
ret_4 = [rf_r r_r];
d=4;
for rr = 1:length(gamma)
    [Q, g] = obj(1,beta,C_ND,rr,ret_4,1);
    eval(['Q_',num2str(d),'(rr,1)','=Q',';']);
end   
% First step
gamma_hat_1 = fminsearch(@(x)obj(eye(2,2),beta,C_ND,x,ret_4,1),0);
for j=1:size(ret_4,2)
    g1(:,j)   = beta.*(C_ND(2:end,1)./C_ND(1:end-1,1)).^(-gamma_hat_1)...
                    .*(1+ret_4(:,j))-1;
end
W21       = inv(longrunW(g1',length(ret_4),40));
% Second step
gamma_hat_2 = fminsearch(@(x)obj(W21,beta,C_ND,x,ret_4,1),gamma_hat_1);
for j=1:size(ret_4,2)
    g2(:,j)   = beta.*(C_ND(2:end,1)./C_ND(1:end-1,1)).^(-gamma_hat_2)...
                    .*(1+ret_4(:,j))-1;
end
W2       = inv(longrunW(g2',length(ret_4),40));
for j=1:size(ret_4,2)
    G(j,1)  = mean(-beta.*(C_ND(2:end,1)./C_ND(1:end-1,1)).^(-gamma_hat_2)...
            .*(1+ret_4(:,j)).*log(C_ND(2:end,1)./C_ND(1:end-1,1))); 
end
AVar = inv(G'*W2*G);
for rr = 1:length(gamma)
    [Q, g] = obj(W21,beta,C_ND,rr,ret_4,1);
    eval(['Qopt_',num2str(d),'(rr,1)','=Q',';']);
end  
subplot(2,3,d)
plot(gamma,eval(['Q_',num2str(d)]));
hold on
title(['Part',num2str(d)]);
xlim([0 ngamma]);
ylabel('Q');xlabel('\gamma');
hold on
gamma_CI = [-1.96*sqrt(AVar(1,1)/length(ret_4))+gamma_hat_2 ...
                1.96*sqrt(AVar(1,1)/length(ret_4))+gamma_hat_2];
[value, err] = obj(W2,beta,C_ND,gamma_hat_2,ret_4,1);
J = 2*length(ret_4)*value;
pvalue_J = 1-chi2cdf(J,1); 

eval(['gamma_CI_',num2str(d),'=gamma_CI',';']);
eval(['gamma_hat_1_',num2str(d),'=gamma_hat_1',';']);
eval(['gamma_hat_2_',num2str(d),'=gamma_hat_2',';']);
eval(['err_',num2str(d),'=err',';']);
eval(['J_',num2str(d),'=J',';']);
eval(['pvalue_J_',num2str(d),'=pvalue_J',';']);


% Part 5
ret_5 = [rf_r exr(1:end-1,1)];
d=5;
for rr = 1:length(gamma)
    [Q, g] = objmod(1,beta,C_ND,rr,ret_5);
    eval(['Q_',num2str(d),'(rr,1)','=Q',';']);
end   

% First step
gamma_hat_1 = fminsearch(@(x)objmod(eye(2,2),beta,C_ND,x,ret_5),0);
g1(:,1) = (beta.*(C_ND(2:end,1)./C_ND(1:end-1,1)).^(-gamma_hat_1)...
                    .*(1+ret_5(:,1)))-1;
g1(:,2) = beta.*(C_ND(2:end,1)./C_ND(1:end-1,1)).^(-gamma_hat_1)...
                    .*(ret_5(:,2));
W21       = inv(longrunW(g1',length(ret_5),40));
% Second step
gamma_hat_2 = fminsearch(@(x)objmod(W21,beta,C_ND,x,ret_5),gamma_hat_1);
g2(:,1) = (beta.*(C_ND(2:end,1)./C_ND(1:end-1,1)).^(-gamma_hat_2)...
                    .*(1+ret_5(:,1)))-1;
g2(:,2) = beta.*(C_ND(2:end,1)./C_ND(1:end-1,1)).^(-gamma_hat_2)...
                    .*(ret_5(:,2));
W2       = inv(longrunW(g2',length(ret_5),40));
G(1,1)  = mean(-beta.*(C_ND(2:end,1)./C_ND(1:end-1,1)).^(-gamma_hat_2)...
            .*(1+ret_5(:,1)).*log(C_ND(2:end,1)./C_ND(1:end-1,1))); 
G(2,1)  = mean(-beta.*(C_ND(2:end,1)./C_ND(1:end-1,1)).^(-gamma_hat_2)...
            .*(ret_5(:,2)).*log(C_ND(2:end,1)./C_ND(1:end-1,1))); 
AVar = inv(G'*W2*G);
for rr = 1:length(gamma)
    [Q, g] = objmod(W21,beta,C_ND,rr,ret_5);
    eval(['Qopt_',num2str(d),'(rr,1)','=Q',';']);
end  
subplot(2,3,d)
plot(gamma,eval(['Q_',num2str(d)]));
hold on
title(['Part',num2str(d)]);
xlim([0 ngamma]);
ylabel('Q');xlabel('\gamma');
hold on
gamma_CI = [-1.96*sqrt(AVar(1,1)/length(ret_5))+gamma_hat_2 ...
                1.96*sqrt(AVar(1,1)/length(ret_5))+gamma_hat_2];
[value, err] = objmod(W2,beta,C_ND,gamma_hat_2,ret_5);
J = 2*length(ret_5)*value;
pvalue_J = 1-chi2cdf(J,1); 

eval(['gamma_CI_',num2str(d),'=gamma_CI',';']);
eval(['gamma_hat_1_',num2str(d),'=gamma_hat_1',';']);
eval(['gamma_hat_2_',num2str(d),'=gamma_hat_2',';']);
eval(['err_',num2str(d),'=err',';']);
eval(['J_',num2str(d),'=J',';']);
eval(['pvalue_J_',num2str(d),'=pvalue_J',';']);

% Part 6
ret_6 = [rf_r r_r s1_r s2_r s3_r b1_r b2_r b3_r];
d=6;
for rr = 1:length(gamma)
    [Q, g] = obj(1,beta,C_ND,rr,ret_6,1);
    eval(['Q_',num2str(d),'(rr,1)','=Q',';']);
end   
% First step
gamma_hat_1 = fminsearch(@(x)obj(eye(8,8),beta,C_ND,x,ret_6,1),0);
for j=1:size(ret_6,2)
    g1(:,j)   = beta.*(C_ND(2:end,1)./C_ND(1:end-1,1)).^(-gamma_hat_1)...
                    .*(1+ret_6(:,j))-1;
end
W21       = inv(longrunW(g1',length(ret_5),40));
% Second step
gamma_hat_2 = fminsearch(@(x)obj(W21,beta,C_ND,x,ret_6,1),gamma_hat_1);
for j=1:size(ret_6,2)
    g2(:,j)   = beta.*(C_ND(2:end,1)./C_ND(1:end-1,1)).^(-gamma_hat_2)...
                    .*(1+ret_6(:,j))-1;
end
W2       = inv(longrunW(g2',length(ret_6),40));
for j=1:size(ret_6,2)
    G(j,1)  = mean(-beta.*(C_ND(2:end,1)./C_ND(1:end-1,1)).^(-gamma_hat_2)...
            .*(1+ret_6(:,j)).*log(C_ND(2:end,1)./C_ND(1:end-1,1))); 
end
AVar = inv(G'*W2*G);
for rr = 1:length(gamma)
    [Q, g] = obj(W21,beta,C_ND,rr,ret_6,1);
    eval(['Qopt_',num2str(d),'(rr,1)','=Q',';']);
end  
subplot(2,3,d)
plot(gamma,eval(['Q_',num2str(d)]));
hold on
title(['Part',num2str(d)]);
xlim([0 ngamma]);
ylabel('Q');xlabel('\gamma');
hold on
gamma_CI = [-1.96*sqrt(AVar(1,1)/length(ret_6))+gamma_hat_2 ...
                1.96*sqrt(AVar(1,1)/length(ret_6))+gamma_hat_2];
[value, err] = obj(W2,beta,C_ND,gamma_hat_2,ret_6,1);
J = 2*length(ret_6)*value;
pvalue_J = 1-chi2cdf(J,7); 

eval(['gamma_CI_',num2str(d),'=gamma_CI',';']);
eval(['gamma_hat_1_',num2str(d),'=gamma_hat_1',';']);
eval(['gamma_hat_2_',num2str(d),'=gamma_hat_2',';']);
eval(['err_',num2str(d),'=err',';']);
eval(['J_',num2str(d),'=J',';']);
eval(['pvalue_J_',num2str(d),'=pvalue_J',';']);

saveas(figure(5),'figure5.png')

figure(6)
subplot(2,2,1)
plot(gamma,Qopt_4,'--');
xlim([0 ngamma]);
ylabel('Q');xlabel('\gamma');
title('Part4(optimal weight)');
subplot(2,2,2)
plot(gamma,Qopt_5,'--');
xlim([0 ngamma]);
ylabel('Q');xlabel('\gamma');
title('Part5(optimal weight)');
subplot(2,2,3)
plot(gamma,Qopt_6,'--');
xlim([0 ngamma]);
ylabel('Q');xlabel('\gamma');
title('Part6(optimal weight)');
saveas(figure(6),'figure6.png')

gammas_2  = table(gamma_hat_1_4,gamma_hat_2_4,gamma_hat_1_5,...
            gamma_hat_2_5,gamma_hat_1_6,gamma_hat_2_6)
err_4;
err_5;
err_6
CIs     = table(gamma_CI_1,gamma_CI_2,gamma_CI_3,gamma_CI_4,gamma_CI_5...
                ,gamma_CI_6)
HansenSargan = table(J_4,J_5,J_6,pvalue_J_4,pvalue_J_5,pvalue_J_6)

%% Question 4

% Part a

sigma = 1.5/100;
gamma = 2;
phi   = 0.87;
g     = 1.89/100;
beta  = 0.89;

sbar  = sigma*sqrt(gamma/(1-phi));
lsbar = log(sbar);
smax  = lsbar + (1-sbar^2)/2;

s(1,1) = lsbar;

for j=1:length(C_ND)-1
    if s(j,1)<smax || s(j,1)==smax
        lambda(j,1) = (sbar)^(-1)*sqrt(1-2*(s(j,1)-lsbar))-1;
    else
        lambda(j,1) = 0;
    end
    s(j+1,1) = (1-phi)*lsbar+phi*s(j,1)+lambda(j,1)*(log(C_ND(j+1,1)/C_ND(j,1))-g);
end
S      = exp(s);
SDF_CC = beta*((S(2:end,1)./S(1:end-1,1)).*(C_ND(2:end,1)./C_ND(1:end-1,1))).^(-gamma);

figure(7);
subplot(2,2,1)
plot(years(22:end-1),S(22:end-1));title('Surplus Ratio - Post WW2');
xlim([years(22,1) years(end-1,1)]);
subplot(2,2,2)
plot(years(22:end-1),SDF_CC(22:end));title('SDF - Post WW2');
xlim([years(22,1) years(end-1,1)]);
subplot(2,2,3)
plot(years(1:end-1),S(1:end-1),'-.');title('Surplus Ratio');
xlim([years(1,1) years(end-1,1)]);
subplot(2,2,4)
plot(years(1:end-1),SDF_CC(1:end),'-.');title('SDF');
xlim([years(1,1) years(end-1,1)]);
saveas(figure(7),'figure7.png');


% Part b

summary=table(mean(SDF_CC(22:end)),mean(S(22:end)),std(SDF_CC(22:end)),...
            std(S(22:end)),sacf(SDF_CC(22:end),1),sacf(S(22:end),1))
close all
C1=cov(S(22:end),exr(22:end));C1SDF=cov(SDF_CC(22:end),exr(22:end-1));
C2=cov(S(22:end-1),rf_r(22:end));C2SDF=cov(SDF_CC(22:end),rf_r(22:end));
C3=cov(S(22:end-1),r_r(22:end));C3SDF=cov(SDF_CC(22:end),r_r(22:end));
C4=cov(S(22:end-1),s1_r(22:end));C4SDF=cov(SDF_CC(22:end),s1_r(22:end));
C5=cov(S(22:end-1),s2_r(22:end));C5SDF=cov(SDF_CC(22:end),s2_r(22:end));
C6=cov(S(22:end-1),s3_r(22:end));C6SDF=cov(SDF_CC(22:end),s3_r(22:end));
C7=cov(S(22:end-1),b1_r(22:end));C7SDF=cov(SDF_CC(22:end),b1_r(22:end));
C8=cov(S(22:end-1),b2_r(22:end));C8SDF=cov(SDF_CC(22:end),b2_r(22:end));
C9=cov(S(22:end-1),b3_r(22:end));C9SDF=cov(SDF_CC(22:end),b3_r(22:end));

covS=table(C1(1,2),C2(1,2),C3(1,2),C4(1,2),C5(1,2),C6(1,2),C7(1,2),C8(1,2),C9(1,2))
covSDF=table(C1SDF(1,2),C2SDF(1,2),C3SDF(1,2),C4SDF(1,2),...
             C5SDF(1,2),C6SDF(1,2),C7SDF(1,2),C8SDF(1,2),C9SDF(1,2))

% Computing the Betas
%All sample
portfolio_s = [ r_r s1_r s2_r s3_r b1_r b2_r b3_r];
port_mean_s = [ mean(r_r-rf_r) mean(s1_r-rf_r),...
            mean(s2_r-rf_r),mean(s3_r-rf_r),...
            mean(b1_r-rf_r),mean(b2_r-rf_r),mean(b3_r-rf_r)]';
Y1 = [ones(length(SDF_CC),1) SDF_CC];
Y2 = [ones(length(SDF_CC),1) log(SDF_CC)];
for j=1:7
    beta_s(j,:) = regress(portfolio_s(:,j),Y1)';
    beta_sdf(j,:)   = regress(portfolio_s(:,j),Y2)';
end

figure(8)
subplot(2,2,3)
scatter(beta_s(:,2),port_mean_s,'filled');title('SDF - All Sample');
ylabel('Mean Excess Return');xlabel('\beta_{SDF}');lsline;
subplot(2,2,2)
scatter(beta_sdf(:,2),port_mean_s,'filled');title('log(SDF) - All Sample');
ylabel('Mean Excess Return');xlabel('\beta_{log(SDF)}');lsline;refline(-1,0);
hold on
%Post WW2
portfolio_s = [ r_r(22:end) s1_r(22:end) s2_r(22:end)...
                s3_r(22:end) b1_r(22:end) b2_r(22:end) b3_r(22:end)];
port_mean_s = [mean(r_r(22:end)-rf_r(22:end))...
            mean(s1_r(22:end)-rf_r(22:end)),mean(s2_r(22:end)-rf_r(22:end)),...
            mean(s3_r(22:end)-rf_r(22:end)),mean(b1_r(22:end)-rf_r(22:end)),...
            mean(b2_r(22:end)-rf_r(22:end)),mean(b3_r(22:end)-rf_r(22:end))]';
Y1 = [ones(length(SDF_CC(22:end)),1) SDF_CC(22:end)];
Y2 = [ones(length(SDF_CC(22:end)),1) log(SDF_CC(22:end))];
for j=1:7
    beta_s(j,:) = regress(portfolio_s(:,j),Y1)';
    beta_sdf(j,:)   = regress(portfolio_s(:,j),Y2)';
end

subplot(2,2,4)
scatter(beta_s(:,2),port_mean_s,'filled');title('SDF - PostWW2');
ylabel('Mean Excess Return');xlabel('\beta_{SDF}');lsline;
subplot(2,2,1)
scatter(beta_sdf(:,2),port_mean_s,'filled');title('log(SDF)  - PostWW2');
ylabel('Mean Excess Return');xlabel('\beta_{log(SDF)}');lsline;refline(-1,0);
saveas(figure(8),'figure8.png');

% Part c
post=22;
er1 = log((1+r_r(post+1:end)));
er2 = er1(2:end,1)+er1(1:end-1,1);
er3 = er1(3:end,1)+er2(1:end-1,1);
er4 = er1(4:end,1)+er3(1:end-1,1); 
er5 = er1(5:end,1)+er4(1:end-1,1); 
er6 = er1(6:end,1)+er5(1:end-1,1); 
er7 = er1(7:end,1)+er6(1:end-1,1); 
for k=1:7
    X  = [ones(length(S(22:end-k-1,1)),1) log(S(22:end-k-1,1))];
    y  = eval(['er',num2str(k)]);
    [b,bint,r,rint,stats]=regress(y,X);
    eval(['b_',num2str(k),'=b(2,1)',';']);
    eval(['bint_',num2str(k),'=bint',';']);
    eval(['R2_',num2str(k),'=stats(1,1)',';']);
end
betas=table(b_1,b_2,b_3,b_4,b_5,b_6,b_7)
R2=table(R2_1,R2_2,R2_3,R2_4,R2_5,R2_6,R2_7)



