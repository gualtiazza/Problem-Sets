%% Problem Set 4 code - Gualtiero Azzalini %%

%% Exercise 1

addpath('/Users/Gualtiero/Dropbox/A NYU/Metrics/crump/psets/4')
addpath('/Users/Gualtiero/Documents/MATLAB/Add-Ons/mitbook/utilities_tom')

clear; clc; close all;

% Set Parameters 

beta0   = 1;    beta1   = 1;    beta2   = 1;
rho_values     = [0 0.5 0.9];
T       = 200;
N       = 10000;
k       = 3;

% Set restrictions 

R_A    =    [0 0;1 0;0 1];
c_A    =    [1 1]';
         
R_B    =    [0 1 1]';
c_B    =    2;

h=1;
figure1=figure;
for j=1:3
    rho=rho_values(1,j);
% Set model
e       = randn(T,N);
x1      = sqrt(12)*rand(T,N);
x2      = sqrt(1-rho^2)*sqrt(12)*rand(T,N) + rho*x1;
y       = beta0 + beta1*x1 + beta2*x2 + e;
i       = ones(T,1);
a1      = [1 0 0]';
a2      = [0 0 1]';

% Containers
beta_OLS      = NaN(k,N);
beta_A        = NaN(k,N);
beta_B        = NaN(k,N);

for n = 1:N

    X = [i, x1(:,n), x2(:,n)];
    Y = y(:,n);

    beta_OLS(:,n)   = X \ Y; 
    beta_A(:,n)     = beta_OLS(:,n)-inv(X'*X)*R_A*inv(R_A'*inv(X'*X)*R_A)*(R_A'*beta_OLS(:,n)-c_A);
    beta_B(:,n)     = beta_OLS(:,n)-inv(X'*X)*R_B*inv(R_B'*inv(X'*X)*R_B)*(R_B'*beta_OLS(:,n)-c_B);
    

end
    figure1
    subplot(3,2,h)
    histogram(a1'*beta_OLS,50,'FaceAlpha',0.5);
    ylabel(rho);
    title('$$ \hat{\beta}_0 $$','interpreter','latex');
    hold on
    histogram(a1'*beta_B,50,'FaceAlpha',0.5);
    hold on
    histogram(a1'*beta_A,50,'FaceAlpha',0.5);
    
    subplot(3,2,h+1)
    histogram(a2'*beta_OLS,50,'FaceAlpha',0.5);
    title('$$ \hat{\beta}_2 $$','interpreter','latex');
    hold on
    histogram(a2'*beta_B,50,'FaceAlpha',0.5);
    hold on
    %histogram(a2'*beta_A,50,'FaceAlpha',0.5);
    %hold on
    h=h+2;
    
end
subplot(3,2,1)
h=legend('Unconstrained', '{\beta}_1+{\beta}_2=2','{\beta}_1={\beta}_2=1','interpreter','latex');
set(h,'FontSize',7)
legend('boxoff');
legend('Location','northwest');

saveas(figure1,'figure1.png')
