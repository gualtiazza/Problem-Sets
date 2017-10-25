%% Problem Set 5 code - Gualtiero Azzalini %%

%% Exercise 1

addpath('/Users/Gualtiero/Dropbox/A NYU/Metrics/crump/psets/5')
addpath('/Users/Gualtiero/Documents/MATLAB/Add-Ons/mitbook/utilities_tom')

clear; clc; close all;

% Set Parameters 

beta              = 1; 
sigma_eu          = 0.5;
R_values          = [0.5 0.1 0.01];
N                 = 200;
S                 = 10000;
e                 = zeros(N,S);
u                 = zeros(N,S);
z                 = zeros(N,S);
qtiles            = [0.25 0.50 0.75];
qtiles_OLS        = zeros(4,3);
qtiles_2SLS       = zeros(4,3);
qtiles_OLS(1,:)   = qtiles;
qtiles_2SLS(1,:)  = qtiles;
h=1;
figure1=figure;

% Simulations
for j=1:3
    R             = R_values(1,j);
    gamma         = sqrt(R/(1-R));
    beta_OLS      = NaN(1,S);   
    beta_2SLS     = NaN(1,S);
    mu            = [0 0];
    sigma         = [1 sigma_eu;sigma_eu 1];
    for k=1:S
    r       = mvnrnd(mu,sigma,N);
    w       = randn(N,1);
    z(:,k)  = w(:,1);
    e(:,k)  = r(:,1);
    u(:,k)  = r(:,2);
    x(:,k)  = z(:,k)*gamma+u(:,k);
    y(:,k)  = beta*x(:,k)+e(:,k);
    
    X = x(:,k);
    Y = y(:,k);
    Z = z(:,k);
    beta_OLS(:,k)   = X \ Y; 
    beta_2SLS(:,k)  = inv(Z'*X)*Z'*Y;
    end    

    figure1;
    subplot(3,2,h)
    histogram(beta_OLS',50,'FaceAlpha',0.5);
    ylabel(R);
    title('$$ \hat{\beta}_{OLS} $$','interpreter','latex');
    subplot(3,2,h+1)
    histogram(beta_2SLS',50,'FaceAlpha',0.5,'FaceColor',[1 0.6 0.3]);
    title('$$ \hat{\beta}_{2SLS} $$','interpreter','latex');
    h=h+2;  
    
    % Quantiles
    qtiles_OLS(j+1,:)   = quantile(beta_OLS,qtiles);
    qtiles_2SLS(j+1,:)  = quantile(beta_2SLS,qtiles);
end
saveas(figure1,'figure1.png')
Rows = {'Quantile';'R^2=0.5';'R^2=0.1';'R^2=0.01'}; 
T = table(qtiles_OLS,qtiles_2SLS,'RowNames',Rows)
