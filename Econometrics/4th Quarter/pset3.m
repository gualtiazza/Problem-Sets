%% Econometrics - Part II/Prof.Cogley %%
%  Gualtiero Azzalini
%  Problem Set 3

addpath('/Users/Gualtiero/Dropbox/A NYU/Metrics/cogley/psets/3')
addpath('/Users/Gualtiero/Documents/MATLAB/Add-Ons/mfe_toolbox/utility')
addpath('/Users/Gualtiero/Documents/MATLAB/Add-Ons/mfe_toolbox/timeseries')

clear;

%% Preliminaries
sigma   = 1; 
kappa   = .2;
beta    = .99;
phi_pi  = 1.25;
phi_x   = .25;
rho_i   = .5;
sigma_i = 2;
rho_g   = .75;
rho_u   = .75;
sigma_g = 1;
sigma_u = 1;
sigma_m = .1;
N       = 1000;

%% Write the model as a system of expectational equations
A   = [1                    0                   sigma       -1      0;
       -kappa               1                   0           0       -1;
       -phi_x*(1-rho_i)     -phi_pi*(1-rho_i)   1           0       0;
       0                    0                   0           1       0;
       0                    0                   0           0       1];
   
B   = [1    sigma   0   0   0;
       0    beta    0   0   0;
       0    0       0   0   0;
       0    0       0   0   0;
       0    0       0   0   0];
   
C   = [0    0   0       0       0;
       0    0   0       0       0;
       0    0   rho_i   0       0;
       0    0   0       rho_g   0;
       0    0   0       0       rho_u];  
   
D   = [0    0   0         0         0;
       0    0   0         0         0;
       0    0   sigma_i   0         0;
       0    0   0         sigma_g   0;
       0    0   0         0         sigma_u]; 
   
%% Solve the system and express the solution in state-space form

K = [zeros(5,5)     eye(5,5);
     -C             A       ];
 
L = [eye(5,5)       zeros(5,5);
    zeros(5,5)      B        ];

% QZ decomposition

[S,T,Q,Z] = qz(K,L);
[S,T,Q,Z] = ordqz(S,T,Q,Z,'udi');  

Z11 = Z(1:5,1:5); Z21 = Z(6:10,1:5);

% Matrices of the state equation

F = real(Z21*inv(Z11));
G = inv(A-B*F)*D;

% Martices of measurement equation

H   = [1 0 0 0 0;
       0 1 0 0 0;
       0 0 1 0 0];
   
J   = [sigma_m 0 0;
      0 sigma_m 0;
      0 0 sigma_m]; 

%% Use state-space form to simulate time series

% Generate the random components in the state and measurement equations
eta = [zeros(2,N);randn(3,N)];
v   = randn(3,N);

% Generate the series for x, pi, i

Xt_1 = zeros(5,1);

for i=1:N
    Xt(:,i) = F*Xt_1 + G*eta(:,i);
    Xt_1 = Xt(:,i);
    Yt(:,i) = H*Xt(:,i) + J*v(:,i);
end

figure(1)
subplot(1,3,1);plot(Yt(1,:));title('x_{t}^{obs}');pbaspect([1 1 1]);
subplot(1,3,2);plot(Yt(2,:));title('\pi_{t}^{obs}');pbaspect([1 1 1]);
subplot(1,3,3);plot(Yt(3,:));title('i_{t}^{obs}');pbaspect([1 1 1]);
saveas(figure(1),'Yt_sim.png');

%% Given the data evaluate the log-likelihood 
loglike = kalman(F,G,H,J,Yt,size(Xt,1));

% %% Verify that true parameters are near the peak
figure(2)
step=20;
e=.2;
sigmak = linspace(sigma-e,sigma+e,step);
for i=1:length(sigmak)
loglike_i(i)=pset3function(sigmak(i),.2,.99,1.25,.25,.5,2,.75,.75,1,1,.1,N,eta,v,Xt,Yt);
end
subplot(4,3,1)
plot(sigmak,loglike_i);title('\sigma');vline(sigma);
xlim([sigma-e sigma+e]);
hold on

kappak = linspace(kappa-e,kappa+e,step);
for i=1:length(kappak)
loglike_i(i)=pset3function(1,kappak(i),.99,1.25,.25,.5,2,.75,.75,1,1,.1,N,eta,v,Xt,Yt);
end
subplot(4,3,2)
plot(kappak,loglike_i);title('\kappa');vline(kappa);
xlim([kappa-e kappa+e]);
hold on

betak = linspace(beta-e,beta+e,step);
for i=1:length(betak)
loglike_i(i)=pset3function(1,.2,betak(i),1.25,.25,.5,2,.75,.75,1,1,.1,N,eta,v,Xt,Yt);
end
subplot(4,3,3)
plot(betak,loglike_i);title('\beta');vline(beta);
xlim([beta-e beta+e]);
hold on

phi_pik = linspace(phi_pi-e,phi_pi+e,step);
for i=1:length(phi_pik)
loglike_i(i)=pset3function(1,.2,.99,phi_pik(i),.25,.5,2,.75,.75,1,1,.1,N,eta,v,Xt,Yt);
end
subplot(4,3,4)
plot(phi_pik,loglike_i);title('\phi_{\pi}');vline(phi_pi);
xlim([phi_pi-e phi_pi+e]);
hold on

phi_xk = linspace(phi_x-e,phi_x+e,step);
for i=1:length(phi_xk)
loglike_i(i)=pset3function(1,.2,.99,1.25,phi_xk(i),.5,2,.75,.75,1,1,.1,N,eta,v,Xt,Yt);
end
subplot(4,3,5)
plot(phi_xk,loglike_i);title('\phi_{x}');vline(phi_x);
xlim([phi_x-e phi_x+e]);
hold on

rho_ik = linspace(rho_i-e,rho_i+e,step);
for i=1:length(rho_ik)
loglike_i(i)=pset3function(1,.2,.99,1.25,0.25,rho_ik(i),2,.75,.75,1,1,.1,N,eta,v,Xt,Yt);
end
subplot(4,3,6)
plot(rho_ik,loglike_i);title('\rho_{i}');vline(rho_i);
xlim([rho_i-e rho_i+e]);
hold on

sigma_ik = linspace(sigma_i-e,sigma_i+e,step);
for i=1:length(sigma_ik)
loglike_i(i)=pset3function(1,.2,.99,1.25,0.25,0.5,sigma_ik(i),.75,.75,1,1,.1,N,eta,v,Xt,Yt);
end
subplot(4,3,7)
plot(sigma_ik,loglike_i);title('\sigma_{i}');vline(sigma_i);
xlim([sigma_i-e sigma_i+e]);
hold on

rho_gk = linspace(rho_g-e,rho_g+e,step);
for i=1:length(rho_gk)
loglike_i(i)=pset3function(1,.2,.99,1.25,0.25,0.5,2,rho_gk(i),.75,1,1,.1,N,eta,v,Xt,Yt);
end
subplot(4,3,8)
plot(rho_gk,loglike_i);title('\rho_{g}');vline(rho_g);
xlim([rho_g-e rho_g+e]);
hold on

sigma_gk = linspace(sigma_g-e,sigma_g+e,step);
for i=1:length(sigma_gk)
loglike_i(i)=pset3function(1,.2,.99,1.25,0.25,0.5,2,.75,0.75,sigma_gk(i),1,.1,N,eta,v,Xt,Yt);
end
subplot(4,3,9)
plot(sigma_gk,loglike_i);title('\sigma_{g}');vline(sigma_g);
xlim([sigma_g-e sigma_g+e]);
hold on

rho_uk = linspace(rho_u-e,rho_u+e,step);
for i=1:length(rho_uk)
loglike_i(i)=pset3function(1,.2,.99,1.25,0.25,0.5,2,0.75,rho_uk(i),1,1,.1,N,eta,v,Xt,Yt);
end
subplot(4,3,10)
plot(rho_uk,loglike_i);title('\rho_{u}');vline(rho_u);
xlim([rho_u-e rho_u+e]);
hold on

sigma_uk = linspace(sigma_u-e,sigma_u+e,step);
for i=1:length(sigma_uk)
loglike_i(i)=pset3function(1,.2,.99,1.25,0.25,0.5,2,.75,.75,1,sigma_uk(i),.1,N,eta,v,Xt,Yt);
end
subplot(4,3,11)
plot(sigma_uk,loglike_i);title('\sigma_{u}');vline(sigma_u);
xlim([sigma_u-e sigma_u+e]);
saveas(figure(2),'Log-likelihood as function of each parameter.png')





