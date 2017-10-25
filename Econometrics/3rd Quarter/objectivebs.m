%% Set the objective function for bootstrap (re-centered)

% This program computes the objective function for Euler Equation conditions

%%INPUTS%% 
% x = (delta,gamma) parameters to be estimated
% W = weighting matrix (KxK) 
% z = observables at t 
% T=obs number 
% c_1= log(C_t+1/C_t)
% c= log(C_t/C_t-1), r_1= log(R_t+1)
% g_center=it is the average of g computed at two-step GMM 

%%OUTPUTS%%
% Sdf = stochastic discount factor (1xT)
% Euler = set of euler conditions equal to zero (K*T)
% g = sample average of each Euler (Kx1)
% Qn = objective function to minimize (1x1)
% Euler_bs = Euler recentered
% default output is Qn only 

function Qn = objectivebs(x,g_center, W, z, T, c_1, r_1);
euler=(x(1)*exp(-x(2)*c_1+r_1)-1).*z;
euler_bs=euler-g_center;
g=[sum(euler_bs')/T]';
Qn = 0.5*(g'*W*g);
end