%% Euler conditions Bootstrap

% This program computes Euler conditions given parameters

%%INPUTS%% 
% x = (delta,gamma) estimated parameters  
% z = observables at t 
% T=obs number 
% c_1= log(C_t+1/C_t)
% c= log(C_t/C_t-1), r_1= log(R_t+1)
% g_center=it is the average of g computed at two-step GMM 


%%OUTPUTS%%
% sdf = stochastic discount factor (1xT)
% euler = set of euler conditions equal to zero (K*T) - default only ee
% Euler_bs = Euler recentered

function euler_bs = eebs(x,g_center, z, T, c_1, r_1);
euler=(x(1)*exp(-x(2)*c_1+r_1)-1).*z;
euler_bs=euler-g_center;
end