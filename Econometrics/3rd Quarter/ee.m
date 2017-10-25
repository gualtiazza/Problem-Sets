%% Euler conditions 

% This program computes Euler conditions given parameters

%%INPUTS%% 
% x = (delta,gamma) estimated parameters  
% z = observables at t 
% T=obs number 
% c_1= log(C_t+1/C_t)
% c= log(C_t/C_t-1), r_1= log(R_t+1)

%%OUTPUTS%%
% sdf = stochastic discount factor (1xT)
% ee = set of euler conditions equal to zero (K*T) - default only ee

function ee = ee(x, z, T, c_1, r_1, r, c);
for i=1:T
sdf(i)=x(1)*exp(-x(2)*c_1(i)+r_1(i))-1;
ee(:,i)=kron(sdf(:,i),z(:,i));
end
end