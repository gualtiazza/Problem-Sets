%% Computation of long-run variance matrix

%%INPUTS%%
% data = matrix of dependent data (KxT)
% T = obs number 
% lags = number of lags to include 

%%OUTPUTS%%
% W_hat = estimated long-run variance (KxK)

function W_hat = longrunW(data, T, lags);
dataT=data';
W_hat=data*data'/T;
for j=1:lags   
gammaj=(data(:,1+j:T)*dataT(1:(T-j),:))/T;
gammaj_sq=gammaj+gammaj';
W_hat=W_hat+(1-j/(lags+1))*gammaj_sq;
end
end