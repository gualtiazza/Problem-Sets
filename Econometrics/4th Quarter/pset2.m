%% Econometrics - Part II/Prof.Cogley %%
%  Gualtiero Azzalini
%  Problem Set 2

addpath('/Users/Gualtiero/Dropbox/A NYU/Metrics/cogley/psets/2')
addpath('/Users/Gualtiero/Documents/MATLAB/Add-Ons/mfe_toolbox/utility')
addpath('/Users/Gualtiero/Documents/MATLAB/Add-Ons/mfe_toolbox/timeseries')

clear;

%% Exercise 4

load ps2.mat;

% Estimate 2nd order VAR - equivalent to OLS eq by eq in this case %

beta_dlny = regress(dlny(3:end),[dlny(2:end-1) dlny(1:end-2) ...
            Zic(2:end-1) Zic(1:end-2)]); 
beta_Zic  = regress(Zic(3:end),[dlny(2:end-1) dlny(1:end-2) ...
            Zic(2:end-1) Zic(1:end-2)]);  
         
% Collect the estimated companion matrices

F1        = [beta_dlny(1,1) beta_dlny(3,1);beta_Zic(1,1) beta_Zic(3,1)];
F2        = [beta_dlny(2,1) beta_dlny(4,1);beta_Zic(2,1) beta_Zic(4,1)];
        
% Beveridge-Nelson output gap %

betaBN_dlny = regress(dlny(3:end),[dlny(2:end-1)-mean(dlny) dlny(1:end-2)-mean(dlny) ...
            Zic(2:end-1) Zic(1:end-2)]); 
betaBN_Zic  = regress(Zic(3:end),[dlny(2:end-1)-mean(dlny) dlny(1:end-2)-mean(dlny) ...
            Zic(2:end-1) Zic(1:end-2)]);  
         
% Collect the estimated companion matrices

F1_BN     = [betaBN_dlny(1,1) betaBN_dlny(3,1);betaBN_Zic(1,1) betaBN_Zic(3,1)];
F2_BN     = [betaBN_dlny(2,1) betaBN_dlny(4,1);betaBN_Zic(2,1) betaBN_Zic(4,1)];

A         = [F1_BN F2_BN;eye(2,2) zeros(2,2)];
e         = [1 0 0 0];

Zt  = [(dlny(3:end)-mean(dlny))';(Zic(3:end))';...
      (dlny(2:end-1)-mean(dlny))';(Zic(2:end-1))'];

for i=1:size(Zt,2)
    
    cycle(i,1) = - e*inv(eye(4,4)-A)*A*Zt(:,i);
    
end
date=date';

figure(1)
plot(date(3:end),cycle,'k');
axis([1950 2016 -0.08 0.08])
title('Beveridge-Nelson Cyclical Component');
legend('B-N Cycle');
saveas(figure(1),'BNgap.png');

% Asymptotics SE report range of end-of-sample in +-1 std

B       = e*inv(eye(4,4)-A)*A;
B_prime = gradient(B);
Zt_1    = [(dlny(2:end-1)-mean(dlny))';(Zic(2:end-1))';...
          (dlny(1:end-2)-mean(dlny))';(Zic(1:end-2))'];
eps     = Zt-A*Zt_1;
V       = length(eps)^(-1)*eps*eps';
AVar    = B_prime*V*B_prime';
Asd     = AVar.^(1/2);

for i=1:size(Zt,2)
    
    cycle_l(i,1) = -e*inv(eye(4,4)-A)*A*Zt(:,i)-1.65*Asd;
    cycle_h(i,1) = -e*inv(eye(4,4)-A)*A*Zt(:,i)+1.65*Asd;
    
end

cycleCI = [cycle_l cycle_h];

figure(4)
plot(date(212:end),cycle(210:end),'k');
hold on
plot(date(212:end),cycleCI(210:end,:),':k');
axis([2000 2016 -0.08 0.08])
title('Beveridge-Nelson Cyclical Component');
legend('B-N Cycle','\pm1 St.Dev CI');
saveas(figure(4),'BNgapCI.png');
