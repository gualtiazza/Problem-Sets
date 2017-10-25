%% Econometrics II - Part II %%
%  Problem Set 1
%  Gualtiero Azzalini

addpath('/Users/Gualtiero/Dropbox/A NYU/Metrics/cogley/psets/1')
addpath('/Users/Gualtiero/Documents/MATLAB/Add-Ons/mfe_toolbox/utility')
addpath('/Users/Gualtiero/Documents/MATLAB/Add-Ons/mfe_toolbox/timeseries')

clear all;

%% Exercise 4

rho1    = [1.25; 1.7];
rho2    = [-0.3; -0.8];
dt      = 1000;
w       = linspace(0,pi,dt)';
i       = sqrt(-1);

for j=1:2
        for p = 1:length(w)
            fxx(p,1)   = 1/(2*pi*(1+rho1(j,1)^(2)+rho2(j,1)^(2)+2*(rho1(j,1)*rho2(j,1)-rho1(j,1))*cos(w(p,1))-2*rho2(j,1)*cos(2*w(p,1))));
        end
    gamma = ifft(fxx);
    acf = gamma/gamma(1,1);   
    eval(['gamma_',num2str(j),'=gamma']);
    eval(['fxx_',num2str(j),'=fxx']);
    eval(['acf_',num2str(j),'=acf']);
end

figure(1)
subplot(2,2,1); plot(w/(2*pi),fxx_1);title('\rho_{1}=1.25,\rho_{2}=-0.3');ylabel('Spectrum');xlabel('Cycles per quarter');
subplot(2,2,3); plot(0:1:99,acf_1(1:100));ylabel('Autocorrelation Function');xlabel('Lag');
subplot(2,2,2); plot(w/(2*pi),fxx_2);title('\rho_{1}=1.7,\rho_{2}=-0.8');ylabel('Spectrum');xlabel('Cycles per quarter');
subplot(2,2,4); plot(0:1:99,acf_2(1:100));ylabel('Autocorrelation Function');xlabel('Lag');
saveas(figure(1),'exercise4.png')

%% Exercise 5
clear
rho   = 0.9;
sigma = 1;
T     = 100;
N     = 1000;
nlag  = 4;

% True Variance
[autocorr, gamma0] = acf(rho,0,99);
acv                = autocorr.*gamma0;
fxx                = abs(fft(acv/(pi)));
TrueVar            = 2*pi*fxx(1)/T;



rng(12345);
epsilon = randn(T,N);
for n=1:N
    x(1,n) = epsilon(1,n);
    for t=2:T
        x(t,n) = rho*x(t-1,n)+epsilon(t,n);
    end
end

xhat = mean(x);

% Newey-West
for n=1:N
    xhat_m        = x(:,n)-xhat(:,n);
    nw    = (xhat_m'*xhat_m)/T;
    for pp=1:nlag
        w     = (nlag+1-pp)./(nlag+1);
        gamma = (T^(-1))*((xhat_m(1+pp:T))'*(xhat_m(1:T-pp)));
        gamma_p  = gamma+gamma';
        nw = nw+w*gamma_p;
    end
    AsyNW(n,1)=nw/sqrt(T);

% Smoothed Periodogram

    width     = 5;
    dZ        = fft(x(:,n),2*T);
    I         = inv(pi*T)*abs(dZ(1:T,:)).^2; 
    I_smooth  = conv2(I,ones(width,1)/width,'same');  
    AsySP(:,n)= 2*pi*I_smooth(1)/T;

% Prewhitening & Recoloring

    width       = 51;
    r           = x(1:end-1,n) \ x(2:end,n);
    u           = x(2:end,n)- r*x(1:end-1,n);
    dZ1         = fft(u,2*T);
    I1          = inv(pi*T)*abs(dZ1(1:T,:)).^2;
    I1_smooth   = conv2(I1,ones(width,1)/width,'same');
    fxx_PW      = abs((1 - r*exp(-1i*(1:100)'*(pi/100))).^(-1).*I1_smooth ...
        .*(1 - r*exp(1i*(1:100)'*(pi/100))).^(-1)); 
    AsyPW(:,n)  = 2*pi*fxx_PW(1)/T;
end


figure(2)
subplot(3,1,1)
histogram(AsyNW,'BinWidth',0.2,'FaceColor','green');title('Newey-West');xlim([0 8]), ylim([0 300])
set(vline([TrueVar, mean(AsyNW)],{'r--','b-.'},{'True','Sample'}),'LineWidth',2)
subplot(3,1,2)
histogram(AsySP,'BinWidth',0.2,'FaceColor','green');title('Smoothed Periodogram');xlim([0 8]), ylim([0 300])
set(vline([TrueVar, mean(AsySP)],{'r--','b-.'},{'True','Sample'}),'LineWidth',2)
subplot(3,1,3)
histogram(AsyPW,'BinWidth',0.2,'FaceColor','green');title('Prewhitening & Recoloring');xlim([0 8]), ylim([0 300])
set(vline([TrueVar, mean(AsyPW)],{'r--','b-.'},{'True','Sample'}),'LineWidth',2)
saveas(figure(2),'exercise5.png')




    














