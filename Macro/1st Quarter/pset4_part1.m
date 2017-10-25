%%% PROBLEM SET 4 - Gualtiero Azzalini %%%

%% PROBLEM 1 %%
clear all
%   Set Parameters

    lambda = 0.5;
    delta  = 0.05;
    beta   = 0.80;  

%   Question 1.6

%   Solving for the roots

a  =   (1-beta)^(2)*lambda;
b  =   2/beta*(beta^(2)*(lambda+delta-1)-beta*(lambda+delta-2)-1);
c  =   lambda;

v01    =    (-b-sqrt(b^(2)-4*a*c))/(2*a);
v02    =    (-b+sqrt(b^(2)-4*a*c))/(2*a);

wbar1  =    (1-beta*(1-delta))*v01-delta*beta*v01;
wbar2  =    (1-beta*(1-delta))*v02-delta*beta*v02;

%   Choose v01 because the other gives wbar higher than 1
%   Compute Vn+1 as function of Vn

coef    =   beta/(2*(1-beta*(1-delta)));
btilda  =   2*(lambda*beta*delta+(1-lambda)*(1-beta*(1-delta)));
I       =   1000;
v0      =   zeros(I,1);
v1      =   zeros(I,1);
v0      =   linspace(0,10,I);
v1      =   coef*(a*v0.^(2)+v0*btilda+c); 
v_star  =   coef*(a*v01^(2)+v01*btilda+c);

%   Figure 1 (Vn+1 as function of Vn)
figure1  =   figure; plot(v0,v1); hold on; plot(v0,v0,'-.'); hold on; plot(v01,v_star,'ko');
legend({'V_{n+1}(0)','45° line', 'v^{*}'},'FontSize',11);legend('Location','northwest');
xlabel('V_{n}(0)'); ylabel('V_{n+1}(0)'); title('V_{n+1}(0) as function of V_{n}(0)','FontSize',14);
saveas(figure1,'Figure1.png');

%   Determine v0_high and v0_low of employed

v0_low  =    0/(1-beta*(1-delta));
v0_high =    1/(1-beta*(1-delta));
v_high  =    v0_high;
v_low   =    v0_low;

%   Convergence v0_low and v0_high

tol=0.000001;its=1;maxits=500;dif_high=100;dif_low=100; %Set additional parameters for the loop

while dif_low>tol & dif_high>tol & its<maxits    
v1_low      =   coef*(a*v0_low^(2)+v0_low*btilda+c);
v1_high     =   coef*(a*v0_high^(2)+v0_high*btilda+c);
dif_low     =   abs(v1_low-v0_low);
dif_high    =   abs(v1_high-v0_high);
v0_low      =   v1_low;
v0_high     =   v1_high;
its         =   its+1;
end

%%  PROBLEM 3
clear
%   Set Parameters
delta=0.05;beta=0.98;alfa=1/3;gamma=2;A=1;its=0;maxits=500;tol=0.01;
dif=100;pointsk=500;

%   Steady State
kss    =    ((1-beta*(1-delta))/(alfa*beta*A))^(1/(alfa-1));

%   Generate grid for capital

kmin   =   1/2*kss;
kmax   =   3/2*kss;
k      =   linspace(kmin,kmax,pointsk)';
I      =   length(k);
%   Generate grid for productivity shock 

rho=0.95;sigma=0.02;N=21;T=6500;sigmaz=1;

for j=1:N
    
    p=(1+rho)/2;
    q=p;
    psi=(N-1)^(1/2)*sigmaz;
    s0=floor((N+1)/2);
    
    [z P]=arprocess(p,q,psi,N);z=z';
    %[X_approx,s]=markov(P,T+1,s0,Y);X_approx=X_approx';
end

%   Generate the F matrix
l=zeros(I,N);
for h=1:I
for j=1:N
    for i=1:I
    l(i,j)=k(h,1);
    end 
end
F=zeros(I,N);
for j=1:N
    for i=1:I
    if l(i,j)<=((1-delta)*k(h,1)+A*exp(z(j,1))*k(h,1)^(alfa)) & l(i,j)>=0
        F(i,j)=(((1-delta)*k(h,1)+A*exp(z(j,1))*k(h,1)^(alfa)-l(i,j))^(1-gamma)-1)/(1-gamma);
    else
        F(i,j)=log(0);
    end
    end
end
v0=ones(I,N);v1=zeros(I,N);
v_exp=zeros(I,N);
v_exp=v0*P;
value=zeros(I,N);

%value function iteration
while dif>tol & its<maxits
value=F+beta*v_exp;
    [M K]=max(value');
    v1=M;
    ind=K;
dif=max(abs(v1-v0));
v0=v1;
its=its+1; 
end
v_final(:,h)=M;
end








