%% interpolation
clear all
%set parameters
tic
delta=0.05;beta=0.98;alfa=1/3;gamma=2;A=1;its=0;maxits=500;tol=0.01;
dif=100;pointsk=700;pointsc=700;
%compute the steady state
kss=((1-beta*(1-delta))/(alfa*beta*A))^(1/(alfa-1));
css=(1-delta)*kss+A*kss^(alfa)-kss;iss=A*kss^(alfa)-css;
%create the grid for capital
kmin=1/2*kss;kmax=3/2*kss;
k=linspace(kmin,kmax,pointsk)';
%create the grid for consumption
cmin=1;cmax=5;c=linspace(cmin,cmax,pointsc)';
%other settings
I=length(k);J=length(c);v0=linspace(1,100,I);v1=zeros(I,1);
ind=zeros(I,1);pol=zeros(I,1);
%construct G
G=zeros(J,I);
for i=1:I
    for j=1:J    
        G(j,i)=(1-delta)*k(i,1)+A*k(i,1)^(alfa)-c(j,1);
    end
end
%construct F
F=zeros(J,I);
for i=1:I
    for j=1:J  
        if G(j,i)>=0 | G(j,i)<=(1-delta)*k(i,1)+A*k(i,1)^(alfa) 
        F(j,i)=((c(j,1))^(1-gamma)-1)/(1-gamma);
        else
        F(j,i)=log(0);
        end
    end
end
value=zeros(J,I);
%value function iteration
while dif>tol & its<maxits
value=F+beta*interp1(k,v0,G,'linear');
    [M K]=max(value);
    v1=M';
    ind=K';
dif=max(abs(v1-v0));
v0=v1;
its=its+1; 
end
for i=1:I
pol(i,1)=G(ind(i,1),i);
c(i,1)=(1-delta)*k(i,1)+A*k(i,1)^(alfa)-pol(i,1);
end
toc
%plots
interp=figure;
subplot(2,2,1);
plot(k,v1,'b');
xlim([kmin kmax]);
ylabel('V(k)');
xlabel('k');
title('Value Function');
subplot(2,2,2);
plot(k,pol,'r-');
hold on;
plot(k,k,'k:');
hold on
plot(kss,kss,'ok');
ylim([kmin kmax]);
xlim([kmin kmax]);
ylabel('g(k)');
xlabel('k');
title('Next Period Capital');
subplot(2,2,3);
plot(k,c,'m');
xlim([kmin kmax]);
ylabel('c(k)');
xlabel('k');
title('Consumption');
saveas(interp,'interpolation.png');