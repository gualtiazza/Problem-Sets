%Problem set 3
%% part 1.1 -elementary iteration%% 
clear all
%set parameters
tic
delta=0.05;
beta=0.98;
alfa=1/3;
gamma=2;
A=1;
points=2000;
%compute the steady state
kss=((1-beta*(1-delta))/(alfa*beta*A))^(1/(alfa-1));
css=(1-delta)*kss+A*kss^(alfa)-kss;
iss=A*kss^(alfa)-css;
%create the grid
kmin=0.98*kss;
kmax=1.05*kss;
k=linspace(kmin,kmax,points)';
%other settings
N=length(k);
c=zeros(N,1);
l=zeros(N,1);   %this is k'
v0=zeros(N,1);
v1=zeros(N,1);
v=zeros(N,1);
F=zeros(N,1);
its=0;
maxits=500;
tol=0.01;
dif=100;
pol=zeros(N,1);    %this stores the indexes for the policy function

for i=1:N
    l(i,1)=k(i,1);
end
while dif>tol & its<maxits
for j=1:N
for i=1:N
    k0=k(j,1);
    if l(i,1)<=((1-delta)*k0+A*k0^(alfa)) & l(i,1)>=0
        F(i,1)=(((1-delta)*k0+A*k0^(alfa)-l(i,1))^(1-gamma)-1)/(1-gamma);
    else
        F(i,1)=log(0);
    end
    v(i,1)=F(i,1)+beta*v0(i,1);
    [M,I]=max(v);
    v1(j,1)=M;
    pol(j,1)=l(I,1);
end
end
dif=max(abs(v1-v0));
v0=v1;
its=its+1;
end
for i=1:N-1
    c(i,1)=(1-delta)*k(i,1)+A*k(i,1)^(alfa)-pol(i,1);
    c(N,1)=NaN;
end
toc
time=toc;
%plots
elementary=figure;
subplot(2,2,1);
plot(k,v1,'b');
xlim([kmin kmax]);
ylabel('V(k)');
xlabel('k');
title('Value Function');
subplot(2,2,2);
plot(k,pol,'r-');
hold on;
plot(k,k,'k-.');
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
saveas(elementary,'elementary.png');





