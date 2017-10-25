%% part 1.2-vectorization 
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
l=zeros(N,N);   %this is k'
v0=zeros(N,N);
v1=zeros(N,N);
v=zeros(N,N);
F=zeros(N,N);
its=0;
maxits=500;
tol=0.01;
dif=100;
pol=zeros(N,1);    %this stores the indexes for the policy function

for j=1:N
    for i=1:N
    l(i,j)=k(i,1);
    end 
end
for j=1:N
    for i=1:N
    if l(i,j)<=((1-delta)*k(j,1)+A*k(j,1)^(alfa)) & l(i,j)>=0
        F(i,j)=(((1-delta)*k(j,1)+A*k(j,1)^(alfa)-l(i,j))^(1-gamma)-1)/(1-gamma);
    else
        F(i,j)=log(0);
    end
    end
end
while dif>tol & its<maxits
for j=1:N
    for i=1:N
    v0(i,j)=v0(i,1);
    end
end
    [M,I]=max(F+beta*v0);
    v1=M';
    ind=I';
for j=1:N    
    pol(j,1)=l(ind(j,1),j);
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
%plots
vectorization=figure;
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
saveas(vectorization,'vectorization098&105.png');