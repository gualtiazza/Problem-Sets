%% uncertainty
clear all
%set parameters
tic
delta=0.05;beta=0.98;alfa=1/3;gamma=2;A=1;its=0;maxits=500;tol=0.01;
dif=100;pointsk=1000;pointsc=1000;epsilon=0.4;
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
%construct G_minus
% G_minus=zeros(J,I);
% for i=1:I
%     for j=1:J    
%         G_minus(j,i)=((1-delta)*k(i,1)+A*k(i,1)^(alfa)-c(j,1))*(1-epsilon);
%     end
% end
G_plus=zeros(J,I);
for i=1:I
    for j=1:J    
        G_plus(j,i)=((1-delta)*k(i,1)+A*k(i,1)^(alfa)-c(j,1));
    end
end
%construct F
F=zeros(J,I);
for i=1:I
    for j=1:J  
        %if G_plus(j,i)>=0 | G_plus(j,i)<=(1-delta)*k(i,1)+A*k(i,1)^(alfa) 
        F(j,i)=((c(j,1))^(1-gamma)-1)/(1-gamma);
        %else
        %F(j,i)=log(0);
        %end
    end
end
value=zeros(J,I);
%value function iteration
while dif>tol & its<maxits
value=F+beta*(1/2*interp1(k,v0,(1-epsilon)*G_plus,'linear','extrap')+1/2*interp1(k,v0,(1+epsilon)*G_plus,'linear','extrap'));
    [M K]=max(value);
    v1=M';
    ind=K';
dif=max(abs(v1-v0));
v0=v1;
its=its+1; 
end
for i=1:I
pol(i,1)=(G_plus(ind(i,1),i));
c(i,1)=(1-delta)*k(i,1)+A*k(i,1)^(alfa)-pol(i,1);
end
pol_base=pol;
c_base=c;
v1_base=v1;

% %% epsilon=0.1
% delta=0.05;beta=0.98;alfa=1/3;gamma=2;A=1;its=0;maxits=500;tol=0.01;
% dif=100;pointsk=500;pointsc=500;epsilon=0.1;
% %compute the steady state
% kss=((1-beta*(1-delta))/(alfa*beta*A))^(1/(alfa-1));
% css=(1-delta)*kss+A*kss^(alfa)-kss;iss=A*kss^(alfa)-css;
% %create the grid for capital
% kmin=1/2*kss;kmax=3/2*kss;
% k=linspace(kmin,kmax,pointsk)';
% %create the grid for consumption
% cmin=1;cmax=5;c=linspace(cmin,cmax,pointsc)';
% %other settings
% I=length(k);J=length(c);v0=linspace(1,100,I);v1=zeros(I,1);
% ind=zeros(I,1);pol=zeros(I,1);
% %construct G_minus
% G_minus=zeros(J,I);
% for i=1:I
%     for j=1:J    
%         G_minus(j,i)=((1-delta)*k(i,1)+A*k(i,1)^(alfa)-c(j,1))*(1-epsilon);
%     end
% end
% G_plus=zeros(J,I);
% for i=1:I
%     for j=1:J    
%         G_plus(j,i)=((1-delta)*k(i,1)+A*k(i,1)^(alfa)-c(j,1))*(1+epsilon);
%     end
% end
% %construct F
% F=zeros(J,I);
% for i=1:I
%     for j=1:J  
%         if (G_minus(j,i)+G_plus(j,i))/2>=0 | (G_minus(j,i)+G_plus(j,i))/2<=(1-delta)*k(i,1)+A*k(i,1)^(alfa) 
%         F(j,i)=((c(j,1))^(1-gamma)-1)/(1-gamma);
%         else
%         F(j,i)=log(0);
%         end
%     end
% end
% value=zeros(J,I);
% %value function iteration
% while dif>tol & its<maxits
% value=F+beta*(1/2*interp1(k,v0,G_plus,'linear','extrap')+1/2*interp1(k,v0,G_minus,'linear','extrap'));
%     [M K]=max(value);
%     v1=M';
%     ind=K';
% dif=max(abs(v1-v0));
% v0=v1;
% its=its+1; 
% end
% for i=1:I
% pol(i,1)=(G_plus(ind(i,1),i)+G_minus(ind(i,1),i))/2;
% c(i,1)=(1-delta)*k(i,1)+A*k(i,1)^(alfa)-pol(i,1);
% end
% pol_01=pol;
% c_01=c;
% v1_01=v1;
% 
% %% epsilon=0.25
% delta=0.05;beta=0.98;alfa=1/3;gamma=2;A=1;its=0;maxits=500;tol=0.01;
% dif=100;pointsk=500;pointsc=500;epsilon=0.25;
% %compute the steady state
% kss=((1-beta*(1-delta))/(alfa*beta*A))^(1/(alfa-1));
% css=(1-delta)*kss+A*kss^(alfa)-kss;iss=A*kss^(alfa)-css;
% %create the grid for capital
% kmin=1/2*kss;kmax=3/2*kss;
% k=linspace(kmin,kmax,pointsk)';
% %create the grid for consumption
% cmin=1;cmax=5;c=linspace(cmin,cmax,pointsc)';
% %other settings
% I=length(k);J=length(c);v0=linspace(1,100,I);v1=zeros(I,1);
% ind=zeros(I,1);pol=zeros(I,1);
% %construct G_minus
% G_minus=zeros(J,I);
% for i=1:I
%     for j=1:J    
%         G_minus(j,i)=((1-delta)*k(i,1)+A*k(i,1)^(alfa)-c(j,1))*(1-epsilon);
%     end
% end
% G_plus=zeros(J,I);
% for i=1:I
%     for j=1:J    
%         G_plus(j,i)=((1-delta)*k(i,1)+A*k(i,1)^(alfa)-c(j,1))*(1+epsilon);
%     end
% end
% %construct F
% F=zeros(J,I);
% for i=1:I
%     for j=1:J  
%         if (G_minus(j,i)+G_plus(j,i))/2>=0 | (G_minus(j,i)+G_plus(j,i))/2<=(1-delta)*k(i,1)+A*k(i,1)^(alfa) 
%         F(j,i)=((c(j,1))^(1-gamma)-1)/(1-gamma);
%         else
%         F(j,i)=log(0);
%         end
%     end
% end
% value=zeros(J,I);
% %value function iteration
% while dif>tol & its<maxits
% value=F+beta*(1/2*interp1(k,v0,G_plus,'linear','extrap')+1/2*interp1(k,v0,G_minus,'linear','extrap'));
%     [M K]=max(value);
%     v1=M';
%     ind=K';
% dif=max(abs(v1-v0));
% v0=v1;
% its=its+1; 
% end
% for i=1:I
% pol(i,1)=(G_plus(ind(i,1),i)+G_minus(ind(i,1),i))/2;
% c(i,1)=(1-delta)*k(i,1)+A*k(i,1)^(alfa)-pol(i,1);
% end
% pol_25=pol;
% c_25=c;
% v1_25=v1;
% 
% %% epsilon=0.4
% delta=0.05;beta=0.98;alfa=1/3;gamma=2;A=1;its=0;maxits=500;tol=0.01;
% dif=100;pointsk=500;pointsc=500;epsilon=0.4;
% %compute the steady state
% kss=((1-beta*(1-delta))/(alfa*beta*A))^(1/(alfa-1));
% css=(1-delta)*kss+A*kss^(alfa)-kss;iss=A*kss^(alfa)-css;
% %create the grid for capital
% kmin=1/2*kss;kmax=3/2*kss;
% k=linspace(kmin,kmax,pointsk)';
% %create the grid for consumption
% cmin=1;cmax=5;c=linspace(cmin,cmax,pointsc)';
% %other settings
% I=length(k);J=length(c);v0=linspace(1,100,I);v1=zeros(I,1);
% ind=zeros(I,1);pol=zeros(I,1);
% %construct G_minus
% G_minus=zeros(J,I);
% for i=1:I
%     for j=1:J    
%         G_minus(j,i)=((1-delta)*k(i,1)+A*k(i,1)^(alfa)-c(j,1))*(1-epsilon);
%     end
% end
% G_plus=zeros(J,I);
% for i=1:I
%     for j=1:J    
%         G_plus(j,i)=((1-delta)*k(i,1)+A*k(i,1)^(alfa)-c(j,1))*(1+epsilon);
%     end
% end
% %construct F
% F=zeros(J,I);
% for i=1:I
%     for j=1:J  
%         if (G_minus(j,i)+G_plus(j,i))/2>=0 | (G_minus(j,i)+G_plus(j,i))/2<=(1-delta)*k(i,1)+A*k(i,1)^(alfa) 
%         F(j,i)=((c(j,1))^(1-gamma)-1)/(1-gamma);
%         else
%         F(j,i)=log(0);
%         end
%     end
% end
% value=zeros(J,I);
% %value function iteration
% while dif>tol & its<maxits
% value=F+beta*(1/2*interp1(k,v0,G_plus,'linear','extrap')+1/2*interp1(k,v0,G_minus,'linear','extrap'));
%     [M K]=max(value);
%     v1=M';
%     ind=K';
% dif=max(abs(v1-v0));
% v0=v1;
% its=its+1; 
% end
% for i=1:I
% pol(i,1)=(G_plus(ind(i,1),i)+G_minus(ind(i,1),i))/2;
% c(i,1)=(1-delta)*k(i,1)+A*k(i,1)^(alfa)-pol(i,1);
% end
% pol_04=pol;
% c_04=c;
% v1_04=v1;
% 
% 
% toc
%plots
uncert=figure;
subplot(2,2,1);
plot(k,v1_base,'b');
% hold on
% plot(k,v1_01,'b-.');
% hold on
% plot(k,v1_25,'b--');
% hold on
% plot(k,v1_04,'b:');
xlim([kmin kmax]);
ylabel('V(k)');
xlabel('k');
%legend('0.01','0.1','0.25','0.4')
legend('Location','northwest')
title('Value Function');
subplot(2,2,2);
plot(k,pol_base,'r');
% hold on;
% plot(k,pol_01,'r-.');
% hold on;
% plot(k,pol_25,'r--');
% hold on;
% plot(k,pol_04,'r:');
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
plot(k,c_base,'m');
% hold on
% plot(k,c_01,'m-.');
% hold on
% plot(k,c_25,'m--');
% hold on
% plot(k,c_04,'m:');
% hold on
ylim([0 3]);
xlim([kmin kmax]);
ylabel('c(k)');
xlabel('k');
title('Consumption');
legend('0.01','0.1','0.25','0.4')
saveas(uncert,'uncertainty.png');