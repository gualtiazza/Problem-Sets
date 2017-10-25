%Problem Set 2
%set parameters
alfa=1/3;
beta=0.98;
gamma=2;
delta=0.05;
kss=10.3;

%compute partial derivatives in linearization

Fyx=gamma*(kss^alfa-delta*kss)^(-(gamma+1))*((1-delta)+alfa*kss^(alfa-1));
Fyy=-gamma*(kss^alfa-delta*kss)^(-(gamma+1));
Fxy=Fyx;
Fxx=(kss^alfa-delta*kss)^(-(gamma+1))*(alfa-1)*alfa*kss^(alfa-2)-gamma*(kss^alfa-delta*kss)^(-(gamma+1))*((1-delta)+alfa*kss^(alfa-1))^2;

%plot evolution of kt from both states
figure1 = figure
t=1:1:200
k1 = 5.15*0.9621.^t+10.3;
k2 = 15.45*0.9621.^t+10.3;
k3 = 10.3*t.^0
plot(t,k1,'r',t,k2,'b',t,k3,'g')
ylim([0 30])
saveas(figure1,'kt1.pdf')

%plot evolution of ct from both states 
figure2 = figure
c1=zeros(1,200);
c2=zeros(1,200);
for i=1:199
c1(i)=(1-delta)*k1(i)+k1(i).^(alfa)-k1(i+1)
c2(i)=(1-delta)*k2(i)+k2(i).^(alfa)-k2(i+1)
end
c1(200)=NaN;
c2(200)=NaN;
plot(t,c1,'m',t,c2,'c')
ylim([0 3])
saveas(figure2,'ct2.pdf')

%plot the phase diagram
figure3 = figure
plot(k1,c1,'or')
hold on
plot(k2,c2,'ob')
xlim([8 26])
ylim([1.5 2.5])
saveas(figure3,'phasefirstpart.pdf')

%shooting
K1=zeros(400,200);
for j=1:400
    K1(j,1)=1/2*kss;
end
    K1(1,2)=1/2*kss;
for j=1:400 
    for i=3:200
    K1(j,i)=(1-delta)*K1(j,i-1)+K1(j,i-1)^(alfa)-(beta*((1-delta)+alfa*K1(j,i-1)^(alfa-1)))^(1/gamma)*((1-delta)*K1(j,i-2)+K1(j,i-2)^(alfa)-K1(j,i-1));
    if K1(j,i)<0
        break
    else
    end
    end
    if abs(K1(j,200)-kss)<0.001
    break
    else
    end
    if K1(j,200)-kss>0
    K1(j+1,2)=K1(j,2)-1/j;
    elseif K1(j,200)-kss<0
    K1(j+1,2)=K1(j,2)+1/j;
    end
end

C1=zeros(400,199);
for j=1:400
    C1(j,1)=NaN;
    C1(j,200)=NaN;
    for i=2:199
        C1(j,i)=(1-delta)*K1(j,i)+K1(j,i).^(alfa)-K1(j,i+1);
    end
end

K2=zeros(400,200);
for j=1:400
    K2(j,1)=3/2*kss;
end
    K2(1,2)=3/2*kss;
for j=1:400 
    for i=3:200 
    K2(j,i)=(1-delta)*K2(j,i-1)+K2(j,i-1)^(alfa)-(beta*((1-delta)+alfa*K2(j,i-1)^(alfa-1)))^(1/gamma)*((1-delta)*K2(j,i-2)+K2(j,i-2)^(alfa)-K2(j,i-1));
    if K2(j,i)<0
        break
    else
    end
    end
    if abs(K2(j,200)-kss)<0.001
    break
    else
    end
    if K2(j,200)-kss>0
    K2(j+1,2)=K2(j,2)-1/j;
    elseif K2(j,200)-kss<0
    K2(j+1,2)=K2(j,2)+1/j;
    end
end
C2=zeros(400,200);
for j=1:400
    C2(j,1)=NaN;
    C2(j,200)=NaN;
    for i=2:199
        C2(j,i)=(1-delta)*K2(j,i)+K2(j,i).^(alfa)-K2(j,i+1);
    end
end

capital=figure
M=200
plot(1:M,K1(1,:))
hold on
plot(1:M,K1(100,:))
hold on
plot(1:M,K1(200,:))
hold on
plot(1:M,K1(304,:))
capital2=figure
plot(1:M,K2(1,:))
hold on
plot(1:M,K2(100,:))
hold on
plot(1:M,K2(200,:))
hold on
plot(1:M,K2(314,:))
ylim([0 90])
saveas(capital,'capital.pdf')
saveas(capital2,'capital2.pdf')

consumption=figure
M=200
plot(1:M,C1(1,:))
hold on
plot(1:M,C1(100,:))
hold on
plot(1:M,C1(200,:))
hold on
plot(1:M,C1(304,:))
consumption2=figure
plot(1:M,C2(1,:))
hold on
plot(1:M,C2(100,:))
hold on
plot(1:M,C2(200,:))
hold on
plot(1:M,C2(314,:))
saveas(consumption,'cons.pdf')
saveas(consumption2,'cons2.pdf')

phase=figure
plot(K1(1,:),C1(1,:),'o')
hold on
plot(K1(100,:),C1(100,:),'o')
hold on
plot(K1(200,:),C1(200,:),'o')
hold on
plot(K2(1,:),C2(1,:),'o')
hold on
plot(K2(100,:),C2(100,:),'o')
hold on
plot(K2(200,:),C2(200,:),'o')
saveas(phase,'phase.pdf')

phase2=figure
plot(K1(304,:),C1(304,:),'og')
hold on
plot(K2(314,:),C2(314,:),'om')
saveas(phase2,'phase2.pdf')

phaseboth=figure
plot(k1,c1,'or')
hold on
plot(K1(304,:),C1(304,:),'ob')
saveas(phaseboth,'phaseboth.pdf')
phaseboth2=figure
plot(k2,c1,'or')
hold on
plot(K2(314,:),C2(314,:),'ob')
saveas(phaseboth2,'phaseboth2.pdf')
