var y l k a c i;
varexo e;

parameters beta g alpha delta phi rho sigmaa;

beta   = 0.984;
g      = 0.004;
alpha  = 1/3;
delta  = 0.025;
phi    = 1;
rho    = 0.95;
sigmaa = 0.0072;

model;
exp(y) = exp(a)^(1-alpha)*exp(k(-1))^(alpha)*exp(l)^(1-alpha);
(1-alpha)*exp(k(-1))^(alpha)*exp(l)^(-alpha)*exp(a)^(1-alpha) = exp(c)*exp(l)^(phi);
exp(c)^(-1) = beta*(1+g)^(-1)*(exp(c(+1))^(-1)*(alpha*exp(y(+1))/exp(k)+1-delta));
exp(k) = (1+g)^(-1)*(exp(y)+(1-delta)*exp(k(-1))-exp(c));
a = rho*a(-1)+e;
exp(i)=exp(y)-exp(c);
end;

initval;
k =log(18.3591);
y =log(2.4964);
c =log(1.9640);
l =log(0.9205);
a = 0;
end;

shocks;
var e=sigmaa^2;
end;

steady;

stoch_simul(periods=1000,irf=100,order=1);
