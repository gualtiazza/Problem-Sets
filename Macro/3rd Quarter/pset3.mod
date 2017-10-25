var y y_star y_tilda pi rn r r_star g a;
varexo ea eg;

parameters share_g lambda kappa phi beta theta phi_pi rhoa rhog;

share_g = 0.2;
phi = 1;
beta = 0.99;
theta = 0.9;
phi_pi =1.1;
kappa = phi+(1-share_g)^(-1);
lambda = (1-theta)*(1-beta*theta)/theta;
rhoa = 0.95;
rhog = 0.95;

model;
y = y(+1)+share_g*(g-g(+1))-(1-share_g)*r(+1);
r = rn-pi(+1); 
y_star = (1+phi)*(kappa)^(-1)*a + share_g*(1-share_g)^(-1)*(kappa)^(-1)*g;
r_star = (1-share_g)^(-1)*(kappa)^(-1)*(1+phi)*(a(+1)-a)-share_g*(1-share_g)^(-1)*(1-(1-share_g)^(-1)*(kappa)^(-1))*(g(+1)-g);
pi = lambda*kappa*y_tilda+beta*pi(+1);
y_tilda = y-y_star;
rn = phi_pi*pi;
a = rhoa*a(-1)+ea;
g = rhog*g(-1)+eg;
end;

initval;
y = 0; 
y_star = 0; 
y_tilda = 0; 
pi = 0;
rn = 0;
r = 0;
r_star = 0;
a=0;
g=0;
end;

shocks;
var ea=0.05^2;
end;

steady;

stoch_simul(irf=100,order=1);