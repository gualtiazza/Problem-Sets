var y c rnom pi a l mu p m v r y_star r_star gap;
varexo em ea;

parameters phi am theta beta phipi rhoa rhom;

phi = 1;
am = 0.001;
theta = 0.9;
beta = 0.99;
phipi = 1.1;
rhoa = 0.95;
rhom = 0.75;

model;
y = c;
c = c(+1)-rnom+pi(+1);
phi*l+c = y-l-mu;
y = a+l;
p-p(-1) = pi;
pi = ((theta-1)*(1-theta*beta)*mu/theta)+beta*theta*pi(+1);
rnom = phipi*pi+v;
m-p = c-(rnom/((1/beta)-1));
a = rhoa*a(-1)+ea;
v = rhom*v(-1)+em;
r = rnom-pi(+1);
y_star = a;
r_star = a(+1)-a;
gap = -mu/(1+phi);
end;

initval;
y = 0;
c = 0;
rnom = 0;
pi = 0;
l = 0;
mu = 0;
p = 0;
m = 0;
r = 0;
a = 0;
v = 0;
y_star = 0;
r_star = 0;
end;

shocks;
var em=0.01^2;
end;

steady;

stoch_simul(irf=20,order=1);