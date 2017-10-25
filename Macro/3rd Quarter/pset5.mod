var y c inv l rn pi k q z mu x v;
varexo ex ev;

parameters alpha phi little_c theta beta phi_pi mu_ss delta rho_v rho_x r_ss share_c share_inv y_over_k ni omega;

alpha=1/3; 
phi=1; 
little_c=2; 
theta=0.9; 
beta=0.99; 
phi_pi=1.1; 
mu_ss=0.1; 
delta=0.025; 
rho_v=0.75; 
rho_x=0.95;
y_over_k=(beta^(-1)+delta-1)/alpha;
share_inv=delta*(y_over_k)^(-1);
share_c=1-share_inv;
ni=((alpha*y_over_k)*(1+mu_ss)^(-1)+1-delta)^(-1);
omega=(1-theta)*(1-theta*beta)*(theta)^(-1);
r_ss=(beta)^(-1)-1;

model;
//AD
y = share_c*c + share_inv*inv;
c = -(rn-pi(+1))+c(+1);
inv-k(-1) = (delta*little_c)^(-1)*q + x;
q = (r_ss)^(-1)*((1-ni)*z(+1)-(rn-pi(+1))+ni*q(+1)); 
z = y-k(-1)-mu;

//AS
y = alpha*k + (1-alpha)*l;
y-l = mu + phi*l + c;
pi = omega*(-mu) + beta*pi(+1);
k = delta*inv + (1-delta)*k(-1);

//MP
rn = phi_pi*pi+v;

//SHOCKS
x = rho_x*x(-1)+ex;
v = rho_v*v(-1)+ev;
end;

initval;
y = 0; 
c = 0; 
inv = 0; 
l = 0;
rn = 0;
pi = 0;
k = 0;
q=0;
z=0;
mu=0;
x=0;
v=0;
end;

shocks;
var ev=0.25^2;
end;

steady;

stoch_simul(irf=50,order=1);