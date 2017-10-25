%% Problem Set 6 Code - Gualtiero Azzalini %%

%% Problem 7.1

addpath('/Users/Gualtiero/Dropbox/A NYU/Macro/psets/6')
addpath('/Users/Gualtiero/Documents/MATLAB/Add-Ons/mitbook/utilities_tom')
clear
% Set Parameters

A0 = 100; A1 = 0.05; beta = 0.95; d = 10; H0 = 95.5; H1 = 0.95;

R = [0 -A0/2 0;
     -A0/2 d/2 A1/2;
     0 A1/2 0];
 
Q = d/2;

H = [0 -d/2 0]';

A = [1 0 0;
     0 0 0;
     H0 0 H1];

B = [0 1 0]';

[F,P] = olrp(beta,A,B,R,Q,H);

%% Problem 7.2

% Set Parameters
H0_a=zeros(3,1);H1_a=zeros(3,1); F_a=zeros(3,3);
H0_a(1,1) = 94.0888           ; H1_a(1,1) = 0.9211;
H0_a(2,1) = 93.22             ; H1_a(2,1) = 0.9433;
H0_a(3,1) = 95.08187459215024 ; H1_a(3,1) = 0.95245906270392;
for i=1:3
    A = [1 0 0;
             0 0 0;
             H0_a(i,1) 0 H1_a(i,1)];
         
    [F,P] = olrp(beta,A,B,R,Q,H);
    F_a(i,:) = -F(:) 
end    

%% Problem 7.3 - Max Welfare

R = [A1/2 -A0/2;
     -A0/2 0];
 
Q = d/2;

H = [0 0]';

A = [1 0;
     0 1];

B = [1 0]';

[F,P] = olrp(beta,A,B,R,Q,H);

% Extension part

s0 = -F(1,2); s1 = 1-F(1,1);
Y_bar = s0/(1-s1); Y0 = 0.9*Y_bar; Y1 = 1.1*Y_bar;

% Compute the trajectories
periods = 200;
Y_up = zeros(1,periods); Y_down = zeros(1,periods); 
Y_up(1,1) = Y1; Y_down(1,1) = Y0; 

for i=2:periods
Y_up(1,i) = s0 + s1*Y_up(1,i-1);
Y_down(1,i) = s0 + s1*Y_down(1,i-1);
end

exercise73 = figure

plot(1:periods,Y_up,'k');hold on; plot(1:periods,Y_down,'k');
hold on
ylim([Y0 Y1]); xlim([1,periods]);
hold on
legend({'Y_{1}','Y_{0}'},'FontSize',11);legend('Location','southeast');
xlabel('Y_{t}'); ylabel('Y_{t+1}'); title('Evolution of Y from Y_{0} and Y_{1} - Competition','FontSize',14);
saveas(exercise73,'Figure1.png');

%% Problem 7.4 - Monopoly

R = [0 -A0/2;
     -A0/2 +A1+d/2];
 
Q = d/2;

H = [0 -d/2]';

A = [1 0;
     0 0];

B = [0 1]';

[F,P] = olrp(beta,A,B,R,Q,H);

% Extension part

s0_m = -F(1,1); s1_m = -F(1,2);
Y_bar_m = s0_m/(1-s1_m); Y0_m = 0.9*Y_bar_m; Y1_m = 1.1*Y_bar_m;

% Compute the trajectories
periods = 200;
Y_up_m = zeros(1,periods); Y_down_m = zeros(1,periods); 
Y_up_m(1,1) = Y1; Y_down_m(1,1) = Y0; 

for i=2:periods
Y_up_m(1,i) = s0_m + s1_m*Y_up_m(1,i-1);
Y_down_m(1,i) = s0_m + s1_m*Y_down_m(1,i-1);
end

exercise74 = figure

plot(1:periods,Y_up_m,'k-.');hold on; plot(1:periods,Y_down_m,'k-.');
hold on
ylim([Y0_m Y1]); xlim([1,periods]);
hold on
legend({'Y_{1}','Y_{0}'},'FontSize',11);legend('Location','northeast');
xlabel('Y_{t}'); ylabel('Y_{t+1}'); title('Evolution of Y from Y_{0} and Y_{1} - Monopoly','FontSize',14);
saveas(exercise74,'Figure2.png');

%% Problem 7.5 - Duopoly  

R1 = [0 -A0/2 0;
      -A0/2 A1+d/2 A1/2;
      0 A1/2 0];
R2 = [0 0 -A0/2;
      0 0 A1/2;
      -A0/2 A1/2 A1+d/2];
Q1 = d/2;
Q2 = Q1;
H1 = [0 -d/2 0]';
H2 = [0 0 -d/2]';
S1 = 0;
S2 = 0;
M1 = 0;
M2 = 0;

A = [1 0 0;
     0 0 0;
     0 0 0];
B1 = [0 1 0]';
B2 = [0 0 1]';

[F1,F2,ao] = nash(A,B1,B2,R1,R2,H1,H2,Q1,Q2,S1,S2,M1,M2);

s0_d = -(F1(1,1)+F2(1,1)); s1_d = -(F1(1,2)+F1(1,3)+F2(1,2)+F2(1,3))/2;
Y_up_d = zeros(1,periods); Y_down_d = zeros(1,periods); 
Y_up_d(1,1) = Y1; Y_down_d(1,1) = Y0; 
Y_bar_d = s0_d/(1-s1_d); Y0_d = 0.9*Y_bar_d; Y1_d = 1.1*Y_bar_d

for i=2:periods
Y_up_d(1,i) = s0_d + s1_d*Y_up_d(1,i-1);
Y_down_d(1,i) = s0_d + s1_d*Y_down_d(1,i-1);
end

exercise75 = figure
plot(1:periods,Y_up_d,'k--');hold on; plot(1:periods,Y_down_d,'k--');
hold on
ylim([Y0_d Y1]); xlim([1,periods]);
hold on
legend({'Y_{1}','Y_{0}'},'FontSize',11);legend('Location','northeast');
xlabel('Y_{t}'); ylabel('Y_{t+1}'); title('Evolution of Y from Y_{0} and Y_{1} - Duopoly','FontSize',14);
saveas(exercise75,'Figure3.png');

%% Plot together results from monopoly, perfect competition and duopoly

monopoly_comp_duop = figure
plot(1:periods,Y_up,'k');hold on; plot(1:periods,Y_down,'k');
hold on
plot(1:periods,Y_up_m,'k-.');hold on; plot(1:periods,Y_down_m,'k-.');
hold on
plot(1:periods,Y_up_d,'k--');hold on; plot(1:periods,Y_down_d,'k--');
hold on
ylim([Y0_m Y1]); xlim([1,periods]);
hold on
legend({'Y_{1} - Competition','Y_{0} - Competition','Y_{1} - Monopoly','Y_{0} - Monopoly','Y_{1} - Duoopoly','Y_{0} - Duopoly'},'FontSize',11);legend('Location','southeastoutside');
xlabel('Y_{t}'); ylabel('Y_{t+1}'); title('Evolution of Y - Monopoly, Duopoly, Competition','FontSize',14);
saveas(monopoly_comp_duop,'Figure4.png');


