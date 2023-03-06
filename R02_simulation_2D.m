%% 

clc;
clear;
close all;

addpath generated_function\;
addpath subfunction\;
load param_2D.mat;

%dt = eps^0.2;   % sampling period
tspan = [0 60];
reltol = 10^(-7);   abstol = 10^(-7);   % tolerances for ode45


%% Set Up 
r0 = l0*1.01;
theta0 = -85 *pi/180;

%% Simulation
% x = [q; q_dot]; q = [x; y];
% p = [r; theta];

p0 = [r0; theta0];
q0 = [p0(1)*cos(p0(2));...
      p0(1)*sin(p0(2))];
p_dot_initial = [0.01; 0.01];
x0 = [q0; func_J_2D(q0) \ p_dot_initial]; %initial condition

ode_option = odeset('reltol', reltol, 'abstol', abstol);
[t, X] = ode45(@func_f_2D, tspan, x0, ode_option);



%% Figure
figure()
plot(X(:,1), X(:,2));
xlabel('$X [m]$', 'interpreter', 'latex', 'fontsize', 12');
ylabel('$Y [m]$', 'interpreter', 'latex', 'fontsize', 12');
grid on; grid minor;
axis equal;
title('2D elastic pendulum trajectory',...
    'interpreter', 'latex', 'fontsize', 12')

save('database\chaotic_data_2D.mat', 'X');
disp('Done.')