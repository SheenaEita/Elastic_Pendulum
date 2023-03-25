%% 

clc;
clear;
close all;

addpath generated_function\;
addpath subfunction\;
load param.mat;

dt = eps^0.2;   % sampling period
reltol = eps^0.6;   abstol = eps^0.8;   % tolerances for ode45


%% Set Up 
r0 = l0;
theta0 = 175 *pi/180;
phi0 = 5 *pi/180;

%%----------------------------- D I V I D E R -----------------------------

%% Simulation
% x = [q; q_dot]; q = [x; y; z];
% p = [r; theta; phi];

% p0 = [r0; theta0; phi0];
% q0 = [p0(1)*sin(p0(2))*cos(p0(3));...
%     p0(1)*sin(p0(2))*sin(p0(3));...
%     p0(1)*cos(p0(2))];
% p_dot_initial = [0.01; 0; -0.01];
% x0 = [q0; func_J(q0) \ p_dot_initial]; %initial condition

q0_1 = [0.04; 0; 0.08-l0];
x0_1 = [q0_1; 0; 0.03427; 0];

q0_2 = [0.041; 0; 0.08-l0];
x0_2 = [q0_2; 0; 0.03427; 0];

opt = odeset('mass', @func_M, 'reltol', reltol, 'abstol', abstol);
[t_1, X_1] = ode45(@func_f, 0:dt:25, x0_1, opt);

[t_2, X_2] = ode45(@func_f, 0:dt:25, x0_2, opt);


%% Figure
figure()
plot3(X_1(:,1), X_1(:,2), X_1(:,3));
hold on;
plot3(X_2(:,1), X_2(:,2), X_2(:,3));
hold on;
xlabel('$X [m]$', 'interpreter', 'latex', 'fontsize', 12');
ylabel('$Y [m]$', 'interpreter', 'latex', 'fontsize', 12');
zlabel('$Z [m]$', 'interpreter', 'latex', 'fontsize', 12');
grid on; grid minor;
axis equal;
title('Two trajectories starting 1mm apart',...
    'interpreter', 'latex', 'fontsize', 12')

figure()
plot(X_1(:,1), X_1(:,3));
hold on;
plot(X_2(:,1), X_2(:,3));
hold on;
xlabel('$X [m]$', 'interpreter', 'latex', 'fontsize', 12');
ylabel('$Z [m]$', 'interpreter', 'latex', 'fontsize', 12');
grid on; grid minor;
axis equal;
title('Two trajectories starting 1mm apart',...
    'interpreter', 'latex', 'fontsize', 12')

save('database\chaotic_data.mat', 'X_1', 'X_2')