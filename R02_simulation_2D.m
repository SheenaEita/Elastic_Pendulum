%% 

clc;
clear;
close all;

addpath generated_function\;
addpath subfunction\;
load param_2D.mat;

dt = eps^0.2;   % sampling period
reltol = eps^0.6;   abstol = eps^0.8;   % tolerances for ode45


%% Set Up 
% Chaos
% r0 = l0*0.01;
% theta0 = -3 *pi/180;

% Regular
r0 = l0*1.01;
theta0 = -85 *pi/180;
%%----------------------------- D I V I D E R -----------------------------

%% Simulation
% x = [q; q_dot]; q = [x; y];
% p = [r; theta];

p0 = [r0; theta0];
q0 = [p0(1)*cos(p0(2));...
      p0(1)*sin(p0(2))];
p_dot_initial = [0.01; 0.01];
x0 = [q0; func_J_2D(q0) \ p_dot_initial]; %initial condition

% q0_1 = [0.04; 0; 0.08-l0];
% x0_1 = [q0_1; 0; 0.03427; 0];
% 
% q0_2 = [0.041; 0; 0.08-l0];
% x0_2 = [q0_2; 0; 0.03427; 0];

opt = odeset('mass', @func_M_2D, 'reltol', reltol, 'abstol', abstol);
[t, X] = ode45(@func_f_2D, 0:dt:120, x0, opt);



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

%% Animation
figure()
pause(10)
for i = 1:8000
    plot([0, X(i, 1)], [0, X(i, 2)], 'color', '#AE2012', 'LineWidth', 1);
    hold on;
    plot(0, 0, '.');
    hold on;
    plot(X(i, 1), X(i, 2), 'o');
    hold on;
    plot(X(1:i, 1), X(1:i, 2));
    xlabel('$X [m]$', 'interpreter', 'latex', 'fontsize', 12');
    ylabel('$Y [m]$', 'interpreter', 'latex', 'fontsize', 12');
    grid on; grid minor;
    axis equal;
    axis([-2 2 -1.5 0.2]);
    title('2D elastic pendulum trajectory',...
        'interpreter', 'latex', 'fontsize', 12');
    pause(0.01)
    hold off;
end