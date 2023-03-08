

clc;
clear;
close all;

addpath gen_function\;
load param.mat;

tspan = [0 30];
reltol = 10^(-7);
abstol = 10^(-7);
v = 0.2;
alpha = -135 * pi / 180;

r0 = l0 * 1.01;
theta0 = -85 * pi / 180;

%initial state
initial_q = [r0; theta0];
initial_p = [v*cos(alpha); v*sin(alpha)];
initial_cond = [initial_q; func_J(initial_q) \ initial_p]; %initial condition [r_0, theta_0 ]

ode_option = odeset('RelTol', reltol, 'AbsTol', abstol);

[t, X] = ode45(@func_f, tspan, initial_cond, ode_option);


Y = NaN(length(t(:, 1)), 4);
for i = 1:length(t(:, 1))
    Y(i, 1) = X(i, 1)*cos(X(i, 2));
    Y(i, 2) = X(i, 1)*sin(X(i, 2));
    Y(i, 3) = X(i, 3)*cos(X(i, 2)) - X(i, 1)*sin(X(i, 2))*X(i, 4);
    Y(i, 4) = X(i, 3)*sin(X(i, 2)) + X(i, 1)*cos(X(i, 2))*X(i, 4);
end


%% Figure
figure()
plot(Y(:, 1), Y(:, 2));
xlabel('$X [m]$', 'interpreter', 'latex', 'fontsize', 12);
ylabel('$Y [m]$', 'interpreter', 'latex', 'fontsize', 12);
grid on; grid minor;
axis equal;
ylim([-1 1]);
ylim([-2 0.5]);
title('2D elastic pendulum trajectory',...
    'interpreter', 'latex', 'fontsize', 12);

save('database\elastic_pendulum_simulation.mat', "t", "X", "Y");
disp('Done.');
