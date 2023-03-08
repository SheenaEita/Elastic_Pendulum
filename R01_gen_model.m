clc;
clear;
close all;

m = 1; % mass [kg]
g = 9.81; % gravity [m/s^2]
k = 36; % elastic constant [N/m]
l0 = 1; % initial length [l]

save("gen_function\param.mat","l0","k","g","m");

%% System coordinates
x = sym('x', [4, 1], 'real');
y = sym('y', [4, 1], 'real');
q = x(1:2); q_dot = x(3:4);
    % q = [r; theta]
    % q(1) = r; q(2) = theta;
    % q_dot(1) = r_dot; q_dot(2) = theta_dot;
p = [q(1) * cos(q(2));...
     q(1) * sin(q(2))];

J = jacobian(p, q);
p_dot = J * q_dot;
y(1:2) = p; y(3:4) = p_dot;


%% Lagrangian
T = 0.5*m*(p_dot(1)*p_dot(1) + p_dot(2)*p_dot(2));
% T = 1/2 * m * (x_dot^2 + y_dot^2);

V = m*g*p(2) + 0.5*k*(q(1)-l0)*(q(1)-l0);
% V = mgz + 1/2 * k * (r-l0)^2 ;

L = T - V;

%% Generate and simplify ODE and jacobian
f = [q_dot; jacobian(L, q)' - jacobian(jacobian(L, q_dot), q)*q_dot];
%state space : d/dt[q; q_dot] = [q_dot; f] 

% use matlabFunction to simplify and save the function
syms t real;
matlabFunction(simplify(J), 'file', 'gen_function\func_J', 'vars', {q});
matlabFunction(simplify(f), 'file', 'gen_function\func_f', 'vars', {t, x});

disp('Done.');