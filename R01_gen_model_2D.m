clc;
clear;
close all;

%%  Save Parameters for the Swinging Spring Model
m = 1;
g = 9.8;
k = 4*pi*pi;
l0 = 1;

save('generated_function\param_2D.mat', 'm', 'l0', 'k', 'g');

%%  Define Coordinate
x = sym('x', [4, 1], 'real');
q = x(1:2); q_dot = x(3:4); 
% q = [x; y]
    % q(1) = x; q(2) = y;
    % q_dot(1) = x_dot; q_dot(2) = y_dot;

r = sqrt( q(1)*q(1) + q(2)*q(2) );
p = [r; atan2(q(2), q(1))]; 
% p = [r; theta; phi]
    % r = x^2 + y^2
    % theta = arctan(y/x)

%把卡式座標轉到極座標
J = jacobian(p, q);


%% Calculate Lagrangian
% T = 1/2 * m * (x_dot^2 + y_dot^2 + z_dot^2)
T = 0.5*m* (q_dot(1)*q_dot(1) + q_dot(2)*q_dot(2));

% V = mgz + 1/2 * k * (r-l0)^2 ; r = sqrt(x^2 + y^2 + z^2);
V = m*g*q(2) + 0.5 * k * (r-l0)*(r-l0);

L = T - V;


%% Simplify Equation
%state space : d/dt[q; q_dot] = [q_dot; f] 
%%% find f %%%
M = blkdiag(eye(2), jacobian(jacobian(L, q_dot), q_dot));
f = [q_dot; jacobian(L, q)' - jacobian(jacobian(L, q_dot), q) *q_dot];

%use matlabFunction to simplify and save the equation
syms t real
matlabFunction(simplify(M), 'file', 'generated_function\func_M_2D', 'vars', {t, x});
matlabFunction(simplify(J), 'file', 'generated_function\func_J_2D', 'vars', {q});
matlabFunction(simplify(f), 'file', 'generated_function\func_f_2D', 'vars', {t, x});

disp('Done.');