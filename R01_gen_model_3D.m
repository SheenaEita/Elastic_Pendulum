clc;
clear;
close all;

%%  Save Parameters for the Swinging Spring Model
m = 1;
g = 9.81;
k = 4*pi*pi;
l0 = 1;

save('generated_function\param.mat', 'm', 'l0', 'k', 'g');

%%  Define Coordinate
x = sym('x', [6, 1], 'real');
q = x(1:3); q_dot = x(4:6); 
% q = [x; y; z]
    % q(1) = x; q(2) = y; q(3) = z
    % q_dot(1) = x_dot; q_dot(2) = y_dot; q_dot(3) = z_dot

r = sqrt( q(1)*q(1) + q(2)*q(2) + q(3)*q(3) );
p = [r; acos(q(3)/r); atan2(q(2), q(3))]; 
% p = [r; theta; phi]
    % r = x^2 + y^2 + z^2
    % theta = arccos(z/r)
    % phi = arctan2(y, z)

J = jacobian(p, q);


%% Calculate Lagrangian
% T = 1/2 * m * (x_dot^2 + y_dot^2 + z_dot^2)
T = 0.5*m* (q_dot(1)*q_dot(1) + q_dot(2)*q_dot(2) + q_dot(3)*q_dot(3));

% V = mgz + 1/2 * k * (r-l0)^2 ; r = sqrt(x^2 + y^2 + z^2);
V = m*g*q(3) + 0.5 * k * (r-l0)*(r-l0);

L = T - V;


%% Simplify Equation
%state space : d/dt[q; q_dot] = [q_dot; f] 
%%% find f %%%
M = blkdiag(eye(3), jacobian(jacobian(L, q_dot), q_dot));
f = [q_dot; jacobian(L, q)' - jacobian(jacobian(L, q_dot), q) *q_dot];

%use matlabFunction to simplify and save the equation
syms t real
matlabFunction(simplify(M), 'file', 'generated_function\func_M', 'vars', {t, x});
matlabFunction(simplify(J), 'file', 'generated_function\func_J', 'vars', {q});
matlabFunction(simplify(f), 'file', 'generated_function\func_f', 'vars', {t, x});

disp('Done.');