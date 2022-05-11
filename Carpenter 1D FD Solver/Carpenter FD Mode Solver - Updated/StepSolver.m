% Ralf Mouthaan
% University of Cambridge
% September 2019
%
% Finite difference mode solver for step index fiber.
% Uses Joel Carpenter's finite difference engine for the moment.

clc; clear variables; close all;

%% User-defined paramters

a = 50e-6; % fibre diameter
dr = 0.1e-6; % step size
r = 0:dr:3*a; % radius vector
n = zeros(size(r)); % refractive index vector
n(r<=a/2) = 1.45; % Core refractive index
n(r>a/2) = 1.44; % Cladding refractive index
lambda = 1000e-9;
L = 2;

%% Calculation

[beta, R, tau, A, r] = RadialModeSolver2(n, r, L, lambda, 1);

plot(r*1e6, R);
xlabel('r (\mum)')
xlim([0 50])
legend

disp(beta);