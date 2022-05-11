% Ralf Mouthaan
% University of Cambridge
% July 2021
% 
% Script to plot losses for ARROW waveeguides as per Archambault paper.
% * Archambault, Black, Lacroix & Bures "Loss Calculations for Antiresonant
% Waveguides", JLWT, 1993.

clc; clear variables; close all;

%% User-Defined

lambda = 633e-9;
l = 0; m = 1; % Mode indices we're solving for.

% Geometry corresponding to Arrambault Fig. 2
n = [1 1.5 1];
rho = [5*lambda (5 + 2.46)*lambda 10*lambda];

% Our experimental geometry?
%n = [1 1.45 1];
%rho = [15e-6 (15e-6 + 200e-9) 30e-6];

%%

N = length(n);
U_infty = besselzero(l,m,1);
u = U_infty./rho;
U = u*rho(1);
alpha2 = lambda * U_infty * U(1)/pi/n(1)/rho(1)^2/U(2)
alpha3_min = alpha2 * U(3)/U(2)