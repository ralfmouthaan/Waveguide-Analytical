% Ralf Mouthaan
% University of Cambridge
% August 2021
%
% Modelling of a planar waveguide. Based on Snydder + Love Section 12-3.
% Note: Newton-Raphson doesn't quite converge for higher order values

clc; clear variables; close all;

%%

width = 25e-6;
a = width/2;
lambda = 633e-9;
n_co =  1.457; % Silica
n_cl =  1.4403;
k = 2*pi/lambda;
NA = sqrt(n_co^2 - n_cl^2); 
V = a*k*NA;
Delta = (1 - n_cl^2/n_co^2)/2;
bolDebug = false;

%% Waveguide data
disp('PLANAR WAVEGUIDE');
disp(['width = ' num2str(width*1e6) 'um']);
disp(['Wavelength = ' num2str(lambda*1e9) 'nm']);
disp(['Core refractive index = ' num2str(n_co)]);
disp(['Cladding refractive index = ' num2str(n_cl)]);
disp(['Numerical Aperture = ' num2str(NA)]);
disp(['Profile height parameter = ' num2str(Delta)]);
disp(['V-number = ' num2str(V)]);
disp(['Number of modes = ' num2str(4*V/pi)]);
disp(['Beta min = ' num2str(k*n_cl) '; Beta max = ' num2str(k*n_co)]);

if n_co < n_cl
    disp('WARNING: n_co must be > n_cl');
    return
end

if Delta > 1
    disp('WARNING: Waveguide is not weakly guiding');
    return
end

%% Even Eigenvalue equation

% First, determine eigenequation values for range of U values.
U = 0:0.001:V;
W = sqrt(V^2 - U.^2);

% Plot function
if bolDebug == true
    figure;
    plot(U, abs(W - U.*tan(U)));
    ylim([-30 30]);
end

% Get approximate zeros of eigenequation;
[~,idxs] = findpeaks(-abs(W - U.*tan(U)));
U = U(idxs);
W = W(idxs);

% Improve estimates of U with Newton-Raphson
for i = 1:5
    Uapprox = U(i);
    Wapprox = W(i);
    for j = 1:100
        f = Uapprox*tan(Uapprox) - Wapprox;
        dfdU = Uapprox/Wapprox + tan(Uapprox) + Uapprox*sec(Uapprox)^2;
        Uapprox = Uapprox - f/dfdU;
        Wapprox = sqrt(V^2 - Uapprox^2);
    end
    U(i) = Uapprox;
    W(i) = Wapprox;
end

if bolDebug == true
    for i = 1:5
        xline(U(i), ':');
    end
end

fprintf('EVEN MODES\n');
for i = 1:3
    fprintf('U = %0.5f; W = %0.5f\n', U(i), W(i));
end

W - U.*tan(U)

%% Plot Even Modes

x = -2:0.01:2;

figure;
for i = 1:3
    E(abs(x) < 1) = cos(U(i)*x(abs(x) < 1))/cos(U(i));
    E(abs(x) >= 1) = exp(-W(i)*abs(x(abs(x) >= 1)) + W(i));
    plot(x, E);
    hold on
end

xline(-1, ':');
xline(1, ':');

%% Odd Eigenvalue Equation

% First, determine eigenequation values for range of U values.
U = 0:0.001:V;
W = sqrt(V^2 - U.^2);

% Plot function
if bolDebug == true
    figure;
    plot(U, abs(W + U.*cot(U)));
    ylim([-30 30]);
end

% Get approximate zeros of eigenequation;
[~,idxs] = findpeaks(-abs(W + U.*cot(U)));
U = U(idxs);
W = W(idxs);

% Improve estimates of U with Newton-Raphson
for i = 1:5
    Uapprox = U(i);
    Wapprox = W(i);
    for j = 1:100
        f = Uapprox*cot(Uapprox) + Wapprox;
        dfdU = -Uapprox/Wapprox + cot(Uapprox) - Uapprox*csc(Uapprox)^2;
        Uapprox = Uapprox - f/dfdU;
        Wapprox = sqrt(V^2 - Uapprox^2);
    end
    U(i) = Uapprox;
    W(i) = Wapprox;
end

if bolDebug == true
    for i = 1:5
        xline(U(i), ':');
    end
end

fprintf('ODD MODES\n');
for i = 1:3
    fprintf('U = %0.5f; W = %0.5f\n', U(i), W(i));
end

%% Plot Odd Modes

x = -2:0.01:2;

figure;
for i = 1:3
    E(abs(x) < 1) = sin(U(i)*x(abs(x) < 1))/sin(U(i));
    E(abs(x) >= 1) = x(abs(x) >= 1)./abs(x(abs(x) >= 1)).*exp(-W(i)*abs(x(abs(x) >= 1)))/exp(-W(i));
    plot(x, E);
    hold on
end

xline(-1, ':');
xline(1, ':');
