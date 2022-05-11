% Ralf Mouthaan
% University of Cambridge
% August 2019
%
% Code to calculate spatial LP modes and associated propagation constants
% for a step index fiber.
% See:
%   * Snyder & Love Ch. 13 + Ch. 14
%   * Gloge, Applied Optics, 1971.

function StepIndex

clc; clear variables; close all;
global l V;
bolDebug = true;

%% User-Defined Parameters + Derived Parameters

diameter = 25e-6;
a = diameter/2;
lambda = 633e-9;
n_co =  1.457; % Silica
n_cl =  1.4403;
k = 2*pi/lambda;
NA = sqrt(n_co^2 - n_cl^2); 
V = a*k*NA;
Delta = (1 - n_cl^2/n_co^2)/2;
l = 2;
m = 2;

%% Waveguide data
disp('STEP INDEX CYLINDRICAL WAVEGUIDE');
disp(['Diameter = ' num2str(diameter*1e6) 'um']);
disp(['Wavelength = ' num2str(lambda*1e9) 'nm']);
disp(['Core refractive index = ' num2str(n_co)]);
disp(['Cladding refractive index = ' num2str(n_cl)]);
disp(['Numerical Aperture = ' num2str(NA)]);
disp(['Profile height parameter = ' num2str(Delta)]);
disp(['V-number = ' num2str(V)]);
disp(['Number of modes = ' num2str(V^2/2)]);
disp(['LP(' num2str(l) ',' num2str(m) ') mode']);
disp(['Beta min = ' num2str(k*n_cl) '; Beta max = ' num2str(k*n_co)]);

if n_co < n_cl
    disp('WARNING: n_co must be > n_cl');
    return
end

if l < 0
    disp('WARNING: l must be >= 0');
    return
end

if m < 1
    disp('WARNING: m must be >= 1');
    return
end

if Delta > 1
    disp('WARNING: Waveguide is not weakly guiding');
    return
end

%% Eigenvalue equation

% The eigenvalue equation corresponds to the boundary matching conditions
% between the core and the cladding. The values U that satisfy the
% eigenvalue equation correspond to propagating modes.

% The eigenvalue equation f is evaluated at a range of U-values between 0
% and V.
arrU = 0:0.001:V;
funEigEq = @(x) EigenvalueEquation(x);

if bolDebug == true
    figure;
    plot(arrU, EigenvalueEquation(arrU), 'rx');
    hold on
    fplot(funEigEq, [0 V]);
    ylim([-100 100])
end

% The findpeaks() function is used on arrEigEq to find the approximate
% position of the roots
[~,idxs] = findpeaks(-abs(EigenvalueEquation(arrU)));

if m > length(idxs)
    disp('WARNING: This mode cannot be supported in the waveguide');
    return
end

% fzero() is used to get a better approximation for the root of interest
% U = fzero(funEigEq, arrU(idxs(m)));
U = arrU(idxs(m));

% Associated values are then calculated.
W = sqrt(V^2 - U^2);
beta = V/a/(2*Delta)^(1/2)*(1-2*Delta*U^2/V^2)^(1/2);

% Print derived parameters
disp(['U-number = ' num2str(U)]);
disp(['W-number = ' num2str(W)]);
disp(['Propagation Constant = ' num2str(beta) '/m']);

if beta < k*n_cl || beta > k*n_co
    disp('WARNING: This mode is not supported in the waveguide');
    return
end

%% Plot

x = linspace(-2*a,2*a,500);
[x_mesh, y_mesh] = meshgrid(x, x.');
normradius = sqrt(x_mesh.^2 + y_mesh.^2)/a;
theta = atan2(y_mesh, x_mesh);

F = zeros(size(normradius));
F(normradius <= 1) = besselj(l, U.*normradius(normradius<=1))./besselj(l, U);
F(normradius >  1) = besselk(l, W.*normradius(normradius> 1))./besselk(l, W);

F = F.*cos(l*theta);

figure;
imagesc(x*1e6, x*1e6, abs(F));
xlabel('\mum');
ylabel('\mum');
colormap hot
axis square
title({['|E| - LP(' num2str(l) ',' num2str(m) ')']; ...
    ['\phi = ' num2str(diameter*1e6) '\mum; NA = ' num2str(NA) ...
    '; \lambda = ' num2str(lambda*1e9) 'nm; \beta = ' ...
    num2str(beta*1e-6) '\mum^{-1}']})

figure;
imagesc(x*1e6, x*1e6, angle(F));
xlabel('\mum');
ylabel('\mum');
axis square
title({['\angleE - LP(' num2str(l) ',' num2str(m) ')']; ...
    ['\phi = ' num2str(diameter*1e6) '\mum; NA = ' num2str(NA) ...
    '; \lambda = ' num2str(lambda*1e9) 'nm; \beta = ' ...
    num2str(beta*1e-6) '\mum^{-1}']})


end
function RetVal = EigenvalueEquation(U)

    % Eigenvalue equation as given in Snyder & Love. Is equivalent to the
    % Eigenvalue equation given in Gloge, who just used a different
    % recursive bessel function relationship for the derivation.

    global l V;

    W = sqrt(V^2 - U.^2);
    RetVal = U.*besselj(l+1,U)./besselj(l,U) - ...
        W.*besselk(l+1,W)./besselk(l,W);

end