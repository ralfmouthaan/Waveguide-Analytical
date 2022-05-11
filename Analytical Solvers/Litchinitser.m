% Ralf Mouthaan
% University of Cambridge
% April 2021
% 
% Model of an antiresosonant nodeless photonic crystal fibre. I am using
% the equations defined in Jan Heck's thesis to predict where the
% antiresonant and resonant conditions occur. The waveguide should guide at
% antiresonant wavelenths.

clc; clear variables; close all;

d = linspace(100e-9, 1000e-9, 1000); % Struct thickness
n_silica = 1.45;
%n_core = 1.358; % Pentane
n_core = 1.33;


figure;
for m = 0:10
    
    lambda_ares = 4*d/(2*m+1)*sqrt(n_silica^2 - n_core^2);
    plot(d*1e9, lambda_ares*1e9, 'k-');
    
    hold on
    
    if m >= 1
        lambda_res = 4*d/2/m*sqrt(n_silica^2 - n_core^2);
        plot(d*1e9, lambda_res*1e9, 'k--');
    end
    
end

xlabel('Struct Thickness (nm)');
ylabel('(anti-)resonant wavelength (nm)');

title(['n_{core} = ' num2str(n_core) '; n_{silica} = ' num2str(n_silica)]);

xline(381, 'r');
yline(425, 'r');

ylim([200 2500]);