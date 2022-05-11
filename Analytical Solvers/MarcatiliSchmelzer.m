% Ralf Mouthaan
% University of Cambridge
% June 2021
%
% Script to calculate propagation constants using Marcatili-Schmelzer model

clc; clear variables; close all;

%% User-Defined

% Marcatili-Schmelzer Example
lambda = 1e-6; 
a = 1e-3;
v = 1.5;

% Erlangen ARPCF
lambda = 633e-9;
a = 25e-6/2;
v = 1.45;

%% Build up array of propagation constants

k0 = 2*pi/lambda;
arrIdent = {};
arrBeta = [];

idx = 0;

for n = 0:5
    for m = 1:5
        
        u_nm = besselzero(n-1,m,1);
        u_nm = u_nm(m);
        
        beta_real = 2*pi/lambda * (1-1/2*(u_nm * lambda/2/pi/a)^2);
        
        if n == 0
        
            % TE Mode
            idx = idx + 1;
            beta_imag = (u_nm/2/pi)^2 * lambda^2/a^3 / sqrt(v^2 - 1);
            strIdent{idx} = ['TE' num2str(n) num2str(m)];
            arrBeta(idx) = beta_real + 1i*beta_imag;

            % TM Mode
            idx = idx + 1;
            beta_imag = (u_nm/2/pi)^2 * lambda^2/a^3 *v^2/ sqrt(v^2 - 1);
            strIdent{idx} = ['TM' num2str(n) num2str(m)];
            arrBeta(idx) = beta_real + 1i*beta_imag;
        
        end
        
        % EH Mode
        idx = idx + 1;
        beta_imag = (u_nm/2/pi)^2 * lambda^2/a^3 /2*(v^2+1)/ sqrt(v^2 - 1);
        strIdent{idx} = ['EH' num2str(n) num2str(m)];
        arrBeta(idx) = beta_real + 1i*beta_imag;
        
    end
end

%% Sort

[~,idx] = sort(real(arrBeta), 'descend');
arrBeta = arrBeta(idx);
strIdent = strIdent(idx);

%% Display

for i = 1:length(arrBeta)
    fprintf('%s; Effective n = %0.10f; Loss = %0.3f dB/m\n', strIdent{i}, real(arrBeta(i))/k0, imag(arrBeta(i))*20/log(10));
end