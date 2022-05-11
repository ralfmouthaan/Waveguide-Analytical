clear all;
close all;

%MMF properties
% alpha parameter of refractive index profile
% 2 is parabolic graded-index, infinite is step-index
alfa = 1.835;
%alfa = 2.0;
% peak germanium doping concentration of core
% cladding is pure silica
doping =10;
% core radius in metres
core_r = 50e-6/2;
% centre wavelength is metres
LAMBDA = 1550e-9;
% number of wavelengths to solve
% at least 3 to get group delay, at least 5 to get chromatic dispersion
lambdaCount = 5;
% wavelength step
dLambda = 0.1e-9;
% wavelength vector to simulate
lambdas = LAMBDA+dLambda.*linspace(-1,1,lambdaCount);
% find the index of the centre wavelength
[minV lambda0Idx] = (min(abs(lambdas-LAMBDA)));

% number of radial points to simulate
r_points = 5001;
% how many core radii away from the centre to simulate
extent = 10;
% radial vector
r = (linspace(0,1,r_points)).*core_r.*extent;
%radial vector in micron
r_um = r*1e6;

%maximum mode-group order to include
maxMG =9;

%maximum possible L and M values for LP(l,m) modes.
l_start = 0;
l_stop = maxMG;
m_start = 1;
m_stop = maxMG;


%Constants
c = 299792458;
u0 = 4.*pi.*10^-7;
e0 = (1./((c.^2).*u0));

% main loop
for lambdaIdx=1:lambdaCount
    lambda = lambdas(lambdaIdx);
    fprintf('%10.10fnm\n',lambda*1e9);
    
    k0 = 2.*pi./lambda;
    f = c/lambda;
    w = 2.*pi.*f;
    
    modeCount = ((l_stop-l_start)+1)*((m_stop-m_start)+1);
    L = l_start:l_stop;
    M = m_start:m_stop;
    
    % fibre properties
    fiber.doping = doping;
    fiber.lambda = lambda;
    fiber.a = core_r;
    fiber.alpha = alfa;
    fiber.cladDefect = 0;
    fiber.function = 'alpha_nSell3';
    
    % refractive index profile
    n_r = feval(fiber.function,r,fiber);
    fiber.nr = n_r;
    fiber.r = r;
    
    % refractive index of core and cladding
    n_clad = min(min(n_r));
    n_core = max(max(n_r));
    % fibre refractive index contrast parameter
    delta = (n_core.^2-n_clad.^2)./(2.*n_core.^2);
    
    % number of simulation points in the core and cladding
    pointsCore = sum(n_r>n_clad);
    pointsClad = sum(n_r<=n_clad);
    
    if (lambdaIdx==1)
        fprintf('Points [Core] : %10.0f\n',pointsCore);
        fprintf('Points [Clad] : %10.0f\n',pointsClad);
    end
    
    % maximum & minimum propagation constants
    beta_min = n_clad*k0;
    beta_max = n_core*k0;
    
    Er = zeros(modeCount,length(r));
    beta = zeros(modeCount,1);
    tau = zeros(modeCount,1);
    
    idx=1;
    disp('	Solving Eigenfunctions...');
    % Solve the eigenfunction for each value of L
    for i=1:length(L)
        [beta_,F_l,tau_,A,r_] = RadialModeSolver2(n_r, r, L(i), lambda,0);
        for l=1:length(M)
            if ((M(l))<=length(beta_))
                %Normalise
                F_l(:,M(l)) = F_l(:,M(l))./sqrt(sum(sum(abs(F_l(:,M(l))).^2)));
                if (F_l(1,M(l))<0)
                    F_l(:,M(l)) = -F_l(:,M(l));
                end
                Er(idx,1:length(F_l)) = squeeze(F_l(:,M(l)));
                
                beta(idx) = beta_(M(l));
                tau(idx)  = tau_(M(l));
            else
                beta(idx) = 0;
                tau(idx) = -1;
            end
            idx = idx+1;
        end
        fprintf('		%10.1f percent complete\n', 100.0.*i./length(L));
    end
    
    % Check which modes are actually bound
    isBound_L = zeros(modeCount,1);
    isBound_M = zeros(modeCount,1);
    isBound = zeros(modeCount,1);
    
    totalModes=0;
    beta = abs(beta);
    k=1;
    for i=1:length(L)
        for j=1:length(M)
            MG = 2.*M(j)+L(i)-1;
            [maxPos maxIdx] = max(abs(Er(k,:)));
            if (beta(k)>beta_max || beta(k)<beta_min || r(maxIdx)>core_r || MG>maxMG)
                if (beta(k)==0)
                    isBound(k)=-1;
                else
                    isBound(k)=0;
                end
                isBound_L(k) = nan;
                isBound_M(k) = nan;
            else
                isBound_L(k) = L(i);
                isBound_M(k) = M(j);
                isBound(k) = 1;
                totalModes=totalModes+1;
            end
            k=k+1;
        end
    end
    fprintf('%10.0f total modes\n',totalModes);
    
    MODES = 1:modeCount;
    modeIdxes = MODES(~isnan(isBound_L) & ~isnan(isBound_M));
    clear MODES;
    
    % discard unbound modes
    L = isBound_L(modeIdxes);
    M = isBound_M(modeIdxes);
    beta = beta(modeIdxes);
    tau = tau(modeIdxes);
    Er = Er(modeIdxes,:);
    
    % sort by mode-group rather than by L
    MG = 2.*M+L-1;
    [mg mgIdxes] = sort(MG);
    L = L(mgIdxes);
    M = M(mgIdxes);
    Er = (Er(mgIdxes,:));
    beta = (beta(mgIdxes));
    tau = (tau(mgIdxes));
    
    % initial array of all BETA, TAU and W values;
    if (lambdaIdx==1)
        names = cell(1,length(L));
        BETA = zeros(lambdaCount,length(beta));
        TAU = zeros(lambdaCount,length(tau));
        W = zeros(1,lambdaCount);
        
        for i=1:length(L)
            names(i) = cellstr(sprintf('%1.1i,%1.1i',L(i),M(i)));
        end
    end
    
    BETA(lambdaIdx,:) = beta;
    TAU(lambdaIdx,:) = tau;
    W(lambdaIdx) = w;
    
    % If it's the centre wavelength, print out some info, particularly
    % about the fundamental mode
    if (lambdaIdx==lambda0Idx)
        
        figure(2);
        plot(r_um,n_r);
        title('Refractive Index Used');
        xlim([0 2.*core_r].*1e6);
        ylabel('n_r');
        xlabel('\mum');
        grid on;
        
        beta_MAX = beta_max;
        beta_MIN = beta_min;
        
        E_r = Er;
        %Fundamental mode properties
        LP01 = Er(1,:);
        LP01_MAX = max(max(abs(LP01)));
        
        LP01_e2 = LP01_MAX/exp(1);
        tempL = length(Er(1,:));
        MFD = 0;
        for idx=1:tempL-1
            if (abs(Er(1,idx))>LP01_e2 && abs(Er(1,idx+1))<LP01_e2)
                MFD = (r(idx)+r(idx+1))/2*10^6;
            end
        end
        %Peterman II
        dr = r(2)-r(1);
        dLP01_dr = gradient(LP01)./dr;
        LP01_INTENSITY = abs(LP01).^2;
        numeratorr = sum(LP01_INTENSITY.*r.*dr);
        denominatorr = sum((dLP01_dr.^2).*r.*dr);
        MFD_Petermann2 = 2.*sqrt(2).*(numeratorr./denominatorr).^0.5;
        
        
        NA = sqrt(n_core.^2-n_clad.^2);
        theta_c = acos(n_clad./n_core);
        
        divergence = (4.*lambda)./((MFD.*10^-6).*2.*pi);
        
        fprintf('Numerical Aperture (refractive index)  : %10.10f\n',NA);
        fprintf('Numerical Aperture (divergence of LP01): %10.10f\n',sin(divergence/2));
        fprintf('Complement of critical angle (rad): %10.10f\n',theta_c);
        
        fprintf('Divergence %10.10f (deg)\n',divergence.*180./pi);
        
        fprintf('Core Diameter (um)      : %10.10f\n',2*core_r*1e6);
        fprintf('MFD (1/e)  (um): %10.10f\n',2*MFD)
        fprintf('MFD (P2)   (um): %10.10f\n',MFD_Petermann2*1e6)
    end
end

% Propagation constant
figure(1);
subplot(2,2,1);
plot(1:length(beta),BETA(lambda0Idx,:),'-x');
xlabel('LP_l_,_m');
ylabel('\beta (1/m)');
set(gca,'XTick',1:length(beta));
set(gca,'XTickLabel',names);
xtickangle(45)
grid on;
title('Propagation Constants');

%figure;

if (lambdaCount>=3)
    
    % Differential Mode Delay (derivative of beta)
    dBeta = zeros(lambdaCount,totalModes);
    for modeIdx=1:totalModes
        dBeta(:,modeIdx) = gradient(BETA(:,modeIdx));
    end
    dOmega = gradient(W);
    
    for lambdaIdx=1:lambdaCount
        TAU(lambdaIdx,:) = dBeta(lambdaIdx,:)./dOmega(lambdaIdx);
    end
    tau = TAU(lambda0Idx,:);
    
    subplot(2,2,2);
    plot(1:length(tau),(tau-tau(1)).*1e12,'-x');
    xlabel('LP_l_,_m');
    ylabel('\delta\tau (ps/m)');
    set(gca,'XTick',1:length(tau));
    set(gca,'XTickLabel',names);
    title(sprintf('LP_0_1 v_g = %3.3fc\n',(1./tau(1))./c));
    grid on;
    title('Differential Mode Delay');
    
    % Chromatic dispersion of each mode (derivative of the group delay)
    if (lambdaCount>=5)
        dTAU = zeros(lambdaCount,totalModes);
        dLambda_nm = gradient(lambdas.*1e9);
        
        for modeIdx=1:totalModes
            dTAU(:,modeIdx) = gradient(TAU(:,modeIdx).*1e12);
        end
        
        for lambdaIdx=1:lambdaCount
            dTAU(lambdaIdx,:) = dTAU(lambdaIdx,:)./dLambda_nm(lambdaIdx);
        end
        
        dtau = dTAU(lambda0Idx,:);
        
        subplot(2,2,3);
        plot(1:length(dtau),(dtau).*1000,'-*');
        xlabel('LP_l_,_m');
        ylabel('\delta\tau (ps/km/nm)');
        set(gca,'XTick',1:length(dtau));
        set(gca,'XTickLabel',names);
        title('Chromatic Dispersion');
        grid on;
        
    end
end

% subplot(2,2,3);
% plot(r_um,abs(LP01)./LP01_MAX);
% hold on;
% plot([MFD MFD],[0 1],'-r');
% plot([1 1].*(1e6.*MFD_Petermann2./2),[0 1],'-black');
% hold off;
% xlim([0 3*MFD]);

subplot(2,2,4);
plot(r_um,E_r);
xlim([0 2*core_r]*1e6);
xlabel('\mum');
ylabel('Er(r)');
ylim(max(max(abs(Er))).*1.1.*[-1 1]);
grid on;
title('Mode Profiles');