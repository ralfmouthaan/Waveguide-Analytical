function n=alpha_nSell3(r,fiber)

% index profile
% inputs r -> radial positions for which n is desired
% fibers -> parameters of fiber
% fiber.n1 -> core index
% fiber.alpha -> fiber alpha
% fiber.a -> fiber core radius
% fiber.delta -> fiber contrast
% i.e. delta = ~(n core-n clad)/n core
% fiber.function -> ’alpha n’
%
% n^2 = nc^2*(1-2*delta*(r/a)^2

% use default params if none is specified (use 50-um Corning MMF)
if (nargin<2)
    fiber.n1 = 1.47;
    fiber.alpha = 2;
    fiber.a = 25e-6;
    fiber.delta = 0.012;
    fiber.y = 0.1;
    fiber.length = 1000;
    fiber.function = 'alpha_n';
end
% if no arguments, return template for fiber (use 50-um Corning MMF)
if (nargin==0)
    n=fiber;
    return
end
    a = fiber.a;
    alpha = fiber.alpha;
    if (length(alpha)==1)
        alphas = ones(1,length(r)).*alpha;
    else
        %if (length(alpha)==2)
        %    alphas = alpha(1)+(linspace(0,1,length(r))).*(alpha(2)-alpha(1));
       % else
            %x = linspace(0,fiber.a,length(alpha));
            x = fiber.xpos;
            y = alpha;
            
          %  x
         %   y
           % fo = fit(x',y','cubicinterp');
            fo = fit(x',y','linearinterp');
          %  fo = fit(x',y','poly1');
            alphas = fo(r);
            figure(123);
            plot(r.*1e6,alphas);
            xlim([0 fiber.a].*1e6);
       % end
    end
    doping = fiber.doping;
    lambda = fiber.lambda;
    cladDefect = fiber.cladDefect(1);
    if (length(fiber.cladDefect)>1)
        coreCladMismatch = fiber.cladDefect(2);
    else
        coreCladMismatch = 0;
    end
    defect = cladDefect;
    [minV minIdx] = min(abs((a-cladDefect/2)-r));
    
    cladDefectStart = r(minIdx);
    [minV maxIdx] = min(abs((a+cladDefect/2)-r));
    cladDefectStop = r(maxIdx);
    %'D:\Documents\PhD\Papers\Joel\My Papers\CLEO PacificRim\PD\Data\06cRESULT.mat'
    n = ones(size(r));
if (fiber.alpha(1) < 1000)
    r = abs(r);
    alpha = fiber.alpha;
    if (alpha>0)   
        for i=1:length(r)
            n(i) = Sellmeier2(lambda,doping.*(1-(abs(r(i))./(a+coreCladMismatch)).^alphas(i)));
        end
        n(r>a) = Sellmeier2(lambda,0);
        if (abs(defect)>0)
        defectStartDoping = doping.*(1-(abs(r(minIdx))./(a+coreCladMismatch)).^alphas(minIdx));
        defectStopDoping = 0;
        
        for i=minIdx:length(r)
            exponent = -2.*((r(i)-r(minIdx))./(r(maxIdx)-r(minIdx))).^1;
            n(i)=Sellmeier2(lambda,defectStartDoping+((defectStopDoping-defectStartDoping).*(1-exp(exponent))));
         %  n(i)=Sellmeier2(lambda,0);
        end
        end
    else
      %  disp('2');
        alpha = -alpha;
        %r = min(r,fiber.a); % for r larger than core radius,
    % fix it to a so that proper value of n is returned
        for i=1:length(r)
            n(i) = Sellmeier2(lambda,doping.*(1-(abs(r(i))./a).^alphas(i)));
        end
    end
else
    n1 = Sellmeier2(lambda,0);
     n2 = Sellmeier2(lambda,doping);
    n = ones(size(r)).*n2;
  % disp('3');
    idx = find(abs(r) >= fiber.a);
   
    n(idx) = n1;
    %n(idx) = fiber.n1 * sqrt(abs(1-2*fiber.delta));
    
    %for i=1:length(r)
    %    n(i) = Sellmeier2(lambda,doping.*(1-(abs(r(i))./a).^alpha));
    %end
end

