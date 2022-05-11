function [beta, R, tau, A, r] = RadialModeSolver2(n_r, r, L, lambda, boundOnly)

% compute solution to radial scalar wave equation
% D^2 R(r) + 1/r D R(r) + [ (k0 * n(r))^2 - (L/r)^2 ] R(r) = beta^2 R(r)
%
% with the following approximations
%
%          R(i+1) - R(i-1)
% D R(r) = ---------------
%          r(i+1) - r(i-1)
%
%
%
%             R(i+1) - R(i)      R(i) - R(i-1)
%             ------------   -   -------------
%             r(i+1) - r(i)      r(i) - r(i-1)
% D^2 R(r) = ---------------------------------
%                    r(i+1) - r(i-1)
%                    ---------------
%                           2
%
% thus turning into eigenvalue problem
%
% inputs:
% fiber -> fiber parameters (see function stored in fiber.function)
% fiber.function -> string storing function used for evaluation of
% fiber index profile (required)
%
% dr max -> maximum allowed step size for radial mesh
% L -> azimuthal mode number: sin(L*phi); cos(L*phi)
% lambda -> freespace wavelength


options.disp = 0;
options.maxit = 30000;
options.isreal = false;
c = 299792458;

% define delta's given that r can be arbitrary step sizes
% delta_0 -> r(i) - r(i-1)
% delta_1 -> r(i+1) - r(i)
% delta_2 -> r(i+1) - r(i-1) = delta 0 + delta 1

delta_0 = r(2:end) - r(1:end-1);
delta_1 = r(3:end) - r(2:end-1);
delta_1(end+1) = delta_1(end); % fix end point of delta 1
delta_2 = delta_0 + delta_1;

% remove r=0 point for subsequent computations
r(1) = [];
n_r(1) = [];

% set slope of R at r=0
dR0 = mod(L,2);

% compute freespace wave vector
k0 = 2*pi/lambda;
N = length(r);

% compute diagonals of scalar wave equation matrix
%
% [diagC digaR 0 ] * R = beta * R
% [diagL diagC diagR]
% [ 0 diagL diagC]
%
% i.e. R(i+1) * diagR(i) + R(i) diagC(i) + R(i-1) * diagL(i) = beta * R(i)

K2 = (k0 * n_r).^2 - (L./r).^2;

diagC = -2./(delta_1 .* delta_0) + K2; % center diagonal
diagL = (2./delta_0 - 1./r) ./ delta_2; % left diagonal
diagR = (2./delta_1 + 1./r) ./ delta_2; % right diagonal
% for end points, the boundary value must be folded into elements of matrix
% for i=0 (i.e. r(i=0) = 0) the slope is either 0 or 1
% for slope of 1, R(i) = 0
% for slope of 0, R(i) = R(i+1) = R(i-1)
% thus eigen eq. is
% R(i+1) * diagR(i) + R(i) diagC(i) + R(i) * diagL(i) = beta * R(i)
% so for i=1, the diagL(i) term is folded into diagC(i) if R(i) can be non-zero

diagC(1) = diagC(1) + ~dR0 * diagL(1);

% similarly for i=end (last point), the boundary value must be folded into diagC
% however, the slope of R(i) is forced to be 0 and R(i) = 0
% thus eigen eq. is
% 0 * diagR(i) + R(i) diagC(i) + R(i) * diagL(i) = beta * R(i)
% and diagR(end) folds into diagC(end) trivially

diagC(end) = diagC(end) + 0 * diagR(end);

% create sparse matrix with diagonals
A = spdiags(diagC',0,N,N);
A = spdiags(diagL',1,A); % insert into right side but later transpose matrix
A = spdiags(diagR',-1,A); % insert into left side but later transpose matrix

A = A'; % transpose matrix

% estimate number of solutions to wave equation
beta_min = k0*n_r(end); % compute minimum possible beta
kt = sqrt(K2-beta_min^2); % compute maximum tranverse wave number
pos = find((r>0) & (kt.^2 > 0)); % find where kt is real
dr = diff(r);
%pos = pos(1:end-1);
% integrate kt to estimate # of nodes in highest order mode (i.e. mode number)
% integrate only where kt^2 > 0 (and in r>0 domain)
%size(kt)
%size(dr)
%pos
M = ceil( kt(pos) * dr(pos)' / pi);
M = M+2; % padM to allow for error in estimate

% find M eigenvalues (beta^2) closest to maximum possible of beta^2
nn = k0.*n_r(2);
%[R,beta2] = eigs(A,M,max(k0*n_r)^2,options);
[R,beta2] = eigs(A,M,nn.^2,options);
%[R,beta2] = eigs(A,M*60,'LM',options);
%beta = diag(sqrt(real(beta2)));
beta = diag(sqrt((beta2)));

%[R,beta2] = eigs(A,M,'LM',options);
%beta = diag(sqrt(real(beta2)));

% sort for unique beta (remove degenerate values) and keep guiding modes (beta > k min)
if (boundOnly==1) 
    idx = find(beta > beta_min);
	beta = beta(idx);
	R = R(:,idx);
end
if isempty(beta)
    R = [];
    tau = [];
    return
end

% [dmp,idx] = max(abs(R),[],1);
% scale = diag(R(idx,1:length(idx)))';
% scale = repmat(scale,[size(R,1),1]);
% R = R ./ scale;

% replace smallest value of r with 0
% r(1) = 0;
% if (L~=0)
%     R(1,:) = 0;
% end
[beta,IDX] = sort(beta,1,'descend');

R = R(:,IDX);
tau = zeros(1,length(beta));

for j=1:length(beta)
    tau(j) = (((c.*beta(j)./k0).*sum(trapz(r,r.*R(:,j)'.^2))./sum(trapz(r,r.*(n_r.*R(:,j)').^2))).^-1);
end

% optionally plot eigenfunction
if (0)
    clf; hold on;
    iidx = round(linspace(1,length(r),101));
    plot(r(iidx)/lambda,R(iidx,:),'k.-');
    pause(0.1);
end