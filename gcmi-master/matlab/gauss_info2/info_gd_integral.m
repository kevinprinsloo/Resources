function [I, LLR] = info_gd_integral(x, y, Ym, demeaned, dollr, mcN)
% out = info_gg(x, y, Ym, demeaned, dollr)
% Calculate Gaussian - Discrete information
% x - gaussian input (Ntrl x Nvar)
% y - discrete integer input (Ntrl in [0, Ym-1]
% demeaned = 1 => x, y already have zero mean
% dollr = 1 => also calculate likelihood ratio test statistic
% mcN => number of samples for monte-carlo integration of unconditional
% mixture entropy 
% mcN=0 => use data as MC samples [likelihood]
if isvector(x)
    x = x(:);
end
if isvector(y)
    y = y(:);
else
    error('info_gd: only univariate discrete variable supported');
end

Ntrl = size(x,1);
Nvar = size(x,2);

if size(y,1) ~= Ntrl
    error('info_gd: number of trials do not match');
end

%TODO: default values for options (or put in wrapper
% convenience interface)

% build gaussian mixture model
gmm.G = Ym;
gmm.N = Nvar;
gmm.C = zeros(Nvar, Nvar, Ym);
gmm.chC = zeros(Nvar, Nvar, Ym);
gmm.m = zeros(Nvar, Ym);
Ntrl_y = zeros(1,Ym);

% TODO: move this loop to fortran single loop to 
% avoid having to extract conditional repsonses
Hcond = zeros(1,Ym);
LLRcond = zeros(1,Ym);

c = 0.5*log(2*pi) + 0.5;
for yi=1:Ym
    idx = y==(yi-1);
    tmp = x(idx,:);
    thsNtrl = size(tmp,1);
    gmm.w(yi) = thsNtrl ./ Ntrl;
    gmm.m(:,yi) = sum(tmp, 1)/thsNtrl;
    tmp = bsxfun(@minus, tmp, gmm.m(:,yi)');
    cx = (tmp'*tmp) / (thsNtrl - 1);
    chC = chol(cx);
    gmm.C(:,:,yi) = cx;
    gmm.chC(:,:,yi) = chC;
    Hcond(yi) = sum(log(diag(chC))) + c*Nvar;
    if dollr
        LLRcond(yi) = sum(norm_innerv(tmp, chC)) ... 
                        - thsNtrl*(Hcond(yi) - 0.5*Nvar - log(gmm.w(yi)));
    end
end

% single gaussian model for LLR
if dollr
    if demeaned
        gunc.m = zeros(1, Nvar);
        tmp = x;
    else
        gunc.m = sum(x,1)/Ntrl;
        tmp = bsxfun(@minus,x,gunc.m);
    end
    gunc.cx = (tmp'*tmp) / (Ntrl-1);
    chC = chol(gunc.cx);
    LLRunc = sum(norm_innerv(tmp, chC)) ...
        - Ntrl*(sum(log(diag(chC))) + (c-0.5)*Nvar);% ...
        %+ sum(log(gmm.w(y+1)));
    LLRunc = LLRunc + Ntrl*sum(gmm.w.*log(gmm.w));
    LLR = 2*(sum(LLRcond) - LLRunc);
end

% MC numerical integral of gmm model
if mcN == 0
    mcN = Ntrl;
    mcsamp = x;
else
    % sample!!norm_innerv(dx', gmm.chC(:,:,i))
    mcsamp = GMMsample(gmm, mcN)';
end

% calculate gmm likelihood of data
lik = zeros(mcN,1);
for i=1:gmm.G
    dx = bsxfun(@minus, mcsamp, gmm.m(:,i));
    % log-likelihood for this mixture
    ll = norm_innerv(dx, gmm.chC(:,:,i)) ...
            - (Hcond(i) - 0.5*Nvar);
    lik = lik + exp(log(gmm.w(i))+ll);
end
Hunc = -sum(log(lik)) / mcN;

I = Hunc - sum(gmm.w .* Hcond);

I = I / log(2);

%
%

function w = norm_innerv(x, chC)
% sum of normalised innovations
% m = (chC')\x;
m = x/chC;
w = -0.5 * sum(m.*m,2);

function s = GMMsample(gmm, N)
% GMMSAMPLE Generate samples from a Gaussian mixture.
% IN        gmm     The mixture to sample from.
%           N       Number of samples. (optional)
% RET       s       Samples. D-by-N matrix, where D is dimesion of <gmm>
%
% AUTHOR    Marco Huber, 2008, marco.huber@ieee.org

D = gmm.N;

% Determine <N> uniformly distributed numbers in [0;1] to sample from each
% mixture component
prob = rand(N, 1);% mcN=-1 => use single gaussian instead of mixture integral
s = zeros(D, N);

% Cumulates the weights (probabilities) of the components
sum_prob = 0;
% Inceases for each drawn sample.
sampled = 0;
for i = 1:gmm.G
    select = sum( (prob >= sum_prob) & (prob < sum_prob + gmm.w(i)) );
    if sampled+select ~= sampled
        s(:,sampled+1:sampled+select) = ...
            randn_N(gmm.m(:,i), gmm.chC(:,:,i), select, D);
    end
    
    sum_prob = sum_prob + gmm.w(i);
    sampled  = sampled + select;
end

function r = randn_N(m, chC, N, D)
%[eve, eva] = eig(C);
%r = m*ones(1, N) + eve*sqrt(eva)*randn(D, N);
r = m*ones(1, N) + chC'*randn(D,N);

