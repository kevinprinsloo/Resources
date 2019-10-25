function [I, LLR] = info_td(x, y, Ym, demeaned, dollr)
% out = info_td(x, y, Ym, demeaned, dollr)
% Calculate t-distri - Discrete information
% x - t-distributed input (Ntrl x Nvar)
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
    error('info_td: number of trials do not match');
end

%TODO: default values for options (or put in wrapper
% convenience interface)

% fit t-dist mixture model
[grp_mu grp_S nu] = fitt_commonnu(x,y,Ym);

w = zeros(1,Ym);
Ntrl_y = zeros(Ym,1);
Hcond = zeros(1,Ym);

p2 = Nvar / 2;
nup2 = (nu+Nvar)/2;
nu2 = nu/2;
Hnu = log( ((nu*pi).^p2) * beta(p2, nu2) ) - gammaln(p2) ...
               + nup2*(psi(nup2)-psi(nu2));
for yi=1:Ym
    thsNtrl = sum(y==(yi-1));
    Ntrl_y(yi) = thsNtrl;
    w(yi) = thsNtrl ./ Ntrl;
    
    chC = chol(grp_S(:,:,yi));
    Hcond(yi) = sum(log(diag(chC))) + Hnu;
end

% single dist for unconditional P(x)
% use nu from conditional fits
[mu S] = fitt_fixnu(x, nu);
chC = chol(S);
Hunc = sum(log(diag(chC))) + Hnu;

% apply bias corrections
ln2 = log(2);
% psiterms = psi((Ntrl - (1:Nvar))/2) / 2;
% dterm = (ln2 - log(Ntrl-1)) / 2;
% Hunc = (Hunc - Nvar*dterm - sum(psiterms));

% Hcond = Hcond - gaussian_bias(Ntrl_y, Nvar)';
I = Hunc - sum(w .* Hcond);
Hunc;
Hcond;
I = I / ln2;

%
%

function bias = gaussian_bias(Nt,L)
% GAUSSIAN_BIAS Computes bias of gaussian entropy estimates
%   Copyright (C) 2010 Cesare Magri
%   Version: 1.0.0
Nt = Nt-1;
NtMat = Nt(:, ones(L,1));
pvec = 0:L-1;
NtMat = NtMat - pvec(ones(length(Nt),1), :);
bias = (sum(psi(NtMat./2), 2) - L.*log(Nt(:,1)./2)) / 2;
