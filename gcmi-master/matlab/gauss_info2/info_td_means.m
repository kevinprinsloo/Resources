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
    error('info_gd: number of trials do not match');
end

%TODO: default values for options (or put in wrapper
% convenience interface)

% build t-dist mixture model
% fix common S and nu
[grp_mu S nu] = fitt_commonsnu(x,y,Ym);
% w = zeros(1,Ym);
% Ntrl_y = zeros(Ym,1);
% Hcond = zeros(1,Ym);

% for yi=1:Ym
%     idx = y==(yi-1);
%     thsNtrl = sum(idx);
%     Ntrl_y(yi) = thsNtrl;
%     w(yi) = thsNtrl ./ Ntrl;
    
%     tmp = bsxfun(@minus, tmp, sum(tmp, 1)/thsNtrl);
%     x_group(idx,:) = tmp;

    
%     chC = chol(Pd.sigma);
%     p2 = Nvar / 2;
%     nu = Pd.nu;
%     nup2 = (nu+Nvar)/2;
%     nu2 = nu/2;
%     Hcond(yi) = sum(log(diag(chC))) ...
%                + log( ((nu*pi).^p2) * beta(p2, nu2) ) - gammaln(p2) ...
%                + nup2*(psi(nup2)-psi(nu2));
% end

% conditional casee
% Pd = fitdist(x_group, 'tlocationscale');
% chC = chol(Pd.sigma);
% p2 = Nvar / 2;
% nu = Pd.nu;
% nup2 = (nu+Nvar)/2;
% nu2 = nu/2;
% Hcond = sum(log(diag(chC))) ...
%                + log( ((nu*pi).^p2) * beta(p2, nu2) ) - gammaln(p2) ...
%                + nup2*(psi(nup2)-psi(nu2));
% Hcond_nu = log( ((nu*pi).^p2) * beta(p2, nu2) ) - gammaln(p2) ...
%                + nup2*(psi(nup2)-psi(nu2));

% [mu, S, n] = fitt(x_group);
chC = chol(S);
p2 = Nvar / 2;
nup2 = (nu+Nvar)/2;
nu2 = nu/2;
Hnu = log( ((nu*pi).^p2) * beta(p2, nu2) ) - gammaln(p2) ...
               + nup2*(psi(nup2)-psi(nu2));
% entropy the same for each group
Hcond = sum(log(diag(chC))) + Hnu;
           
% single dist for unconditional P(x)
% Pd = fitdist(x, 'tlocationscale');
% chC = chol(Pd.sigma);
% p2 = Nvar / 2;
% nu = Pd.nu;
% nup2 = (nu+Nvar)/2;
% nu2 = nu/2;
% Hunc = sum(log(diag(chC))) ...
%                + log( ((nu*pi).^p2) * beta(p2, nu2) ) - gammaln(p2) ...
%                + nup2*(psi(nup2)-psi(nu2));
% Hunc_nu = log( ((nu*pi).^p2) * beta(p2, nu2) ) - gammaln(p2) ...
%                + nup2*(psi(nup2)-psi(nu2));

% single dist for unconditional P(x) - with same nu from conditional
[c, S] = fitt_fixnu(x,nu);
chC = chol(S);
Hunc = sum(log(diag(chC))) + Hnu;


% apply bias corrections
ln2 = log(2);
% psiterms = psi((Ntrl - (1:Nvar))/2) / 2;
% dterm = (ln2 - log(Ntrl-1)) / 2;
% Hunc = (Hunc - Nvar*dterm - sum(psiterms));

% Hcond = Hcond - gaussian_bias(Ntrl_y, Nvar)';
I = Hunc - Hcond;
I = I / ln2;

