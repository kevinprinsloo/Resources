function [I, LLR] = info_gd(x, y, Ym, demeaned, dollr)
% TEST FIXING PRECISION ACROSS CONDITIONS TO MAKE IT JUST A TEST 
% FOR DIFFERENCE IN MEAN
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
Ntrl_y = zeros(Ym,1);

% TODO: move this loop to fortran single loop to 
% avoid having to extract conditional repsonses
Hcond = zeros(1,Ym);
LLRcond = zeros(1,Ym);

c = 0.5*log(2*pi) + 0.5;
x_group = zeros(size(x)); % x group demeaned
for yi=1:Ym
    idx = y==(yi-1);
    tmp = x(idx,:);
    thsNtrl = size(tmp,1);
    Ntrl_y(yi) = thsNtrl;
    gmm.w(yi) = thsNtrl ./ Ntrl;
    gmm.m(:,yi) = sum(tmp, 1)/thsNtrl;
    tmp = bsxfun(@minus, tmp, gmm.m(:,yi)');
    x_group(idx,:) = tmp;
%     cx = (tmp'*tmp) / (thsNtrl - 1);
%     chC = chol(cx);
%     gmm.C(:,:,yi) = cx;
%     gmm.chC(:,:,yi) = chC;
%     Hcond(yi) = sum(log(diag(chC))) + c*Nvar;
%     if dollr
%         LLRcond(yi) = sum(norm_innerv(tmp, chC)) ... 
%                         - thsNtrl*(Hcond(yi) - 0.5*Nvar - log(gmm.w(yi)));
%     end
end

% single gaussian for unconditional P(x)
if demeaned
    gunc.m = zeros(1, Nvar);
    tmp = x;
else
    gunc.m = sum(x,1)/Ntrl;
    tmp = bsxfun(@minus,x,gunc.m);
end
tmp = bsxfun(@minus,x,gunc.m);
gunc.cx = (tmp'*tmp) / (Ntrl-1);
% gunc.cx = 1.253*mad(tmp,0);
chC = chol(gunc.cx);
Hunc = sum(log(diag(chC))) + c*Nvar;

% ?? dof = Ym, since 1 mean from each group
% ^^ leads to negative offset
gcon.cx = (x_group'*x_group) / (Ntrl-1);
% gcon.cx = 1.253*mad(x_group,0);
chC = chol(gcon.cx);
Hcond = sum(log(diag(chC))) + c*Nvar;

ln2 = log(2);

I = Hunc - Hcond;
I = I / ln2;

%
%

function w = norm_innerv(x, chC)
% sum of normalised innovations
% m = (chC')\x;
m = x/chC;
w = -0.5 * sum(m.*m,2);

function bias = gaussian_bias(Nt,L)
% GAUSSIAN_BIAS Computes bias of gaussian entropy estimates
%   Copyright (C) 2010 Cesare Magri
%   Version: 1.0.0
Nt = Nt-1;
NtMat = Nt(:, ones(L,1));
pvec = 0:L-1;
NtMat = NtMat - pvec(ones(length(Nt),1), :);
bias = (sum(psi(NtMat./2), 2) - L.*log(Nt(:,1)./2)) / 2;
