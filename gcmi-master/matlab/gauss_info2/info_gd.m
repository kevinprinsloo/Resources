function [I, LLR] = info_gd(x, y, Ym, biascorrect, demeaned, dollr)
% out = info_gg(x, y, Ym, demeaned, dollr)
% Calculate Gaussian - Discrete information
% x - gaussian input (Ntrl x Nvar)
% y - discrete integer input (Ntrl in [0, Ym-1]
% biascorrect => apply analytic gaussian bias correction
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
    Ntrl_y(yi) = thsNtrl;
    gmm.w(yi) = thsNtrl ./ Ntrl;
    gmm.m(:,yi) = sum(tmp, 1)/thsNtrl;
    tmp = bsxfun(@minus, tmp, gmm.m(:,yi)');
    cx = (tmp'*tmp) / (thsNtrl - 1);
    chC = mychol(cx);
    gmm.C(:,:,yi) = cx;
    gmm.chC(:,:,yi) = chC;
    Hcond(yi) = sum(log(diag(chC)));% + c*Nvar;
    if dollr
        LLRcond(yi) = sum(norm_innerv(tmp, chC)) ... 
                        - thsNtrl*(Hcond(yi) - 0.5*Nvar - log(gmm.w(yi)));
    end
end

% single gaussian for unconditional P(x)
if demeaned
    gunc.m = zeros(1, Nvar);
    tmp = x;
else
    gunc.m = sum(x,1)/Ntrl;
    tmp = bsxfun(@minus,x,gunc.m);
end
gunc.cx = (tmp'*tmp) / (Ntrl-1);
chC = mychol(gunc.cx);
Hunc = sum(log(diag(chC)));% + c*Nvar;

if dollr
    LLRunc = sum(norm_innerv(tmp, chC)) ...
        - Ntrl*(Hunc - 0.5*Nvar);% ...
        %+ sum(log(gmm.w(y+1)));
    LLRunc = LLRunc + Ntrl*sum(gmm.w.*log(gmm.w));
    LLR = 2*(sum(LLRcond) - LLRunc);
end

% apply bias corrections
ln2 = log(2);
if biascorrect
    vars = 1:Nvar;
    
    psiterms = psi((Ntrl - vars)/2) / 2;
    dterm = (ln2 - log(Ntrl-1)) / 2;
    Hunc = (Hunc - Nvar*dterm - sum(psiterms));
    
    dterm = (ln2 - log(Ntrl_y-1)) / 2;
    psiterms = zeros(1,Ym);
    for vi=vars
        psiterms = psiterms + psi((Ntrl_y-vi)/2);
    end
    Hcond = Hcond - Nvar*dterm - (psiterms/2);
end

% Hunc
% Hcond
I = Hunc - sum(gmm.w .* Hcond);
I = I / ln2;

%
%

function w = norm_innerv(x, chC)
% sum of normalised innovations
% m = (chC')\x;
m = x/chC;
w = -0.5 * sum(m.*m,2);
