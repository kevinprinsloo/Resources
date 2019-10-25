function [I, LLR] = cmi_gdg(x, y, Ym, z, biascorrect, demeaned, dollr)
% out = info_gg(x, y, Ym, demeaned, dollr)
% Calculate Gaussian - Discrete CMI conditioned on Gaussian
% x - gaussian input (Ntrl x Nvar)
% y - discrete integer input (Ntrl in [0, Ym-1]
% z - gaussian input (Ntrl x Nvar)
% biascorrect => apply analytic gaussian bias correction
% demeaned = 1 => x, y already have zero mean
% dollr = 1 => also calculate likelihood ratio test statistic
% mcN => number of samples for monte-carlo integration of unconditional
% mixture entropy 
% mcN=0 => use data as MC samples [likelihood]
if isvector(x)
    x = x(:);
end
if isvector(z)
    z = z(:);
end
if isvector(y)
    y = y(:);
else
    error('info_gd: only univariate discrete variable supported');
end

Ntrl = size(x,1);
Nvarx = size(x,2);
Nvarz = size(z,2);
Nvar = Nvarx + Nvarz;
xz = [x z];

if size(y,1) ~= Ntrl || size(z,1) ~= Ntrl
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
HXZcond = zeros(1,Ym);
HZcond = zeros(1,Ym);
LLRcond = zeros(1,Ym);

c = 0.5*log(2*pi) + 0.5;
zidx = (Nvarx + 1):Nvar;
for yi=1:Ym
    idx = y==(yi-1);
    tmp = xz(idx,:);
    thsNtrl = size(tmp,1);
    Ntrl_y(yi) = thsNtrl;
    gmm.w(yi) = thsNtrl ./ Ntrl;
    gmm.m(:,yi) = sum(tmp, 1)/thsNtrl;
    tmp = bsxfun(@minus, tmp, gmm.m(:,yi)');
    cxz = (tmp'*tmp) / (thsNtrl - 1);
    chC = chol(cxz);
    HXZcond(yi) = sum(log(diag(chC))) + c*Nvar;
    cz = cxz(zidx,zidx);
    chC = chol(cz);
    HZcond(yi) = sum(log(diag(chC))) + c*Nvarz;
    if dollr
%         LLRcond(yi) = sum(norm_innerv(tmp, chC)) ... 
%                         - thsNtrl*(Hcond(yi) - 0.5*Nvar - log(gmm.w(yi)));
    end
end

% single gaussian for unconditional P(x)
if demeaned
    gunc.m = zeros(1, Nvar);
    tmp = xz;
else
    gunc.m = sum(xz,1)/Ntrl;
    tmp = bsxfun(@minus,xz,gunc.m);
end
gunc.cxz = (tmp'*tmp) / (Ntrl-1);
chC = chol(gunc.cxz);
Hunc = sum(log(diag(chC))) + 0.5*Nvar*log(2*pi*exp(1));

% HZ
gz.C = gunc.cxz(zidx,zidx);
gz.chC = chol(gz.C);

HZ = sum(log(diag(gz.chC))) + 0.5*Nvarz*log(2*pi*exp(1));



% likelihood ratio test
if dollr
%     LLRunc = sum(norm_innerv(tmp, chC)) ...
%         - Ntrl*(Hunc - 0.5*Nvar);% ...
%         %+ sum(log(gmm.w(y+1)));
%     LLRunc = LLRunc + Ntrl*sum(gmm.w.*log(gmm.w));
%     LLR = 2*(sum(LLRcond) - LLRunc);
%     LLRZ = sum_norm_innerv(z, gz.chC) - Ntrl*HZ; 
end

% apply bias corrections
ln2 = log(2);
if biascorrect
    vars = 1:Nvar;
    
    % unconditional entropies
    psiterms = psi((Ntrl - vars)/2) / 2;
    dterm = (ln2 - log(Ntrl-1)) / 2;
    Hunc = (Hunc - Nvar*dterm - sum(psiterms));
    HZ = (HZ - Nvarz*dterm - sum(psiterms(1:Nvarz)));
    
    % Y-conditional entropies
    dterm = (ln2 - log(Ntrl_y-1)) / 2;
    psiterms = zeros(1,Ym);
    for vi=1:Nvarz
        psiterms = psiterms + psi((Ntrl_y-vi)/2);
    end
    HZcond = HZcond - Nvarz*dterm - psiterms/2;
    for vi=(Nvarz+1):Nvar
        psiterms = psiterms + psi((Ntrl_y-vi)/2);
    end
    HXZcond = HXZcond - Nvar*dterm - psiterms/2;
end


I = Hunc + sum(gmm.w .* HZcond) - sum(gmm.w .* HXZcond) - HZ;
I = I / ln2;

%
%

function w = norm_innerv(x, chC)
% sum of normalised innovations
% m = (chC')\x;
m = x/chC;
w = -0.5 * sum(m.*m,2);
