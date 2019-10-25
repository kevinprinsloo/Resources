function [I, LLR] = info_gg(x, y, biascorrect, demeaned, dollr)
% out = info_gg(x, y, demeaned, llr)
% Calculate Gaussian - Gaussian information
% demeaned = 1 => x, y already have zero mean
% llr = 1 => also calculate likelihood ratio test statistic
if isvector(x)
    x = x(:);
end
if isvector(y)
    y = y(:);
end
Ntrl = size(x,1);
Nvarx = size(x,2);
Nvary = size(y,2);

if size(y,1) ~= Ntrl
    error('info_gg: number of trials do not match')
end

%TODO: default values for options (or put in wrapper
% convenience interface)
if demeaned
    gx.m = zeros(1, Nvarx);
    gy.m = zeros(1, Nvary);
    gxy.m = zeros(1, Nvarx+Nvary);
else
    gx.m = sum(x,1)/Ntrl;
    gy.m = sum(y,1)/Ntrl;
    gxy.m = [gx.m gy.m];

    x = bsxfun(@minus,x,gx.m);
    y = bsxfun(@minus,y,gy.m);
end

% concatenating actually seems relatively expensive which
% is why I dont split into fit_gg llr_gg info_gg functions
% (also don't want to have to demean twice, or stored duplicate
% demeaned data
% TODO: directly calculate covariance without requiring
% concatenation
xy = [x y];

% do covariance once
C = (xy'*xy) / (Ntrl - 1);
gxy.C = C;
% extract right subregions
gx.C = C(1:Nvarx,1:Nvarx);
ystart = Nvarx + 1;
Nvarxy = Nvarx + Nvary;
gy.C = C(ystart:Nvarxy,ystart:Nvarxy);

% chol needed for both entropy calculation and likelihood
gx.chC = chol(gx.C);
gy.chC = chol(gy.C);
gxy.chC = chol(gxy.C);
% TODO: some kind of blocking here to remove redunancy
% gx.chC = gxy.chC(1:Nvarx,1:Nvarx);
% gy.chC = gxy.chC(ystart:Nvarxy,ystart:Nvarxy);

% normalisations cancel for information
% raw nats entropy so it can be used for likelihood
% normalisations not needed since they cancel 
% (both for entropy and likelihodd)
HX = sum(log(diag(gx.chC))); % + 0.5*Nvarx*log(2*pi*exp(1));
HY = sum(log(diag(gy.chC))); % + 0.5*Nvary*log(2*pi*exp(1));
HXY = sum(log(diag(gxy.chC))); % + 0.5*(Nvarx+Nvary)*log(2*pi*exp(1));

LLR = [];
% likelihood ratio test
if dollr
    LLRX = sum_norm_innerv(x, gx.chC) - Ntrl*HX; 
    LLRY = sum_norm_innerv(y, gy.chC) - Ntrl*HY; 
    LLRXY = sum_norm_innerv(xy, gxy.chC) - Ntrl*HXY; 
    LLR = 2*(LLRXY - LLRX - LLRY);
end

% convert to bits
% bias correct - make this optional?
ln2 = log(2);
psiterms = psi((Ntrl - (1:Nvarxy))/2) / 2;
dterm = (ln2 - log(Ntrl-1)) / 2;

if biascorrect
    HX = (HX - Nvarx*dterm - sum(psiterms(1:Nvarx)));
    HY = (HY - Nvary*dterm - sum(psiterms(1:Nvary)));
    HXY = (HXY - Nvarxy*dterm - sum(psiterms));
end

I = (HX + HY - HXY) / ln2;

function w = sum_norm_innerv(x, chC)
% sum of normalised innovations
%m = (chC')\x;
% avoid having to transpose the data
m = x/chC; 
w = -0.5 * sum(sum(m.*m));
