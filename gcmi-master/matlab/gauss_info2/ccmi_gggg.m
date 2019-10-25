function [I, LLR] = ccmi_gggg(x, y, z, k, biascorrect, demeaned, dollr)
% out = ccmi_gggg(x, y, z, k, biascorrect, demeaned, llr)
% Calculate CMI for three gausian variables, conditioned on a 4th
% demeaned = 1 => x, y, z already have zero mean
% dollr = 1 => also calculate likelihood ratio test statistic
if isvector(x)
    x = x(:);
end
if isvector(y)
    y = y(:);
end
if isvector(z)
    z = z(:);
end
if isvector(k)
    k = k(:);
end
Ntrl = size(x,1);
Nvarx = size(x,2);
Nvary = size(y,2);
Nvarz = size(z,2);
Nvark = size(z,2);

if (size(y,1) ~= Ntrl) || (size(z,1) ~= Ntrl) || (size(k,1) ~= Ntrl)
    error('ccmi_gggg: number of trials do not match')
end

%TODO: default values for options (or put in wrapper
% convenience interface)
if demeaned
    gzk.m = zeros(1, Nvarz+Nvark);
    gxzk.m = zeros(1, Nvarx+Nvarz+Nvark);
    gyzk.m = zeros(1, Nvary+Nvarz+Nvark);
    gxyzk.m = zeros(1, Nvarx+Nvary+Nvarz+Nvark);
else
    gx.m = sum(x,1)/Ntrl;
    gy.m = sum(y,1)/Ntrl;
    gz.m = sum(z,1)/Ntrl;
    gk.m = sum(k,1)/Ntrl;
    
    gzk.m = [gz.m gk.m];
    gxzk.m = [gx.m gz.m gk.m];
    gyzk.m = [gy.m gz.m gk.m];
    gxyzk.m = [gx.m gy.m gz.m gk.m];

    x = bsxfun(@minus,x,gx.m);
    y = bsxfun(@minus,y,gy.m);
    z = bsxfun(@minus,z,gz.m);    
    k = bsxfun(@minus,k,gk.m);
end

% TODO: directly calculate covariance without requiring
% concatenation
xyzk = [x y z k];

% do covariance once
C = (xyzk'*xyzk) / (Ntrl - 1);
gxyzk.C = C;
% extract right subregions
Nvarxyzk = Nvarx + Nvary + Nvarz + Nvark;
Nvarxy = Nvarx + Nvary;
Nvarzk = Nvarz + Nvark;
idx = (Nvarxy+1):Nvarxyzk;
gzk.C = gxyzk.C(idx,idx);

Nvaryzk = Nvary + Nvarzk;
idx = (Nvarx+1):Nvarxyzk;
gyzk.C = gxyzk.C(idx,idx);

Nvarxzk = Nvarx + Nvarzk;
gxzk.C = zeros(Nvarxzk, Nvarxzk);
xidx = 1:Nvarx;
gxzk.C = gxyzk.C(xidx,xidx);
idx = (Nvarxy+1):Nvarxyzk;
idx_t = (Nvarx+1):Nvarxzk;
gxzk.C(xidx, idx_t) = gxyzk.C(xidx, idx);
gxzk.C(idx_t, xidx) = gxyzk.C(idx, xidx);
gxzk.C(idx_t,idx_t) = gxyzk.C(idx, idx);


% chol needed for both entropy calculation and likelihood
gzk.chC = chol(gzk.C);
gxzk.chC = chol(gxzk.C);
gyzk.chC = chol(gyzk.C);
gxyzk.chC = chol(gxyzk.C);

% normalisations cancel for information
% raw nats entropy so it can be used for likelihood
% normalisations not needed since they cancel 
% (both for entropy and likelihodd)
HZK = sum(log(diag(gzk.chC))); % + 0.5*Nvarx*log(2*pi*exp(1));
HXZK = sum(log(diag(gxzk.chC))); % + 0.5*(Nvarx+Nvary)*log(2*pi*exp(1));
HYZK = sum(log(diag(gyzk.chC))); % + 0.5*(Nvarx+Nvary)*log(2*pi*exp(1));
HXYZK = sum(log(diag(gxyzk.chC))); % + 0.5*(Nvarx+Nvary)*log(2*pi*exp(1));

LLR = [];
% likelihood ratio test
if dollr
    LLRZK = sum_norm_innerv([z k], gzk.chC) - Ntrl*HZK; 
    LLRXZK = sum_norm_innerv([x z k], gxzk.chC) - Ntrl*HXZK;
    LLRYZK = sum_norm_innerv([y z k], gyzk.chC) - Ntrl*HYZK;
    LLRXYZK = sum_norm_innerv(xyzk, gxyzk.chC) - Ntrl*HXYZK; 
    LLR = 2*(LLRXYZK + LLRZK - LLRXZK - LLRYZK);
end

% convert to bits
% bias correct - make this optional?
ln2 = log(2);

if biascorrect
    psiterms = psi((Ntrl - (1:Nvarxyzk))/2) / 2;
    dterm = (ln2 - log(Ntrl-1)) / 2;
    HZK = (HZK - Nvarzk*dterm - sum(psiterms(1:Nvarzk)));
    HXZK = (HXZK - Nvarxzk*dterm - sum(psiterms(1:Nvarxzk)));
    HYZK = (HYZK - Nvaryzk*dterm - sum(psiterms(1:Nvaryzk)));
    HXYZK = (HXYZK - Nvarxyzk*dterm - sum(psiterms));
end

I = (HXZK + HYZK - HXYZK - HZK) / ln2;

function w = sum_norm_innerv(x, chC)
% sum of normalised innovations
%m = (chC')\x;
% avoid having to transpose the data
m = x/chC; 
w = -0.5 * sum(sum(m.*m));


