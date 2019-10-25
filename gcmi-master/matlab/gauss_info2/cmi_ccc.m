function [I, LLR] = cmi_ccc(x, y, z, Hmarg, biascorrect, dollr)
% out = cmi_ccc(x, y, z, Hmarg, biascorrect, llr)
% Calculate CMI for three Gaussian copula variables
% (transformed to normal gaussian marginals via empirical CDF)
% Hmarg is precomputed marginal entropy (depends only on number of trials)
% (IN NATS)
% biascorrect = analytic bias correction 
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
Ntrl = size(x,1);
Nvarx = size(x,2);
Nvary = size(y,2);
Nvarz = size(z,2);

if (size(y,1) ~= Ntrl) || (size(z,1) ~= Ntrl)
    error('cmi_ccc: number of trials do not match')
end

% HXZ = ent_g([x z], biascorrect, true);
% HYZ = ent_g([y z], biascorrect, true);
% HXYZ = ent_g([x y z], biascorrect, true);
% 
% gxz.m = zeros(1, Nvarx+Nvarz);
% gyz.m = zeros(1, Nvary+Nvarz);
% gxyz.m = zeros(1, Nvarx+Nvary+Nvarz);


% TODO: directly calculate covariance without requiring
% concatenation
xyz = [x y z];

% do covariance once
C = (xyz'*xyz) / (Ntrl - 1);
gxyz.C = C;
% extract right subregions
Nvaryz = Nvary + Nvarz;
Nvarxyz = Nvarx + Nvaryz;
zidx = (Nvarx + Nvary + 1):Nvarxyz;
idx = (Nvarx + 1):Nvarxyz;
gyz.C = C(idx, idx);

Nvarxz = Nvarx + Nvarz;
gxz.C = zeros(Nvarxz);
xidx = 1:Nvarx;
gxz.C(xidx,xidx) = gxyz.C(xidx,xidx);
zidxxz = (Nvarx+1):Nvarxz;
gxz.C(xidx,zidxxz) = gxyz.C(xidx,zidx);
gxz.C(zidxxz,xidx) = gxyz.C(zidx,xidx);
gxz.C(zidxxz,zidxxz) = gxyz.C(zidx,zidx);


% chol needed for both entropy calculation and likelihood
gxz.chC = chol(gxz.C);
gyz.chC = chol(gyz.C);
gxyz.chC = chol(gxyz.C);

% normalisations cancel for information
% raw nats entropy so it can be used for likelihood
% normalisations not needed since they cancel 
% (both for entropy and likelihodd)
k = 0.5*log(2*pi*exp(1));
HXZ = sum(log(diag(gxz.chC))) + k*Nvarxz;
HYZ = sum(log(diag(gyz.chC))) + k*(Nvaryz);
HXYZ = sum(log(diag(gxyz.chC))) + k*Nvarxyz;

LLR = [];
% likelihood ratio test
if dollr
%     LLRZ = sum_norm_innerv(z, gz.chC) - Ntrl*HZ; 
    LLRXZ = sum_norm_innerv([x z], gxz.chC) - Ntrl*HXZ;
    LLRYZ = sum_norm_innerv([y z], gyz.chC) - Ntrl*HYZ;
    LLRXYZ = sum_norm_innerv(xyz, gxyz.chC) - Ntrl*HXYZ; 
    LLR = 2*(LLRXYZ + LLRZ - LLRXZ - LLRYZ);
end


ln2 = log(2);
% bias correct - make this optional?
if biascorrect
    psiterms = psi((Ntrl - (1:Nvarxyz))/2) / 2;
    dterm = (ln2 - log(Ntrl-1)) / 2;
%     HZ = (HZ - Nvarz*dterm - sum(psiterms(1:Nvarz)));
    HXZ = (HXZ - Nvarxz*dterm - sum(psiterms(1:Nvarxz)));
    HYZ = (HYZ - Nvaryz*dterm - sum(psiterms(1:Nvaryz)));
    HXYZ = (HXYZ - Nvarxyz*dterm - sum(psiterms));
end

I = (HXZ + HYZ - HXYZ - Hmarg) / ln2;


function w = sum_norm_innerv(x, chC)
% sum of normalised innovations
%m = (chC')\x;
% avoid having to transpose the data
m = x/chC; 
w = -0.5 * sum(sum(m.*m));