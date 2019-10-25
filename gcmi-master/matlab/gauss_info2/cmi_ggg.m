function [I, LLR] = cmi_ggg(x, y, z, biascorrect, demeaned, dollr)
% out = cmi_ggg(x, y, z, demeaned, llr)
% Calculate CMI for three gausian variables
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
Ntrl = size(x,1);
Nvarx = size(x,2);
Nvary = size(y,2);
Nvarz = size(z,2);

if (size(y,1) ~= Ntrl) || (size(z,1) ~= Ntrl)
    error('cmi_ggg: number of trials do not match')
end

%TODO: default values for options (or put in wrapper
% convenience interface)
if demeaned
    gz.m = zeros(1, Nvarz);
    gxz.m = zeros(1, Nvarx+Nvarz);
    gyz.m = zeros(1, Nvary+Nvarz);
    gxyz.m = zeros(1, Nvarx+Nvary+Nvarz);
else
    gx.m = sum(x,1)/Ntrl;
    gy.m = sum(y,1)/Ntrl;
    gz.m = sum(z,1)/Ntrl;
    
    gxz.m = [gx.m gz.m];
    gyz.m = [gy.m gz.m];
    gxyz.m = [gx.m gy.m gz.m];

    x = bsxfun(@minus,x,gx.m);
    y = bsxfun(@minus,y,gy.m);
    z = bsxfun(@minus,z,gz.m);    
end

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
gz.C = C(zidx,zidx);

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
gz.chC = mychol(gz.C);
gxz.chC = mychol(gxz.C);
gyz.chC = mychol(gyz.C);
gxyz.chC = mychol(gxyz.C);

% normalisations cancel for information
% raw nats entropy so it can be used for likelihood
% normalisations not needed since they cancel 
% (both for entropy and likelihodd)
HZ = sum(log(diag(gz.chC))); % + 0.5*Nvarx*log(2*pi*exp(1));
HXZ = sum(log(diag(gxz.chC))); % + 0.5*(Nvarx+Nvarz)*log(2*pi*exp(1));
HYZ = sum(log(diag(gyz.chC))); % + 0.5*(Nvary+Nvarz)*log(2*pi*exp(1));
HXYZ = sum(log(diag(gxyz.chC))); % + 0.5*(Nvarx+Nvary+Nvarz)*log(2*pi*exp(1));

LLR = [];
% likelihood ratio test
if dollr
    LLRZ = sum_norm_innerv(z, gz.chC) - Ntrl*HZ; 
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
    HZ = (HZ - Nvarz*dterm - sum(psiterms(1:Nvarz)));
    HXZ = (HXZ - Nvarxz*dterm - sum(psiterms(1:Nvarxz)));
    HYZ = (HYZ - Nvaryz*dterm - sum(psiterms(1:Nvaryz)));
    HXYZ = (HXYZ - Nvarxyz*dterm - sum(psiterms));
end

I = (HXZ + HYZ - HXYZ - HZ) / ln2;

function w = sum_norm_innerv(x, chC)
% sum of normalised innovations
%m = (chC')\x;
% avoid having to transpose the data
m = x/chC; 
w = -0.5 * sum(sum(m.*m));