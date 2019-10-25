function [H, C] = ent_g(x, biascorrect, demeaned)
% out = ent_g(x, biascorrect, demeaned)
% Calculate Gaussian Entropy
% demeaned = 1 => x, already has zero mean
if isvector(x)
    x = x(:);
end
Ntrl = size(x,1);
Nvarx = size(x,2);

%TODO: default values for options (or put in wrapper
% convenience interface)
if demeaned
%     gx.m = zeros(1, Nvarx);
else
    gx.m = sum(x,1)/Ntrl;
    x = bsxfun(@minus,x,gx.m);
end


% covariance
C = (x'*x) / (Ntrl - 1);
gx.C = C;

% chol needed for both entropy calculation and likelihood
gx.chC = chol(gx.C);

HX = sum(log(diag(gx.chC))) + 0.5*Nvarx*log(2*pi*exp(1));

% convert to bits
% bias correct - make this optional?
ln2 = log(2);
psiterms = psi((Ntrl - (1:Nvarx))/2) / 2;
dterm = (ln2 - log(Ntrl-1)) / 2;

if biascorrect
    HX = (HX - Nvarx*dterm - sum(psiterms));
end

H = HX / ln2;
