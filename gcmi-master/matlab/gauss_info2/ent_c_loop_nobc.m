function [H, C] = ent_c_loop_nobc(x)
% out = ent_g(x, biascorrect, demeaned)
% Calculate Gaussian Entropy
% demeaned = 1 => x, already has zero mean
if isvector(x)
    x = x(:);
end
Ntrl = size(x,1);
Nvarx = size(x,2);

% covariance
C = (x'*x) / (Ntrl - 1);
gx.C = C;

% chol needed for both entropy calculation and likelihood
gx.chC = chol(gx.C);

% WARNING - THIS IS ENTROPY WITHOUT NORMALISATION
HX = sum(log(diag(gx.chC))); % + bc.entnorm;

% convert to bits
% bias correct - make this optional?
ln2 = log(2);

HX = (HX);

H = HX / ln2;
