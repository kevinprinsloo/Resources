function [I, LLR] = info_vmd(x, y, Ym)
% out = info_vmd(x, y, Ym, demeaned, dollr)
% Calculate Chi2 - Discrete information
% x - Chi2 input (Ntrl x 1)
% y - discrete integer input (Ntrl in [0, Ym-1]

if isvector(x)
    x = x(:);
else
    error('info_vmd: Only univariate von Mises variables are supported')
end
if isvector(y)
    y = y(:);
else
    error('info_vmd: only univariate discrete variable supported');
end

Ntrl = size(x,1);
Nvar = size(x,2);

if size(y,1) ~= Ntrl
    error('info_gd: number of trials do not match');
end

% build gaussian mixture model
gmm.G = Ym;
gmm.N = Nvar;
gmm.k = zeros(1, Ym);
w = zeros(1,Ym);

Hcond = zeros(1,Ym);

for yi=1:Ym
    idx = y==(yi-1);
    tmp = x(idx,:);
    thsNtrl = size(tmp,1);
    w(yi) = thsNtrl ./ Ntrl;
    k = mean(tmp);
    gmm.k(yi) = k;
    k2 = k/2;
    Hcond(yi) = k2 + log(2*gamma(k2)) + (1-k2)*psi(k2);
end

% unconditional response
k = mean(x);
k2 = k/2;
Hunc = k2 + log(2*gamma(k2)) + (1-k2)*psi(k2);


ln2 = log(2);

I = Hunc - sum(w .* Hcond);
I = I / ln2;