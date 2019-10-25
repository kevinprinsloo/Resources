function [I, LLR] = info_vmd(x, y, Ym)
% out = info_vmd(x, y, Ym, demeaned, dollr)
% Calculate von Mises - Discrete information
% x - von Mises input (Ntrl x 1)
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

%TODO: default values for options (or put in wrapper
% convenience interface)

% build gaussian mixture model
gmm.G = Ym;
gmm.N = Nvar;
gmm.kap = zeros(1, Ym);
w = zeros(1,Ym);

Hcond = zeros(1,Ym);

x_group = zeros(size(x)); % group demeaned
for yi=1:Ym
    idx = y==(yi-1);
    tmp = x(idx,:);
    thsNtrl = size(tmp,1);
    w(yi) = thsNtrl ./ Ntrl;
    x_group(idx) = tmp - circ_mean(tmp);
end

% unconditional response
kap = circ_kappa(x);
I0 = besseli(0, kap);
I1 = besseli(1, kap);
Hunc = (-kap*(I1/I0)) + log(2*pi*I0);

% conditional response
kap = circ_kappa(x_group);
I0 = besseli(0, kap);
I1 = besseli(1, kap);
Hcond = (-kap*(I1/I0)) + log(2*pi*I0);

ln2 = log(2);

I = Hunc - Hcond;
I = I / ln2;

