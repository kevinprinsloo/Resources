function x = BLG_copnorm(x)
% COPNORM(x) - copula normalisation
% x (Ntrls, Nvars) - continuous data (TRIALS FIRST AXIS)
% x out = - "standard normal"-ized data

[~,x] = sort(x, 1);
[~,x] = sort(x, 1);
x = x / (size(x, 1) + 1);
x = -sqrt(2).*erfcinv(2*x);
