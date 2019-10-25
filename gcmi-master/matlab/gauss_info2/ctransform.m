function [ ar ] = ctransform(a)
% Copula-transform array - rank and scale to [0, 1]
if ~isvector(a)
    error('ctransform: only vector inputs supported')
end
a = a(:);
[as ai] = sort(a, 1);
[aa ar] = sort(ai, 1);
%     ar = (ar - 1) / (size(ar, 2) - 1);
ar = ar / (size(ar, 1) + 1);
end

