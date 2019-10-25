function I = info_cc_loop2_bc(x, HX, y, HY, bc)
% out = info_cc_loop_bc(x, y, biascorrect)
% Calculate Gaussian - Gaussian information for copula normalised data
% Hx/Hy =  precomputed marginal entropies IN BITS 
%          (must match biascorrect option)
% biascorrect = structure of parameters from biasterms()
% ** NO NORMALISATION - USE ONLY WITH ENTROPYS FROM ENT_C_LOOP

if isvector(x)
    x = x(:);
end
if isvector(y)
    y = y(:);
end
Ntrl = size(x,1);
Nvarx = size(x,2);
Nvary = size(y,2);
ystart = Nvarx + 1;
Nvarxy = Nvarx + Nvary;

if size(y,1) ~= Ntrl
    error('info_cc: number of trials do not match')
end

xy = [x y];
C = (xy'*xy) / (Ntrl - 1);
chC = chol(C);

% NO NORMALISATION
HXY = sum(log(diag(chC)));% + bc.entnorm;

HXY = (HXY - bc.bias);

I = HX + HY - (HXY / log(2));