function I = info_cc_loop_bc(x, Cx, HX, y, Cy, HY, bc)
% out = info_cc_loop_nobc(x, y, biascorrect)
% Calculate Gaussian - Gaussian information for copula normalised data
% Hx/Hy =  precomputed marginal entropies IN BITS 
%          (NO BIAS CORRECTION)
% Cx/Cy =  precomputed covariance

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

C = zeros(Nvarxy,Nvarxy);

% fill in self-cov
C(1:Nvarx,1:Nvarx) = Cx;
C(ystart:Nvarxy,ystart:Nvarxy) = Cy;

% calc cross-cov for this info calculation
Cxy = (x'*y) / (Ntrl-1);
C(1:Nvarx,ystart:Nvarxy) = Cxy;
C(ystart:Nvarxy,1:Nvarx) = Cxy';
chC = chol(C);

HXY = sum(log(diag(chC))) + bc.entnorm;

% HXY = (HXY - Nvarxy*bc.dterm - bc.sumpsiterms);

I = HX + HY - (HXY / log(2));