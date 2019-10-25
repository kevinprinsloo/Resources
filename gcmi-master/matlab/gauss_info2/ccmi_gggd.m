function [I, LLR, Icond] = ccmi_gggd(x, y, z, k, Km, biascorrect, demeaned, dollr)
% out = ccmi_gggd(x, y, z, k, Km, demeaned, llr)
% Calculate CMI for 3 gausian variables conditioned on a further 
% discrete variable
% demeaned = 1 => x, y already have zero mean
% dollr = 1 => also calculate likelihood ratio test statistic
% Returns
% I - CCMI (expectation over conditional infos)
% LLR - loglikelihood ratio test statistic for ex
% Icond - full vector (Ym) of conditional CMI values
if isvector(x)
    x = x(:);
end
if isvector(y)
    y = y(:);
end
Ntrl = size(x,1);

if (size(y,1) ~= Ntrl) || (size(z,1) ~= Ntrl) || (length(k) ~= Ntrl)
    error('ccmi_gggd: number of trials do not match')
end

Icond = zeros(Km,1);
LLRcond = zeros(Km,1);
Pk = zeros(1,Km);
for ki=1:Km
    idx = k==(ki-1);
    Pk(ki) = sum(idx);
    tx = x(idx,:);
    ty = y(idx,:);
    tz = z(idx,:);
    % must demean again
    [I, LLR] = cmi_ggg(tx,ty,tz,biascorrect,false,dollr);
    Icond(ki) = I;
    if dollr
        LLRcond(ki) = LLR;
    end
end
Pk = Pk ./ Ntrl;

I = Pk*Icond;
if dollr
    LLR = sum(LLRcond);
end
