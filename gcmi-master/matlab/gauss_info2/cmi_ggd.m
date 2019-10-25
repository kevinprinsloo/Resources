function [I, LLR, Icond] = cmi_ggd(x, y, z, Zm, biascorrect, demeaned, dollr)
% out = cmi_ggd(x, y, z, Zm, demeaned, llr)
% Calculate CMI for 2 gausian variables conditioned on a discrete variable
% demeaned = 1 => x, y already have zero mean
% dollr = 1 => also calculate likelihood ratio test statistic
% Returns
% I - CMI (expectation over conditional infos)
% Icond - full vector (Ym) of conditional info values
if isvector(x)
    x = x(:);
end
if isvector(y)
    y = y(:);
end
Ntrl = size(x,1);
Nvarx = size(x,2);
Nvary = size(y,2);

if (size(y,1) ~= Ntrl) || (length(z) ~= Ntrl)
    error('cmi_ggd: number of trials do not match')
end

Icond = zeros(Zm,1);
LLRcond = zeros(Zm,1);
Pz = zeros(1,Zm);
for zi=1:Zm
    idx = z==(zi-1);
    Pz(zi) = sum(idx);
    tx = x(idx,:);
    ty = y(idx,:);
    [I, LLR] = info_gg(tx,ty,biascorrect,demeaned,dollr);
    Icond(zi) = I;
    if dollr
        LLRcond(zi) = LLR;
    end
end
Pz = Pz ./ Ntrl;

I = Pz*Icond;
if dollr
    LLR = sum(LLRcond);
end