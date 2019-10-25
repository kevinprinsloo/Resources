function [I, LLR] = info_perms_gg(x, y, Nperm, biascorrect, demeaned)
% out = info_perms_gg(x, y, Nperm, biascorrect, demeaned)
% Calculate shuffled Gaussian - Gaussian information permutations
% NB, y is shuffled so if inputs have different dimensionality make
% y the smaller one
% demeaned = 1 => x, y already have zero mean
% llr = 1 => also calculate likelihood ratio test statistic
if isvector(x)
    x = x(:);
end
if isvector(y)
    y = y(:);
end
Ntrl = size(x,1);
Nvarx = size(x,2);
Nvary = size(y,2);

if size(y,1) ~= Ntrl
    error('info_gg: number of trials do not match')
end

%TODO: default values for options (or put in wrapper
% convenience interface)
if demeaned
    gx.m = zeros(1, Nvarx);
    gy.m = zeros(1, Nvary);
    gxy.m = zeros(1, Nvarx+Nvary);
else
    gx.m = sum(x,1)/Ntrl;
    gy.m = sum(y,1)/Ntrl;
    gxy.m = [gx.m gy.m];

    x = bsxfun(@minus,x,gx.m);
    y = bsxfun(@minus,y,gy.m);
end


% individual entropies
gx.C = (x'*x) / (Ntrl - 1);
gy.C = (y'*y) / (Ntrl - 1);

gx.chC = chol(gx.C);
gy.chC = chol(gy.C);

% normalisations cancel for information
HX = sum(log(diag(gx.chC))); % + 0.5*Nvarx*log(2*pi*exp(1));
HY = sum(log(diag(gy.chC))); % + 0.5*Nvary*log(2*pi*exp(1));

ln2 = log(2);
Nvarxy = Nvarx + Nvary;
psiterms = psi((Ntrl - (1:Nvarxy))/2) / 2;
dterm = (ln2 - log(Ntrl-1)) / 2;
if biascorrect
    HX = (HX - Nvarx*dterm - sum(psiterms(1:Nvarx)));
    HY = (HY - Nvary*dterm - sum(psiterms(1:Nvary)));
end

% permutations only affect joint entropy
HXY = zeros(1,Nperm);

% avoid repeated allocations from concatenation
ystart = Nvarx + 1;
xy = zeros(Ntrl, Nvarxy);
xy(:,1:Nvarx) = x;
for shi=1:Nperm
    xy(:,ystart:Nvarxy) = y(randperm(Ntrl),:);
    
    % do covariance once
    C = (xy'*xy) / (Ntrl - 1);
    chC = chol(C);

    % normalisations cancel for information
    HXY(shi) = sum(log(diag(chC))); % + 0.5*(Nvarx+Nvary)*log(2*pi*exp(1));
end

if biascorrect
    jointbias = Nvarxy*dterm + sum(psiterms);
    HXY = HXY - jointbias;
end

I = (HX + HY - HXY) / ln2;
