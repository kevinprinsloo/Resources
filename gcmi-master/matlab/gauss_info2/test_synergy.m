
% high number of trials for accurate measures
Ntrl = 10000;

% two random features
f1 = randn(Ntrl,1);
f2 = randn(Ntrl,1);

% a random linear response
a = 0.5;
R = f1 + f2 + a*randn(Ntrl,1);

% calculate information theoretic quantities
cf1 = copnorm(f1);
cf2 = copnorm(f2);
cR = copnorm(R);
IRf1 = info_gg(cf1,cR,true,false,false)
IRf2 = info_gg(cf2,cR,true,false,false)
If1f2 = info_gg(cf1,cf2,true,false,false)
IRf1f2 = info_gg([cf1 cf2], cR, true, false, false)
SYN = IRf1f2 - IRf1 - IRf2

% linear regression for comparison
ds = dataset(f1,f2,R);
mdl = fitlm(ds,'interactions')

%%
% [b,dev,stats] = glmfit([f1 f2 f1.*f2], R, 'normal')
mdl = fitglm(ds,'interactions')