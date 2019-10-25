
tic
for i=1:10000
[gx, gx, gxy] = fit_gg(x,y);
end
toc

%%
Ntrl = 1000;
x = copnorm(randn(Ntrl,2));
y = copnorm(randn(Ntrl,2));


%%
N = 10000;
tic
for i=1:N
    I = info_cc(x,y,true,false);
end
toc

%%
tic
[Hx Cx] = ent_g(x, true, true);
[Hy Cy] = ent_g(y, true, true);
bc = biasterms(size(x,2)+size(y,2), Ntrl);
for i=1:N
    I = info_cc_loop_bc(x, Cx, Hx, y, Cy, Hy, bc);
end
toc