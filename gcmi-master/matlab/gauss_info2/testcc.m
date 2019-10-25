%%
Ntrl = 1000;
Nx = 1000;
xdim = 4;
ydim = 2;
X = randn(Ntrl,xdim,Nx);
Y = randn(Ntrl,ydim);

%%
% clear mex
Nthread = 1;
tic
for ri=1:1000
I = info_cc_slice_nobc_omp(X,xdim,Y,Ntrl,Nthread);
end
toc

%%
tic
Iold = zeros(1,Nx);
for i=1:Nx
    Iold(i) = info_gg(X(:,:,i),Y,false,true,false);
end
Iold;
toc

%%
for ri=1:1000
    B = t*t;
end

%%
Ntrl = 1000;
Nx = 1000;
xdim = 4;
ydim = 2;
X = randn(Ntrl,xdim,Nx);
Y = randn(Ntrl,ydim,Nx);

%%
% clear mex
Nthread = 8;
tic
% for ri=1:1000
I = info_cc_multi_nobc_omp(X,xdim,Y,ydim,Ntrl,Nthread);
% end
toc
%%
tic
Iold = zeros(1,Nx);
for i=1:Nx
    Iold(i) = info_gg(X(:,:,i),Y(:,:,i),false,true,false);
end
Iold;
toc

