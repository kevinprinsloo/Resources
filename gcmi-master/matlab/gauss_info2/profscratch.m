Ntrl = 9000;

x = randn(Ntrl,1);
sx = single(x);

Nrep = 5000;

tic
for i=1:Nrep
    cx = copnorm(x);
end
toc

tic
for i=1:Nrep
    cx = copnorm_c_double(x);
end
toc

%%
Ntrl = 9000;
Nx = 10000;

Nrep = 1;
x = randn(Ntrl,Nx);

% tic
% for ri=1:Nrep
%     cx = copnorm(x);
% end
% toc

tic
for ri=1:Nrep
    cx = zeros(size(x));
    for xi=1:Nx
        cx(:,xi) = copnorm_c_double(x(:,xi));
    end
end
toc

tic
Nthread = 12;
for ri=1:Nrep
    cx2 = copnorm_slice_omp_c_double(x,Nthread);
end
toc

%%
sx = single(x);

tic
for ri=1:Nrep
    cx = zeros(size(x));
    for xi=1:Nx
        cx(:,xi) = copnorm_c_float(sx(:,xi));
    end
end
toc

tic
Nthread = 12;
for ri=1:Nrep
    cx2 = copnorm_slice_omp_c_float(sx,Nthread);
end
toc

%%
Ntrl = 9000;
Nx = 1000;

Nrep = 1;
x = randn(Ntrl,2,Nx);

tic
for ri=1:Nrep
    cx = BLG_copnorm(x);
end
toc

%%
tic
for i=1:Nrep
    cx = copnorm_c_float(sx);
end
toc

%%

Ntrl=9000;
Nx=500;

% 1D
y1=rand(Ntrl,1);
x1=y*ones(1,Nx)+rand(Ntrl,Nx);

%
y2=rand(Ntrl,2);
x2=repmat(y2,[1 1 Nx])+rand([Ntrl 2 Nx]);

%%
Ntrl = 9000;
Nx = 1000;

x = randn(Ntrl,2,Nx);

cx = BLG_copnorm(x);

% just check first variable
max(cx(:,1,1) - copnorm(x(:,1,1)))

