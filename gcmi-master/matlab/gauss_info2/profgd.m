Ntrl = 3000;
Nx = 2;
x = copnorm(randn(Ntrl,Nx));

Cx1 = x'*x
Cx2 = test_dgemm(x,1)
%%
Nrep = 10000;
tic
for ri=1:Nrep
    Cx = x'*x;
end
toc

pause(0.5)
% tic
% for ri=1:Nrep
%     Cx = info_dg_slice_omp(x, 1, int16(zeros(2,Ntrl)),Ntrl,1);
% end
% toc

tic
% for ri=1:Nrep
%     Cx = test_dgemm(x);
% end
Cx = test_dgemm(x,Nrep);
toc

%%
sx = single(x);
tic
for ri=1:Nrep
    Cx = x'*x;
end
toc
tic
for ri=1:Nrep
    Cx = sx'*sx;
end
toc

%%
Ntrl = 1000;
Nx =1;
x = copnorm(randn(Ntrl,Nx));

Nrep = 10000;

tic 
for ri=1:Nrep
    E = sum(log(diag(chol(x'*x))));
end
toc

tic 
for ri=1:Nrep
    [q r] = qr(x,0);
    H = sum(log(abs(diag(r))));
end
toc


%% 
% Ntrl = 17000;
% Nx = 7620;
Ntrl = 500;
Nx = 1000;
Xm = 2;
Ny = 1;

X = int16(randi(Xm,[Ntrl,Nx]));
Y = randn(Ntrl,Ny);
Y(:,1) = Y(:,1) + double(X(:,1))/3;
cY = copnorm(Y)';
Yold = Y;
Xold = X-1;

%%
tic
cY = copnorm(Y)';
Inew = info_dc_slice_omp(X, Xm, cY, Ntrl, 64);
toc
%%


tic
I = zeros(1,Nx);
for i=1%:Nx
%     I(i) = info_cd(Yt, X(:,i), Xm, false);
    I(i) = info_gd(Yold, Xold(:,i), Xm, false, true, false);
end
toc

%% log reg AUC
tic
Ax = zeros(1,Nx);
for xi=1:100
    model = logregFit(Yold, Xold(:,xi));
    [yhat, p] = logregPredict(model, Yold);
    Az(xi) = rocarea_mex(p, Xold(:,xi));
end
toc

%%
tic
for xi=1:1000
    Az_m = rocarea(p, Xold(:,xi));
end
toc
Az_m

tic
for xi=1:1000
    Az_mex = rocarea_mex(p, Xold(:,xi));
end
toc
Az_mex
%%
figure
subplot(121)
hist(Az(2:end),30)
vline(Az(1), 'r-')
title('Az')
subplot(122)
hist(Inew(2:end),30)
vline(Inew(1), 'r-')
title('MI')



    