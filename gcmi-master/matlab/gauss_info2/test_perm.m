

Nperm = 10000;
Ntrl = 1000;

x = randn(Ntrl,1);
y = randn(Ntrl,1);

disp('manual')
tic
Iman = zeros(1,Nperm);
for shi=1:Nperm
    Iman(shi) = info_gg(x,y(randperm(Ntrl)), true, false, false);
end
toc

disp('new func')
tic
Ip = info_perms_gg(x,y,Nperm, true, false);
toc