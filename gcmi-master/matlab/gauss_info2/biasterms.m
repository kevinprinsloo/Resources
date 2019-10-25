function bc = biasterms(Nvar,Ntrl)
% calculate Gaussian entropy bias correction terms for given
% Nvar and Ntrl
bc.psiterms = psi((Ntrl - (1:Nvar))/2) / 2;
bc.dterm = (log(2) - log(Ntrl-1)) / 2;
bc.bias = Nvar*bc.dterm + sum(bc.psiterms);
bc.entnorm = 0.5*(Nvar)*log(2*pi*exp(1));