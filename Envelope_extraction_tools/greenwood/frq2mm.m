function mm = frq2mm(frq)

% Greenwood's function for mapping frequency (frq) to place (mm) on the 
% basilar membrane

% Appropriate for measuring basilar membrane length in mm.
a = 0.06; 
k = 165.4;

mm = (1/a)*log10(frq/k+1);





