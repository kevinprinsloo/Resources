function frq = mm2frq(mm)

% Greenwood's function for mapping place (mm) on the basilar membrane to 
% frequency (frq)

% Appropriate for measuring basilar membrane length in mm.
a = 0.06; 
k = 165.4;

frq = k*(10.^(a*mm)-1);




