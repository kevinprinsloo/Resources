function chC = mychol(C)
% chol with correction for small eigenvals
try
    chC = chol(C);
catch 
    % not positive definite
    % try a small correction
    display('Warning: regularising non-positive definite covariance')
    chC = chol(C + 1e-15*eye(size(C,1)));
end