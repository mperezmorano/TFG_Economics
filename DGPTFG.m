function piesimulate = DGPTFG(betas, M, row, err, t)
%
% 
% DGP process function for my BSc Thesis 
%
% Outputs:
%    - piesimulate: returns the PIE values matrix for OLS, MMA, JMA and HRCp
%            depending on error variance and sample size
%    - Y: nx1 dependent variable generated
%    - X: nxp matrix of regressors generated
%
% Inputs:
%    - betas: 1xp vector of coefficients
%    - err: error variance desired
%    - row: sample size
%    - M: number of replications for PIE index
%    - t: type of estimator (OLS, MMA, JMA and HRCp in that order)

[r,col] = size(betas);
X0 = -1 + 2*rand(row,col);
X0(:,1) = ones(row,1);

e = normrnd(0,err, [row,1]);
Y = X0*(betas)' + e;


piesimulate = PIE_average(Y, X0, M, t);
% piesimselect = PIE_select(Y, X0, M, 1);

