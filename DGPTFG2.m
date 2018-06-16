function timesimulate = DGPTFG2(M, col, row, type)
%
% 
% DGP process function for my BSc Thesis 
%
% Outputs:
%    - piesimulate: returns the computation time of matrix for OLS, MMA, JMA and HRCp
%            depending on the sample size and number of regressors
%    - Y: nx1 dependent variable generated
%    - X: nxp matrix of regressors generated
%
% Inputs:
%    - col: p regressors
%    - row: sample size
%    - M: number of replications for PIE index
%    - t: type of estimator (OLS, MMA, JMA and HRCp in that order)

X = rand(row,col);
X(:,1) = ones(row,1);
Y = rand(row,1);
s=candidate(X,2);
time = cputime;
if type == 1
    betashat = regress(Y,X);
elseif type == 2
    betashat = gma(Y, X, 1, 2);
elseif type == 3
    betashat = gma(Y, X, 2, 2);
elseif type == 4
    betashat = ma(Y, X, s);
end
timesimulate = cputime-time;
% piesimselect = PIE_select(Y, X0, M, 1);
