function [piesimulate, betas, Y, X] = DGPTFG4(row, col, err, M, t)
%
%
%
%
%
syms z
X = -1 + 2*rand(row,col);
%X = rand(row,col);
X = [ones(row,1) X];
Y = rand(row,1);
%Y = normrnd(0.5,1,[row,1]);
Ysym = Y.*z;
betas = inv(X'*X)*X'*Ysym;
Yhat = X*betas;
%Yhat = aicregress(Ysym, X);
eqn1 = sigmasqhat(Ysym, Yhat, X) == err;
solvez = solve(eqn1);
solvez = double(solvez);
z = solvez(1,1);
Y = Y.*z;
%----------------------
piesimulate = PIE_average(Y,X,M,t);