function [PIE_value, P, I_P, I_P_Hat] = PIE_average(Y, X, M, type)
%
% This function computes the PIE index proposed in Yuan and Yang (2005), so
% we can determine whether we should stay with model selection or use
% instead a model averaging approach.
%
% In the function, the aim of returning P, I_P and I_P_Hat is to graph them
% by writing the command " plot(P, I_P, P, I_P_Hat) "
%
% type 1: OLS regression (min)
%   "  2: MMA averaging
%   "  3: JMA averaging
%   "  4: HRCp averaging
%
% X: matrix of dependent variables (you can also include an intercept)
% Y: dependent variable vector
% M: number of replications
%


s=candidate(X,2); % Generate 2^p candidate models 

% Possible Econometric Modelling
if type == 1
    betashat = regress(Y,X);
elseif type == 2
    betashat = gma(Y, X, 1, 2);
elseif type == 3
    betashat = gma(Y, X, 2, 2);
elseif type == 4
    betashat = ma(Y, X, s);
end

n = numel(Y); % count the number of observations

Yhat = X*betashat; % Fitted values vector

Sigmasqhat = sigmasqhat(Y, Yhat, X); % Error variance of original model


%% Core section of the code
for k = 1:21 % Index have to be positive integers ; P is computed then as a function of k
      for j = 1:M
    P = (k-1)*0.05; % given that P = 0:0.05:20
      
    W_pert = normrnd(0, (P^2)*Sigmasqhat,[n,1]); % Perturbation vector 
    Y_pert = Y + W_pert; % Perturbed dataset
    %Betas_Hat_Pert = regress(Y_Pert, X); % Betahat = (X'X)^(-1)X'Y slow procedure
    
  if type == 1
      Betashat_pert = regress(Y_pert,X);
    elseif type == 2
    Betashat_pert = gma(Y_pert, X, 1);
elseif type == 3
    Betashat_pert = gma(Y_pert, X, 2);
elseif type == 4
    Betashat_pert = ma(Y_pert, X, s);
    end
 
 Y_pert_hat = X*Betashat_pert; % fitted values for perturbed dataset

%%
    Sum_P = sum((Y_pert_hat - Yhat).^2); % sum of array elements (perturbed-original estimation)
    Sub_I_P(j) = ((Sum_P/n)^0.5)/(Sigmasqhat^0.5);
    
    Ag_Sub_I_P(k) = sum(Sub_I_P);
    I_P(k) = (1/M)*Ag_Sub_I_P(k);

      end
end
I_P = I_P'; % Transposing into column vector to compute OLS for PIE
P = 0:0.05:1;
P = P'; % Transposing into column vector to compute OLS for PIE
PIE_value = regress (I_P, P); % PIE Index from Yuan and Yang (2005)
I_P_Hat = P*PIE_value; % fitted value of I_P (useful for graph representation)

% figure PIE_plot
% plot(P, I_P, P, I_P_Hat)
