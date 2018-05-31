%%%%%%%%%%%%%%%%%%%
% DGP simulation
% Specifying error variance, size or type of estimator
%%%%%%%%%%%%%%%%%%%
clear all
betas = [0.9 1.5 1.6 1.7 1.5 0.4 0.3 0.2 0.1]; % DGP betas
M = 5; % number of replications
err_set = [0.01 1 2.25 7]; % values for error variance
n_set = [700 200 100 30]; % values for sample size
n_setflip = fliplr(n_set);

%------------------------------------------------
%------------------------------------------------
%------------------------------------------------

% Specify error variance and size
OLS
for k = 1:4
    for h = 1:4
        row = n_setflip(1,k);
        err = err_set(1,h);
        piematrix3(k,h) = DGPTFG(betas, M, row, err, 1); 
    end
end
% MMA
for k = 1:4
    for h = 1:4
        row = n_setflip(1,k);
        err = err_set(1,h);
        piematrix4(k,h) = DGPTFG(betas, M, row, err, 2); 
    end
end
% JMA
for k = 1:4
    for h = 1:4
        row = n_setflip(1,k);
        err = err_set(1,h);
        piematrix5(k,h) = DGPTFG(betas, M, row, err, 3); 
    end
end
% HRCp
for k = 1
    for h = 1:4
        row = n_set(1,k);
        err = err_set(1,h);
        piematrix6(1,h) = DGPTFG(betas, M, row, err, 4); 
    end
end
