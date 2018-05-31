%--------------------------------------
% Alternative simulation
%--------------------------------------
clear all
clc

row_set = [700 200 100 30];
err_set = [0.01 1 2.25 7];
col = 8;
M = 100;

%OLS
for h = 1:4
    for k = 1:4
        row = row_set(1,h);
        err = err_set(1,k);
        [piematrix3(h,k),betas] = DGPTFG4(row, col, err, M, 1);
    end
end

%MMA
for h = 1:4
    for k = 1:4
        row = row_set(1,h);
        err = err_set(1,k);
        piematrix4(h,k) = DGPTFG4(row, col, err, M, 2);
    end
end

%JMA
for h = 1:4
    for k = 1:4
        row = row_set(1,h);
        err = err_set(1,k);
        piematrix5(h,k) = DGPTFG4(row, col, err, M, 3);
    end
end

%HRCp
 h = 1;
 for k = 1:4;
     row = row_set(1,h);
     err = err_set(1,k);
     pievector6(k) = DGPTFG4(row, col, err, M, 4);
 end


