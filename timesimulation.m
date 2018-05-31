%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

row_set = [30 100 400 1000 1500];
col_set = [3 5 8 10 12];
M = 50;
%%

% OLS
for j = 1:5
    for k = 1:5
        row = row_set(1,j);
        col = col_set(1,k);
timesimulate1(j,k) = DGPTFG2(M, col, row, 1);
    end
end

% MMA
for j = 1:5
    for k = 1:5
        row = row_set(1,j);
        col = col_set(1,k);
timesimulate2(j,k) = DGPTFG2(M, col, row, 2);
    end
end

%JMA
for j = 1:5
    for k = 1:5
        row = row_set(1,j);
        col = col_set(1,k);
timesimulate3(j,k) = DGPTFG2(M, col, row, 3);
    end
end

%HRCp
for j = 1:5
    for k = 1:5
        row = row_set(1,j);
        col = col_set(1,k);
timesimulate4(j,k) = DGPTFG2(M, col, row, 4);
    end
end